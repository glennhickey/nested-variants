#!/bin/bash

################################################################################
# bcftools-call.sh
#
# Description:
#   Runs bcftools mpileup + call (via Docker) in parallel to produce a per-sample
#   VCF from a BAM file against a FASTA reference. Uses GNU parallel to split by
#   genomic regions; each region runs in its own Docker container.
#
#   Short vs long read presets: by default bcftools mpileup uses the --config
#   illumina preset for short reads; pass --long-read to switch to pacbio-ccs
#   (HiFi). The preset affects quality/indel handling parameters only.
#
# Usage:
#   bcftools-call.sh --bam <file.bam> --ref <file.fa.gz> --sample <name> \
#                    --out-dir <dir> --out-name <name> \
#                    [--long-read] [--extra-args STR] \
#                    [--region-size N] [--docker IMAGE] \
#                    [--cpus N] [--mem size] [--local]
#
################################################################################

set -e

# Initialize variables
BAM=""
REF=""
SAMPLE=""
OUTPUT_DIR="."
OUTPUT_NAME=""
REGION_SIZE=1000000    # bcftools is less region-sensitive than freebayes; bigger chunks OK
LONG_READ=false
EXTRA_ARGS=""
DOCKER_IMAGE="staphb/bcftools:1.21"

# Resource defaults
CPUS="16"
MEM="64gb"
LOCAL=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --bam)
            BAM="$2"
            shift 2
            ;;
        --ref)
            REF="$2"
            shift 2
            ;;
        --sample)
            SAMPLE="$2"
            shift 2
            ;;
        --out-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --out-name)
            OUTPUT_NAME="$2"
            shift 2
            ;;
        --region-size)
            REGION_SIZE="$2"
            shift 2
            ;;
        --long-read)
            LONG_READ=true
            shift
            ;;
        --extra-args)
            EXTRA_ARGS="$2"
            shift 2
            ;;
        --docker)
            DOCKER_IMAGE="$2"
            shift 2
            ;;
        --cpus)
            CPUS="$2"
            shift 2
            ;;
        --mem)
            MEM="$2"
            shift 2
            ;;
        --local)
            LOCAL=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 --bam <file.bam> --ref <file.fa> --sample <name> --out-dir <dir> --out-name <name> [options]"
            echo ""
            echo "Required Options:"
            echo "  --bam <file>          Sorted BAM file"
            echo "  --ref <file>          FASTA reference file (bgzipped OK)"
            echo "  --sample <name>       Sample name"
            echo "  --out-dir <dir>       Output directory for VCF file"
            echo "  --out-name <name>     Output name for VCF file"
            echo ""
            echo "Caller Options:"
            echo "  --long-read           Use pacbio-ccs preset (default: illumina)"
            echo "  --extra-args <str>    Extra args appended to bcftools mpileup"
            echo "  --region-size <N>     Region chunk size for parallelization (default: 1000000)"
            echo "  --docker <image>      Docker image (default: staphb/bcftools:1.21)"
            echo ""
            echo "Resource Options:"
            echo "  --cpus <N>            CPUs per task (default: 16)"
            echo "  --mem <size>          Memory per job (default: 64gb)"
            echo "  --local               Unused (for symmetry with other callers)"
            exit 0
            ;;
        *)
            echo "Error: Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Validate required parameters
if [ -z "$BAM" ]; then
    echo "Error: --bam is required"
    exit 1
fi

if [ -z "$REF" ]; then
    echo "Error: --ref is required"
    exit 1
fi

if [ -z "$SAMPLE" ]; then
    echo "Error: --sample is required"
    exit 1
fi

if [ -z "$OUTPUT_NAME" ]; then
    echo "Error: --out-name is required"
    exit 1
fi

# Ensure OUTPUT_NAME has .vcf.gz extension
if [[ ! "$OUTPUT_NAME" =~ \.vcf\.gz$ ]]; then
    OUTPUT_NAME="${OUTPUT_NAME}.vcf.gz"
fi

set -x

mkdir -p "$OUTPUT_DIR"

VCF="${OUTPUT_DIR}/${OUTPUT_NAME}"

# bcftools mpileup accepts bgzipped FASTA but needs both .fai and .gzi indexes
if [[ "$REF" == *.gz ]]; then
    if [ ! -f "${REF}.fai" ] || [ ! -f "${REF}.gzi" ]; then
        samtools faidx "$REF"
    fi
else
    if [ ! -f "${REF}.fai" ]; then
        samtools faidx "$REF"
    fi
fi

# Generate regions from FAI, filtered to contigs present in the BAM.
# Same filtering as freebayes.sh: HPRC augref has 283k contigs, most empty.
FAI="${REF}.fai"
BAM_CONTIGS="${OUTPUT_DIR}/${OUTPUT_NAME%.vcf.gz}.bam-contigs.txt"
samtools idxstats "$BAM" | awk '$3 > 0 {print $1}' > "$BAM_CONTIGS"
echo "BAM has reads on $(wc -l < "$BAM_CONTIGS") of $(wc -l < "$FAI") contigs"

REGIONS_FILE="${OUTPUT_DIR}/${OUTPUT_NAME%.vcf.gz}.regions.txt"
awk -v size="$REGION_SIZE" 'NR==FNR{keep[$1]=1;next} ($1 in keep) {
    chrom = $1; len = $2; pos = 0
    while (pos < len) {
        end = pos + size
        if (end > len) end = len
        print chrom ":" (pos+1) "-" end
        pos = end
    }
}' "$BAM_CONTIGS" "$FAI" > "$REGIONS_FILE"
rm -f "$BAM_CONTIGS"

echo "Generated $(wc -l < "$REGIONS_FILE") regions (${REGION_SIZE}bp chunks)"

# Pick preset and FORMAT annotations based on read type.
# pacbio-ccs preset (bcftools 1.20+) tunes indel handling for HiFi.
# illumina preset uses indels-cns indel caller.
if [ "$LONG_READ" = "true" ]; then
    CONFIG_PRESET="pacbio-ccs"
    FORMAT_ANNOTS="FORMAT/AD,FORMAT/DP"
else
    CONFIG_PRESET="illumina"
    FORMAT_ANNOTS="FORMAT/AD,FORMAT/DP,FORMAT/SP"
fi

# Absolute paths for Docker bind mounts; dedup parent dirs so we don't
# pass redundant -v flags when REF/BAM live in the same directory.
REF_ABS="$(realpath "$REF")"
BAM_ABS="$(realpath "$BAM")"
declare -A MOUNT_SET
MOUNT_SET["$(dirname "$REF_ABS")"]=1
MOUNT_SET["$(dirname "$BAM_ABS")"]=1
MOUNT_FLAGS=""
for d in "${!MOUNT_SET[@]}"; do
    MOUNT_FLAGS="$MOUNT_FLAGS -v $d:$d"
done

echo "Running bcftools mpileup|call via Docker: $DOCKER_IMAGE (preset: $CONFIG_PRESET)"

# Run bcftools mpileup | call in parallel over regions (one docker container per
# region). Each region emits a full VCF (with header); awk keeps first header
# only and sets FILTER=PASS on all records (bcftools call emits "." by default).
/usr/bin/time -v cat "$REGIONS_FILE" \
  | parallel -k -j "$CPUS" \
      "docker run --rm --user $(id -u):$(id -g) $MOUNT_FLAGS $DOCKER_IMAGE bash -c \"bcftools mpileup -Ou --config $CONFIG_PRESET -a $FORMAT_ANNOTS -f '$REF_ABS' -r {} '$BAM_ABS' $EXTRA_ARGS | bcftools call -mv\"" \
  | awk 'BEGIN{OFS="\t"; p=1} /^#/{if(p)print; if(/^#CHROM/)p=0; next} {if($7==".") $7="PASS"; print}' \
  | bcftools reheader -s <(echo "$SAMPLE") \
  | bgzip > "$VCF"
tabix -p vcf "$VCF"
