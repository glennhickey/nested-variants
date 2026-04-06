#!/bin/bash

################################################################################
# freebayes.sh
#
# Description:
#   Runs FreeBayes in parallel to call variants from a BAM file against a
#   FASTA reference.  Uses GNU parallel to split by genomic regions.
#
# Usage:
#   freebayes.sh --bam <file.bam> --ref <file.fa.gz> --sample <name> \
#                --out-dir <dir> --out-name <name> \
#                [--region-size N] \
#                [--cpus N] [--mem size] [--time HH:MM:SS] [--partition name] [--local]
#
################################################################################

set -e

# Initialize variables
BAM=""
REF=""
SAMPLE=""
OUTPUT_DIR="."
OUTPUT_NAME=""
REGION_SIZE=100000

# SLURM resource defaults
CPUS="16"
MEM="128gb"
TIME="16:00:00"
PARTITION="long"
JOB_NAME="freebayes"
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
        --cpus)
            CPUS="$2"
            shift 2
            ;;
        --mem)
            MEM="$2"
            shift 2
            ;;
        --time)
            TIME="$2"
            shift 2
            ;;
        --partition)
            PARTITION="$2"
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
            echo "FreeBayes Options:"
            echo "  --region-size <N>     Region chunk size for parallelization (default: 100000)"
            echo ""
            echo "Execution Options:"
            echo "  --local               Run commands locally instead of via SLURM"
            echo ""
            echo "SLURM Resource Options (optional, with defaults):"
            echo "  --cpus <N>            CPUs per task (default: 16)"
            echo "  --mem <size>          Memory per job (default: 128gb)"
            echo "  --time <time>         Wall clock limit (default: 16:00:00)"
            echo "  --partition <p>       SLURM partition/queue (default: long)"
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

# FreeBayes requires uncompressed FASTA; decompress bgzipped ref if needed
if [[ "$REF" == *.gz ]]; then
    REF_PLAIN="${OUTPUT_DIR}/${OUTPUT_NAME%.vcf.gz}.ref.fa"
    echo "Decompressing reference to $REF_PLAIN"
    gunzip -c "$REF" > "$REF_PLAIN"
    samtools faidx "$REF_PLAIN"
    REF="$REF_PLAIN"
fi

# Generate regions from FAI
FAI="${REF}.fai"
if [ ! -f "$FAI" ]; then
    samtools faidx "$REF"
fi

REGIONS_FILE="${OUTPUT_DIR}/${OUTPUT_NAME%.vcf.gz}.regions.txt"
awk -v size="$REGION_SIZE" '{
    chrom = $1; len = $2; pos = 0
    while (pos < len) {
        end = pos + size
        if (end > len) end = len
        print chrom ":" pos "-" end
        pos = end
    }
}' "$FAI" > "$REGIONS_FILE"

echo "Generated $(wc -l < "$REGIONS_FILE") regions (${REGION_SIZE}bp chunks)"

# Run freebayes in parallel over regions, keep first header only,
# set FILTER=PASS on all records (FreeBayes outputs "." by default),
# then bgzip and index.
/usr/bin/time -v cat "$REGIONS_FILE" \
  | parallel -k -j "$CPUS" \
      freebayes -f "$REF" "$BAM" --region {} \
  | awk 'BEGIN{OFS="\t"; p=1} /^#/{if(p)print; if(/^#CHROM/)p=0; next} {if($7==".") $7="PASS"; print}' \
  | bgzip > "$VCF"
tabix -p vcf "$VCF"
