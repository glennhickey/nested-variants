#!/bin/bash

################################################################################
# pangenie.sh
#
# Description:
#   Runs PanGenie via Docker to genotype variants from FASTQ reads against
#   a panel VCF and reference FASTA.  All inputs are decompressed to scratch
#   since PanGenie requires uncompressed files.
#
# Usage:
#   pangenie.sh --reads <reads.idx> --ref <file.fa.gz> --vcf <panel.vcf.gz> \
#               --sample <name> --out-dir <dir> --out-name <name> \
#               [--docker IMAGE] [--cpus N] [--mem size] [--local]
#
################################################################################

set -e

# Initialize variables
READS=""
REF=""
VCF=""
SAMPLE=""
OUTPUT_DIR="."
OUTPUT_NAME=""
DOCKER_IMAGE="mgibio/pangenie:v4.2.1-bookworm"
JELLYFISH_SIZE="3000000000"

# SLURM resource defaults
CPUS="16"
MEM="128gb"
TIME="16:00:00"
PARTITION="long"
JOB_NAME="pangenie"
LOCAL=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --reads)
            READS="$2"
            shift 2
            ;;
        --ref)
            REF="$2"
            shift 2
            ;;
        --vcf)
            VCF="$2"
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
        --docker)
            DOCKER_IMAGE="$2"
            shift 2
            ;;
        --jellyfish-size)
            JELLYFISH_SIZE="$2"
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
            echo "Usage: $0 --reads <reads.idx> --ref <file.fa> --vcf <panel.vcf> --sample <name> --out-dir <dir> --out-name <name> [options]"
            echo ""
            echo "Required Options:"
            echo "  --reads <file>        Reads index file (one FASTQ path per line)"
            echo "  --ref <file>          FASTA reference (bgzipped OK, will be decompressed)"
            echo "  --vcf <file>          Panel VCF (bgzipped OK, will be decompressed)"
            echo "  --sample <name>       Sample name for output VCF"
            echo "  --out-dir <dir>       Output directory"
            echo "  --out-name <name>     Output VCF filename"
            echo ""
            echo "PanGenie Options:"
            echo "  --docker <image>      Docker image (default: mgibio/pangenie:v4.2.1-bookworm)"
            echo ""
            echo "Execution Options:"
            echo "  --local               Run commands locally instead of via SLURM"
            echo "  --cpus <N>            CPUs per task (default: 16)"
            echo "  --mem <size>          Memory per job (default: 128gb)"
            exit 0
            ;;
        *)
            echo "Error: Unknown option: $1"
            exit 1
            ;;
    esac
done

# Validate required parameters
for var in READS REF VCF SAMPLE OUTPUT_NAME; do
    if [ -z "${!var}" ]; then
        echo "Error: --$(echo $var | tr '[:upper:]' '[:lower:]') is required"
        exit 1
    fi
done

# Ensure OUTPUT_NAME has .vcf.gz extension
if [[ ! "$OUTPUT_NAME" =~ \.vcf\.gz$ ]]; then
    OUTPUT_NAME="${OUTPUT_NAME}.vcf.gz"
fi

set -x

mkdir -p "$OUTPUT_DIR"
OUT_VCF="${OUTPUT_DIR}/${OUTPUT_NAME}"

# Create scratch directory for uncompressed inputs
WORK_TMPDIR=$(mktemp -d "${TMPDIR:-${OUTPUT_DIR}}/pangenie.${SAMPLE}.XXXXXX")
trap '[ -n "${WORK_TMPDIR}" ] && rm -rf "${WORK_TMPDIR}"' EXIT

echo "PanGenie scratch: ${WORK_TMPDIR}"

# Decompress reference if needed
if [[ "$REF" == *.gz ]]; then
    echo "Decompressing reference"
    gunzip -c "$REF" > "${WORK_TMPDIR}/ref.fa"
else
    cp "$REF" "${WORK_TMPDIR}/ref.fa"
fi

# Decompress panel VCF if needed
if [[ "$VCF" == *.gz ]]; then
    echo "Decompressing panel VCF"
    gunzip -c "$VCF" > "${WORK_TMPDIR}/panel.vcf"
else
    cp "$VCF" "${WORK_TMPDIR}/panel.vcf"
fi

# Concatenate all FASTQ files from reads index into one uncompressed file
echo "Preparing reads"
> "${WORK_TMPDIR}/reads.fq"
while IFS= read -r fq || [ -n "$fq" ]; do
    case "$fq" in
        *.gz)
            gunzip -c "$fq" >> "${WORK_TMPDIR}/reads.fq"
            ;;
        *)
            cat "$fq" >> "${WORK_TMPDIR}/reads.fq"
            ;;
    esac
done < "$READS"

echo "Running PanGenie via Docker"
/usr/bin/time -v docker run \
    --user "$(id -u):$(id -g)" \
    -v "${WORK_TMPDIR}":"${WORK_TMPDIR}" \
    "${DOCKER_IMAGE}" \
    PanGenie \
    -i "${WORK_TMPDIR}/reads.fq" \
    -r "${WORK_TMPDIR}/ref.fa" \
    -v "${WORK_TMPDIR}/panel.vcf" \
    -o "${WORK_TMPDIR}/out" \
    -s "${SAMPLE}" \
    -t "${CPUS}" \
    -j "${CPUS}" \
    -e "${JELLYFISH_SIZE}"

# PanGenie outputs {prefix}_genotyping.vcf
PG_OUT="${WORK_TMPDIR}/out_genotyping.vcf"
if [ ! -f "$PG_OUT" ]; then
    echo "Error: PanGenie output not found: $PG_OUT"
    ls "${WORK_TMPDIR}"/out* 2>/dev/null
    exit 1
fi

# Set FILTER=PASS (PanGenie outputs "." by default), bgzip and index
awk 'BEGIN{OFS="\t"} /^#/{print;next} {if($7==".") $7="PASS"; print}' "$PG_OUT" \
    | bgzip > "$OUT_VCF"
tabix -p vcf "$OUT_VCF"

echo "PanGenie complete: $OUT_VCF"
