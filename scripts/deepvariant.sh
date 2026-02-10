#!/bin/bash

################################################################################
# deepvariant.sh
#
# Description:
#   Runs DeepVariant via Docker to call variants from a BAM file against a
#   FASTA reference.
#
# Usage:
#   deepvariant.sh --bam <file.bam> --ref <file.fa> --sample <name> \
#                  --out-dir <dir> --out-name <name> \
#                  [--dv-version VER] \
#                  [--cpus N] [--mem size] [--time HH:MM:SS] [--partition name] [--local]
#
################################################################################

set -e

# Initialize variables
BAM=""
REF=""
SAMPLE=""
OUTPUT_DIR="."
OUTPUT_NAME=""
DV_VERSION="1.9.0"

# SLURM resource defaults
CPUS="16"
MEM="128gb"
TIME="16:00:00"
PARTITION="long"
JOB_NAME="deepvariant"
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
        --dv-version)
            DV_VERSION="$2"
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
            echo "  --bam <file>          Sorted BAM file (from surject)"
            echo "  --ref <file>          FASTA reference file (with .fai index)"
            echo "  --sample <name>       Sample name for VCF output"
            echo "  --out-dir <dir>       Output directory for VCF file"
            echo "  --out-name <name>     Output name for VCF file"
            echo ""
            echo "DeepVariant Options:"
            echo "  --dv-version <ver>    DeepVariant Docker image version (default: 1.9.0)"
            echo ""
            echo "Execution Options:"
            echo "  --local               Run commands locally instead of via SLURM"
            echo ""
            echo "SLURM Resource Options (optional, with defaults):"
            echo "  --cpus <N>            CPUs per task (default: 16)"
            echo "  --mem <size>          Memory per job (default: 128gb)"
            echo "  --time <time>         Wall clock limit (default: 16:00:00)"
            echo "  --partition <p>       SLURM partition/queue (default: long)"
            echo ""
            echo "Examples:"
            echo "  # Basic usage"
            echo "  $0 --bam sample.bam --ref augref.fa --sample HG002 --out-dir ./output --out-name HG002.deepvariant.vcf.gz"
            echo ""
            echo "  # Run locally with specific version"
            echo "  $0 --bam sample.bam --ref augref.fa --sample HG002 --out-dir ./output --out-name HG002.deepvariant.vcf.gz --dv-version 1.9.0 --local"
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

# Resolve all paths to absolute for Docker bind mounts and container args
REF_ABS="$(cd "$(dirname "$REF")" && pwd)/$(basename "$REF")"
BAM_ABS="$(cd "$(dirname "$BAM")" && pwd)/$(basename "$BAM")"
OUT_ABS="$(cd "$OUTPUT_DIR" && pwd)"
VCF_ABS="${OUT_ABS}/${OUTPUT_NAME}"

# Build the command to run
CMD="/usr/bin/time -v docker run \
  --user \"$(id -u):$(id -g)\" \
  -v \"$(dirname "${REF_ABS}")\":\"$(dirname "${REF_ABS}")\" \
  -v \"$(dirname "${BAM_ABS}")\":\"$(dirname "${BAM_ABS}")\" \
  -v \"${OUT_ABS}\":\"${OUT_ABS}\" \
  google/deepvariant:${DV_VERSION} \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=\"${REF_ABS}\" \
  --reads=\"${BAM_ABS}\" \
  --output_vcf=\"${VCF_ABS}\" \
  --num_shards=${CPUS} \
  --sample_name=\"${SAMPLE}\" \
  --make_examples_extra_args=\"min_mapping_quality=0,keep_legacy_allele_counter_behavior=true,normalize_reads=true\""

if $LOCAL; then
    # Run locally
    bash -c "$CMD"
else
    # Submit SLURM job with resource requirements
    sbatch -W \
        --job-name="${JOB_NAME}" \
        --partition="${PARTITION}" \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task="${CPUS}" \
        --mem="${MEM}" \
        --time="${TIME}" \
        --output=/dev/null \
        --error="${OUTPUT_DIR}/${OUTPUT_NAME%.vcf.gz}.deepvariant.log" \
        --wrap="$CMD"
fi
