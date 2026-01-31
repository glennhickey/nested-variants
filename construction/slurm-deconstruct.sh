#!/bin/bash

################################################################################
# slurm-deconstruct.sh
#
# Description:
#   Runs vg deconstruct to generate a VCF file from a GBZ (variation graph)
#   file. Extracts variant calls relative to a specified reference path prefix.
#
# Usage:
#   slurm-deconstruct.sh --gbz <file.gbz> --ref <ref> --out-dir <dir> --out-name <name> \
#                        [--cpus N] [--mem size] [--time HH:MM:SS] [--partition name] [--local]
#
################################################################################

set -e

# Initialize variables
GBZ=""
REF=""
OUTPUT_DIR="."
OUTPUT_NAME=""
CLUSTER=""
STAR_ALLELE=false

# SLURM resource defaults
CPUS="20"
MEM="200gb"
TIME="16:00:00"
PARTITION="long"
JOB_NAME="deconstruct"
LOCAL=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --gbz)
            GBZ="$2"
            shift 2
            ;;
        --ref)
            REF="$2"
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
        --cluster)
            CLUSTER="$2"
            shift 2
            ;;
        --star-allele)
            STAR_ALLELE=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 --gbz <file.gbz> --ref <ref> --out-dir <dir> --out-name <name> [options]"
            echo ""
            echo "Required Options:"
            echo "  --gbz <file>          GBZ file to deconstruct"
            echo "  --ref <ref>           Reference path prefix (passed to -P)"
            echo "  --out-dir <dir>       Output directory for VCF file"
            echo "  --out-name <name>     Output name for VCF file"
            echo ""
            echo "Deconstruct Options:"
            echo "  --cluster <F>         Cluster traversals with Jaccard >= F (passed to -L)"
            echo "  --star-allele         Use *-alleles for spanning haplotypes (passed to -R)"
            echo ""
            echo "Execution Options:"
            echo "  --local               Run commands locally instead of via SLURM"
            echo ""
            echo "SLURM Resource Options (optional, with defaults):"
            echo "  --cpus <N>            CPUs per task (default: 20)"
            echo "  --mem <size>          Memory per job (default: 200gb)"
            echo "  --time <time>         Wall clock limit (default: 16:00:00)"
            echo "  --partition <p>       SLURM partition/queue (default: long)"
            echo ""
            echo "Examples:"
            echo "  # Basic usage"
            echo "  $0 --gbz graph.gbz --ref GRCh38 --out-dir ./output --out-name result.vcf.gz"
            echo ""
            echo "  # Run locally"
            echo "  $0 --gbz graph.gbz --ref GRCh38 --out-dir ./output --out-name result.vcf.gz --local"
            echo ""
            echo "  # With custom SLURM resources"
            echo "  $0 --gbz graph.gbz --ref GRCh38 --out-dir ./output --out-name result.vcf.gz --cpus 16 --mem 100gb"
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
if [ -z "$GBZ" ]; then
    echo "Error: --gbz is required"
    exit 1
fi

if [ -z "$REF" ]; then
    echo "Error: --ref is required"
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

# Build the command to run
DECONSTRUCT_OPTS="-P ${REF} -a -t ${CPUS}"
if [ -n "$CLUSTER" ]; then
    DECONSTRUCT_OPTS="${DECONSTRUCT_OPTS} -L ${CLUSTER}"
fi
if $STAR_ALLELE; then
    DECONSTRUCT_OPTS="${DECONSTRUCT_OPTS} -R"
fi
CMD="vg deconstruct \"$GBZ\" ${DECONSTRUCT_OPTS} | bgzip > \"${VCF}\" && tabix -fp vcf \"${VCF}\""

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
        --wrap="$CMD"
fi
