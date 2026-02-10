#!/bin/bash

################################################################################
# fasta.sh
#
# Description:
#   Extracts augmented reference paths from a GBZ file as a FASTA, then indexes
#   it with samtools faidx.
#
# Usage:
#   fasta.sh --gbz <file.gbz> --ref <ref> \
#            --out-dir <dir> --out-name <name> \
#            [--cpus N] [--mem size] [--time HH:MM:SS] [--partition name] [--local]
#
################################################################################

set -e

# Initialize variables
GBZ=""
REF=""
OUTPUT_DIR="."
OUTPUT_NAME=""

# SLURM resource defaults
CPUS="16"
MEM="128gb"
TIME="16:00:00"
PARTITION="long"
JOB_NAME="fasta"
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
        -h|--help)
            echo "Usage: $0 --gbz <file.gbz> --ref <ref> --out-dir <dir> --out-name <name> [options]"
            echo ""
            echo "Required Options:"
            echo "  --gbz <file>          GBZ graph file (with augmented reference paths)"
            echo "  --ref <ref>           Augmented reference sample name (e.g. augref_CHM13)"
            echo "  --out-dir <dir>       Output directory for FASTA file"
            echo "  --out-name <name>     Output name for FASTA file"
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
            echo "  $0 --gbz graph.gbz --ref augref_CHM13 --out-dir ./output --out-name graph.fa.gz"
            echo ""
            echo "  # Run locally"
            echo "  $0 --gbz graph.gbz --ref augref_CHM13 --out-dir ./output --out-name graph.fa.gz --local"
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

# Ensure OUTPUT_NAME has .fa.gz extension
if [[ ! "$OUTPUT_NAME" =~ \.fa\.gz$ ]]; then
    OUTPUT_NAME="${OUTPUT_NAME}.fa.gz"
fi

set -x

mkdir -p "$OUTPUT_DIR"

FASTA="${OUTPUT_DIR}/${OUTPUT_NAME}"

# Build the command to run
# vg paths: extract augmented reference paths as FASTA
# bgzip: block-gzip for indexed random access
# samtools faidx: create .fai + .gzi indexes
CMD="/usr/bin/time -v bash -c 'vg paths -x \"${GBZ}\" -S \"${REF}\" -F -t ${CPUS} | bgzip -@ ${CPUS} > \"${FASTA}\"' && \
/usr/bin/time -v samtools faidx \"${FASTA}\""

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
        --error="${OUTPUT_DIR}/${OUTPUT_NAME%.fa.gz}.fasta.log" \
        --wrap="$CMD"
fi
