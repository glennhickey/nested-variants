#!/bin/bash

################################################################################
# haplotypes.sh
#
# Description:
#   Builds a haplotype subsampling index (.hapl) from a GBZ file.  This index
#   is used by vg giraffe for haplotype-aware read mapping.
#
#   Steps:
#     1. vg index -j  → distance index (--no-nested-distance to reduce memory)
#     2. vg gbwt  -r  → r-index
#     3. vg haplotypes -H → .hapl index
#
# Usage:
#   haplotypes.sh --gbz <file.gbz> --ref <ref> \
#                 --out-dir <dir> --out-name <name> \
#                 [--cpus N] [--mem size] [--time HH:MM:SS] [--partition name] [--local]
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
JOB_NAME="haplotypes"
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
            echo "  --gbz <file>          GBZ graph file"
            echo "  --ref <ref>           Reference name (original, e.g. GRCh38, CHM13, S288C)"
            echo "  --out-dir <dir>       Output directory for .hapl file"
            echo "  --out-name <name>     Output name for .hapl file"
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
            echo "  $0 --gbz graph.gbz --ref GRCh38 --out-dir ./output --out-name graph.hapl"
            echo ""
            echo "  # Run locally"
            echo "  $0 --gbz graph.gbz --ref GRCh38 --out-dir ./output --out-name graph.hapl --local"
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

# Ensure OUTPUT_NAME has .hapl extension
if [[ ! "$OUTPUT_NAME" =~ \.hapl$ ]]; then
    OUTPUT_NAME="${OUTPUT_NAME}.hapl"
fi

set -x

mkdir -p "$OUTPUT_DIR"

HAPL="${OUTPUT_DIR}/${OUTPUT_NAME}"
DIST="${OUTPUT_DIR}/${OUTPUT_NAME%.hapl}.dist.tmp"
RI="${OUTPUT_DIR}/${OUTPUT_NAME%.hapl}.ri.tmp"

# Build the command to run
# 1. Distance index (--no-nested-distance to reduce memory)
# 2. R-index from GBWT
# 3. Haplotype index from distance + r-index
# 4. Clean up intermediate files
CMD="/usr/bin/time -v vg index -t ${CPUS} -j \"${DIST}\" \"${GBZ}\" --no-nested-distance -P ${REF} && \
/usr/bin/time -v vg gbwt --num-threads ${CPUS} -r \"${RI}\" -Z \"${GBZ}\" && \
/usr/bin/time -v vg haplotypes -t ${CPUS} -H \"${HAPL}\" -d \"${DIST}\" -r \"${RI}\" \"${GBZ}\" && \
rm -f \"${DIST}\" \"${RI}\""

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
        --error="${OUTPUT_DIR}/${OUTPUT_NAME%.hapl}.haplotypes.log" \
        --wrap="$CMD"
fi
