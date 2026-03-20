#!/bin/bash

################################################################################
# slurm-call.sh
#
# Description:
#   Runs vg pack followed by vg call to call variants from alignments.
#   Produces a VCF file.
#
# Usage:
#   slurm-call.sh --gbz <file.gbz> --gam <file.gam> --ref <ref> --sample <name> \
#                 --out-dir <dir> --out-name <name> \
#                 [--cpus N] [--mem size] [--time HH:MM:SS] [--partition name] [--local]
#
################################################################################

set -eo pipefail

# Initialize variables
GBZ=""
GAM=""
REF=""
SAMPLE=""
OUTPUT_DIR="."
OUTPUT_NAME=""
# SLURM resource defaults
CPUS="16"
MEM="128gb"
TIME="16:00:00"
PARTITION="long"
JOB_NAME="call"
LOCAL=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --gbz)
            GBZ="$2"
            shift 2
            ;;
        --gam)
            GAM="$2"
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
        --filter-contigs)
            echo "Warning: --filter-contigs is deprecated (filtering moved to separate rule)" >&2
            shift 2
            ;;
        --local)
            LOCAL=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 --gbz <file.gbz> --gam <file.gam> --ref <ref> --sample <name> --out-dir <dir> --out-name <name> [options]"
            echo ""
            echo "Required Options:"
            echo "  --gbz <file>          GBZ graph file"
            echo "  --gam <file>          GAM alignment file"
            echo "  --ref <ref>           Reference path prefix (passed to -P)"
            echo "  --sample <name>       Sample name (passed to -s)"
            echo "  --out-dir <dir>       Output directory for VCF file"
            echo "  --out-name <name>     Output name for VCF file"
            echo ""
            echo "Deprecated Options:"
            echo "  --filter-contigs <file>  (ignored, filtering moved to separate Snakemake rule)"
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
            echo "  $0 --gbz graph.gbz --gam sample.gam --ref GRCh38 --sample NA12878 --out-dir ./output --out-name sample.vcf.gz"
            echo ""
            echo "  # Run locally"
            echo "  $0 --gbz graph.gbz --gam sample.gam --ref GRCh38 --sample NA12878 --out-dir ./output --out-name sample.vcf.gz --local"
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

if [ -z "$GAM" ]; then
    echo "Error: --gam is required"
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
PACK="${OUTPUT_DIR}/${OUTPUT_NAME%.vcf.gz}.pack"

# Build the command to run
CMD="set -eo pipefail; /usr/bin/time -v vg pack -x \"${GBZ}\" -g \"${GAM}\" -Q 5 -o \"${PACK}\" -t ${CPUS} && \\
/usr/bin/time -v vg call \"${GBZ}\" -k \"${PACK}\" -z -a -A -S ${REF} -s ${SAMPLE} -t ${CPUS} | bgzip > \"${VCF}\" && \\
tabix -fp vcf \"${VCF}\""

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
        --error="${OUTPUT_DIR}/${OUTPUT_NAME%.vcf.gz}.call.log" \
        --wrap="$CMD"
fi
