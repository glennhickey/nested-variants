#!/bin/bash

################################################################################
# surject.sh
#
# Description:
#   Runs vg surject to project GAM alignments onto reference paths, producing
#   a sorted BAM file.
#
# Usage:
#   surject.sh --gbz <file.gbz> --gam <file.gam> --ref <ref> --sample <name> \
#              --out-dir <dir> --out-name <name> \
#              [--cpus N] [--mem size] [--time HH:MM:SS] [--partition name] [--local]
#
################################################################################

set -e

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
JOB_NAME="surject"
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
        --local)
            LOCAL=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 --gbz <file.gbz> --gam <file.gam> --ref <ref> --sample <name> --out-dir <dir> --out-name <name> [options]"
            echo ""
            echo "Required Options:"
            echo "  --gbz <file>          GBZ graph file (with augmented reference paths)"
            echo "  --gam <file>          GAM alignment file"
            echo "  --ref <ref>           Augmented reference name (passed to -n)"
            echo "  --sample <name>       Sample name (passed to -N)"
            echo "  --out-dir <dir>       Output directory for BAM file"
            echo "  --out-name <name>     Output name for BAM file"
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
            echo "  $0 --gbz graph.gbz --gam sample.gam --ref augref_CHM13 --sample HG002 --out-dir ./output --out-name sample.bam"
            echo ""
            echo "  # Run locally"
            echo "  $0 --gbz graph.gbz --gam sample.gam --ref augref_CHM13 --sample HG002 --out-dir ./output --out-name sample.bam --local"
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

# Ensure OUTPUT_NAME has .bam extension
if [[ ! "$OUTPUT_NAME" =~ \.bam$ ]]; then
    OUTPUT_NAME="${OUTPUT_NAME}.bam"
fi

set -x

mkdir -p "$OUTPUT_DIR"

BAM="${OUTPUT_DIR}/${OUTPUT_NAME}"

# Build the command to run
# Use TMPDIR if set (node-local scratch on clusters), otherwise fall back to output dir
# Stage GBZ to local scratch for fast random I/O
# Use local scratch for samtools sort temp files and intermediate BAM
GBZ_BASE=$(basename "$GBZ")
BAM_BASE=$(basename "$BAM")
CMD="WORK_TMPDIR=\${TMPDIR:-${OUTPUT_DIR}} && \\
echo \"Staging GBZ to \${WORK_TMPDIR}\" && \\
cp \"${GBZ}\" \"\${WORK_TMPDIR}/${GBZ_BASE}\" && \\
/usr/bin/time -v vg surject -x \"\${WORK_TMPDIR}/${GBZ_BASE}\" -n ${REF} -N ${SAMPLE} -i -b -t ${CPUS} \"${GAM}\" | \\
/usr/bin/time -v samtools sort -@ ${CPUS} -T \"\${WORK_TMPDIR}/sort_${SAMPLE}\" -o \"\${WORK_TMPDIR}/${BAM_BASE}\" && \\
mv \"\${WORK_TMPDIR}/${BAM_BASE}\" \"${BAM}\" && \\
/usr/bin/time -v samtools index -@ ${CPUS} \"${BAM}\" && \\
rm -f \"\${WORK_TMPDIR}/${GBZ_BASE}\""

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
        --error="${OUTPUT_DIR}/${OUTPUT_NAME%.bam}.surject.log" \
        --wrap="$CMD"
fi
