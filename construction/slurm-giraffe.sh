#!/bin/bash

################################################################################
# slurm-giraffe.sh
#
# Description:
#   Runs kmc followed by vg giraffe to align reads to a GBZ graph.
#   Produces a GAM file.
#
# Usage:
#   slurm-giraffe.sh --gbz <file.gbz> --hapl <file.hapl> --reads <file.fq.gz> \
#                    --sample <name> --out-dir <dir> --out-name <name> \
#                    [--cpus N] [--mem size] [--time HH:MM:SS] [--partition name] [--local]
#
################################################################################

set -e

# Initialize variables
GBZ=""
HAPL=""
READS=""
SAMPLE=""
OUTPUT_DIR="."
OUTPUT_NAME=""

# SLURM resource defaults
CPUS="16"
MEM="128gb"
TIME="16:00:00"
PARTITION="long"
JOB_NAME="giraffe"
LOCAL=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --gbz)
            GBZ="$2"
            shift 2
            ;;
        --hapl)
            HAPL="$2"
            shift 2
            ;;
        --reads)
            READS="$2"
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
            echo "Usage: $0 --gbz <file.gbz> --hapl <file.hapl> --reads <file.fq.gz> --sample <name> --out-dir <dir> --out-name <name> [options]"
            echo ""
            echo "Required Options:"
            echo "  --gbz <file>          GBZ graph file"
            echo "  --hapl <file>         Haplotype index file (.hapl)"
            echo "  --reads <file>        Input reads index (file containing list of fastq paths)"
            echo "  --sample <name>       Sample name"
            echo "  --out-dir <dir>       Output directory for GAM file"
            echo "  --out-name <name>     Output name for GAM file"
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
            echo "  $0 --gbz graph.gbz --hapl graph.hapl --reads sample.fq.gz --sample NA12878 --out-dir ./output --out-name sample.gam"
            echo ""
            echo "  # Run locally"
            echo "  $0 --gbz graph.gbz --hapl graph.hapl --reads sample.fq.gz --sample NA12878 --out-dir ./output --out-name sample.gam --local"
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

if [ -z "$HAPL" ]; then
    echo "Error: --hapl is required"
    exit 1
fi

if [ -z "$READS" ]; then
    echo "Error: --reads is required"
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

# Ensure OUTPUT_NAME has .gam extension
if [[ ! "$OUTPUT_NAME" =~ \.gam$ ]]; then
    OUTPUT_NAME="${OUTPUT_NAME}.gam"
fi

set -x

mkdir -p "$OUTPUT_DIR"

GAM="${OUTPUT_DIR}/${OUTPUT_NAME}"

# Extract numeric value from MEM for kmc (e.g., "128gb" -> "128")
MEM_NUM=$(echo "$MEM" | sed 's/[^0-9]//g')

# Read the two fastq paths from the index file and build -f arguments for giraffe
FASTQ_ARGS=""
while read -r fq; do
    FASTQ_ARGS="${FASTQ_ARGS} -f \"${fq}\""
done < "${READS}"

# Build the command to run
# Use TMPDIR if set, otherwise use output directory for temp files
CMD_TMPDIR="\${TMPDIR:-${OUTPUT_DIR}}"
CMD="WORK_TMPDIR=${CMD_TMPDIR} && \\
kmc -k29 -m${MEM_NUM} -okff -t${CPUS} -hp \"@${READS}\" \"\${WORK_TMPDIR}/${SAMPLE}\" \"\${WORK_TMPDIR}\" && \\
/usr/bin/time -v vg giraffe -p -t ${CPUS} -Z \"${GBZ}\" --haplotype-name \"${HAPL}\" --kff-name \"\${WORK_TMPDIR}/${SAMPLE}.kff\" \\
    --index-basename \"${OUTPUT_DIR}/${OUTPUT_NAME%.gam}.${SAMPLE}\" -N ${SAMPLE} ${FASTQ_ARGS} > \"${GAM}\" && \\
rm -f \"\${WORK_TMPDIR}/${SAMPLE}.kff\" \"\${WORK_TMPDIR}/${SAMPLE}.kff.kmc_pre\" \"\${WORK_TMPDIR}/${SAMPLE}.kff.kmc_suf\""

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
        --error="${OUTPUT_DIR}/${OUTPUT_NAME%.gam}.giraffe.log" \
        --wrap="$CMD"
fi
