#!/bin/bash

################################################################################
# slurm-paths.sh
#
# Description:
#   Runs vg paths in parallel on SLURM to compute augmented reference paths
#   from VG (variation graph) files. This script processes multiple VG files
#   (typically one per chromosome) concurrently using SLURM job scheduling,
#   computes augmented reference paths relative to a specified reference,
#   and outputs GFA format files.
#
# Usage:
#   slurm-paths.sh --vg <file1.vg> [--vg <file2.vg> ...] \
#                  --ref <ref> --out-dir <dir> --out-name <name> \
#                  [--min-augref-len N] [--cpus N] [--mem size] [--time HH:MM:SS] [--partition name]
#
################################################################################

set -e

# Initialize variables
VG_FILES=()
REF=""
OUTPUT_DIR="."
OUTPUT_NAME=""
MIN_AUGREF_LEN="50"

# SLURM resource defaults
CPUS="20"
MEM="200gb"
TIME="16:00:00"
PARTITION="long"
JOB_NAME="paths"
LOCAL=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --vg)
            # Expand glob patterns and add all matching files
            # Enable extended globbing for patterns like !(*.d9).vg
            shopt -s nullglob extglob
            for file in $2; do
                VG_FILES+=("$file")
            done
            shopt -u nullglob extglob
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
        --min-augref-len)
            MIN_AUGREF_LEN="$2"
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
            echo "Usage: $0 --vg <file.vg|pattern> [--vg <file.vg|pattern> ...] --ref <ref> --out-dir <dir> --out-name <name> [options]"
            echo ""
            echo "Required Options:"
            echo "  --vg <file>           VG file or glob pattern (e.g., chr*.vg) to process"
            echo "                        Can be specified multiple times"
            echo "  --ref <ref>           Reference name (passed to -Q)"
            echo "  --out-dir <dir>       Output directory for final merged GFA and GBZ files"
            echo "  --out-name <name>     Output name prefix (produces .gbz, .augref-segs.tsv)"
            echo ""
            echo "Augmented Reference Options:"
            echo "  --min-augref-len <N>  Minimum augref fragment length (default: 50)"
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
            echo "  $0 --vg 'chr*.vg' --ref GRCh38 --out-dir ./output --out-name merged"
            echo ""
            echo "  # With custom min-augref-len"
            echo "  $0 --vg 'chr*.vg' --ref GRCh38 --out-dir ./output --out-name merged --min-augref-len 100"
            echo ""
            echo "  # With custom SLURM resources"
            echo "  $0 --vg 'chr*.vg' --ref GRCh38 --out-dir ./output --out-name merged --cpus 16 --mem 100gb --time 8:00:00"
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
if [ ${#VG_FILES[@]} -eq 0 ]; then
    echo "Error: At least one --vg file must be specified"
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

# Strip common extensions to get a clean base name
OUTPUT_NAME="${OUTPUT_NAME%.gbz}"
OUTPUT_NAME="${OUTPUT_NAME%.gfa.gz}"
OUTPUT_NAME="${OUTPUT_NAME%.gfa}"

# Compute augmented reference sample name
AUGREF_SAMPLE="augref_${REF}"

set -x

mkdir -p "$OUTPUT_DIR"

# Create a temporary working directory for intermediate files
WORK_DIR="${OUTPUT_DIR}/tmp_paths_$$"
mkdir -p "$WORK_DIR"

# Track output files in order
GFA_FILES=()

# Count non-EBV VG files to split threads across parallel jobs
NUM_JOBS=0
for VG in "${VG_FILES[@]}"; do
    BASE=$(basename "$VG")
    [[ $BASE != "chrEBV.vg" ]] && ((NUM_JOBS++)) || true
done
THREADS_PER_JOB=$(( CPUS / (NUM_JOBS > 0 ? NUM_JOBS : 1) ))
(( THREADS_PER_JOB < 1 )) && THREADS_PER_JOB=1

# Process each VG file
PIDS=()
PID_NAMES=()
for VG in "${VG_FILES[@]}"; do
    BASE=$(basename "$VG")
    if [[ $BASE != "chrEBV.vg" ]]; then
        BASE=${BASE%.vg}
        GFA="${WORK_DIR}/${BASE}.augref.gfa"
        GFA_FILES+=("$GFA")

        SEGS="${WORK_DIR}/${BASE}.augref-segs.tsv"

        # Build the command to run (threads split across parallel jobs)
        CMD="/usr/bin/time -v vg paths -x \"$VG\" -Q ${REF} --compute-augref --min-augref-len ${MIN_AUGREF_LEN} --augref-sample ${AUGREF_SAMPLE} --augref-segs \"${SEGS}\" -t ${THREADS_PER_JOB} | /usr/bin/time -v vg convert -f - > \"${GFA}\""

        if $LOCAL; then
            # Run locally in background
            bash -c "$CMD" &
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
                --error="${WORK_DIR}/${OUTPUT_NAME}.${BASE}.log" \
                --wrap="$CMD" &
        fi
        PIDS+=($!)
        PID_NAMES+=("$BASE")
    fi
done

# Wait for each job individually and check exit codes
FAILED=()
for i in "${!PIDS[@]}"; do
    if ! wait "${PIDS[$i]}"; then
        FAILED+=("${PID_NAMES[$i]}")
        echo "ERROR: ${PID_NAMES[$i]} failed (PID ${PIDS[$i]})" >&2
    fi
done
if [ ${#FAILED[@]} -gt 0 ]; then
    echo "ERROR: ${#FAILED[@]} job(s) failed: ${FAILED[*]}" >&2
    exit 1
fi

# Merge per-chromosome augref segment tables into one file
SEGS_MERGED="${OUTPUT_DIR}/${OUTPUT_NAME}.augref-segs.tsv"
FIRST=true
for SEGS in "$WORK_DIR"/*.augref-segs.tsv; do
    if [ -f "$SEGS" ]; then
        if $FIRST; then
            cat "$SEGS"
            FIRST=false
        else
            # Skip header line from subsequent files
            tail -n +2 "$SEGS"
        fi
    fi
done > "$SEGS_MERGED"

# Verify all expected GFA files exist before merging
MISSING=()
for GFA in "${GFA_FILES[@]}"; do
    if [ ! -s "$GFA" ]; then
        MISSING+=("$(basename "$GFA")")
    fi
done
if [ ${#MISSING[@]} -gt 0 ]; then
    echo "ERROR: ${#MISSING[@]} GFA file(s) missing or empty: ${MISSING[*]}" >&2
    exit 1
fi

# Merge per-chromosome GFAs into single file (vg gbwt needs a real file)
WORK_TMPDIR="${TMPDIR:-${OUTPUT_DIR}}"
MERGED_GFA="${WORK_TMPDIR}/${OUTPUT_NAME}.gfa"
FIRST=true
for GFA in "${GFA_FILES[@]}"; do
    if $FIRST; then
        cat "$GFA"
        FIRST=false
    else
        tail -n +2 "$GFA"
    fi
done > "$MERGED_GFA"

# Clean up temporary working directory (per-chromosome GFAs no longer needed)
rm -rf "$WORK_DIR"

# Convert merged GFA to GBZ format
GBZ_OUTPUT="${OUTPUT_DIR}/${OUTPUT_NAME}.gbz"

if $LOCAL; then
    /usr/bin/time -v vg gbwt -G "${MERGED_GFA}" --gbz-format -g "${GBZ_OUTPUT}"
else
    CMD="/usr/bin/time -v vg gbwt -G \"${MERGED_GFA}\" --gbz-format -g \"${GBZ_OUTPUT}\""
    sbatch -W \
        --job-name="gbz" \
        --partition="${PARTITION}" \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task="${CPUS}" \
        --mem="${MEM}" \
        --time="${TIME}" \
        --output=/dev/null \
        --error="${OUTPUT_DIR}/${OUTPUT_NAME}.gbz.log" \
        --wrap="$CMD"
fi

# Clean up merged GFA
rm -f "${MERGED_GFA}"
