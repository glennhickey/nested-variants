#!/bin/bash

################################################################################
# augref-depth.sh
#
# Description:
#   Compute per-segment haplotype depth on augref _alt paths only.
#   For every node in every augref _alt path, counts how many graph paths
#   (reference + haplotype senses) traverse it. Output: one row per _alt path
#   in the format produced by `vg depth -P ... -b BIG`:
#
#     <path_name>  <bin_start>  <bin_end>  <mean_depth>  <stddev>
#
#   Implementation: enumerate main chromosomes in the graph, then run
#   `vg depth -P augref_<REF>#0#<chrom>_ ...` once per chrom (in parallel
#   when --parallel > 1). The trailing underscore anchors the prefix so it
#   matches `chr1_N_alt` but not `chr1` or `chr10` — this skips the
#   (~250 Mb each) main reference paths entirely, which is what makes this
#   tractable on HPRC scale.
#
#   Each vg depth invocation loads the GBZ, so total RAM ≈ (per-job GBZ
#   footprint) × --parallel. Size --parallel to fit your memory allocation.
#
# Requirements:
#   - A vg build whose `vg depth -P` honours haplotype-sense paths on GBZ
#     (older builds silently ignore haplotype traversals and you'll get
#     depth ≈ 1 for every segment).
#   - GNU parallel on PATH when --parallel > 1.
#
# Usage:
#   augref-depth.sh --gbz <graph.gbz> --ref <REF_NAME> --out <out.tsv>
#                   [--vg <vg-binary>] [--threads N] [--parallel N]
#                   [--bin-size N]
#
#   --ref      reference sample name (e.g. CHM13, GRCh38, S288C).
#   --threads  TOTAL threads across all parallel jobs; each vg depth is
#              given max(1, threads/parallel) threads.
#   --parallel Number of chromosomes to process concurrently (default 1).
################################################################################

set -euo pipefail

GBZ=""
REF=""
OUT=""
VG="vg"
THREADS=8
PARALLEL=1
BIN_SIZE=1000000000

while [[ $# -gt 0 ]]; do
    case $1 in
        --gbz) GBZ="$2"; shift 2 ;;
        --ref) REF="$2"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        --vg) VG="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --parallel) PARALLEL="$2"; shift 2 ;;
        --bin-size) BIN_SIZE="$2"; shift 2 ;;
        -h|--help)
            sed -n '3,/^######/p' "$0" | sed 's/^# //; s/^#//'
            exit 0 ;;
        *) echo "Error: unknown option: $1" >&2; exit 1 ;;
    esac
done

for required in GBZ REF OUT; do
    if [ -z "${!required}" ]; then
        echo "Error: --${required,,} is required" >&2
        exit 1
    fi
done

AUGREF_PREFIX="augref_${REF}#0#"

# Threads per vg depth invocation. vg depth on HPRC doesn't parallelize well
# internally (seen pinning at ~1 CPU), so giving each job 1 thread and running
# more in parallel is usually better, but we respect the supplied budget.
THREADS_PER_JOB=$(( THREADS / PARALLEL ))
if [ "$THREADS_PER_JOB" -lt 1 ]; then
    THREADS_PER_JOB=1
fi

# 1. Enumerate main chromosome names (reference-sense augref paths without
#    the _alt suffix).
echo "[augref-depth] listing main-chromosome paths from $GBZ" >&2
CHROMS=$("$VG" paths -x "$GBZ" --list -R \
    | grep -E "^${AUGREF_PREFIX}" \
    | grep -v "_alt\$" \
    | awk -F'#0#' '{print $2}')

if [ -z "$CHROMS" ]; then
    echo "Error: no reference paths found with prefix ${AUGREF_PREFIX} (non-_alt)" >&2
    exit 1
fi

N_CHROMS=$(echo "$CHROMS" | wc -l)
echo "[augref-depth] found $N_CHROMS main chromosomes; parallel=$PARALLEL, threads-per-job=$THREADS_PER_JOB" >&2

# 2. Per-chromosome vg depth.
if [ "$PARALLEL" -le 1 ]; then
    # Serial loop, append as we go.
    : > "$OUT"
    for chrom in $CHROMS; do
        echo "[augref-depth] $chrom" >&2
        "$VG" depth \
            -P "${AUGREF_PREFIX}${chrom}_" \
            -b "$BIN_SIZE" \
            -t "$THREADS_PER_JOB" \
            "$GBZ" >> "$OUT"
    done
else
    # Parallel: each chrom writes to its own temp TSV; concatenate in chrom
    # order afterward so output ordering is deterministic.
    WORK_TMPDIR=$(mktemp -d "${TMPDIR:-/tmp}/augref-depth.XXXXXX")
    trap 'rm -rf "$WORK_TMPDIR"' EXIT

    # Export everything the worker function needs.
    export VG AUGREF_PREFIX BIN_SIZE THREADS_PER_JOB GBZ WORK_TMPDIR
    do_chrom() {
        local chrom="$1"
        echo "[augref-depth] $chrom start" >&2
        "$VG" depth \
            -P "${AUGREF_PREFIX}${chrom}_" \
            -b "$BIN_SIZE" \
            -t "$THREADS_PER_JOB" \
            "$GBZ" > "$WORK_TMPDIR/${chrom}.tsv"
        echo "[augref-depth] $chrom done ($(wc -l < "$WORK_TMPDIR/${chrom}.tsv") rows)" >&2
    }
    export -f do_chrom

    # Feed chroms to GNU parallel on stdin.
    echo "$CHROMS" | parallel --will-cite -j "$PARALLEL" --line-buffer do_chrom {}

    # Concatenate per-chrom outputs in the order of the original chrom list.
    : > "$OUT"
    for chrom in $CHROMS; do
        f="$WORK_TMPDIR/${chrom}.tsv"
        if [ -f "$f" ]; then
            cat "$f" >> "$OUT"
        else
            echo "Warning: missing output for $chrom ($f)" >&2
        fi
    done
fi

N_ROWS=$(wc -l < "$OUT")
echo "[augref-depth] done: $N_ROWS rows written to $OUT" >&2
