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
#   Implementation: enumerate main chromosomes in the graph, then loop and run
#   `vg depth -P augref_<REF>#0#<chrom>_ ...` once per chrom. The trailing
#   underscore anchors the prefix so it matches `chr1_N_alt` but not `chr1` or
#   `chr10` — this skips the (~250 Mb each) main reference paths entirely,
#   which is what makes this tractable on HPRC scale.
#
# Requirements:
#   - A vg build whose `vg depth -P` honours haplotype-sense paths on GBZ
#     (older builds silently ignore haplotype traversals and you'll get
#     depth ≈ 1 for every segment).
#
# Usage:
#   augref-depth.sh --gbz <graph.gbz> --ref <REF_NAME> --out <out.tsv>
#                   [--vg <vg-binary>] [--threads N]
#
#   --ref is the reference sample name used as augref_<REF> prefix
#   (e.g. CHM13, GRCh38, S288C).
################################################################################

set -euo pipefail

GBZ=""
REF=""
OUT=""
VG="vg"
THREADS=8
BIN_SIZE=1000000000

while [[ $# -gt 0 ]]; do
    case $1 in
        --gbz) GBZ="$2"; shift 2 ;;
        --ref) REF="$2"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        --vg) VG="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
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
echo "[augref-depth] found $N_CHROMS main chromosomes; scanning _alt paths per chrom" >&2

# 2. Per-chromosome vg depth, appending to the combined output.
: > "$OUT"
for chrom in $CHROMS; do
    echo "[augref-depth] $chrom" >&2
    "$VG" depth \
        -P "${AUGREF_PREFIX}${chrom}_" \
        -b "$BIN_SIZE" \
        -t "$THREADS" \
        "$GBZ" >> "$OUT"
done

N_ROWS=$(wc -l < "$OUT")
echo "[augref-depth] done: $N_ROWS rows written to $OUT" >&2
