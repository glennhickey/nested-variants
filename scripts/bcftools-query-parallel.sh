#!/bin/bash

################################################################################
# bcftools-query-parallel.sh
#
# Parallel wrapper around `bcftools view | bcftools +fill-tags | bcftools query`.
# Shards the input VCF by main chromosome (grouping `_alt` contigs under their
# parent chrom) and runs one bcftools pipeline per shard concurrently. Emits the
# concatenated TSV output to stdout.
#
# Motivation: bcftools query is single-threaded. On a 26 GB bgzipped deconstruct
# VCF a single pass takes hours; parallelising by main chrom cuts wall clock by
# ~N on an N-core node. The HPRC graph has ~25 main chromosomes plus hundreds of
# thousands of short `_alt` contigs — grouping the alts under their main chrom
# keeps the shard count tractable while still giving each shard uniform work.
#
# Shard strategy:
#   - Contigs matching `<prefix>#0#chr<something>[_<N>_alt]` are grouped by the
#     `chr<something>` main chromosome (trailing `_N_alt` stripped).
#   - Contigs that don't match (no `#0#`) become their own shard.
#   - One BED file per shard with `contig<TAB>0<TAB>length` lines (from the VCF
#     header `##contig` records). bcftools view -R reads it once per shard.
#
# Usage:
#   bcftools-query-parallel.sh --vcf <vcf.gz> --format <bcftools-query-format>
#                              [--filter pass]      (apply `bcftools view -f PASS`)
#                              [--norm]             (apply `bcftools norm -m-` first)
#                              [--fill-tags STR]    (default: "-t AF")
#                              [--parallel N]       (default: 16)
#                              [--no-fill-tags]     (skip +fill-tags entirely)
#
# Output: concatenated bcftools query TSV on stdout (one line per record across
# all shards). Order follows the VCF header's contig declaration order. Stderr
# is used for progress; bcftools stderr is silenced.
################################################################################

set -euo pipefail

VCF=""
FORMAT=""
FILTER=""
NORM=false
FILL_TAGS="-t AF"
USE_FILL_TAGS=true
PARALLEL=16

while [[ $# -gt 0 ]]; do
    case $1 in
        --vcf)         VCF="$2"; shift 2 ;;
        --format)      FORMAT="$2"; shift 2 ;;
        --filter)      FILTER="$2"; shift 2 ;;
        --norm)        NORM=true; shift ;;
        --fill-tags)   FILL_TAGS="$2"; shift 2 ;;
        --no-fill-tags) USE_FILL_TAGS=false; shift ;;
        --parallel)    PARALLEL="$2"; shift 2 ;;
        -h|--help)     sed -n '3,/^######/p' "$0" | sed 's/^# //; s/^#//'; exit 0 ;;
        *) echo "Error: unknown option: $1" >&2; exit 1 ;;
    esac
done

[ -n "$VCF" ]    || { echo "Error: --vcf is required" >&2; exit 1; }
[ -n "$FORMAT" ] || { echo "Error: --format is required" >&2; exit 1; }
[ -r "$VCF" ]    || { echo "Error: cannot read VCF: $VCF" >&2; exit 1; }

WORK=$(mktemp -d "${TMPDIR:-/tmp}/bcfqp.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

# 1. Enumerate contigs from VCF header, group by main chrom.
#    augref_<REF>#0#<chrom>[_<N>_alt]  →  main-chrom = <chrom>
#    anything else: main-chrom = the contig itself
echo "[bcftools-query-parallel] enumerating contigs from $VCF" >&2
bcftools view -h "$VCF" 2>/dev/null \
  | awk '/^##contig=<ID=/' \
  | sed -E 's/^##contig=<ID=([^,>]+)(,length=([0-9]+))?.*$/\1\t\3/' \
  | awk -F'\t' -v OFS='\t' '
      {
        contig = $1
        length_bp = ($2 == "" ? 0 : $2)
        n = split(contig, a, "#0#")
        if (n < 2) {
          main = contig
        } else {
          tail = a[2]
          if (match(tail, /_[0-9]+_alt$/)) {
            main = substr(tail, 1, RSTART - 1)
          } else {
            main = tail
          }
        }
        # Record original header order so we can preserve it across shards.
        print NR, main, contig, length_bp
      }' > "$WORK/contig_map.tsv"

N_CONTIGS=$(wc -l < "$WORK/contig_map.tsv")
if [ "$N_CONTIGS" -eq 0 ]; then
    echo "Error: no ##contig records in VCF header" >&2
    exit 1
fi

# Shard order = first header-position of each main chrom. Preserves chromosome
# order for downstream consumers that care (most aggregations don't).
SHARD_ORDER=$(awk -F'\t' '!seen[$2]++ {print $2}' "$WORK/contig_map.tsv")
N_SHARDS=$(echo "$SHARD_ORDER" | wc -l)

# 2. Write per-shard BEDs: contig\t0\tlength for every contig in the shard.
mkdir -p "$WORK/beds"
awk -F'\t' -v OFS='\t' -v dir="$WORK/beds" '
  $4 > 0 { print $3, 0, $4 >> dir"/"$2".bed" }
' "$WORK/contig_map.tsv"

echo "[bcftools-query-parallel] $N_CONTIGS contigs → $N_SHARDS shards, parallel=$PARALLEL" >&2

# 3. Build shard-level pipeline.
FILTER_STAGE=""
if [ "$FILTER" = "pass" ]; then
    FILTER_STAGE="| bcftools view -f PASS 2>/dev/null"
fi

FILL_STAGE=""
if $USE_FILL_TAGS; then
    FILL_STAGE="| bcftools +fill-tags - -- $FILL_TAGS 2>/dev/null"
fi

if $NORM; then
    READ_STAGE='bcftools norm -m- -R "$BED" "$VCF" 2>/dev/null | bcftools view -c1 2>/dev/null'
else
    READ_STAGE='bcftools view -c1 -R "$BED" "$VCF" 2>/dev/null'
fi

# 4. Worker function: one shard → stdout (stdout captured by parallel).
process_shard() {
    local shard="$1"
    local BED="$WORK/beds/$shard.bed"
    if [ ! -s "$BED" ]; then
        return 0
    fi
    eval "$READ_STAGE $FILTER_STAGE $FILL_STAGE | bcftools query -f \"\$FORMAT\" 2>/dev/null"
}
export -f process_shard
export VCF FORMAT WORK READ_STAGE FILTER_STAGE FILL_STAGE

# 5. Parallel dispatch. `--line-buffer` streams each completed line as it
#    arrives without serialising shards behind each other (vs `-k`, which
#    made workers wait to drain stdout in FAI order — effectively single-
#    threaded under uneven shard difficulty). Downstream consumers (R
#    fread into a groupby aggregation, aggregator scripts) don't care
#    about line order.
if [ "$PARALLEL" -le 1 ]; then
    for shard in $SHARD_ORDER; do
        process_shard "$shard"
    done
else
    echo "$SHARD_ORDER" | parallel --will-cite -j "$PARALLEL" --line-buffer process_shard {}
fi

echo "[bcftools-query-parallel] done" >&2
