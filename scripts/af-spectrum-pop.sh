#!/bin/bash
# af-spectrum-pop.sh — ad-hoc population-AF spectrum from a deconstruct VCF.
#
# Bypasses the default vcf-stats compute path's AF (which is bcftools
# +fill-tags's AC/AN where AN excludes missing genotypes, and so spikes at
# 0.5/1.0 for off-reference variants whose denominator is "haplotypes that
# carry the underlying _alt segment", not the full panel).
#
# Pipeline (re-shards by chromosome via scripts/bcftools-query-parallel.sh):
#   bcftools view <vcf>
#     | bcftools +setGT -t . -n 0       # missing → ref/ref
#     | bcftools +fill-tags -t AF        # AF = AC / (2 * N_samples_in_VCF)
#     | bcftools query -f '<chr>\t<af>\t<ref_context>'
# Then R re-bins and emits {prefix}.af-spectrum-pop.{tsv,png,pdf}.
#
# Usage:
#   bash scripts/af-spectrum-pop.sh <vcf> <out_prefix> [--af-step 0.01] [--parallel 32]
# Output: {out_prefix}.af-spectrum-pop.{tsv,png,pdf}

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: $0 <vcf> <out_prefix> [--af-step S] [--parallel N]"
    exit 1
fi

VCF="$1"
PREFIX="$2"
shift 2
AF_STEP=0.01
PAR=32
while [ $# -gt 0 ]; do
    case $1 in
        --af-step)  AF_STEP="$2"; shift 2 ;;
        --parallel) PAR="$2";     shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WRAPPER="$SCRIPT_DIR/bcftools-query-parallel.sh"
REBIN_R="$SCRIPT_DIR/af-spectrum-pop-bin.R"

AC_TSV="${PREFIX}.af-spectrum-pop.ac.tsv"

# Pull (CHROM, ref_context_marker via the _N_alt suffix, AF) from every
# record. AF here is computed AFTER setGT, so missing → ref and AN = 2N.
echo "[af-spectrum-pop] dumping per-record AF (this is the slow step)" >&2
bash "$WRAPPER" \
    --vcf "$VCF" \
    --format '%CHROM\t%INFO/AF\n' \
    --fill-tags "-t AF" \
    --missing-as-ref \
    --parallel "$PAR" \
  > "$AC_TSV"

echo "[af-spectrum-pop] dumped $(wc -l < "$AC_TSV") records to $AC_TSV" >&2
echo "[af-spectrum-pop] re-binning at af_step=$AF_STEP" >&2

Rscript "$REBIN_R" "$AC_TSV" "$PREFIX" --af-step "$AF_STEP" --pdf
