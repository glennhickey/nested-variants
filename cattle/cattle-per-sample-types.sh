#!/bin/bash
################################################################################
# cattle-per-sample-types.sh
#
# Description:
#   Extract per-sample variant counts from the cattle graph VCF, split by
#   variant type (SNP / MNP / Insertion / Deletion / SV Insertion / SV Deletion)
#   and by ref_context (On-reference on the main chromosome / Off-reference on
#   a _alt contig). Emits a TSV compatible with
#   scripts/per-sample-types-slide.R so the same plotting script that renders
#   the HPRC slide can also render the cattle slide.
#
#   Multi-allelic sites are split with `bcftools norm -m-` so each output
#   record is biallelic, keeping the SNP/Indel/SV classification unambiguous.
#
# Usage:
#   scripts/cattle/cattle-per-sample-types.sh <in.vcf.gz> <out.tsv>
#
# Output columns: sample \t variant_type \t ref_context \t count
################################################################################

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: $0 <in.vcf.gz> <out.tsv>" >&2
    exit 1
fi

VCF="$1"
OUT="$2"

# Grab sample list once, so awk can tag counts by sample name.
SAMPLES=$(bcftools query -l "$VCF" | paste -sd,)
if [ -z "$SAMPLES" ]; then
    echo "Error: no samples found in $VCF" >&2
    exit 1
fi
echo "[cattle-per-sample-types] samples: $SAMPLES" >&2

# Threshold at which small indels flip to SV Insertion / SV Deletion,
# matching the HPRC pipeline's default (50 bp).
SV_THRESH=50

bcftools norm -m- -Ou "$VCF" 2>/dev/null \
  | bcftools query -f '%CHROM\t%REF\t%ALT[\t%GT]\n' \
  | awk -v samples="$SAMPLES" -v sv_thresh="$SV_THRESH" '
BEGIN {
    n_samp = split(samples, s, ",")
}
{
    chrom = $1; ref = $2; alt = $3
    reflen = length(ref); altlen = length(alt)

    if (chrom ~ /_alt$/) ctx = "Off-reference"; else ctx = "On-reference"

    if (reflen == 1 && altlen == 1) {
        type = "SNP"
    } else if (reflen == altlen) {
        type = "MNP"
    } else {
        diff = altlen - reflen
        if (diff > 0) {
            type = (diff >= sv_thresh) ? "SV Insertion" : "Insertion"
        } else {
            type = (-diff >= sv_thresh) ? "SV Deletion" : "Deletion"
        }
    }

    for (i = 1; i <= n_samp; i++) {
        gt = $(3 + i)
        # Sample "has" the variant if any allele in the GT is the alt (1).
        # After `bcftools norm -m-` each record is biallelic, so any "1" in
        # the GT (accounting for phased / unphased / haploid) counts.
        if (gt ~ /1/) {
            key = s[i] SUBSEP type SUBSEP ctx
            counts[key]++
        }
    }
}
END {
    print "sample\tvariant_type\tref_context\tcount"
    for (k in counts) {
        n = split(k, parts, SUBSEP)
        print parts[1] "\t" parts[2] "\t" parts[3] "\t" counts[k]
    }
}' > "$OUT"

echo "[cattle-per-sample-types] wrote $OUT ($(wc -l < "$OUT") rows)" >&2
