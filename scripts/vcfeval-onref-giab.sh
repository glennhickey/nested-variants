#!/bin/bash
# vcfeval-onref-giab.sh — Stratify on-ref vcfeval TP/FP/FN by GIAB regions
#
# Usage: vcfeval-onref-giab.sh <tp.vcf.gz> <fp.vcf.gz> <fn.vcf.gz> \
#          <giab_beds> <giab_names> <output.tsv>
#
# giab_beds: comma-separated GIAB BED files
# giab_names: comma-separated display names (same order)
#
# Output: TSV with columns: giab_region, category, variant_type, count

set -euo pipefail

# Args: <tp_vcfs> <fp_vcfs> <fn_vcfs> <giab_beds> <giab_names> <output.tsv>
# VCF args are comma-separated lists of per-sample VCFs
TP_VCFS="$1"
FP_VCFS="$2"
FN_VCFS="$3"
GIAB_BEDS="$4"
GIAB_NAMES="$5"
OUTPUT="$6"

IFS=',' read -ra BEDS <<< "$GIAB_BEDS"
IFS=',' read -ra NAMES <<< "$GIAB_NAMES"

# Pantree-paper 3-way partition: Easy / Segdup / Hard. The Hard bucket is
# computed as "not in Easy AND not in Segdup" — absorbs the GIAB "Other
# Difficult" input plus any position off every BED. Build a combined
# Easy∪Segdup BED once so each category×filter pass can subtract it.
EASY_BED=""
SEGDUP_BED=""
for i in "${!BEDS[@]}"; do
    case "${NAMES[$i]}" in
        Easy)   EASY_BED="${BEDS[$i]}" ;;
        Segdup) SEGDUP_BED="${BEDS[$i]}" ;;
    esac
done
TMPD=$(mktemp -d)
trap 'rm -rf "$TMPD"' EXIT
UNION_BED="$TMPD/easy_or_segdup.bed"
if [ -n "$EASY_BED" ] && [ -n "$SEGDUP_BED" ]; then
    LC_ALL=C sort -m -k1,1 -k2,2n "$EASY_BED" "$SEGDUP_BED" \
      | bedtools merge > "$UNION_BED"
elif [ -n "$EASY_BED" ]; then
    cp "$EASY_BED" "$UNION_BED"
elif [ -n "$SEGDUP_BED" ]; then
    cp "$SEGDUP_BED" "$UNION_BED"
else
    : > "$UNION_BED"
fi

echo -e "giab_region\tcategory\tvariant_type\tcount" > "$OUTPUT"

emit_variants() {
    # Emit on-ref variants from per-sample VCFs as a sorted BED of
    # (chrom, pos-1, pos, variant_type) — used by both the per-BED
    # intersect pass and the Hard pass.
    local vcfs="$1"; local out="$2"
    IFS=',' read -ra VCF_ARR <<< "$vcfs"
    { bcftools concat -a "${VCF_ARR[@]}" 2>/dev/null | bcftools view -H 2>/dev/null \
      | awk -F'\t' '$1 !~ /_alt$/ {
            r=length($4); a=length($5);
            if(r==1 && a==1) t="SNP"; else t="Indel";
            chrom=$1; sub(/^augref_/, "", chrom);
            print chrom"\t"$2-1"\t"$2"\t"t
        }' \
      > "$out"; } || true
}

for CAT_LABEL in "TP:$TP_VCFS" "FP:$FP_VCFS" "FN:$FN_VCFS"; do
    CAT="${CAT_LABEL%%:*}"
    VCFS="${CAT_LABEL#*:}"
    VARS_BED="$TMPD/${CAT,,}.bed"
    emit_variants "$VCFS" "$VARS_BED"

    # Per-BED intersects (Easy, Segdup; "Other Difficult" gets absorbed
    # into Hard below, so skip it here if present).
    for i in "${!BEDS[@]}"; do
        NAME="${NAMES[$i]}"
        case "$NAME" in
            Easy|Segdup) ;;
            *)           continue ;;
        esac
        BED="${BEDS[$i]}"
        { bedtools intersect -a "$VARS_BED" -b "$BED" -u 2>/dev/null \
          | cut -f4 | sort | uniq -c \
          | awk -v region="$NAME" -v cat="$CAT" '{print region"\t"cat"\t"$2"\t"$1}' \
          >> "$OUTPUT"; } || true
    done

    # Hard = variants NOT in (Easy ∪ Segdup). bedtools intersect -v emits
    # the complement.
    { bedtools intersect -a "$VARS_BED" -b "$UNION_BED" -v 2>/dev/null \
      | cut -f4 | sort | uniq -c \
      | awk -v cat="$CAT" '{print "Hard\t"cat"\t"$2"\t"$1}' \
      >> "$OUTPUT"; } || true
done

echo "On-ref GIAB stratification written to $OUTPUT" >&2
