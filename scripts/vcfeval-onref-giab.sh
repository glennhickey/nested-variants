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

echo -e "giab_region\tcategory\tvariant_type\tcount" > "$OUTPUT"

for i in "${!BEDS[@]}"; do
    BED="${BEDS[$i]}"
    NAME="${NAMES[$i]}"

    for CAT_LABEL in "TP:$TP_VCFS" "FP:$FP_VCFS" "FN:$FN_VCFS"; do
        CAT="${CAT_LABEL%%:*}"
        VCFS="${CAT_LABEL#*:}"

        # Concatenate per-sample VCFs, extract on-ref variants, intersect with GIAB
        IFS=',' read -ra VCF_ARR <<< "$VCFS"
        { bcftools concat -a "${VCF_ARR[@]}" 2>/dev/null | bcftools view -H 2>/dev/null \
          | awk -F'\t' '$1 !~ /_alt$/ {
                r=length($4); a=length($5);
                if(r==1 && a==1) t="SNP"; else t="Indel";
                print $1"\t"$2-1"\t"$2"\t"t
            }' \
          | bedtools intersect -a - -b "$BED" -u 2>/dev/null \
          | cut -f4 | sort | uniq -c \
          | awk -v region="$NAME" -v cat="$CAT" '{print region"\t"cat"\t"$2"\t"$1}' \
          >> "$OUTPUT"; } || true
    done
done

echo "On-ref GIAB stratification written to $OUTPUT" >&2
