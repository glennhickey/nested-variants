#!/usr/bin/env bash
# Subset a per-chromosome .vg to CHM13 + requested haplotypes and produce
# an augref-unified GFA suitable for Bandage.
#
# Optionally restrict to a CHM13 path range via vg snarls + vg chunk -S.
#
# Usage:
#   subset-cartoon.sh <chrom.vg> <out_prefix> [--range PATH:START-END] \
#       <sample.hap> [<sample.hap>...]
#
# Example (whole chromosome):
#   ./subset-cartoon.sh chr12.vg chr12.trio HG00097.1 HG00099.2
#
# Example (CHM13 range):
#   ./subset-cartoon.sh chr12.vg chr12.win1 \
#       --range CHM13#0#CHM13.0:574914-575776 HG00253.1 NA20806.1
#
# Outputs:
#   <prefix>.paths     PanSN path list fed to vg paths -p
#   <prefix>.gfa       subsetted graph (walks, CHM-rooted)
#   <prefix>.aug.gfa   augref-unified graph for Bandage
#   <prefix>.aug.tsv   augref segment table

set -euo pipefail

RANGE=""
if [ $# -ge 4 ] && [ "$3" = "--range" ]; then
    RANGE=$4
    VG=$1
    PREFIX=$2
    shift 4
elif [ $# -ge 3 ] && [ "$1" != "--help" ]; then
    VG=$1
    PREFIX=$2
    shift 2
else
    echo "Usage: $0 <chrom.vg> <out_prefix> [--range PATH:START-END] <sample.hap> [<sample.hap>...]" >&2
    echo "CHM13 is always included." >&2
    exit 1
fi

{
    echo "CHM13#0#CHM13.0"
    for hap in "$@"; do
        [[ "$hap" == CHM13* ]] && continue
        sample=${hap%.*}
        h=${hap##*.}
        echo "${sample}#${h}#${hap}#0"
    done
} > "${PREFIX}.paths"

# Subset to requested paths, clip CHM-edge tips, unchop linear runs, emit GFA
vg paths -v "$VG" -p "${PREFIX}.paths" -r \
    | vg clip -d 1 -P CHM - \
    | vg mod -u - \
    | vg convert -fW - > "${PREFIX}.sub.gfa"

# Rewrite sample P-line names to PanSN phase-0 form (so vg accepts them
# as reference paths) and declare them in the H-line RS:Z tag.
# HG00253#1#HG00253.1#0  ->  HG00253.1#0#HG00253.1
# CHM13#0#CHM13.0        ->  unchanged (already phase 0, 3 fields)
RS_SAMPLES="CHM13"
for hap in "$@"; do
    [[ "$hap" == CHM13* ]] && continue
    RS_SAMPLES="${RS_SAMPLES} ${hap}"
done

awk -v rs="$RS_SAMPLES" 'BEGIN{FS=OFS="\t"}
  $1=="H" {
    found=0
    for (i=1; i<=NF; i++) if ($i ~ /^RS:Z:/) { $i="RS:Z:" rs; found=1 }
    if (!found) { $(NF+1) = "RS:Z:" rs }
    print; next
  }
  $1=="P" {
    n = split($2, parts, "#")
    if (n >= 4) { $2 = parts[3] "#0#" parts[3] }
    print; next
  }
  { print }
' "${PREFIX}.sub.gfa" > "${PREFIX}.sub.gfa.tmp" && mv "${PREFIX}.sub.gfa.tmp" "${PREFIX}.sub.gfa"

if [ -n "$RANGE" ]; then
    vg convert -g "${PREFIX}.sub.gfa" -p > "${PREFIX}.sub.vg"
    vg snarls "${PREFIX}.sub.vg" > "${PREFIX}.snarls"
    vg chunk -x "${PREFIX}.sub.vg" -p "$RANGE" -S "${PREFIX}.snarls" \
        | vg convert -fW - > "${PREFIX}.gfa"
    rm -f "${PREFIX}.sub.vg"
else
    mv "${PREFIX}.sub.gfa" "${PREFIX}.gfa"
fi

vg paths -x "${PREFIX}.gfa" -u -Q CHM13 -N augref_CHM13 \
    --augref-segs "${PREFIX}.aug.tsv" \
    | vg convert -fW - > "${PREFIX}.aug.gfa"

rm -f "${PREFIX}.sub.gfa"
