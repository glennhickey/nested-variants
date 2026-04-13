#!/usr/bin/env bash
# Subset a per-chromosome .vg to CHM13 + requested haplotypes and produce
# an augref-unified GFA suitable for Bandage.
#
# Usage:
#   subset-cartoon.sh <chrom.vg> <out_prefix> <sample.hap> [<sample.hap>...]
#
# Example:
#   ./subset-cartoon.sh ../chr12.vg chr12.trio HG00097.1 HG00099.2
#
# Outputs:
#   <prefix>.paths     PanSN path list fed to vg paths -p
#   <prefix>.gfa       subsetted graph (walks, CHM-rooted)
#   <prefix>.aug.gfa   augref-unified graph for Bandage
#   <prefix>.aug.tsv   augref segment table

set -euo pipefail

if [ $# -lt 3 ]; then
    echo "Usage: $0 <chrom.vg> <out_prefix> <sample.hap> [<sample.hap>...]" >&2
    echo "CHM13 is always included." >&2
    exit 1
fi

VG=$1
PREFIX=$2
shift 2

{
    echo "CHM13#0#CHM13.0"
    for hap in "$@"; do
        [[ "$hap" == CHM13* ]] && continue
        sample=${hap%.*}
        h=${hap##*.}
        echo "${sample}#${h}#${hap}#0"
    done
} > "${PREFIX}.paths"

vg paths -v "$VG" -p "${PREFIX}.paths" -r \
    | vg clip -d 1 -P CHM - \
    | vg mod -u - \
    | vg convert -fW - > "${PREFIX}.gfa"

vg paths -x "${PREFIX}.gfa" -u -Q CHM13 -N augref_CHM13 \
    --augref-segs "${PREFIX}.aug.tsv" \
    | vg convert -fW - > "${PREFIX}.aug.gfa"
