#!/bin/bash
# run-v2.1-chm13-paper.sh — paper-figure run for HPRC v2.1 MC CHM13-eval.
#
# Mirrors the README's "Run CHM13 pipeline" command and replaces it with a
# single-line wrapper so iteration on cluster runs is paste-free.  The
# heavy compute lives under output/v2.1-chm13/ from previous runs; this
# script keeps the same out_dir / out_name so finished outputs are reused
# and only the gaps fill in.
#
# Scope: graph_only + short-read mapping + vg call + FreeBayes + bcftools
# (the existing enabled set) + call-vs-FB concordance + pantree comparison.
# Disabled: long reads (no longread_samples=...), DeepVariant
# (enable_deepvariant=false), PanGenie (enable_pangenie=false).
#
# Usage:
#   bash scripts/run-v2.1-chm13-paper.sh            # live run
#   bash scripts/run-v2.1-chm13-paper.sh --dry-run  # dry-run

set -euo pipefail

DRY=""
if [ "${1:-}" = "--dry-run" ] || [ "${1:-}" = "-n" ]; then
    DRY="-n"
fi

ANNOT=data/hprc-v2-annotations
GIAB=data/giab-reads
GRAPH_DIR=/private/groups/hprc/hprc-graphs/hprc-v2.1-dec23/hprc-v2.1-mc-chm13-eval

snakemake --profile profiles/slurm $DRY all \
    --rerun-incomplete \
    --default-resources \
        slurm_partition=high_priority \
        runtime=960 \
        mem_mb=32000 \
        tmpdir=/data/tmp \
    --config \
        ref=CHM13 \
        vg="$GRAPH_DIR/hprc-v2.1-mc-chm13-eval.chroms/!(*.d*).vg" \
        out_dir=output/v2.1-chm13 \
        out_name=hprc-v2.1-mc-chm13.nested \
        min_surject_len=1000 \
        refgaps_bed="$GRAPH_DIR/hprc-v2.1-mc-chm13-eval.refgaps.bed" \
        paths_mem_gb=1024 \
        annot_genes=$ANNOT/hprc-v2-genes-grch38-chm13.bed \
        annot_repeats=$ANNOT/hprc-v2-rm-grch38-chm13.bed \
        annot_segdups=$ANNOT/hprc-v2-sd-grch38-chm13.bed \
        annot_censat=$ANNOT/hprc-v2-censat-grch38-chm13.bed \
        annot_pclai=$ANNOT/hprc-v2-pclai-grch38-chm13.bed \
        giab_strat=$ANNOT/hprc-v2-giab \
        pantree_vcf=/private/home/ghickey/dev/work/pantree/CHM13-464.MCv2.0.noY.vcf.gz \
        enable_deepvariant=false \
        enable_pangenie=false \
        "samples={HG001: $GIAB/HG001.novaseq.pcr-free.gs.paths, HG002: $GIAB/HG002.novaseq.pcr-free.gs.paths, HG003: $GIAB/HG003.novaseq.pcr-free.gs.paths, HG004: $GIAB/HG004.novaseq.pcr-free.gs.paths, HG005: $GIAB/HG005.novaseq.pcr-free.gs.paths, HG006: $GIAB/HG006.novaseq.pcr-free.gs.paths, HG007: $GIAB/HG007.novaseq.pcr-free.gs.paths, NA12891: $GIAB/NA12891.novaseq.pcr-free.gs.paths, NA12892: $GIAB/NA12892.novaseq.pcr-free.gs.paths}" \
        figure_titles=false \
        emit_pdf=true
