#!/bin/bash
# run-v2.0-grch38-paper.sh — paper-figure graph_only run for HPRC v2.0 MC GRCh38.
# Builds figures 1, 2, and 6 (incl. pantree-compare-3panel) on high_priority
# SLURM partition with titles suppressed and PDFs emitted.
#
# Usage:
#   bash scripts/run-v2.0-grch38-paper.sh            # live run
#   bash scripts/run-v2.0-grch38-paper.sh --dry-run  # dry-run

set -euo pipefail

DRY=""
if [ "${1:-}" = "--dry-run" ] || [ "${1:-}" = "-n" ]; then
    DRY="-n"
fi

snakemake --profile profiles/slurm $DRY graph_only \
    --default-resources \
        slurm_partition=high_priority \
        runtime=960 \
        mem_mb=32000 \
        tmpdir=/data/tmp \
    --config \
        ref=GRCh38 \
        vg='/private/home/ghickey/dev/work/hprc-v2.0-feb28/hprc-v2.0-mc-grch38/hprc-v2.0-mc-grch38.chroms/chr*.vg' \
        out_dir=output/v2.0-grch38 \
        out_name=hprc-v2.0-mc-grch38.nested \
        min_augref_len=50 \
        annot_genes=data/hprc-v2-annotations/hprc-v2-genes-grch38-chm13.bed \
        annot_repeats=data/hprc-v2-annotations/hprc-v2-rm-grch38-chm13.bed \
        annot_segdups=data/hprc-v2-annotations/hprc-v2-sd-grch38-chm13.bed \
        annot_censat=data/hprc-v2-annotations/hprc-v2-censat-grch38-chm13.bed \
        annot_pclai=data/hprc-v2-annotations/hprc-v2-pclai-grch38-chm13.bed \
        giab_strat=data/hprc-v2-annotations/hprc-v2-giab \
        refgaps_bed=/private/home/ghickey/dev/work/hprc-v2.0-feb28/hprc-v2.0-mc-grch38/hprc-v2.0-mc-grch38.refgaps.bed \
        pantree_vcf=/private/home/ghickey/dev/work/pantree/GRCh38-464.MCv2.0.noY.vcf.gz \
        figure_titles=false \
        emit_pdf=true
