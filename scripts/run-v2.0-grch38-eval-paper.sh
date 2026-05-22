#!/bin/bash
# run-v2.0-grch38-eval-paper.sh — paper-figure run for HPRC v2.0 MC GRCh38-EVAL.
#
# Same scope as run-v2.0-grch38-paper.sh (figures 1, 2, 6) PLUS short-read
# mapping + vg call genotyping + FreeBayes + call-vs-FB concordance:
#   3.call-summary.png        — vg call genotyping summary
#   4b.freebayes-summary.png  — FreeBayes summary + call-vs-FB vcfeval
#   5c.coverage-summary.png   — short-read coverage
#   5d.mapq-summary.png       — short-read MAPQ
#   10.freebayes-concordance-summary.png        — call-vs-FB chromsplit
#   10b.freebayes-concordance-onref-summary.png — call-vs-FB on-ref chromsplit
#
# Disabled: long-read samples, DeepVariant, PanGenie, bcftools caller.
#
# Usage:
#   bash scripts/run-v2.0-grch38-eval-paper.sh            # live run
#   bash scripts/run-v2.0-grch38-eval-paper.sh --dry-run  # dry-run

set -euo pipefail

DRY=""
if [ "${1:-}" = "--dry-run" ] || [ "${1:-}" = "-n" ]; then
    DRY="-n"
fi

GIAB=data/giab-reads

snakemake --profile profiles/slurm $DRY all \
    --default-resources \
        slurm_partition=high_priority \
        runtime=960 \
        mem_mb=32000 \
        tmpdir=/data/tmp \
    --config \
        ref=GRCh38 \
        vg='/private/home/ghickey/dev/work/hprc-v2.0-feb28/hprc-v2.0-mc-grch38-eval/hprc-v2.0-mc-grch38-eval.chroms/chr*.vg' \
        out_dir=output/v2.0-grch38-eval \
        out_name=hprc-v2.0-mc-grch38-eval.nested \
        min_augref_len=50 \
        annot_genes=data/hprc-v2-annotations/hprc-v2-genes-grch38-chm13.bed \
        annot_repeats=data/hprc-v2-annotations/hprc-v2-rm-grch38-chm13.bed \
        annot_segdups=data/hprc-v2-annotations/hprc-v2-sd-grch38-chm13.bed \
        annot_censat=data/hprc-v2-annotations/hprc-v2-censat-grch38-chm13.bed \
        annot_pclai=data/hprc-v2-annotations/hprc-v2-pclai-grch38-chm13.bed \
        giab_strat=data/hprc-v2-annotations/hprc-v2-giab \
        refgaps_bed=/private/home/ghickey/dev/work/hprc-v2.0-feb28/hprc-v2.0-mc-grch38-eval/hprc-v2.0-mc-grch38-eval.refgaps.bed \
        pantree_vcf=/private/home/ghickey/dev/work/pantree/GRCh38-464.MCv2.0.noY.vcf.gz \
        enable_deepvariant=false \
        enable_pangenie=false \
        enable_bcftools=false \
        "samples={HG001: $GIAB/HG001.novaseq.pcr-free.gs.paths, HG002: $GIAB/HG002.novaseq.pcr-free.gs.paths, HG003: $GIAB/HG003.novaseq.pcr-free.gs.paths, HG004: $GIAB/HG004.novaseq.pcr-free.gs.paths, HG005: $GIAB/HG005.novaseq.pcr-free.gs.paths, HG006: $GIAB/HG006.novaseq.pcr-free.gs.paths, HG007: $GIAB/HG007.novaseq.pcr-free.gs.paths, NA12891: $GIAB/NA12891.novaseq.pcr-free.gs.paths, NA12892: $GIAB/NA12892.novaseq.pcr-free.gs.paths}" \
        figure_titles=false \
        emit_pdf=true
