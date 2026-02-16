# Yeast chrI Test

End-to-end test of the nested-variants pipeline on *S. cerevisiae* chromosome I
with two samples (SK1 and YPS128), exercising multi-sample processing and
VCF merging.

## Prerequisites

- vg (>= 1.56)
- samtools / bgzip / tabix
- bcftools (for `split_vcf`, `merge_call_vcfs`, `merge_dv_vcfs`)
- Docker (for `deepvariant`, optional)
- R + ggplot2, data.table, dplyr, scales (for `plots` / `call_plots` / `dv_plots` / `annotation_plots`, optional)
- bedtools (for `annotation_intersect`, optional)
- snakemake (>= 8)

All commands below are run from the **repository root** (`nested-variants/`).

```bash
source venv-nested-variants/bin/activate
```

## 1. Simulate paired-end reads

Simulate 10 000 read pairs (150 bp reads, ~500 bp fragments) from the SK1 and
YPS128 haplotypes directly from the VG file, then deinterleave into separate
R1/R2 files.

```bash
for SAMPLE in SK1 YPS128; do
  # Simulate interleaved paired-end FASTQ
  vg sim -x yeast-test/chrI.vg \
    -n 10000 -l 150 -p 500 -v 50 -q -m "$SAMPLE" -t 4 \
    > yeast-test/${SAMPLE}.interleaved.fq

  # Deinterleave into R1/R2 (every 8 lines = one read pair)
  paste - - - - - - - - < yeast-test/${SAMPLE}.interleaved.fq \
    | tee >(cut -f 1-4 | tr '\t' '\n' | gzip > yeast-test/${SAMPLE}_1.fq.gz) \
    | cut -f 5-8 | tr '\t' '\n' | gzip > yeast-test/${SAMPLE}_2.fq.gz

  rm yeast-test/${SAMPLE}.interleaved.fq

  # Create reads index (one FASTQ path per line, absolute paths)
  printf '%s\n' \
    "$(pwd)/yeast-test/${SAMPLE}_1.fq.gz" \
    "$(pwd)/yeast-test/${SAMPLE}_2.fq.gz" \
    > yeast-test/${SAMPLE}.reads.idx
done
```

## 2. Run the pipeline

```bash
# Graph construction + analysis (no genotyping)
snakemake --cores 4 graph_only \
  --config vg=yeast-test/chrI.vg ref=S288C \
           out_dir=yeast-test/output out_name=chrI.nested \
           mem_gb=4

# Graph construction + annotation overlap analysis
snakemake --cores 4 graph_only \
  --config vg=yeast-test/chrI.vg ref=S288C \
           out_dir=yeast-test/output out_name=chrI.nested \
           mem_gb=4 \
           annot_genes=yeast-test/fake-genes.bed \
           annot_repeats=yeast-test/fake-repeats.bed \
           annot_segdups=yeast-test/fake-segdups.bed \
           annot_censat=yeast-test/fake-censat.bed

# Full pipeline: genotype + DeepVariant + merge for both samples
snakemake --cores 4 all \
  --config vg=yeast-test/chrI.vg ref=S288C \
           out_dir=yeast-test/output out_name=chrI.nested \
           'samples={SK1: yeast-test/SK1.reads.idx, YPS128: yeast-test/YPS128.reads.idx}' \
           mem_gb=4 \
           annot_genes=yeast-test/fake-genes.bed \
           annot_repeats=yeast-test/fake-repeats.bed \
           annot_segdups=yeast-test/fake-segdups.bed \
           annot_censat=yeast-test/fake-censat.bed
```

The fake annotation BED files (`fake-genes.bed`, `fake-repeats.bed`,
`fake-segdups.bed`, `fake-censat.bed`) contain synthetic intervals placed to
overlap real segment coordinates, useful for testing the annotation overlap
pipeline. The repeats file uses 6-column BED with RepeatMasker-style classes
in column 6.

## Expected output

```
yeast-test/output/
#
# --- Graph construction & deconstruct ---
#
├── chrI.nested.gbz                        # augmented reference graph
├── chrI.nested.augref-segs.tsv            # augref segment table
├── chrI.nested.hapl                       # haplotype index
├── chrI.nested.fa.gz                      # augmented reference FASTA
├── chrI.nested.vcf.gz                     # deconstructed VCF
├── chrI.nested.augref-length-hist.png     # augref segment length histogram
├── chrI.nested.offref.png                 # offref density ideogram
├── chrI.nested.sites.vcf-stats.tsv        # deconstruct site-level stats
├── chrI.nested.sites.variant-types.png
├── chrI.nested.sites.size-dist.png
├── chrI.nested.sites.af-spectrum.png
├── chrI.nested.variants.vcf-stats.tsv     # deconstruct variant-level stats
├── chrI.nested.variants.variant-types.png
├── chrI.nested.variants.size-dist.png
├── chrI.nested.variants.af-spectrum.png
#
# --- Annotation overlap (when annot_* configured) ---
#
├── chrI.nested.annot-per-segment.tsv      # per-segment annotation overlap table
├── chrI.nested.annot-summary.png          # annotation overlap bar chart
├── chrI.nested.annot-scatter.png          # length vs overlap scatter
├── chrI.nested.annot-cooccur.png          # annotation co-occurrence heatmap
├── chrI.nested.annot-repeats.png          # repeat class breakdown (if repeats)
├── chrI.nested.annot-stats.tsv            # annotation overlap summary stats
├── chrI.nested.annot-snp-counts.all.png   # deconstruct SNP count by annotation
├── chrI.nested.annot-snp-tstv.all.png     # deconstruct Ts/Tv by annotation
#
# --- Per-sample call (SK1 shown; YPS128 identical) ---
#
├── SK1.vcf.gz                             # genotyped VCF (vg call)
├── SK1.call-offref.png                    # call off-reference density
├── SK1.call.{sites,variants}.{all,pass}.vcf-stats.tsv
├── SK1.call.{sites,variants}.{all,pass}.variant-types.png
├── SK1.call.{sites,variants}.{all,pass}.size-dist.png
├── SK1.annot-snp-counts.{all,pass}.png    # call SNP count by annotation
├── SK1.annot-snp-tstv.{all,pass}.png      # call Ts/Tv by annotation
#
# --- Per-sample DeepVariant (SK1 shown; YPS128 identical) ---
#
├── SK1.deepvariant.vcf.gz                 # DeepVariant VCF
├── SK1.dv-offref.png                      # DV off-reference density
├── SK1.dv.{sites,variants}.{all,pass}.vcf-stats.tsv
├── SK1.dv.{sites,variants}.{all,pass}.variant-types.png
├── SK1.dv.{sites,variants}.{all,pass}.size-dist.png
├── SK1.deepvariant.annot-snp-counts.{all,pass}.png
├── SK1.deepvariant.annot-snp-tstv.{all,pass}.png
#
# --- Merged call ---
#
├── merged.call.vcf.gz
├── merged.call-offref.png
├── merged.call.{sites,variants}.{all,pass}.vcf-stats.tsv
├── merged.call.{sites,variants}.{all,pass}.variant-types.png
├── merged.call.{sites,variants}.{all,pass}.size-dist.png
├── merged.call.{sites,variants}.{all,pass}.af-spectrum.png
├── merged.call.annot-snp-counts.{all,pass}.png
├── merged.call.annot-snp-tstv.{all,pass}.png
#
# --- Merged DeepVariant ---
#
├── merged.deepvariant.vcf.gz
├── merged.dv-offref.png
├── merged.dv.{sites,variants}.{all,pass}.vcf-stats.tsv
├── merged.dv.{sites,variants}.{all,pass}.variant-types.png
├── merged.dv.{sites,variants}.{all,pass}.size-dist.png
├── merged.dv.{sites,variants}.{all,pass}.af-spectrum.png
├── merged.deepvariant.annot-snp-counts.{all,pass}.png
└── merged.deepvariant.annot-snp-tstv.{all,pass}.png
```

## Clean up

```bash
rm -rf yeast-test/output yeast-test/*.fq.gz yeast-test/*.reads.idx
```
