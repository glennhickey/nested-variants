# Yeast chrI Test

End-to-end test of the nested-variants pipeline on *S. cerevisiae* chromosome I
with two samples (SK1 and YPS128), exercising multi-sample processing and
VCF merging.

## Prerequisites

- vg (>= 1.56)
- samtools / bgzip / tabix
- bcftools (for `split_vcf`, `merge_call_vcfs`, `merge_dv_vcfs`)
- Docker (for `deepvariant`, optional)
- rtg (RTG Tools, for `vcfeval`)
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

# Simulate HiFi-like long reads (15 kb, single-end) from SK1
vg sim -x yeast-test/chrI.vg -n 1000 -l 15000 -q -m SK1 -t 4 \
  | gzip > yeast-test/SK1-hifi.fq.gz
echo "$(pwd)/yeast-test/SK1-hifi.fq.gz" > yeast-test/SK1-hifi.reads.idx
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

# Graph construction + pantree comparison (optional)
snakemake --cores 4 graph_only \
  --config vg=yeast-test/chrI.vg ref=S288C \
           out_dir=yeast-test/output out_name=chrI.nested \
           mem_gb=4 pantree_vcf=yeast-test/pantree-chrI.vcf

# Full pipeline: genotype + DeepVariant + merge + vcfeval + long reads
snakemake --cores 4 all \
  --config vg=yeast-test/chrI.vg ref=S288C \
           out_dir=yeast-test/output out_name=chrI.nested \
           'samples={SK1: yeast-test/SK1.reads.idx, YPS128: yeast-test/YPS128.reads.idx}' \
           'longread_samples={SK1-hifi: yeast-test/SK1-hifi.reads.idx}' \
           mem_gb=4 min_vcfeval_len=1000 min_surject_len=1000 \
           vcfeval_cpus=4 vcfeval_mem_gb=4 \
           annot_genes=yeast-test/fake-genes.bed \
           annot_repeats=yeast-test/fake-repeats.bed \
           annot_segdups=yeast-test/fake-segdups.bed \
           annot_censat=yeast-test/fake-censat.bed \
           pantree_vcf=yeast-test/pantree-chrI.vcf \
           giab_strat=yeast-test/fake-giab
```

The fake annotation BED files (`fake-genes.bed`, `fake-repeats.bed`,
`fake-segdups.bed`, `fake-censat.bed`) contain synthetic intervals placed to
overlap real segment coordinates, useful for testing the annotation overlap
pipeline. The repeats file uses 6-column BED with RepeatMasker-style classes
in column 6.

**Note on FreeBayes versions:** the default `freebayes_extra_args` is
`--max-coverage 500` (added to suppress pileup-blowup on HPRC-scale repeats).
FreeBayes ≥ 1.3 implements this correctly; FreeBayes **v1.0.2 has an inverted
`--max-coverage` sense** that drops *everything* on low-coverage data. If
you're testing locally with an older freebayes, override the config to bypass
the flag:

```bash
snakemake --cores 4 all --config ... freebayes_extra_args=''
```

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
├── chrI.nested.segment-polymorphism.tsv   # per-segment polymorphism table
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
├── merged.deepvariant.annot-snp-tstv.{all,pass}.png
#
# --- vcfeval: call vs DeepVariant comparison (per-sample + merged) ---
#
├── vcfeval/SK1/tp.vcf.gz                    # true positives (call perspective)
├── vcfeval/SK1/tp-baseline.vcf.gz           # true positives (baseline perspective)
├── vcfeval/SK1/fp.vcf.gz                    # false positives
├── vcfeval/SK1/fn.vcf.gz                    # false negatives
├── vcfeval/SK1/summary.txt                  # precision/recall summary
├── vcfeval/YPS128/...                       # same structure
├── merged.call-vs-dv.vcfeval-compare.tsv    # aggregated comparison table
├── merged.call-vs-dv.vcfeval-compare.png    # comparison bar chart
#
# --- Long-read sample (SK1-hifi, when longread_samples configured) ---
#
├── SK1-hifi.gam                               # long-read GAM (giraffe -b hifi)
├── SK1-hifi.bam                               # surjected BAM (-D long)
├── SK1-hifi.vcf.gz                            # genotyped VCF (vg call)
├── SK1-hifi.call-offref.png                   # call off-reference density
├── SK1-hifi.call.{sites,variants}.{all,pass}.vcf-stats.tsv
├── SK1-hifi.call.{sites,variants}.{all,pass}.variant-types.png
├── SK1-hifi.contig-depth.tsv                  # pack depth per contig
├── SK1-hifi.bam-depth.tsv                     # BAM depth per contig
├── SK1-hifi.gam-mapq.tsv                      # GAM MAPQ distribution
├── SK1-hifi.bam-mapq.tsv                      # BAM MAPQ distribution
├── merged.longread.call.vcf.gz                # merged long-read call VCF
├── merged.longread.call.sites.pass.vcf-stats.tsv
├── merged.longread.call.sites.pass.variant-types.png
├── merged.longread.call.sites.pass.call-summary-panel.png
├── 7.call-summary-longread.png                # long-read call summary figure
├── 8.deepvariant-summary-longread.png         # long-read DeepVariant summary
├── 9.concordance-summary-longread.png         # long-read concordance (off-ref)
├── 9b.concordance-onref-summary-longread.png  # long-read concordance (on-ref)
├── 9c.coverage-summary-longread.png           # long-read coverage summary
├── 9d.mapq-summary-longread.png               # long-read MAPQ summary
#
# --- Pantree comparison (when pantree_vcf configured) ---
#
├── pantree.records.tsv                      # extracted pantree records
├── pantree.variant-types.png                # pantree standalone variant types
├── pantree.density.png                      # pantree on-reference density ideogram
├── chrI.nested.records.tsv                  # our deconstruct records (for comparison)
├── chrI.nested.pantree-types.png            # side-by-side variant type counts
├── chrI.nested.pantree-types-pct.png        # side-by-side variant type percentages
├── chrI.nested.pantree-size-dist.png        # overlaid indel size distributions
├── chrI.nested.pantree-af.png               # overlaid AF spectra
└── chrI.nested.pantree-compare.tsv          # summary comparison table
```

## Clean up

```bash
rm -rf yeast-test/output yeast-test/*.fq.gz yeast-test/*.reads.idx
```
