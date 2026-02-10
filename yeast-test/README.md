# Yeast chrI Test

End-to-end test of the nested-variants pipeline on *S. cerevisiae* chromosome I
with two samples (SK1 and YPS128), exercising multi-sample processing and
VCF merging.

## Prerequisites

- vg (>= 1.56)
- samtools / bgzip / tabix
- bcftools (for `split_vcf`, `merge_call_vcfs`, `merge_dv_vcfs`)
- Docker (for `deepvariant`, optional)
- R + ggplot2, data.table, dplyr, scales (for `plots` / `call_plots` / `dv_plots`, optional)
- snakemake (>= 8)

All commands below are run from the **repository root** (`nested-variants/`).

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

# Full pipeline: genotype + DeepVariant + merge for both samples
snakemake --cores 4 all \
  --config vg=yeast-test/chrI.vg ref=S288C \
           out_dir=yeast-test/output out_name=chrI.nested \
           'samples={SK1: yeast-test/SK1.reads.idx, YPS128: yeast-test/YPS128.reads.idx}' \
           mem_gb=4
```

## Expected output

```
yeast-test/output/
├── chrI.nested.gbz                   # augmented reference graph
├── chrI.nested.augref-segs.tsv       # augref segment table
├── chrI.nested.gfa.gz                # augmented GFA
├── chrI.nested.hapl                  # haplotype index
├── chrI.nested.vcf.gz                # deconstructed VCF
├── chrI.nested.onref.vcf.gz          # on-reference variants
├── chrI.nested.nestedref.vcf.gz      # nested-reference variants
├── chrI.nested.offref.vcf.gz         # off-reference variants
├── chrI.nested.augref-length-hist.png # augref segment length histogram
├── chrI.nested.offref.png            # offref density ideogram
├── chrI.nested.fa.gz                 # augmented reference FASTA (bgzipped)
├── chrI.nested.fa.gz.fai             # FASTA index
├── chrI.nested.fa.gz.gzi             # bgzip index
├── SK1.gam                           # SK1 read alignments
├── SK1.bam                           # SK1 surjected BAM
├── SK1.bam.bai                       # SK1 BAM index
├── SK1.pack                          # SK1 coverage pileup
├── SK1.vcf.gz                        # SK1 genotyped VCF (vg call)
├── SK1.deepvariant.vcf.gz            # SK1 DeepVariant VCF
├── SK1.call-offref.png               # SK1 call off-reference density
├── SK1.dv-offref.png                 # SK1 DeepVariant off-reference density
├── YPS128.gam                        # YPS128 read alignments
├── YPS128.bam                        # YPS128 surjected BAM
├── YPS128.bam.bai                    # YPS128 BAM index
├── YPS128.pack                       # YPS128 coverage pileup
├── YPS128.vcf.gz                     # YPS128 genotyped VCF (vg call)
├── YPS128.deepvariant.vcf.gz         # YPS128 DeepVariant VCF
├── YPS128.call-offref.png            # YPS128 call off-reference density
├── YPS128.dv-offref.png              # YPS128 DeepVariant off-reference density
├── merged.call.vcf.gz                # merged call VCFs (bcftools merge)
├── merged.deepvariant.vcf.gz         # merged DeepVariant VCFs (bcftools merge)
├── merged.call-offref.png            # merged call density
└── merged.dv-offref.png              # merged DV density
```

## Clean up

```bash
rm -rf yeast-test/output yeast-test/*.fq.gz yeast-test/*.reads.idx
```
