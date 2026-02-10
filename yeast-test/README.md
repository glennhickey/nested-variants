# Yeast chrI Test

End-to-end test of the nested-variants pipeline on *S. cerevisiae* chromosome I.

## Prerequisites

- vg (>= 1.56)
- samtools / bgzip / tabix
- bcftools (for `make split-vcf`)
- Docker (for `make deepvariant`, optional)
- R + ggplot2, data.table, dplyr, scales (for `make plots` / `make call-plots` / `make dv-plots`, optional)

All commands below are run from the **repository root** (`nested-variants/`).

## 1. Simulate paired-end reads

Simulate 10 000 read pairs (150 bp reads, ~500 bp fragments) from the SK1
haplotypes directly from the VG file, then deinterleave into separate R1/R2
files.

```bash
# Simulate interleaved paired-end FASTQ from SK1 haplotypes
vg sim -x yeast-test/chrI.vg \
  -n 10000 -l 150 -p 500 -v 50 -q -m SK1 -t 4\
  > yeast-test/sim.interleaved.fq

# Deinterleave into R1/R2 (every 8 lines = one read pair)
paste - - - - - - - - < yeast-test/sim.interleaved.fq \
  | tee >(cut -f 1-4 | tr '\t' '\n' | gzip > yeast-test/sim_1.fq.gz) \
  | cut -f 5-8 | tr '\t' '\n' | gzip > yeast-test/sim_2.fq.gz

rm yeast-test/sim.interleaved.fq

# Create reads index (one FASTQ path per line, absolute paths)
printf '%s\n' \
  "$(pwd)/yeast-test/sim_1.fq.gz" \
  "$(pwd)/yeast-test/sim_2.fq.gz" \
  > yeast-test/sim.reads.idx
```

## 2. Run the pipeline

```bash
make paths haplotypes deconstruct split-vcf \
     length-hist plots \
     genotype surject fasta deepvariant call-plots dv-plots \
  VG=yeast-test/chrI.vg \
  REF=S288C \
  OUT_DIR=yeast-test/output \
  OUT_NAME=chrI.nested \
  SAMPLE=SK1 \
  READS=yeast-test/sim.reads.idx \
  HAPL=yeast-test/output/chrI.nested.hapl \
  CPUS=4 MEM=4gb
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
├── SK1.gam                     # read alignments
├── SK1.bam                     # surjected BAM
├── SK1.bam.bai                 # BAM index
├── SK1.pack                    # coverage pileup
├── SK1.vcf.gz                  # genotyped VCF (vg call)
├── SK1.deepvariant.vcf.gz      # DeepVariant VCF
├── SK1.call-offref.png         # call off-reference density
└── SK1.dv-offref.png           # DeepVariant off-reference density
```

## Clean up

```bash
rm -rf yeast-test/output yeast-test/sim*.fq.gz yeast-test/sim.reads.idx
```
