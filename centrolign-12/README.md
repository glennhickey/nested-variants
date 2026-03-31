# centrolign-12

Augref + deconstruct + call + DeepVariant analysis of chr12 centromere
pangenome graphs. Uses scripts from `../scripts/` and a local copy of
`vcf-stats.R` (with a workaround for a `data.table::fread` segfault on
wide genotype tables >2 GB).

## Prerequisites

```bash
source ../venv-nested-variants/bin/activate
```

Requires: vg, bcftools, samtools, bgzip, tabix, R (ggplot2, data.table,
scales, gridExtra), Docker (for DeepVariant).

## Input graphs

Two VG graphs are available:

- **chr12.vg** — Full centrolign graph (373 haplotypes, 357 MB)
- **chr12.subsample150.100kb_flanks.fixed.vg** — Subsampled graph
  (150 haplotypes, 100 kb flanks). Built from
  `chr12.subsample150.100kb_flanks.gfa` with path name fixes applied
  via `fix-path-names.sed` to match the PanSN format of `chr12.vg`.

Both graphs use CHM13 as the reference.

### Rebuilding the fixed VG from GFA

```bash
sed -f fix-path-names.sed chr12.subsample150.100kb_flanks.gfa \
  | vg convert -g - > chr12.subsample150.100kb_flanks.fixed.vg
```

## Running

### Local

```bash
snakemake --cores 8 -np   # dry run
snakemake --cores 8        # full run (default: chr12.vg -> output/)
```

### Cluster (SLURM)

```bash
snakemake --profile ../profiles/slurm \
  --config vg=chr12.vg out=output name=chr12
```

### Subsampled graph

```bash
# Rebuild the fixed VG from GFA (if not already present):
sed -f fix-path-names.sed chr12.subsample150.100kb_flanks.gfa \
  | vg convert -g - > chr12.subsample150.100kb_flanks.fixed.vg

# Local:
snakemake --cores 8 \
  --config vg=chr12.subsample150.100kb_flanks.fixed.vg \
           out=output-flanks name=chr12-flanks use_hapl=False

# Cluster:
snakemake --profile ../profiles/slurm \
  --config vg=chr12.subsample150.100kb_flanks.fixed.vg \
           out=output-flanks name=chr12-flanks use_hapl=False
```

`use_hapl=False` skips the haplotype index (vg haplotypes fails on
graphs with loops) and runs giraffe in basic mode.

## Output figures

| Figure | Description |
|--------|-------------|
| `1.augref-summary.png` | Segment lengths (cumulative + log-log) + centromere density |
| `2.deconstruct-summary.png` | Variant types, size distribution, AF spectrum, per-sample |
| `3.call-summary.png` | vg call genotyping (PASS): variant types, per-sample, summary |
| `4.deepvariant-summary.png` | DeepVariant (PASS): variant types, per-sample |

## Samples

Three HiFi long-read samples are used for genotyping:

- HG00099.1 (`real_HG00099.1.chr12.hifi.fastq`)
- HG01167.1 (`real_HG01167.1.chr12.hifi.fastq`)
- NA18943.2 (`real_NA18943.2.chr12.hifi.fastq`)

## Notes

- Per-sample call VCFs are PASS-filtered **before** merging to avoid
  cross-sample filter contamination (see commit `2153063` in the main repo).
- The centromere density plot uses a local R script (`centromere-density.R`)
  since `chrom-density-tsv.R` expects full-chromosome coordinates.
- DeepVariant uses the PACBIO model for HiFi reads.
