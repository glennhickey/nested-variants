# centrolign-all

Multi-chromosome centromere augref analysis. Runs the augref +
deconstruct pipeline on each chromosome independently, then aggregates
results into cross-chromosome summary figures.

## Prerequisites

```bash
source ../venv-nested-variants/bin/activate
```

## Input

Place `*.gbz` files (one per chromosome) in this directory. Chromosomes
are auto-discovered from GBZ filenames (excluding `*.giraffe.gbz`).

## Running

### Local (subset)

```bash
snakemake --cores 8 --config 'chroms=chr12,chr17'
```

### Cluster (all chromosomes)

```bash
snakemake --profile ../profiles/slurm
```

### Cluster (subset)

```bash
snakemake --profile ../profiles/slurm --config 'chroms=chr12,chr17'
```

Output goes to `output-all/`.

## Output figures

| Figure | Description |
|--------|-------------|
| `1.augref-summary.png` | Segment lengths (cumulative + log-log) + density ideogram |
| `2.deconstruct-summary.png` | Aggregate variant types, size distribution, AF spectrum, per-sample |
| `3.deconstruct-by-chrom.png` | SNP, MNP/Indel, SV counts + AF spectrum by chromosome |
| `4.per-sample-by-chrom.png` | Per-sample SNP, MNP/Indel, SV box plots by chromosome |
