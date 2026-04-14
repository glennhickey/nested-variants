# centrolign-all

Multi-chromosome centromere augref analysis. Runs the augref +
deconstruct pipeline on each chromosome independently, then aggregates
results into cross-chromosome summary figures. Optionally maps reads
and calls variants with vg call and DeepVariant.

## Prerequisites

```bash
source ../venv-nested-variants/bin/activate
```

## Input

Place `*.gbz` files (one per chromosome) in this directory. Chromosomes
are auto-discovered from GBZ filenames (excluding `*.giraffe.gbz`).

For mapping/calling, place HiFi read files in `reads/` named as
`{chrom}.{sample}.{hap}.real.fastq.gz`. Both haplotypes (.1 and .2) for
a sample are concatenated before mapping. Samples are auto-discovered.

## Running

### Local

```bash
# Subset of chromosomes:
snakemake --cores 8 --config 'chroms=chr12,chr17' mem_gb=14

# All chromosomes (if enough memory):
snakemake --cores 8
```

### Cluster (SLURM)

```bash
snakemake --profile ../profiles/slurm \
  --default-resources slurm_partition=long
```

Subset:

```bash
snakemake --profile ../profiles/slurm \
  --default-resources slurm_partition=long \
  --config 'chroms=chr12,chr17'
```

If a chromosome has loops that prevent haplotype index construction:

```bash
snakemake --profile ../profiles/slurm \
  --default-resources slurm_partition=long \
  --config 'skip_hapl=chr6'
```

Cap samples per chromosome (useful for a quick end-to-end run that hits
every panel before committing to the full cohort). Samples that appear
in the most chromosomes are preferred, so the cross-chromosome
intersection used for aggregate panels stays large:

```bash
snakemake --profile ../profiles/slurm \
  --default-resources slurm_partition=long \
  --config max_samples_per_chrom=10
```

Emit `.svg` alongside every `.png` (for Illustrator/Inkscape retouching
of individual panels — composite summary figures stay PNG since they're
stitched bitmap):

```bash
snakemake --profile ../profiles/slurm \
  --default-resources slurm_partition=long \
  --config emit_svg=true
```

Output goes to `../output/centrolign-all/`.

## Output figures

### Deconstruct (panels 1-4, no reads required)

| Figure | Description |
|--------|-------------|
| `1.augref-summary.png` | Segment lengths (cumulative + log-log) + density ideogram |
| `2.deconstruct-summary.png` | Aggregate variant types, size distribution, AF spectrum, per-sample |
| `3.deconstruct-by-chrom.png` | SNP, MNP/Indel, SV counts + AF spectrum by chromosome |
| `4.per-sample-by-chrom.png` | Deconstruct per-sample counts by chromosome (off-ref + on-ref) |

### Mapping and calling (panels 5-8b, require reads)

| Figure | Description |
|--------|-------------|
| `5.call-summary.png` | Aggregate vg call variant types, size dist, AF, per-sample |
| `6.call-by-chrom.png` | vg call SNP/MNP-Indel/SV counts + AF by chromosome |
| `6b.call-per-sample-by-chrom.png` | vg call per-sample counts by chromosome (off-ref + on-ref, 3×2 grid) |
| `6c.call-per-sample-by-chrom-combined.png` | vg call per-sample counts by chromosome, off/on-ref overlaid in 3-panel stack (same layout as 4b) |
| `7.deepvariant-summary.png` | Aggregate DeepVariant variant types, size dist, AF, per-sample |
| `8.deepvariant-by-chrom.png` | DeepVariant SNP/MNP-Indel/SV counts + AF by chromosome |
| `8b.dv-per-sample-by-chrom.png` | DeepVariant per-sample counts by chromosome (off-ref + on-ref) |

## Bandage cartoons (cartoons/)

Scripts for building small CHM13 + few-sample graphs suitable for loading
into BandageNG with per-augref colouring.

| Script | Purpose |
|--------|---------|
| `cartoons/subset-cartoon.sh` | Subset a `.vg` to CHM13 + chosen haps, optionally clip to a CHM13 range via snarls + `vg chunk`, emit augref-unified GFA |
| `cartoons/path-similarity.py` | Jaccard similarity of paths (by node set), used to pick haplotype pairs when no cenhap TSV is available |
| `cartoons/augref-colors.py` | Generate a BandageNG CSV (`Name,Colour`) colouring the CHM13 backbone + top-N longest `_alt` contigs with distinct rainbow colours |

### Example: chr12 trio cartoon (HG00099.1 + HG01891.1)

HG00099.1 and HG01891.1 were chosen via `path-similarity.py` for being
mutually similar (Jaccard ≈ 0.95) while ~44 % similar to CHM13 — enough
nested variation to be interesting but still anchored to the reference.

```bash
cd cartoons

# 1. Extract path-only GFA from chr12.vg for similarity ranking
vg convert -fW ../chr12.vg | awk '$1=="P"' > chr12.paths.gfa

# 2. Rank pairs with moderate CHM13 similarity and high mutual similarity
python3 path-similarity.py chr12.paths.gfa CHM13

# 3. Subset to CHM13 + HG00099.1 + HG01891.1, clip to CHM13:1-1300000
./subset-cartoon.sh ../chr12.vg chr12.trio5.win \
    --range CHM13#0#CHM13.0:1-1300000 \
    HG00099.1 HG01891.1

# 4. Generate BandageNG colour CSV (CHM13 backbone black, top-20 alts rainbow)
python3 augref-colors.py chr12.trio5.win.aug.gfa chr12.trio5.win.colours.csv \
    --ref-color '#000000' --rest-color '#FAD7D7'
```

Load `chr12.trio5.win.aug.gfa` in BandageNG, then *File → Load CSV data*
on `chr12.trio5.win.colours.csv` and switch the *Colour* dropdown to
**Custom colours**.
