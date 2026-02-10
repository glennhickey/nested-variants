# nested-variants

Discover and analyze off-reference (nested) variants in pangenome graphs built with [vg](https://github.com/vgteam/vg) and [Minigraph-Cactus](https://github.com/ComparativeGenomicsToolkit/cactus). The pipeline computes augmented reference paths from VG files, deconstructs the graph into a nested VCF, splits variants by reference context, and produces chromosome-density ideogram plots.

## Prerequisites

| Tool | Required for | Install |
|------|-------------|---------|
| **vg** (>= 1.56) | `paths`, `deconstruct`, `genotype`, `fasta` | [vg releases](https://github.com/vgteam/vg/releases) |
| **samtools** | `fasta`, `surject` | `apt install samtools` / `conda install samtools` |
| **docker** | `deepvariant` | [Docker install](https://docs.docker.com/get-docker/) |
| **bcftools** | `split_vcf`, `merge_call_vcfs`, `merge_dv_vcfs` | `apt install bcftools` / `conda install bcftools` |
| **bgzip / tabix** (htslib) | VCF compression & indexing | `apt install tabix` / `conda install htslib` |
| **snakemake** (>= 8) | Pipeline orchestration | `pip install snakemake` / `conda install snakemake` |
| Rscript | `plots` (optional) | `apt install r-base` |
| kmc | `genotype` (optional) | [kmc releases](https://github.com/refresh-bio/KMC) |
| shellcheck | testing | `apt install shellcheck` |

## Quick Start

```bash
# 1. Copy or symlink chr20 test data into data/
mkdir -p data
cp ~/dev/work/test-altpaths-chr20/chr20.vg data/

# 2. Run the pipeline (local mode, chr20 defaults)
snakemake --cores 8 graph_only
# This runs: paths → deconstruct → split_vcf + plots
```

## Pipeline Stages

### 1. Augmented Reference Paths (`paths`)

Runs `vg paths` on each input VG file to compute augmented reference paths, then merges and converts to GBZ format.

- **Input:** VG file(s) (`data/chr20.vg`)
- **Output:** `output/chr20.nested.gbz`
- **Script:** `scripts/paths.sh`

### 2. Deconstruct (`deconstruct`)

Runs `vg deconstruct` on the GBZ to produce a nested VCF with on-reference and off-reference variant calls.

- **Input:** GBZ from stage 1
- **Output:** `output/chr20.nested.vcf.gz` (+ .tbi index)
- **Script:** `scripts/deconstruct.sh`

### 3. Split VCF (`split_vcf`)

Splits the nested VCF into three categories based on reference context:

| Category | File suffix | Description |
|----------|-----------|-------------|
| On-reference | `.onref.vcf.gz` | Variants on contigs starting with the reference prefix |
| Nested-reference | `.nestedref.vcf.gz` | On-reference variants at nesting level > 0 |
| Off-reference | `.offref.vcf.gz` | Variants on non-reference contigs |

- **Script:** `scripts/split-ref.sh`

### 4. Genotype (optional, `genotype_all`)

Aligns reads to the graph with `vg giraffe` and calls variants with `vg call`. Requires `samples` to be configured (see [Configuration](#configuration)). The pipeline-built GBZ and `.hapl` index are used automatically for read mapping.

- Each sample entry maps a sample name to its reads index file (a text file listing FASTQ paths, one per line). Paths can be local files or remote URLs (`gs://`, `http://`, `https://`); remote files are automatically downloaded to node-local scratch before mapping.

```bash
snakemake --cores 8 genotype_all \
  --config 'samples={HG002: data/HG002.reads.idx}'
```

Example reads index file:
```
/local/data/HG002.R1.fastq.gz
/local/data/HG002.R2.fastq.gz
```

Or with remote URLs:
```
gs://deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG002.novaseq.pcr-free.30x.R1.fastq.gz
gs://deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG002.novaseq.pcr-free.30x.R2.fastq.gz
```

### 5. Haplotype Index (`haplotypes`)

Builds a `.hapl` index from the augmented GBZ for haplotype-aware read mapping with giraffe. Runs `vg index` (distance index), `vg gbwt` (r-index), and `vg haplotypes` in sequence; intermediate files are cleaned up automatically. This is built automatically when genotyping is requested.

- **Input:** Augmented GBZ from `paths`
- **Output:** `output/<out_name>.hapl`
- **Script:** `scripts/haplotypes.sh`

### 6. Surject (optional, `surject`)

Projects GAM alignments onto the augmented reference paths to produce a coordinate-sorted BAM file with index. Uses the augmented GBZ (from `paths`) so that reads are placed on the nested reference contigs.

- **Input:** GAM from giraffe, augmented GBZ from `paths`
- **Output:** `output/<sample>.bam` (+ `.bam.bai` index)
- **Script:** `scripts/surject.sh`

### 7. FASTA Extraction (optional, `fasta`)

Extracts augmented reference paths from the GBZ as a bgzipped FASTA file and creates `.fai` and `.gzi` indexes. This is needed as the reference for DeepVariant.

- **Input:** Augmented GBZ from `paths`
- **Output:** `output/<out_name>.fa.gz` (+ `.fa.gz.fai` and `.fa.gz.gzi` indexes)
- **Script:** `scripts/fasta.sh`

### 8. DeepVariant (optional, `deepvariant_all`)

Runs [DeepVariant](https://github.com/google/deepvariant) via Docker to call variants from the surjected BAM against the augmented reference FASTA. Requires `samples` and Docker.

- **Input:** BAM from `surject`, FASTA from `fasta`
- **Output:** `output/<sample>.deepvariant.vcf.gz`
- **Script:** `scripts/deepvariant.sh`

### 9. Full Pipeline (`all`)

Runs genotyping and DeepVariant for all configured samples, then merges the per-sample VCFs with `bcftools merge` and produces density plots on the merged VCFs.

```bash
snakemake --cores 8 all \
  --config 'samples={HG002: data/HG002.reads.idx, NA12878: data/NA12878.reads.idx}'
```

## Configuration

All settings live in `config.yaml` (committed defaults) and can be overridden in `config.local.yaml` (gitignored) or on the command line with `--config`.

| Variable | Default | Description |
|----------|---------|-------------|
| `ref` | `GRCh38` | Reference name for vg (`GRCh38`, `CHM13`, etc.) |
| `vg` | `data/chr20.vg` | Input VG file(s) |
| `out_dir` | `output` | Output directory |
| `out_name` | `chr20.nested` | Output filename prefix |
| `min_augref_len` | `50` | Minimum augref fragment length |
| `dv_version` | `1.9.0` | DeepVariant Docker image version |
| `refgaps_bed` | *(empty)* | BED file for reference gap overlay on plots |
| `scale_type` | `log1p` | Scale type for density plots |
| `samples` | `{}` | Map of sample name → reads index file path |
| `samples_tsv` | *(unset)* | TSV file with `sample` and `reads_index` columns |
| `cpus` | `8` | Global CPU fallback for all rules |
| `mem_gb` | `200` | Global memory (GB) fallback for all rules |
| `runtime_min` | `960` | Global runtime (minutes) fallback for all rules |

Each rule has built-in defaults that are used when neither `{rule}_cpus`/`{rule}_mem_gb` nor the global `cpus`/`mem_gb` is set:

| Rule | CPUs | Memory (GB) | Notes |
|------|------|-------------|-------|
| `paths` | 128 | 512 | Graph construction |
| `deconstruct` | 128 | 512 | VCF extraction from graph |
| `haplotypes` | 128 | 512 | Haplotype index construction |
| `giraffe` | 128 | 512 | Read mapping |
| `surject` | 128 | 512 | GAM → BAM projection |
| `call` | 128 | 512 | Variant calling |
| `deepvariant` | 128 | 512 | Deep learning variant calling |
| `fasta` | 128 | 512 | Reference extraction |

Override per-rule: `--config giraffe_cpus=64 giraffe_mem_gb=256`

See `config.yaml` for the full list.

To use a local config file, create `config.local.yaml` and pass it with `--configfile`:

```bash
snakemake --cores 8 --configfile config.local.yaml graph_only
```

## Testing

```bash
bash test/test-pipeline.sh   # shellcheck + --help flag tests
```

The test suite runs [shellcheck](https://www.shellcheck.net/) on all shell scripts and verifies that each script's `--help` flag exits cleanly. Full pipeline tests require `vg` and test data.

For a small end-to-end test using *S. cerevisiae* chromosome I, see [yeast-test/README.md](yeast-test/README.md).

## Cluster Usage

Snakemake handles SLURM scheduling natively. All shell scripts are always called with `--local`; Snakemake submits each rule as a separate SLURM job.

**Local run:**
```bash
snakemake --cores 8 graph_only
```

**SLURM run (direct):**
```bash
snakemake --executor slurm \
  --default-resources mem_mb=200000 runtime=960 \
  --jobs 50 graph_only
```

**SLURM run (profile):**
```bash
snakemake --profile profiles/slurm graph_only
```

For SLURM execution, install the executor plugin:
```bash
pip install snakemake-executor-plugin-slurm
```

### Node-local scratch

The `giraffe` and `surject` steps do heavy random I/O on the GBZ file. On clusters with slow shared filesystems, these scripts automatically stage the GBZ to node-local scratch (`$TMPDIR`) before processing, and `samtools sort` temp files are also written there. If your cluster's SLURM prolog sets `$TMPDIR` to node-local storage, this works automatically. Otherwise, set `tmpdir` in the SLURM profile:

```yaml
# profiles/slurm/config.yaml
default-resources:
  tmpdir: /scratch/$USER
```

If `$TMPDIR` is not set, the scripts fall back to the output directory.

### Running multiple graphs

Override config values on the command line to run different inputs into separate output directories:

```bash
# HPRC v2.1 CHM13 — full pipeline including genotyping + DeepVariant
GIAB=/private/home/ghickey/dev/work/giab-reads
snakemake --profile profiles/slurm all \
  --config \
    ref=CHM13 \
    vg='/private/groups/hprc/hprc-graphs/hprc-v2.1-dec23/hprc-v2.1-mc-chm13-eval/hprc-v2.1-mc-chm13-eval.chroms/!(*.d*).vg' \
    out_dir=output/v2.1-chm13 \
    out_name=hprc-v2.1-mc-chm13.nested \
    refgaps_bed=data/hprc-v2.1-mc-chm13.refgaps.bed \
    "samples={HG001: $GIAB/HG001.novaseq.pcr-free.gs.paths, HG002: $GIAB/HG002.novaseq.pcr-free.gs.paths, HG003: $GIAB/HG003.novaseq.pcr-free.gs.paths, HG004: $GIAB/HG004.novaseq.pcr-free.gs.paths, HG005: $GIAB/HG005.novaseq.pcr-free.gs.paths, HG006: $GIAB/HG006.novaseq.pcr-free.gs.paths, HG007: $GIAB/HG007.novaseq.pcr-free.gs.paths}"

# HPRC v2.1 GRCh38 — full pipeline including genotyping + DeepVariant
snakemake --profile profiles/slurm all \
  --config \
    ref=GRCh38 \
    vg='/private/groups/hprc/hprc-graphs/hprc-v2.1-dec23/hprc-v2.1-mc-grch38-eval/hprc-v2.1-mc-grch38-eval.chroms/!(*.d*).vg' \
    out_dir=output/v2.1-grch38 \
    out_name=hprc-v2.1-mc-grch38.nested \
    refgaps_bed=data/hprc-v2.1-mc-grch38.refgaps.bed \
    'samples={HG002: data/HG002.reads.idx}'

# HPRC v1.1 CHM13 — without genotyping
snakemake --profile profiles/slurm graph_only \
  --config \
    ref=CHM13 \
    vg='/private/groups/hprc/hprc-graphs/hprc-v1.1-jul4/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.chroms/!(*.d9).vg' \
    out_dir=output/v1.1-chm13 \
    out_name=hprc-v1.1-mc-chm13.nested \
    refgaps_bed=data/hprc-v1.1-mc-chm13.refgaps.bed
```

Each run gets its own output directory with the full set of outputs:
```
output/v2.1-chm13/
├── hprc-v2.1-mc-chm13.nested.gbz               # augmented reference graph
├── hprc-v2.1-mc-chm13.nested.gfa.gz            # augmented GFA
├── hprc-v2.1-mc-chm13.nested.augref-segs.tsv   # augref segment table
├── hprc-v2.1-mc-chm13.nested.vcf.gz            # nested VCF (deconstruct)
├── hprc-v2.1-mc-chm13.nested.onref.vcf.gz      # on-reference variants
├── hprc-v2.1-mc-chm13.nested.nestedref.vcf.gz  # nested-reference variants
├── hprc-v2.1-mc-chm13.nested.offref.vcf.gz     # off-reference variants
├── hprc-v2.1-mc-chm13.nested.augref-length-hist.png  # segment length histogram
├── hprc-v2.1-mc-chm13.nested.offref.png        # off-reference density ideogram
├── hprc-v2.1-mc-chm13.nested.vcf-stats.tsv     # variant statistics (deconstruct)
├── hprc-v2.1-mc-chm13.nested.variant-types.png  # variant type bar chart
├── hprc-v2.1-mc-chm13.nested.size-dist.png     # indel/SV size distribution
├── hprc-v2.1-mc-chm13.nested.af-spectrum.png   # allele frequency spectrum
├── hprc-v2.1-mc-chm13.nested.hapl              # haplotype index (for giraffe)
├── hprc-v2.1-mc-chm13.nested.fa.gz             # augmented reference FASTA (bgzipped)
├── hprc-v2.1-mc-chm13.nested.fa.gz.fai         # FASTA index
├── hprc-v2.1-mc-chm13.nested.fa.gz.gzi         # bgzip index
├── HG002.gam                                    # read alignments (genotype)
├── HG002.bam                                    # surjected alignments (sorted BAM)
├── HG002.bam.bai                                # BAM index
├── HG002.pack                                   # coverage pileup
├── HG002.vcf.gz                                 # genotyped VCF (vg call)
├── HG002.deepvariant.vcf.gz                     # DeepVariant VCF
├── HG002.call-offref.png                        # call off-reference density
├── HG002.dv-offref.png                          # DeepVariant off-reference density
├── merged.call.vcf.gz                           # merged call VCFs (all)
├── merged.deepvariant.vcf.gz                    # merged DeepVariant VCFs (all)
├── merged.call-offref.png                       # merged call density
├── merged.dv-offref.png                         # merged DV density
├── merged.call.af-spectrum.png                  # merged call AF spectrum
└── merged.dv.af-spectrum.png                    # merged DV AF spectrum
```

### Genotyping a sample

To genotype a sample against the augmented reference and plot the results:

```bash
snakemake --cores 8 genotype_all deepvariant_all \
  --config \
    ref=CHM13 \
    out_dir=output/v2.1-chm13 \
    out_name=hprc-v2.1-mc-chm13.nested \
    'samples={NA12878: data/sample.reads.idx}'
```

### Multi-sample batch processing

Configure multiple samples to process them in parallel and merge their VCFs:

```bash
snakemake --profile profiles/slurm all \
  --config \
    ref=CHM13 \
    out_dir=output/v2.1-chm13 \
    out_name=hprc-v2.1-mc-chm13.nested \
    'samples={HG002: data/HG002.reads.idx, NA12878: data/NA12878.reads.idx}'
```

Or use a samples TSV file:
```bash
snakemake --profile profiles/slurm all \
  --config samples_tsv=samples.tsv ...
```

Where `samples.tsv` contains:
```
sample	reads_index
HG002	data/HG002.reads.idx
NA12878	data/NA12878.reads.idx
```

## Repository Structure

```
nested-variants/
├── Snakefile                   Pipeline orchestrator
├── config.yaml                 Default config (chr20 local test)
├── profiles/slurm/config.yaml  SLURM profile template
├── scripts/
│   ├── paths.sh                VG → GBZ (augmented reference paths)
│   ├── deconstruct.sh          GBZ → VCF (vg deconstruct)
│   ├── giraffe.sh              GBZ + reads → GAM (vg giraffe)
│   ├── surject.sh              GAM → sorted BAM (vg surject)
│   ├── haplotypes.sh           GBZ → .hapl index (vg haplotypes)
│   ├── fasta.sh                GBZ → augmented reference FASTA
│   ├── deepvariant.sh          BAM + FASTA → VCF (DeepVariant Docker)
│   ├── call.sh                 GBZ + GAM → VCF (vg call)
│   ├── split-ref.sh            VCF → onref/nestedref/offref VCFs
│   ├── offref-length-hist.R    Size distribution histograms
│   ├── chrom-density-common.R  Shared ideogram plotting code
│   └── chrom-density-segs.R    Off-reference density ideogram (VCF + segments table)
├── annotation/
│   ├── download-hprc-annotations.py
│   └── intersect-annotations.py
├── test/
│   └── test-pipeline.sh        CI test script
├── yeast-test/                 End-to-end test (S. cerevisiae chrI)
├── data/                       Gitignored; put input data here
├── .github/workflows/ci.yml
└── LICENSE
```

## License

MIT License. See [LICENSE](LICENSE).
