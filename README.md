# nested-variants

Discover and analyze off-reference (nested) variants in pangenome graphs built with [vg](https://github.com/vgteam/vg) and [Minigraph-Cactus](https://github.com/ComparativeGenomicsToolkit/cactus). The pipeline computes augmented reference paths from VG files, deconstructs the graph into a nested VCF, splits variants by reference context, and produces chromosome-density ideogram plots.

## Prerequisites

| Tool | Required for | Install |
|------|-------------|---------|
| **vg** (>= 1.56) | `make paths`, `make deconstruct`, `make genotype`, `make fasta` | [vg releases](https://github.com/vgteam/vg/releases) |
| **samtools** | `make fasta`, `make surject` | `apt install samtools` / `conda install samtools` |
| **docker** | `make deepvariant` | [Docker install](https://docs.docker.com/get-docker/) |
| **bcftools** | `make split-vcf` | `apt install bcftools` / `conda install bcftools` |
| **bgzip / tabix** (htslib) | VCF compression & indexing | `apt install tabix` / `conda install htslib` |
| **make** | Pipeline orchestration | Usually pre-installed |
| Rscript | `make plots` (optional) | `apt install r-base` |
| kmc | `make genotype` (optional) | [kmc releases](https://github.com/refresh-bio/KMC) |
| shellcheck | `make test` | `apt install shellcheck` |

## Quick Start

```bash
# 1. Copy or symlink chr20 test data into data/
mkdir -p data
cp ~/dev/work/test-altpaths-chr20/chr20.vg data/

# 2. Run the pipeline (local mode, chr20 defaults)
make paths          # VG → GBZ
make deconstruct    # GBZ → VCF
make split-vcf      # VCF → onref / nestedref / offref
make plots          # offref VCF → ideogram PNG
# or simply:
make all            # runs split-vcf + plots (requires deconstruct output)
```

## Pipeline Stages

### 1. Augmented Reference Paths (`make paths`)

Runs `vg paths` on each input VG file to compute augmented reference paths, then merges and converts to GBZ format.

- **Input:** VG file(s) (`data/chr20.vg`)
- **Output:** `output/chr20.nested.gbz`
- **Script:** `scripts/paths.sh`

### 2. Deconstruct (`make deconstruct`)

Runs `vg deconstruct` on the GBZ to produce a nested VCF with on-reference and off-reference variant calls.

- **Input:** GBZ from stage 1
- **Output:** `output/chr20.nested.vcf.gz` (+ .tbi index)
- **Script:** `scripts/deconstruct.sh`

### 3. Split VCF (`make split-vcf`)

Splits the nested VCF into three categories based on reference context:

| Category | File suffix | Description |
|----------|-----------|-------------|
| On-reference | `.onref.vcf.gz` | Variants on contigs starting with the reference prefix |
| Nested-reference | `.nestedref.vcf.gz` | On-reference variants at nesting level > 0 |
| Off-reference | `.offref.vcf.gz` | Variants on non-reference contigs |

- **Script:** `scripts/split-ref.sh`

### 4. Genotype (optional, `make genotype`)

Aligns reads to the graph with `vg giraffe` and calls variants with `vg call`. Requires `MAP_GBZ`, `READS`, `HAPL`, and `SAMPLE` to be set.

- `MAP_GBZ` — pre-built GBZ for read mapping. The `.hapl` index must match this GBZ (i.e., the original pangenome GBZ distributed with the HPRC release, **not** the augmented-reference GBZ built by `make paths`). Giraffe uses this for alignment; `vg call` then uses the augmented GBZ for variant calling.
- `READS` — a text file listing input FASTQ paths (one per line, typically two lines for paired-end reads). This is user-provided sequencing data.
- `HAPL` — haplotype index file (`.hapl`) for the graph, typically distributed alongside the HPRC pangenome release.
- `SAMPLE` — sample name to embed in the output GAM/VCF.

```bash
make genotype \
  MAP_GBZ=data/hprc-v2.0-mc-chm13.gbz \
  READS=data/HG002.reads.idx \
  HAPL=data/hprc-v2.0-mc-chm13.hapl \
  SAMPLE=HG002
```

### 5. Haplotype Index (optional, `make haplotypes`)

Builds a `.hapl` index from the augmented GBZ for haplotype-aware read mapping with giraffe. Runs `vg index` (distance index), `vg gbwt` (r-index), and `vg haplotypes` in sequence; intermediate files are cleaned up automatically.

- **Input:** Augmented GBZ from `make paths`
- **Output:** `output/<OUT_NAME>.hapl`
- **Script:** `scripts/haplotypes.sh`

If you already have a `.hapl` index (e.g. from an HPRC release), set `HAPL=<path>` directly and skip this step.

### 6. Surject (optional, `make surject`)

Projects GAM alignments onto the augmented reference paths to produce a coordinate-sorted BAM file with index. Uses the augmented GBZ (from `make paths`) so that reads are placed on the nested reference contigs.

- **Input:** GAM from giraffe, augmented GBZ from `make paths`
- **Output:** `output/<SAMPLE>.bam` (+ `.bam.bai` index)
- **Script:** `scripts/surject.sh`

```bash
make surject \
  OUT_DIR=output/v2-chm13 \
  OUT_NAME=hprc-v2.0-mc-chm13.nested.95 \
  SAMPLE=HG002
```

### 7. FASTA Extraction (optional, `make fasta`)

Extracts augmented reference paths from the GBZ as a bgzipped FASTA file and creates `.fai` and `.gzi` indexes. This is needed as the reference for DeepVariant.

- **Input:** Augmented GBZ from `make paths`
- **Output:** `output/<OUT_NAME>.fa.gz` (+ `.fa.gz.fai` and `.fa.gz.gzi` indexes)
- **Script:** `scripts/fasta.sh`

```bash
make fasta \
  OUT_DIR=output/v2-chm13 \
  OUT_NAME=hprc-v2.0-mc-chm13.nested.95
```

### 8. DeepVariant (optional, `make deepvariant`)

Runs [DeepVariant](https://github.com/google/deepvariant) via Docker to call variants from the surjected BAM against the augmented reference FASTA. Requires `SAMPLE` and Docker.

- **Input:** BAM from `make surject`, FASTA from `make fasta`
- **Output:** `output/<SAMPLE>.deepvariant.vcf.gz`
- **Script:** `scripts/deepvariant.sh`

```bash
make deepvariant \
  OUT_DIR=output/v2-chm13 \
  OUT_NAME=hprc-v2.0-mc-chm13.nested.95 \
  SAMPLE=HG002
```

## Configuration

All settings live in `config.mk` (committed defaults) and can be overridden in `config.local.mk` (gitignored) or on the command line.

| Variable | Default | Description |
|----------|---------|-------------|
| `EXEC_MODE` | `local` | `local` or `slurm` |
| `REF` | `GRCh38` | Reference name for vg (`GRCh38`, `CHM13`, etc.) |
| `VG` | `data/chr20.vg` | Input VG file(s) |
| `OUT_DIR` | `output` | Output directory |
| `OUT_NAME` | `chr20.nested` | Output filename prefix |
| `MIN_AUGREF_LEN` | `50` | Minimum augref fragment length |
| `MAP_GBZ` | *(empty)* | Pre-built GBZ for read mapping (must match `.hapl`); falls back to pipeline GBZ if unset |
| `DV_VERSION` | `1.9.0` | DeepVariant Docker image version |
| `REFGAPS_BED` | *(empty)* | BED file for reference gap overlay on plots |

See `config.mk` for the full list.

## Testing

```bash
make test          # shellcheck + --help flag tests
```

The test suite runs [shellcheck](https://www.shellcheck.net/) on all shell scripts and verifies that each script's `--help` flag exits cleanly. Full pipeline tests require `vg` and test data.

For a small end-to-end test using *S. cerevisiae* chromosome I, see [yeast-test/README.md](yeast-test/README.md).

## Cluster Usage

When `EXEC_MODE=slurm`, each pipeline step submits SLURM jobs via `sbatch -W` and waits for completion before proceeding to the next step. SLURM resources (`CPUS`, `MEM`, `TIME`, `PARTITION`) are all configurable.

### Running multiple graphs

Override `VG`, `REF`, `OUT_DIR`, and `OUT_NAME` on the command line to run different inputs into separate output directories. The `VG` variable accepts glob patterns (including bash extended globs like `!(*.d9).vg`).

```bash
# HPRC v2.0 CHM13 — full pipeline including genotyping + DeepVariant
make paths deconstruct genotype surject fasta deepvariant split-vcf plots call-plots \
  EXEC_MODE=slurm \
  REF=CHM13 \
  VG='/path/to/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.chroms/!(*.d9).vg' \
  OUT_DIR=output/v2-chm13 \
  OUT_NAME=hprc-v2.0-mc-chm13.nested.95 \
  MAP_GBZ=/path/to/hprc-v2.0-mc-chm13.gbz \
  READS=data/HG002.reads.idx \
  HAPL=data/hprc-v2.0-mc-chm13.hapl \
  SAMPLE=HG002 \
  REFGAPS_BED=data/hprc-v2.0-mc-chm13.refgaps.bed

# HPRC v2.0 GRCh38 — full pipeline including genotyping + DeepVariant
make paths deconstruct genotype surject fasta deepvariant split-vcf plots call-plots \
  EXEC_MODE=slurm \
  REF=GRCh38 \
  VG='/path/to/hprc-v2.0-mc-grch38/hprc-v2.0-mc-grch38.chroms/!(*.d9).vg' \
  OUT_DIR=output/v2-grch38 \
  OUT_NAME=hprc-v2.0-mc-grch38.nested.95 \
  MAP_GBZ=/path/to/hprc-v2.0-mc-grch38.gbz \
  READS=data/HG002.reads.idx \
  HAPL=data/hprc-v2.0-mc-grch38.hapl \
  SAMPLE=HG002 \
  REFGAPS_BED=data/hprc-v2.0-mc-grch38.refgaps.bed

# HPRC v1.1 CHM13 — without genotyping
make paths deconstruct split-vcf plots \
  EXEC_MODE=slurm \
  REF=CHM13 \
  VG='/path/to/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.chroms/!(*.d9).vg' \
  OUT_DIR=output/v1-chm13 \
  OUT_NAME=hprc-v1.1-mc-chm13.nested.95 \
  REFGAPS_BED=data/hprc-v1.1-mc-chm13.refgaps.bed
```

Each run gets its own output directory with the full set of outputs:
```
output/v2-chm13/
├── hprc-v2.0-mc-chm13.nested.95.gbz              # augmented reference graph
├── hprc-v2.0-mc-chm13.nested.95.augref-segs.tsv   # augref segment table
├── hprc-v2.0-mc-chm13.nested.95.vcf.gz            # nested VCF (deconstruct)
├── hprc-v2.0-mc-chm13.nested.95.onref.vcf.gz      # on-reference variants
├── hprc-v2.0-mc-chm13.nested.95.nestedref.vcf.gz  # nested-reference variants
├── hprc-v2.0-mc-chm13.nested.95.offref.vcf.gz     # off-reference variants
├── hprc-v2.0-mc-chm13.nested.95.offref.png        # deconstruct density ideogram
├── HG002.gam                                       # read alignments (genotype)
├── HG002.bam                                       # surjected alignments (sorted BAM)
├── HG002.bam.bai                                   # BAM index
├── hprc-v2.0-mc-chm13.nested.95.fa.gz               # augmented reference FASTA (bgzipped)
├── hprc-v2.0-mc-chm13.nested.95.fa.gz.fai           # FASTA index
├── hprc-v2.0-mc-chm13.nested.95.fa.gz.gzi           # bgzip index
├── HG002.deepvariant.vcf.gz                        # DeepVariant VCF
├── HG002.pack                                      # coverage pileup
├── HG002.vcf.gz                                    # genotyped VCF (vg call)
└── HG002.call-density.png                          # genotyped density ideogram
```

### Genotyping a sample

To genotype a sample against the augmented reference and plot the results:

```bash
make genotype call-plots \
  EXEC_MODE=slurm \
  REF=CHM13 \
  OUT_DIR=output/v2-chm13 \
  OUT_NAME=hprc-v2.0-mc-chm13.nested.95 \
  MAP_GBZ=/path/to/hprc-v2.0-mc-chm13.gbz \
  READS=data/sample.reads.idx \
  HAPL=data/hprc-v2.0-mc-chm13.hapl \
  SAMPLE=NA12878
```

### SLURM resource tuning

The default SLURM resources can be overridden per run:

```bash
make paths EXEC_MODE=slurm CPUS=32 MEM=400gb TIME=24:00:00 PARTITION=long ...
```

### Using config.local.mk

For repeated use, create a `config.local.mk` with your common settings:

```bash
cp config.local.mk.example config.local.mk
# Edit config.local.mk with your defaults
```

Then override only what changes per run:

```bash
make paths deconstruct split-vcf plots \
  REF=CHM13 \
  VG='/path/to/chm13-chroms/!(*.d9).vg' \
  OUT_DIR=output/v2-chm13 \
  OUT_NAME=hprc-v2.0-mc-chm13.nested.95
```

## Repository Structure

```
nested-variants/
├── Makefile                    Pipeline orchestrator
├── config.mk                   Default config (chr20 local test)
├── config.local.mk.example     Template for cluster config
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
│   ├── chrom-density-tsv.R     Ideogram from TSV nesting files
│   ├── chrom-density-vcf.R     Ideogram from VCF INFO fields
│   └── chrom-density-call.R    Ideogram from genotyped VCF (vg call)
├── annotation/
│   ├── download-hprc-annotations.py
│   └── intersect-annotations.py
├── test/
│   └── test-pipeline.sh        CI test script
├── data/                        Gitignored; put input data here
├── .github/workflows/ci.yml
└── LICENSE
```

## License

MIT License. See [LICENSE](LICENSE).
