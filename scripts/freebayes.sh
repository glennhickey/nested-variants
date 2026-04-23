#!/bin/bash

################################################################################
# freebayes.sh
#
# Description:
#   Runs FreeBayes in parallel to call variants from a BAM file against a
#   FASTA reference.  Uses GNU parallel to split by genomic regions; each
#   region is a separate freebayes process. Runs natively (freebayes binary
#   on PATH) by default. Pass --docker IMAGE to fall back to containerised
#   execution when the binary isn't available.
#
# Usage:
#   freebayes.sh --bam <file.bam> --ref <file.fa.gz> --sample <name> \
#                --out-dir <dir> --out-name <name> \
#                [--region-size N] [--docker IMAGE] \
#                [--cpus N] [--mem size] [--time HH:MM:SS] [--partition name] [--local]
#
################################################################################

set -eo pipefail

# Initialize variables
BAM=""
REF=""
SAMPLE=""
OUTPUT_DIR="."
OUTPUT_NAME=""
REGION_SIZE=1000000   # 1 Mb per shard: ~9.7k regions on HPRC (was 100 kb /
                      # ~97k regions); with parallel -j 96 that's ~100
                      # regions per worker — enough load balancing, far
                      # less per-region startup / file-system overhead.
EXTRA_ARGS=""
DOCKER_IMAGE=""  # empty = run freebayes natively; set to a tag for docker fallback

# SLURM resource defaults
CPUS="16"
MEM="128gb"
TIME="16:00:00"
PARTITION="long"
JOB_NAME="freebayes"
LOCAL=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --bam)
            BAM="$2"
            shift 2
            ;;
        --ref)
            REF="$2"
            shift 2
            ;;
        --sample)
            SAMPLE="$2"
            shift 2
            ;;
        --out-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --out-name)
            OUTPUT_NAME="$2"
            shift 2
            ;;
        --region-size)
            REGION_SIZE="$2"
            shift 2
            ;;
        --extra-args)
            EXTRA_ARGS="$2"
            shift 2
            ;;
        --docker)
            DOCKER_IMAGE="$2"
            shift 2
            ;;
        --cpus)
            CPUS="$2"
            shift 2
            ;;
        --mem)
            MEM="$2"
            shift 2
            ;;
        --time)
            TIME="$2"
            shift 2
            ;;
        --partition)
            PARTITION="$2"
            shift 2
            ;;
        --local)
            LOCAL=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 --bam <file.bam> --ref <file.fa> --sample <name> --out-dir <dir> --out-name <name> [options]"
            echo ""
            echo "Required Options:"
            echo "  --bam <file>          Sorted BAM file"
            echo "  --ref <file>          FASTA reference file (bgzipped OK)"
            echo "  --sample <name>       Sample name"
            echo "  --out-dir <dir>       Output directory for VCF file"
            echo "  --out-name <name>     Output name for VCF file"
            echo ""
            echo "FreeBayes Options:"
            echo "  --region-size <N>     Region chunk size for parallelization (default: 1000000)"
            echo "  --docker <image>      Docker image (default: staphb/freebayes:1.3.7)"
            echo ""
            echo "Execution Options:"
            echo "  --local               Run commands locally instead of via SLURM"
            echo ""
            echo "SLURM Resource Options (optional, with defaults):"
            echo "  --cpus <N>            CPUs per task (default: 16)"
            echo "  --mem <size>          Memory per job (default: 128gb)"
            echo "  --time <time>         Wall clock limit (default: 16:00:00)"
            echo "  --partition <p>       SLURM partition/queue (default: long)"
            exit 0
            ;;
        *)
            echo "Error: Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Validate required parameters
if [ -z "$BAM" ]; then
    echo "Error: --bam is required"
    exit 1
fi

if [ -z "$REF" ]; then
    echo "Error: --ref is required"
    exit 1
fi

if [ -z "$SAMPLE" ]; then
    echo "Error: --sample is required"
    exit 1
fi

if [ -z "$OUTPUT_NAME" ]; then
    echo "Error: --out-name is required"
    exit 1
fi

# Ensure OUTPUT_NAME has .vcf.gz extension
if [[ ! "$OUTPUT_NAME" =~ \.vcf\.gz$ ]]; then
    OUTPUT_NAME="${OUTPUT_NAME}.vcf.gz"
fi

set -x

mkdir -p "$OUTPUT_DIR"

# All intermediate work — input copies, decompressed ref, per-region shards,
# merged VCF — lives on node-local scratch. Only the final VCF and .tbi are
# copied back to $OUTPUT_DIR at the end.
#
# Scratch priority: $TMPDIR (set by SLURM to node-local dir in most setups)
# → /data/tmp (cluster-specific) → $OUTPUT_DIR (fallback, shared FS).
SCRATCH_BASE="${TMPDIR:-}"
if [ -z "$SCRATCH_BASE" ] && [ -d /data/tmp ]; then
    SCRATCH_BASE=/data/tmp
fi
if [ -z "$SCRATCH_BASE" ]; then
    SCRATCH_BASE="$OUTPUT_DIR"
fi
WORK_TMPDIR=$(mktemp -d "${SCRATCH_BASE}/freebayes.${SAMPLE}.XXXXXX")
trap 'rm -rf "$WORK_TMPDIR"' EXIT INT TERM
echo "Scratch dir: $WORK_TMPDIR"

# Final outputs land here.
VCF="${OUTPUT_DIR}/${OUTPUT_NAME}"
# Intermediate outputs are on scratch.
VCF_SCRATCH="${WORK_TMPDIR}/out.vcf.gz"

# Stage inputs to scratch. BAM + BAI are the big ones (tens of GB on HPRC);
# the reference FASTA is ~3 GB bgzipped / ~9 GB decompressed.
echo "Staging BAM + BAI to scratch"
BAM_ABS="$(realpath "$BAM")"
cp "$BAM_ABS" "${WORK_TMPDIR}/input.bam"
if [ -f "${BAM_ABS}.bai" ]; then
    cp "${BAM_ABS}.bai" "${WORK_TMPDIR}/input.bam.bai"
else
    samtools index "${WORK_TMPDIR}/input.bam"
fi
BAM_LOCAL="${WORK_TMPDIR}/input.bam"

echo "Staging reference to scratch"
REF_ABS_SRC="$(realpath "$REF")"
if [[ "$REF_ABS_SRC" == *.gz ]]; then
    gunzip -c "$REF_ABS_SRC" > "${WORK_TMPDIR}/ref.fa"
else
    cp "$REF_ABS_SRC" "${WORK_TMPDIR}/ref.fa"
fi
samtools faidx "${WORK_TMPDIR}/ref.fa"
REF_LOCAL="${WORK_TMPDIR}/ref.fa"
FAI_LOCAL="${REF_LOCAL}.fai"

# Generate regions from FAI, filtered to contigs present in the BAM.
# Without this, HPRC-scale FASTAs with 283k contigs generate hundreds of
# thousands of freebayes invocations on empty contigs.
BAM_CONTIGS="${WORK_TMPDIR}/bam-contigs.txt"
samtools idxstats "$BAM_LOCAL" | awk '$3 > 0 {print $1}' > "$BAM_CONTIGS"
echo "BAM has reads on $(wc -l < "$BAM_CONTIGS") of $(wc -l < "$FAI_LOCAL") contigs"

REGIONS_FILE="${WORK_TMPDIR}/regions.txt"
awk -v size="$REGION_SIZE" 'NR==FNR{keep[$1]=1;next} ($1 in keep) {
    chrom = $1; len = $2; pos = 0
    while (pos < len) {
        end = pos + size
        if (end > len) end = len
        print chrom ":" pos "-" end
        pos = end
    }
}' "$BAM_CONTIGS" "$FAI_LOCAL" > "$REGIONS_FILE"
rm -f "$BAM_CONTIGS"

echo "Generated $(wc -l < "$REGIONS_FILE") regions (${REGION_SIZE}bp chunks)"

if [ -n "$DOCKER_IMAGE" ]; then
    # Docker fallback (legacy path). Mount the scratch dir so the container
    # sees the staged inputs and writes shards back to the same place.
    echo "Running FreeBayes via Docker: $DOCKER_IMAGE"
    FB_CMD="docker run --rm --user $(id -u):$(id -g) -v ${WORK_TMPDIR}:${WORK_TMPDIR} $DOCKER_IMAGE freebayes"
else
    # Native path: freebayes binary on PATH (preferred — avoids ~30k docker
    # starts per sample on HPRC-scale augref).
    if ! command -v freebayes >/dev/null 2>&1; then
        echo "Error: freebayes not on PATH (and --docker not set)" >&2
        exit 1
    fi
    echo "Running FreeBayes natively: $(command -v freebayes)"
    FB_CMD="freebayes"
fi

# Per-region temp-file parallelism: each worker writes its region's VCF to
# its own file on scratch, named with a zero-padded FAI-order index so
# glob-order = final-output order. Avoids `parallel -k`, which forced
# workers to block on stdout until all preceding regions drained.
SHARD_DIR="${WORK_TMPDIR}/shards"
mkdir -p "$SHARD_DIR"

# Number the regions so the shard filenames sort in FAI order.
NUMBERED_REGIONS="${WORK_TMPDIR}/regions.numbered.tsv"
awk '{printf "%07d\t%s\n", NR, $0}' "$REGIONS_FILE" > "$NUMBERED_REGIONS"

N_REGIONS=$(wc -l < "$NUMBERED_REGIONS")
echo "Launching freebayes on $N_REGIONS regions with -j $CPUS workers"

# Probe freebayes once up-front so a bad flag (e.g. --limit-coverage on
# freebayes < 1.3) fails loudly instead of killing 96 workers silently.
if [ -n "$EXTRA_ARGS" ]; then
    echo "Flag probe: $FB_CMD $EXTRA_ARGS --help"
    # shellcheck disable=SC2086
    if ! $FB_CMD $EXTRA_ARGS --help >/dev/null 2>&1; then
        # --help itself always exits 0 on sane CLIs; any non-zero here means
        # the flag parser choked before reaching --help.
        echo "Error: freebayes rejects extra-args '$EXTRA_ARGS' (unknown flag?)." >&2
        echo "       Installed freebayes: $($FB_CMD --version 2>&1 | head -1)" >&2
        exit 1
    fi
fi

/usr/bin/time -v parallel -j "$CPUS" --colsep '\t' --halt now,fail=1 \
    "$FB_CMD -f $REF_LOCAL $BAM_LOCAL --region {2} $EXTRA_ARGS > $SHARD_DIR/region-{1}.vcf" \
  :::: "$NUMBERED_REGIONS"

# Sanity-check shard count. An empty shard dir would silently produce an
# empty 28-byte BGZF through the cat pipeline below (cat-missing-glob is
# not fatal under set -eo pipefail because it's the first element of a
# pipe). Fail loudly instead.
N_SHARDS=$(find "$SHARD_DIR" -maxdepth 1 -name 'region-*.vcf' | wc -l)
if [ "$N_SHARDS" -ne "$N_REGIONS" ]; then
    echo "Error: produced $N_SHARDS region shards, expected $N_REGIONS." >&2
    exit 1
fi

# Concat shards in glob order (= FAI order by construction, because the
# shard filenames are zero-padded). `cat region-*.vcf` would inline all
# 70k+ filenames on the command line on HPRC-scale runs and hit Linux
# ARG_MAX (~2 MB) → `cat: Argument list too long`. Use `find | sort |
# xargs cat` which streams filenames across multiple cat invocations.
# awk keeps first header only and forces FILTER=PASS (FreeBayes emits "."
# by default). Merged VCF is built on scratch first, then mv'd to
# $OUTPUT_DIR — avoids half-written output hitting shared FS on cancel.
find "$SHARD_DIR" -maxdepth 1 -name 'region-*.vcf' -print0 \
  | LC_ALL=C sort -z \
  | xargs -0 cat \
  | awk 'BEGIN{OFS="\t"; p=1} /^#/{if(p)print; if(/^#CHROM/)p=0; next} {if($7==".") $7="PASS"; print}' \
  | bcftools annotate -x FORMAT/DPR \
  | bcftools reheader -s <(echo "$SAMPLE") \
  | bgzip > "$VCF_SCRATCH"
tabix -p vcf "$VCF_SCRATCH"

# Final sanity: bgzipped VCF must contain at least one non-header record.
# A 28-byte output is BGZF EOF only — a dead giveaway that the concat
# pipeline saw nothing. Cheap first check: file size. Second check: try
# to read one record, but disable pipefail for the `head -1` pipeline —
# `head -1` can close its input after the first record, giving bcftools
# SIGPIPE, which under set -o pipefail would abort the whole script.
VCF_SIZE=$(stat -c '%s' "$VCF_SCRATCH")
if [ "$VCF_SIZE" -le 100 ]; then
    echo "Error: produced VCF '$VCF_SCRATCH' is ${VCF_SIZE} bytes — concat saw nothing." >&2
    exit 1
fi
set +o pipefail
N_RECS=$(bcftools view -H "$VCF_SCRATCH" 2>/dev/null | head -1 | wc -l)
set -o pipefail
if [ "$N_RECS" -eq 0 ]; then
    echo "Error: produced VCF '$VCF_SCRATCH' has zero records (size ${VCF_SIZE})." >&2
    exit 1
fi

# Atomic-ish publish to shared output dir. mv on a same-filesystem target
# is atomic; cross-FS it's copy + unlink, which leaves a partial file
# visible for a moment — but the upstream failure-mode guarantees we only
# publish a complete VCF.
echo "Publishing merged VCF to $VCF"
mv "$VCF_SCRATCH" "$VCF"
mv "${VCF_SCRATCH}.tbi" "${VCF}.tbi"
