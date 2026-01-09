#!/bin/bash

################################################################################
# slurm-deconstruct.sh
#
# Description:
#   Runs vg deconstruct in parallel on SLURM to generate VCF files from VG
#   (variation graph) files. This script processes multiple VG files
#   (typically one per chromosome) concurrently using SLURM job scheduling,
#   extracts variant calls relative to a specified reference path, and
#   concatenates the results into a single multi-chromosome VCF file.
##
# Usage:
#   slurm-deconstruct.sh --vg <file1.vg> [--vg <file2.vg> ...] \
#                        --ref <ref> --L <L> --out-dir <dir> --out-name <name> \
#                        [--cpus N] [--mem size] [--time HH:MM:SS] [--partition name]
#
################################################################################

set -e

# Initialize variables
VG_FILES=()
REF=""
L=""
OUTPUT_DIR="."
OUTPUT_NAME=""

# SLURM resource defaults
CPUS="20"
MEM="200gb"
TIME="16:00:00"
PARTITION="long"
JOB_NAME="deconstruct"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --vg)
            # Expand glob patterns and add all matching files
            # Enable extended globbing for patterns like !(*.d9).vg
            shopt -s nullglob extglob
            for file in $2; do
                VG_FILES+=("$file")
            done
            shopt -u nullglob extglob
            shift 2
            ;;
        --ref)
            REF="$2"
            shift 2
            ;;
        --L)
            L="$2"
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
        -h|--help)
            echo "Usage: $0 --vg <file.vg|pattern> [--vg <file.vg|pattern> ...] --ref <ref> --L <L> --out-dir <dir> --out-name <name> [SLURM options]"
            echo ""
            echo "Required Options:"
            echo "  --vg <file>      VG file or glob pattern (e.g., chr*.vg) to process"
            echo "                   Can be specified multiple times"
            echo "  --ref <ref>      Reference name"
            echo "  --L <L>          L parameter for nesting"
            echo "  --out-dir <dir>  Output directory for intermediate files"
            echo "  --out-name <name> Output name for final VCF file"
            echo ""
            echo "SLURM Resource Options (optional, with defaults):"
            echo "  --cpus <N>       CPUs per task (default: 20)"
            echo "  --mem <size>     Memory per job (default: 200gb)"
            echo "  --time <time>    Wall clock limit (default: 16:00:00)"
            echo "  --partition <p>  SLURM partition/queue (default: long)"
            echo ""
            echo "Examples:"
            echo "  # Basic usage"
            echo "  $0 --vg chr*.vg --ref GRCh38 --L 50 --out-dir ./output --out-name result.vcf.gz"
            echo ""
            echo "  # With custom SLURM resources"
            echo "  $0 --vg chr*.vg --ref GRCh38 --L 50 --out-dir ./output --out-name result.vcf.gz --cpus 16 --mem 100gb --time 8:00:00"
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
if [ ${#VG_FILES[@]} -eq 0 ]; then
    echo "Error: At least one --vg file must be specified"
    exit 1
fi

if [ -z "$REF" ]; then
    echo "Error: --ref is required"
    exit 1
fi

if [ -z "$L" ]; then
    echo "Error: --L is required"
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

# Extract base name for FASTA output
OUTPUT_FA="${OUTPUT_NAME%.vcf.gz}.fa"

set -x

mkdir -p "$OUTPUT_DIR"

# Create a temporary working directory for intermediate files
WORK_DIR="${OUTPUT_DIR}/tmp_deconstruct_$$"
mkdir -p "$WORK_DIR"
rm -rf "${WORK_DIR}"/*.vcf*

# Process each VG file
for VG in "${VG_FILES[@]}"; do
    BASE=$(basename "$VG")
    if [[ $BASE != "chrEBV.vg" ]]; then
	BASE=${BASE%.vg}
	VCF="${WORK_DIR}/${BASE}.${REF}.${L}.vcf.gz"
	FASTA="${WORK_DIR}/${BASE}.${REF}.${L}.fa"

	# Build the command to run
	CMD="vg deconstruct \"$VG\" -R -n -L ${L} -P ${REF} -f \"${FASTA}\" -t ${CPUS} | bgzip > \"${VCF}\" && tabix -fp vcf \"${VCF}\""

	# Submit SLURM job with resource requirements (no log files)
	sbatch -W \
	    --job-name="${JOB_NAME}" \
	    --partition="${PARTITION}" \
	    --nodes=1 \
	    --ntasks=1 \
	    --cpus-per-task="${CPUS}" \
	    --mem="${MEM}" \
	    --time="${TIME}" \
	    --output=/dev/null \
	    --wrap="$CMD" &
    fi
done

wait

cat "${WORK_DIR}"/*.fa | bgzip > "${OUTPUT_DIR}/${OUTPUT_FA}.gz"
cat "${WORK_DIR}"/*.fa.nesting.tsv | bgzip > "${OUTPUT_DIR}/${OUTPUT_FA}.nesting.tsv.gz"

# Handle sample normalization for chromosomes that may have missing samples
for VG in "${VG_FILES[@]}"; do
    BASE=$(basename "$VG")
    if [[ $BASE == *"X.vg" || $BASE == *"Y.vg" || $BASE == *"M.vg" || $BASE == *[Oo]ther.vg ]]; then
	BASE=${BASE%.vg}
	VCF="${WORK_DIR}/${BASE}.${REF}.${L}.vcf.gz"
	CHR1_VCF="${WORK_DIR}/chr1.${REF}.${L}.vcf.gz"
	# get the samples from chr1
	bcftools query -l "${CHR1_VCF}" | sort > "${WORK_DIR}/all-samples"
	# get the samples from the vcf
	bcftools query -l "${VCF}"  | sort > "${WORK_DIR}/chrom.samples"
	if [[ $(diff "${WORK_DIR}/all-samples" "${WORK_DIR}/chrom.samples") != 0 ]]; then
	    # samples that are missing from this chromosome
	    comm -32 "${WORK_DIR}/all-samples" "${WORK_DIR}/chrom.samples" | awk '{print $1}' > "${WORK_DIR}/missing-samples"
	    bcftools view -h "${CHR1_VCF}" -S "${WORK_DIR}/missing-samples" | bgzip > "${WORK_DIR}/missing-header.vcf.gz"
	    tabix -fp vcf  "${WORK_DIR}/missing-header.vcf.gz"
	    # add these samples to the header
	    bcftools merge "${VCF}" "${WORK_DIR}/missing-header.vcf.gz" -Oz > "${VCF}.merge.vcf.gz"
	    tabix -fp vcf  "${VCF}.merge.vcf.gz"
	    # finally, sort these samples to be same as chr1
	    bcftools view "${VCF}.merge.vcf.gz" -S "${WORK_DIR}/all-samples" -Oz > "${WORK_DIR}/${BASE}.${REF}.${L}.vcf.gz"
	    tabix -fp vcf "${WORK_DIR}/${BASE}.${REF}.${L}.vcf.gz"
	    # remove temp stuff since we're using wildcards below
	    rm -f  "${VCF}.merge.vcf.gz" "${WORK_DIR}/missing-header.vcf.gz" "${WORK_DIR}/missing-header.vcf.gz.tbi"
	fi
    fi
done

bcftools concat "${WORK_DIR}"/*.vcf.gz | bgzip > "${OUTPUT_DIR}/${OUTPUT_NAME}"
tabix -fp vcf "${OUTPUT_DIR}/${OUTPUT_NAME}"

# Clean up temporary working directory
rm -rf "${WORK_DIR}"
