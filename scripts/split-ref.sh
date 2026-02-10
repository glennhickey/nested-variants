#!/bin/bash

set -euo pipefail

# Function to display help message
show_help() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS] -v VCF_FILE -p PREFIX

Split a VCF file into three categories: on-reference, nested-reference, and off-reference variants.
Outputs are generated in parallel and placed in the current directory.

Augmented reference contigs follow the naming convention from vg's augref system:
  - On-reference contigs:  PREFIX#0#chr1, PREFIX#0#chr2, ...
  - Off-reference contigs: PREFIX#0#chr1_42_alt, PREFIX#0#chr2_7_alt, ...
Off-reference (alt) contigs are identified by the "_alt" suffix (matching vg's is_augref_name()).

Required arguments:
    -v, --vcf FILE      Input VCF file (can be in any directory, must end in .vcf.gz)
    -p, --prefix STR    Augmented reference sample name (e.g., 'augref_CHM13')

Optional arguments:
    -t, --threads NUM   Number of threads for bgzip compression (default: 3)
    -h, --help          Display this help message and exit

Output files (created in current directory):
    <basename>.onref.vcf.gz        Variants on reference contigs (no _alt suffix)
    <basename>.nestedref.vcf.gz    Nested reference variants (on-ref contigs with LV>0)
    <basename>.offref.vcf.gz       Variants on off-reference (alt) contigs (_alt suffix)

Each output file is automatically indexed with tabix.

Examples:
    $(basename "$0") -v /path/to/file.vcf.gz -p augref_CHM13
    $(basename "$0") --vcf data/variants.vcf.gz --prefix augref_GRCh38 --threads 4

EOF
    exit 0
}

# Function to display error message
error_exit() {
    echo "Error: $1" >&2
    echo "Try '$(basename "$0") --help' for more information." >&2
    exit 1
}

# Default values
THREADS=3
VCF=""
PREFIX=""

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -v|--vcf)
            VCF="$2"
            shift 2
            ;;
        -p|--prefix)
            PREFIX="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            ;;
        *)
            error_exit "Unknown option: $1"
            ;;
    esac
done

# Validate required arguments
if [[ -z "$VCF" ]]; then
    error_exit "VCF file is required (use -v or --vcf)"
fi

if [[ -z "$PREFIX" ]]; then
    error_exit "Prefix is required (use -p or --prefix)"
fi

# Validate VCF file exists
if [[ ! -f "$VCF" ]]; then
    error_exit "VCF file does not exist: $VCF"
fi

# Validate VCF file extension
if [[ ! "$VCF" =~ \.vcf\.gz$ ]]; then
    error_exit "VCF file must end with .vcf.gz: $VCF"
fi

# Validate threads is a number
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
    error_exit "Threads must be a positive integer: $THREADS"
fi

# Get basename for output files (remove path and .vcf.gz extension)
BASENAME=$(basename "$VCF" .vcf.gz)

echo "Processing VCF file: $VCF"
echo "Using prefix: $PREFIX"
echo "Output basename: $BASENAME"
echo "Compression threads: $THREADS"
echo "Generating outputs in parallel..."

# Alt contig detection: contigs ending with _{N}_alt (matching vg's is_augref_name())
# The grep pattern matches a tab-separated CHROM field ending in _<digits>_alt
ALT_PATTERN='_[0-9]\+_alt	'

# Function to generate onref output (contigs that do NOT end with _alt)
generate_onref() {
    local vcf=$1
    local pattern=$2
    local threads=$3
    local basename=$4

    local out_name="${basename}.onref.vcf.gz"
    echo "  [onref] Starting: $out_name"

    bcftools view "$vcf" -h | bgzip --threads "$threads" > "$out_name"
    bcftools view "$vcf" -H | { grep -v "${pattern}" || true; } | bgzip --threads "$threads" >> "$out_name"
    tabix -fp vcf "$out_name"

    echo "  [onref] Complete: $out_name"
}

# Function to generate nestedref output (on-ref contigs with LV>0)
generate_nestedref() {
    local vcf=$1
    local pattern=$2
    local threads=$3
    local basename=$4

    local out_name="${basename}.nestedref.vcf.gz"
    echo "  [nestedref] Starting: $out_name"

    bcftools view "$vcf" -h | bgzip --threads "$threads" > "$out_name"
    bcftools view "$vcf" -H -i "LV>0" | { grep -v "${pattern}" || true; } | bgzip --threads "$threads" >> "$out_name"
    tabix -fp vcf "$out_name"

    echo "  [nestedref] Complete: $out_name"
}

# Function to generate offref output (contigs ending with _alt)
generate_offref() {
    local vcf=$1
    local pattern=$2
    local threads=$3
    local basename=$4

    local out_name="${basename}.offref.vcf.gz"
    echo "  [offref] Starting: $out_name"

    bcftools view "$vcf" -h | bgzip --threads "$threads" > "$out_name"
    bcftools view "$vcf" -H | { grep "${pattern}" || true; } | bgzip --threads "$threads" >> "$out_name"
    tabix -fp vcf "$out_name"

    echo "  [offref] Complete: $out_name"
}

# Export functions and variables for parallel execution
export -f generate_onref generate_nestedref generate_offref

# Run all three processes in parallel
generate_onref "$VCF" "$ALT_PATTERN" "$THREADS" "$BASENAME" &
PID1=$!

generate_nestedref "$VCF" "$ALT_PATTERN" "$THREADS" "$BASENAME" &
PID2=$!

generate_offref "$VCF" "$ALT_PATTERN" "$THREADS" "$BASENAME" &
PID3=$!

# Wait for all background processes to complete
wait $PID1 || error_exit "onref generation failed"
wait $PID2 || error_exit "nestedref generation failed"
wait $PID3 || error_exit "offref generation failed"

echo ""
echo "All outputs generated successfully!"
echo "Output files:"
echo "  ${BASENAME}.onref.vcf.gz"
echo "  ${BASENAME}.nestedref.vcf.gz"
echo "  ${BASENAME}.offref.vcf.gz"
