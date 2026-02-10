#!/usr/bin/env Rscript

# Chromosome density plot as ideograms for variants using reference coordinates from VCF INFO field
# Usage: ./chrom-density-vcf.R <input.vcf.gz> <output.png> [title] [min_length] [bed_file] [scale] [--ref REF] [--offref]

# Source shared functions (relative to this script's location)
script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
  grep("--file=", commandArgs(trailingOnly = FALSE))
] |> sub("--file=", "", x = _), mustWork = FALSE))
source(file.path(script_dir, "chrom-density-common.R"))

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Extract --ref flag if present
ref <- NULL
ref_idx <- which(args == "--ref")
if (length(ref_idx) > 0) {
  ref <- args[ref_idx + 1]
  args <- args[-c(ref_idx, ref_idx + 1)]
}

# Extract --offref flag if present (keep only alt contigs)
offref <- FALSE
offref_idx <- which(args == "--offref")
if (length(offref_idx) > 0) {
  offref <- TRUE
  args <- args[-offref_idx]
}

if (length(args) < 2) {
  cat("Usage: ./chrom-density-vcf.R <input.vcf.gz> <output.png> [title] [min_length] [bed_file] [scale] [--ref REF] [--offref]\n")
  cat("Example: ./chrom-density-vcf.R decon.vcf.gz density.png \"Title\" 50 refgaps.bed log1p --ref CHM13 --offref\n")
  cat("\nPlots variant density using reference coordinates from INFO fields (RC, RS, RD, RL).\n")
  cat("  --offref   Keep only variants on alt contigs (CHROM ends with '_alt')\n")
  cat("Scale options: log1p (default), sqrt, log, identity (linear)\n")
  cat("Ref options: CHM13, GRCh38, or omit for auto-detect\n")
  quit(status = 1)
}

input_vcf   <- args[1]
output_file <- args[2]
plot_title  <- if (length(args) >= 3) args[3] else "Chromosome Density of Variants"
min_length  <- if (length(args) >= 4) as.numeric(args[4]) else 50
bed_file    <- if (length(args) >= 5) args[5] else NULL
scale_type  <- if (length(args) >= 6) args[6] else "log1p"

cat("Reading VCF from:", input_vcf, "\n")

# Read VCF body (skip # header lines), extract CHROM, POS, INFO columns
vcf_data <- fread(cmd = paste0("zcat ", input_vcf, " | grep -v '^#'"),
                  header = FALSE, sep = "\t",
                  select = c(1, 2, 8),
                  col.names = c("chrom", "pos", "info"),
                  showProgress = TRUE)

cat("Read", nrow(vcf_data), "VCF variants\n")

# Filter to off-reference (alt) contigs if requested
if (offref) {
  vcf_data <- vcf_data[grepl("_alt$", chrom)]
  cat("Filtered to", nrow(vcf_data), "off-reference (alt contig) variants\n")
  if (nrow(vcf_data) == 0) {
    cat("ERROR: No off-reference variants found\n")
    quit(status = 1)
  }
}

# Extract reference coordinates from INFO field (RC, RS, RD, RL)
cat("Extracting reference coordinates from INFO field...\n")

extract_info <- function(info_string, tag) {
  pattern <- paste0(tag, "=([^;]+)")
  matches <- regmatches(info_string, regexec(pattern, info_string))
  sapply(matches, function(x) if (length(x) >= 2) x[2] else NA_character_)
}

vcf_data[, ref_contig := extract_info(info, "RC")]
vcf_data[, ref_start  := as.integer(extract_info(info, "RS"))]
vcf_data[, ref_end    := as.integer(extract_info(info, "RD"))]
# Use RL if present, otherwise compute from RD - RS
rl_values <- extract_info(vcf_data$info, "RL")
if (all(is.na(rl_values))) {
  vcf_data[, ref_length := ref_end - ref_start]
} else {
  vcf_data[, ref_length := as.integer(rl_values)]
}

# Drop rows missing reference coordinates
vcf_data <- vcf_data[!is.na(ref_contig) & !is.na(ref_start) & !is.na(ref_end)]
cat("Extracted reference coordinates for", nrow(vcf_data), "variants\n")

# Filter by minimum length
cat("Filtering for regions >=", min_length, "bp\n")
vcf_data <- vcf_data[ref_length >= min_length]
cat("After length filter:", nrow(vcf_data), "rows\n")

# Extract chromosome from ref_contig (format: reference#haplotype#chr)
vcf_data[, chromosome := sub(".*#", "", ref_contig)]

# Filter and order chromosomes
cat("Filtering to standard chromosomes\n")
result <- filter_standard_chroms(vcf_data, ref)
vcf_data     <- result$data
chrom_levels <- result$chrom_levels

cat("Filtered to", nrow(vcf_data), "rows on standard chromosomes\n")

if (nrow(vcf_data) == 0) {
  cat("ERROR: No data remaining after filtering\n")
  quit(status = 1)
}

# Get chromosome lengths
chrom_lengths <- get_chrom_lengths(ref)

# Compute density bins
cat("Calculating density with 1 Mb bins\n")
density_data <- compute_density_bins(vcf_data, chrom_lengths, chrom_levels)
cat("Created", nrow(density_data), "density bins\n")

# If no preset chrom_lengths, build from data for plotting
if (is.null(chrom_lengths)) {
  chrom_lengths <- vcf_data[, .(length = max(ref_end, na.rm = TRUE)), by = chromosome]
  chrom_lengths <- as.data.frame(chrom_lengths)
}
chrom_lengths$chromosome <- factor(chrom_lengths$chromosome, levels = chrom_levels)

# Read BED overlay
bed_data <- read_bed_overlay(bed_file, chrom_levels)

# Build ideogram and save
p <- plot_ideogram(density_data, chrom_lengths, bed_data, plot_title, scale_type)
save_and_summarize(p, output_file, as.data.frame(vcf_data))
