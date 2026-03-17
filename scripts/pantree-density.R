#!/usr/bin/env Rscript

# pantree-density.R — Chromosome density ideogram for pantree off-reference variants
#
# Usage: Rscript scripts/pantree-density.R <records.tsv> <output.png> [title] [--ref REF]
#
# Reads a pantree records TSV (from pantree-extract.py), filters to off-reference
# variants, and plots an ideogram of variant density on the base reference.

# Source shared functions
script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
  grep("--file=", commandArgs(trailingOnly = FALSE))
] |> sub("--file=", "", x = _), mustWork = FALSE))
source(file.path(script_dir, "chrom-density-common.R"))

args <- commandArgs(trailingOnly = TRUE)

# Extract --ref flag if present
ref <- NULL
ref_idx <- which(args == "--ref")
if (length(ref_idx) > 0) {
  ref <- args[ref_idx + 1]
  args <- args[-c(ref_idx, ref_idx + 1)]
}

if (length(args) < 2) {
  cat("Usage: Rscript pantree-density.R <records.tsv> <output.png> [title] [--ref REF]\n")
  quit(status = 1)
}

input_tsv   <- args[1]
output_file <- args[2]
plot_title  <- if (length(args) >= 3) args[3] else "Pantree Variant Density"

cat("Reading:", input_tsv, "\n")
dt <- data.table::fread(input_tsv)
cat("Total records:", nrow(dt), "\n")

# Keep only off-reference variants (matches deconstruct density plots)
dt <- dt[ref_context == "Off-reference"]
cat("Off-reference records:", nrow(dt), "\n")

if (nrow(dt) == 0) {
  cat("No off-reference records — creating empty plot.\n")
  file.create(output_file)
  quit(status = 0)
}

# Map columns to what chrom-density-common expects
dt[, chromosome := CHROM]
dt[, ref_start := POS]
dt[, ref_length := size]
dt[ref_length == 0, ref_length := 1]  # SNPs/MNPs: use 1 bp

# Filter and order chromosomes
result <- filter_standard_chroms(dt, ref)
dt           <- result$data
chrom_levels <- result$chrom_levels
cat("Chromosomes:", length(chrom_levels), "\n")

# Chromosome lengths
chrom_lengths <- get_chrom_lengths(ref)
if (is.null(chrom_lengths)) {
  chrom_lengths <- dt[, .(length = max(ref_start + ref_length, na.rm = TRUE)), by = chromosome]
  chrom_lengths <- as.data.frame(chrom_lengths)
}
chrom_lengths$chromosome <- factor(chrom_lengths$chromosome, levels = chrom_levels)

# Density bins
density_data <- compute_density_bins(dt, chrom_lengths, chrom_levels)

# Plot
p <- plot_ideogram(density_data, chrom_lengths, bed_data = NULL, plot_title)
save_and_summarize(p, output_file, as.data.frame(dt))
