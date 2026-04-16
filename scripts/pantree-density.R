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

# Extract --bed flag if present (gap regions to black out)
bed_file <- NULL
bed_idx <- which(args == "--bed")
if (length(bed_idx) > 0) {
  bed_file <- args[bed_idx + 1]
  args <- args[-c(bed_idx, bed_idx + 1)]
}

# Extract annotation track flags (--censat, --segdups, --genes) — same as chrom-density-tsv.R
track_display_names <- c(censat = "CenSat", segdups = "SegDups", genes = "Genes")
annot_track_files <- list()
for (flag in c("--censat", "--segdups", "--genes")) {
  idx <- which(args == flag)
  if (length(idx) > 0) {
    key <- sub("^--", "", flag)
    annot_track_files[[track_display_names[key]]] <- args[idx + 1]
    args <- args[-c(idx, idx + 1)]
  }
}

if (length(args) < 2) {
  cat("Usage: Rscript pantree-density.R <records.tsv> <output.png> [title] [--ref REF] [--bed BED] [--censat F] [--segdups F] [--genes F]\n")
  quit(status = 1)
}

input_tsv   <- args[1]
output_file <- args[2]
plot_title  <- if (length(args) >= 3) args[3] else "Pantree Off-Reference Variant Density"

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

# BED overlay (gap regions) and annotation tracks
bed_data <- read_bed_overlay(bed_file, chrom_levels)
annot_tracks <- lapply(annot_track_files, read_bed_overlay, chrom_levels = chrom_levels)

# Plot
p <- plot_ideogram(density_data, chrom_lengths, bed_data, plot_title,
                   annot_tracks = annot_tracks)
save_and_summarize(p, output_file, as.data.frame(dt))
