#!/usr/bin/env Rscript

# Chromosome density plot as ideograms for nested variants from TSV file
# Usage: ./chrom-density-tsv.R <input.tsv> <output.png> [title] [min_length] [bed_file] [scale] [--ref REF]

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

# Extract --censat flag if present
censat_file <- NULL
censat_idx <- which(args == "--censat")
if (length(censat_idx) > 0) {
  censat_file <- args[censat_idx + 1]
  args <- args[-c(censat_idx, censat_idx + 1)]
}

if (length(args) < 2) {
  cat("Usage: ./chrom-density-tsv.R <input.tsv> <output.png> [title] [min_length] [bed_file] [scale] [--ref REF]\n")
  cat("Example: ./chrom-density-tsv.R nesting.tsv density.png \"Nested Variants\" 50 refgaps.bed log1p --ref CHM13\n")
  cat("\nInput TSV format: Uses columns 5-7 (ref_contig, ref_start, ref_end)\n")
  cat("Scale options: log1p (default), sqrt, log, identity (linear)\n")
  cat("Ref options: CHM13, GRCh38, or omit for auto-detect\n")
  quit(status = 1)
}

input_tsv  <- args[1]
output_file <- args[2]
plot_title <- if (length(args) >= 3) args[3] else "Chromosome Density of Nested Variants"
min_length <- if (length(args) >= 4) as.numeric(args[4]) else 50
bed_file   <- if (length(args) >= 5) args[5] else NULL
scale_type <- if (length(args) >= 6) args[6] else "log1p"

cat("Reading TSV from:", input_tsv, "\n")

# Read TSV file — columns 5-7: ref_contig, ref_start, ref_end
tsv_data <- fread(input_tsv,
                  header = FALSE, sep = "\t",
                  select = c(5, 6, 7),
                  col.names = c("ref_contig", "ref_start", "ref_end"),
                  showProgress = TRUE)

cat("Read", nrow(tsv_data), "rows\n")

# Calculate length and filter
tsv_data[, ref_length := ref_end - ref_start]

cat("Filtering for regions >=", min_length, "bp\n")
tsv_data <- tsv_data[ref_length >= min_length]
cat("After length filter:", nrow(tsv_data), "rows\n")

# Extract chromosome from ref_contig (format: reference#chr)
tsv_data[, chromosome := sub(".*#", "", ref_contig)]

# Filter and order chromosomes
cat("Filtering to standard chromosomes\n")
result <- filter_standard_chroms(tsv_data, ref)
tsv_data     <- result$data
chrom_levels <- result$chrom_levels

cat("Filtered to", nrow(tsv_data), "rows on standard chromosomes\n")

if (nrow(tsv_data) == 0) {
  cat("ERROR: No data remaining after filtering\n")
  quit(status = 1)
}

# Get chromosome lengths
chrom_lengths <- get_chrom_lengths(ref)

# Compute density bins
cat("Calculating density with 1 Mb bins\n")
density_data <- compute_density_bins(tsv_data, chrom_lengths, chrom_levels)
cat("Created", nrow(density_data), "density bins\n")

# If no preset chrom_lengths, build from density data for plotting
if (is.null(chrom_lengths)) {
  chrom_lengths <- tsv_data[, .(length = max(ref_end, na.rm = TRUE)), by = chromosome]
  chrom_lengths <- as.data.frame(chrom_lengths)
}
chrom_lengths$chromosome <- factor(chrom_lengths$chromosome, levels = chrom_levels)

# Read BED overlay and optional censat track
bed_data <- read_bed_overlay(bed_file, chrom_levels)
censat_data <- read_bed_overlay(censat_file, chrom_levels)

# Build ideogram and save
p <- plot_ideogram(density_data, chrom_lengths, bed_data, plot_title, scale_type,
                   censat_data = censat_data)
save_and_summarize(p, output_file, as.data.frame(tsv_data))
