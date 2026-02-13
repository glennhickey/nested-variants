#!/usr/bin/env Rscript

# Chromosome density plot using augref segments table to map variants to reference
# For VCFs that lack RC/RS/RD INFO tags (vg call, DeepVariant).
# Usage: ./chrom-density-segs.R <input.vcf.gz> <segments.tsv> <output.png> [title] [min_length] [bed_file] [scale] [--ref REF] [--offref] [--bin-size N]

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

# Extract --bin-size flag if present
bin_size <- NULL
bin_idx <- which(args == "--bin-size")
if (length(bin_idx) > 0) {
  bin_size <- as.numeric(args[bin_idx + 1])
  args <- args[-c(bin_idx, bin_idx + 1)]
}

if (length(args) < 3) {
  cat("Usage: ./chrom-density-segs.R <input.vcf.gz> <segments.tsv> <output.png> [title] [min_length] [bed_file] [scale] [--ref REF] [--offref] [--bin-size N]\n")
  cat("Example: ./chrom-density-segs.R call.vcf.gz augref-segs.tsv density.png \"Density\" 0 refgaps.bed log1p --ref CHM13 --offref\n")
  cat("\nMaps variants to reference coordinates via the augref segments table.\n")
  cat("Handles both full (augref_REF#0#chr_N_alt) and short (chr_N_alt) CHROM names.\n")
  cat("  --offref       Keep only variants on alt contigs (CHROM ends with '_alt')\n")
  cat("  --bin-size N   Density bin size in bp (default: auto-scale to ~200 bins, max 1 Mb)\n")
  cat("Scale options: log1p (default), sqrt, log, identity (linear)\n")
  cat("Ref options: CHM13, GRCh38, or omit for auto-detect\n")
  quit(status = 1)
}

input_vcf    <- args[1]
segments_tsv <- args[2]
output_file  <- args[3]
plot_title   <- if (length(args) >= 4) args[4] else "Chromosome Density of Variants"
min_length   <- if (length(args) >= 5) as.numeric(args[5]) else 0
bed_file     <- if (length(args) >= 6) args[6] else NULL
scale_type   <- if (length(args) >= 7) args[7] else "log1p"

# Read augref segments table (cols 4-7: augref_path, ref_path, ref_start, ref_end)
cat("Reading segments table from:", segments_tsv, "\n")
segs <- fread(segments_tsv, header = FALSE, sep = "\t",
              select = c(4, 5, 6, 7),
              col.names = c("augref_path", "ref_path", "seg_ref_start", "seg_ref_end"))
segs <- unique(segs, by = "augref_path")
cat("Read", nrow(segs), "unique augref segments\n")

# Read VCF body — only need CHROM and POS
cat("Reading VCF from:", input_vcf, "\n")
vcf_data <- fread(cmd = paste0("bcftools view -c1 -H '", input_vcf, "'"),
                  header = FALSE, sep = "\t",
                  select = c(1, 2),
                  col.names = c("chrom", "pos"),
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

# Join with segments table to get reference coordinates
# VCF CHROM may use full names (augref_REF#0#chr_N_alt, e.g. DeepVariant)
# or short names (chr_N_alt, e.g. vg call) — try both
vcf_data <- merge(vcf_data, segs, by.x = "chrom", by.y = "augref_path")
if (nrow(vcf_data) == 0) {
  # Fall back to matching short names (part after last #)
  cat("No full-name matches; trying short CHROM names\n")
  segs[, augref_short := sub(".*#", "", augref_path)]
  vcf_data <- fread(cmd = paste0("bcftools view -c1 -H '", input_vcf, "'"),
                    header = FALSE, sep = "\t",
                    select = c(1, 2),
                    col.names = c("chrom", "pos"),
                    showProgress = FALSE)
  if (offref) vcf_data <- vcf_data[grepl("_alt$", chrom)]
  vcf_data <- merge(vcf_data, segs, by.x = "chrom", by.y = "augref_short")
}
cat("Mapped", nrow(vcf_data), "variants to reference coordinates\n")

if (nrow(vcf_data) == 0) {
  cat("ERROR: No variants could be mapped to reference coordinates\n")
  quit(status = 1)
}

# Extract chromosome from ref_path (format: reference#haplotype#chr)
vcf_data[, chromosome := sub(".*#", "", ref_path)]

# Use segment ref_start as the reference position for density binning
vcf_data[, ref_start := seg_ref_start]
vcf_data[, ref_length := 1L]

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
cat("Calculating density bins\n")
density_data <- compute_density_bins(vcf_data, chrom_lengths, chrom_levels, bin_size)
cat("Created", nrow(density_data), "density bins\n")

# If no preset chrom_lengths, build from data for plotting
if (is.null(chrom_lengths)) {
  chrom_lengths <- vcf_data[, .(length = max(seg_ref_end, na.rm = TRUE)), by = chromosome]
  chrom_lengths <- as.data.frame(chrom_lengths)
}
chrom_lengths$chromosome <- factor(chrom_lengths$chromosome, levels = chrom_levels)

# Read BED overlay
bed_data <- read_bed_overlay(bed_file, chrom_levels)

# Build ideogram and save
p <- plot_ideogram(density_data, chrom_lengths, bed_data, plot_title, scale_type)
save_and_summarize(p, output_file, as.data.frame(vcf_data))
