#!/usr/bin/env Rscript

# Chromosome density plot as ideograms for variants using reference coordinates from VCF INFO field
# Usage: ./chrom-density-vcf-ideogram.R <input.vcf.gz> <output.png> [title] [min_length]

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(data.table)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: ./chrom-density-vcf-ideogram.R <input.vcf.gz> <output.png> [title] [min_length] [bed_file]\n")
  cat("Example: ./chrom-density-vcf-ideogram.R off-ref-calls/hprc-v2.0-mc-chm13.nested.95.offref.vcf.gz chrom-density-ideogram.png \"Title\" 50 hprc-v2.0-mc-chm13.refgaps.bed\n")
  quit(status = 1)
}

input_vcf <- args[1]
output_file <- args[2]
plot_title <- if (length(args) >= 3) args[3] else "Chromosome Density of Variants"
min_length <- if (length(args) >= 4) as.numeric(args[4]) else 50
bed_file <- if (length(args) >= 5) args[5] else NULL

cat("Reading VCF from:", input_vcf, "\n")

# Read VCF using data.table (skip header lines starting with #)
# We need columns: CHROM (1), POS (2), and INFO (8)
vcf_data <- fread(cmd = paste0("zcat ", input_vcf, " | grep -v '^#'"),
                  header = FALSE, sep = "\t",
                  select = c(1, 2, 8),
                  col.names = c("chrom", "pos", "info"),
                  showProgress = TRUE)

cat("Read", nrow(vcf_data), "VCF variants\n")

# Extract RC, RS, RD (note: user said RE but it's actually RD in the VCF) from INFO field
cat("Extracting reference coordinates from INFO field...\n")

# Function to extract INFO field value
extract_info <- function(info_string, tag) {
  pattern <- paste0(tag, "=([^;]+)")
  matches <- regmatches(info_string, regexec(pattern, info_string))
  sapply(matches, function(x) if(length(x) >= 2) x[2] else NA_character_)
}

vcf_data$ref_contig <- extract_info(vcf_data$info, "RC")
vcf_data$ref_start <- as.integer(extract_info(vcf_data$info, "RS"))
vcf_data$ref_end <- as.integer(extract_info(vcf_data$info, "RD"))
vcf_data$ref_length <- as.integer(extract_info(vcf_data$info, "RL"))

# Remove rows with missing reference coordinates
vcf_data <- vcf_data %>%
  filter(!is.na(ref_contig) & !is.na(ref_start) & !is.na(ref_end))

cat("Extracted reference coordinates for", nrow(vcf_data), "variants\n")

# Filter by minimum length
cat("Filtering for regions >=", min_length, "bp\n")
vcf_data <- vcf_data %>%
  filter(ref_length >= min_length)

cat("After length filter:", nrow(vcf_data), "rows\n")

# Extract chromosome from ref_contig (format: reference#haplotype#chr)
vcf_data$chromosome <- sub(".*#", "", vcf_data$ref_contig)

# Filter to standard chromosomes (chr1-chr22, chrX, chrY, chrM)
cat("Filtering to standard chromosomes\n")
vcf_data <- vcf_data %>%
  filter(grepl("^chr", chromosome)) %>%
  mutate(
    chrom_num = case_when(
      chromosome == "chrX" ~ 23,
      chromosome == "chrY" ~ 24,
      chromosome == "chrM" ~ 25,
      grepl("^chr[0-9]+$", chromosome) ~ as.numeric(sub("chr", "", chromosome)),
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(chrom_num))

# Create chromosome factor with proper ordering
chrom_levels <- paste0("chr", c(1:22, "X", "Y", "M"))
vcf_data$chromosome <- factor(vcf_data$chromosome, levels = chrom_levels)
vcf_data <- vcf_data %>% filter(!is.na(chromosome))

cat("Filtered to", nrow(vcf_data), "rows on standard chromosomes\n")

if (nrow(vcf_data) == 0) {
  cat("ERROR: No data remaining after filtering\n")
  quit(status = 1)
}

# Chromosome lengths (CHM13 telomere-to-telomere assembly)
chrom_lengths <- data.frame(
  chromosome = factor(paste0("chr", c(1:22, "X", "Y", "M")), levels = chrom_levels),
  length = c(
    248387328, 242696752, 201105948, 193574945, 182045439,
    172126628, 160567428, 146259331, 150617247, 134758134,
    135127769, 133324548, 113566686, 101991189, 99753195,
    96330374, 84276897, 80542538, 61707364, 66210255,
    45090682, 51324926, 154259566, 62460029, 16569
  )
)

# Centromeres or other regions will come from BED file if provided

# Calculate density using bins
bin_size <- 1e6  # 1 Mb bins
cat("Calculating density with", bin_size / 1e6, "Mb bins\n")

# Create bins for each chromosome (use data.table for speed)
data_dt <- as.data.table(vcf_data)
chrom_lengths_dt <- as.data.table(chrom_lengths)

# Fast binning with data.table - use ref_start for binning position
data_dt[, bin := floor(ref_start / bin_size)]
density_data <- data_dt[, .(count = .N, total_bp = sum(ref_length)), by = .(chromosome, bin)]

# Keep all bins (don't filter out zeros) to show full chromosome coverage
# Bins with zero will show as background gray

# Join with chromosome lengths
density_data <- merge(density_data, chrom_lengths_dt, by = "chromosome")

# Calculate bin positions
density_data[, `:=`(
  bin_start = bin * bin_size,
  bin_end = pmin((bin + 1) * bin_size, length)
)]

# Convert back to data.frame for ggplot
density_data <- as.data.frame(density_data)

cat("Created", nrow(density_data), "density bins\n")

# Read BED file if provided
bed_data <- NULL
if (!is.null(bed_file) && file.exists(bed_file)) {
  cat("Reading BED file from:", bed_file, "\n")
  bed_data <- fread(bed_file, header = FALSE, sep = "\t",
                    col.names = c("contig", "start", "end"),
                    colClasses = c("character", "integer", "integer"))

  # Extract chromosome from contig (same format as VCF)
  bed_data$chromosome <- sub(".*#", "", bed_data$contig)

  # Filter to standard chromosomes
  bed_data <- bed_data %>%
    filter(chromosome %in% chrom_levels)

  bed_data$chromosome <- factor(bed_data$chromosome, levels = chrom_levels)
  bed_data <- bed_data %>% filter(!is.na(chromosome))

  cat("Read", nrow(bed_data), "BED intervals\n")
}

# Create ideogram-style plot
# Calculate layout parameters
ideogram_height <- 0.6  # height of each ideogram bar
ideogram_spacing <- 1.0 # spacing between ideograms
y_positions <- seq(length(chrom_levels), 1, by = -1) * ideogram_spacing

# Create y position mapping
chrom_y_map <- data.frame(
  chromosome = factor(chrom_levels, levels = chrom_levels),
  y_pos = y_positions
)

# Merge y positions into data
density_data <- merge(density_data, chrom_y_map, by = "chromosome")
chrom_lengths <- merge(chrom_lengths, chrom_y_map, by = "chromosome")
if (!is.null(bed_data) && nrow(bed_data) > 0) {
  bed_data <- merge(bed_data, chrom_y_map, by = "chromosome")
}

# Create the ideogram plot
p <- ggplot() +
  # Add chromosome backgrounds (white with black outline)
  geom_rect(data = chrom_lengths,
            aes(xmin = 0, xmax = length,
                ymin = y_pos - ideogram_height/2,
                ymax = y_pos + ideogram_height/2),
            fill = "white", color = "black", linewidth = 0.5) +
  # Add density as colored tiles (colored by variant count)
  geom_rect(data = density_data,
            aes(xmin = bin_start, xmax = bin_end,
                ymin = y_pos - ideogram_height/2,
                ymax = y_pos + ideogram_height/2,
                fill = count)) +
  # Add BED intervals (e.g., centromeres, reference gaps) as black bars if provided
  {if (!is.null(bed_data) && nrow(bed_data) > 0)
    geom_rect(data = bed_data,
              aes(xmin = start, xmax = end,
                  ymin = y_pos - ideogram_height/2,
                  ymax = y_pos + ideogram_height/2),
              fill = "black", alpha = 0.7)
  } +
  scale_fill_gradientn(colors = c("#FFFACD", "yellow", "orange", "#FF4500", "#FF0000"),
                       trans = "sqrt",
                       name = "Variant count\n(sqrt scale)",
                       labels = comma,
                       na.value = "white") +
  scale_y_continuous(breaks = y_positions,
                     labels = chrom_levels,
                     expand = c(0.02, 0)) +
  scale_x_continuous(labels = unit_format(unit = "Mb", scale = 1e-6),
                     expand = c(0.01, 0)) +
  labs(
    title = plot_title,
    x = "Position",
    y = "Chromosome"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 9, hjust = 1),
    axis.text.x = element_text(size = 8),
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Save the plot
cat("Saving plot to:", output_file, "\n")
ggsave(output_file, p, width = 12, height = 10, dpi = 300, bg = "white")

cat("Done!\n")

# Print summary statistics
cat("\nSummary by chromosome:\n")
summary_stats <- vcf_data %>%
  group_by(chromosome) %>%
  summarise(
    n_variants = n(),
    total_bp = sum(ref_length),
    .groups = "drop"
  ) %>%
  arrange(chromosome)

print(summary_stats)
