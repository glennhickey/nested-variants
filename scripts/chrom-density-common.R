#!/usr/bin/env Rscript

# chrom-density-common.R — Shared functions for chromosome density ideogram plots
#
# Sourced by chrom-density-tsv.R and chrom-density-vcf.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(data.table)
})

#' Get chromosome lengths for a known reference assembly
#'
#' @param ref Reference name: "CHM13", "GRCh38", or NULL/other for auto-detect
#' @return data.frame with columns (chromosome, length), or NULL for auto-detect
get_chrom_lengths <- function(ref = NULL) {
  if (is.null(ref)) return(NULL)

  ref_upper <- toupper(ref)

  if (ref_upper == "CHM13") {
    # CHM13 telomere-to-telomere assembly (T2T-CHM13v2.0)
    chrom_names <- paste0("chr", c(1:22, "X", "Y"))
    data.frame(
      chromosome = factor(chrom_names, levels = chrom_names),
      length = c(
        248387328, 242696752, 201105948, 193574945, 182045439,
        172126628, 160567428, 146259331, 150617247, 134758134,
        135127769, 133324548, 113566686, 101991189, 99753195,
        96330374, 84276897, 80542538, 61707364, 66210255,
        45090682, 51324926, 154259566, 62460029
      )
    )
  } else if (ref_upper == "GRCH38") {
    # GRCh38 primary assembly
    chrom_names <- paste0("chr", c(1:22, "X", "Y"))
    data.frame(
      chromosome = factor(chrom_names, levels = chrom_names),
      length = c(
        248956422, 242193529, 198295559, 190214555, 181538259,
        170805979, 159345973, 145138636, 138394717, 133797422,
        135086622, 133275309, 114364328, 107043718, 101991189,
        90338345, 83257441, 80373285, 58617616, 64444167,
        46709983, 50818468, 156040895, 57227415
      )
    )
  } else {
    NULL
  }
}

#' Filter and order chromosomes from data
#'
#' When ref is a known human assembly, keeps chr1-22,X,Y in standard order.
#' Otherwise, auto-detects chromosome names from the data and sorts naturally.
#'
#' @param data data.table with a "chromosome" column
#' @param ref Reference name or NULL
#' @return list(data, chrom_levels) — filtered data and ordered chromosome level names
filter_standard_chroms <- function(data, ref = NULL) {
  ref_upper <- if (!is.null(ref)) toupper(ref) else NULL

  if (!is.null(ref_upper) && ref_upper %in% c("CHM13", "GRCH38")) {
    # Human reference: keep chr1-22, chrX, chrY
    chrom_levels <- paste0("chr", c(1:22, "X", "Y"))
    data <- data[chromosome %in% chrom_levels]
    data[, chromosome := factor(chromosome, levels = chrom_levels)]
  } else {
    # Auto-detect: use all chromosomes found in data, sort naturally
    chroms <- sort(unique(data$chromosome))
    # Attempt natural sort: numeric parts sorted numerically
    chrom_levels <- chroms[order(nchar(chroms), chroms)]
    data[, chromosome := factor(chromosome, levels = chrom_levels)]
  }

  list(data = data, chrom_levels = chrom_levels)
}

#' Compute density bins for ideogram plotting
#'
#' @param data data.table with chromosome, ref_start, ref_length columns
#' @param chrom_lengths data.frame with chromosome, length columns (or NULL for auto)
#' @param chrom_levels character vector of ordered chromosome names
#' @param bin_size Bin size in bp, or NULL for auto-scale (~200 bins across longest chrom, capped at 1 Mb)
#' @return data.frame ready for plotting with bin_start, bin_end, count, total_bp, length
compute_density_bins <- function(data, chrom_lengths, chrom_levels, bin_size = NULL) {
  data_dt <- as.data.table(data)

  # If no chrom_lengths provided, compute from data
  if (is.null(chrom_lengths)) {
    chrom_lengths <- data_dt[, .(length = max(ref_start + ref_length, na.rm = TRUE)),
                             by = chromosome]
  }
  chrom_lengths_dt <- as.data.table(chrom_lengths)

  # Auto-scale bin size if not specified: ~200 bins across longest chrom, capped at 1 Mb
  if (is.null(bin_size)) {
    max_len <- max(chrom_lengths_dt$length, na.rm = TRUE)
    bin_size <- min(1e6, max(1000, ceiling(max_len / 200)))
    cat("Auto-scaled bin size:", format(bin_size, big.mark = ","), "bp\n")
  }

  # Bin by ref_start position
  data_dt[, bin := floor(ref_start / bin_size)]
  density_data <- data_dt[, .(count = .N, total_bp = sum(ref_length)),
                          by = .(chromosome, bin)]

  # Join with chromosome lengths
  density_data <- merge(density_data, chrom_lengths_dt, by = "chromosome")

  # Calculate bin positions
  density_data[, `:=`(
    bin_start = bin * bin_size,
    bin_end   = pmin((bin + 1) * bin_size, length)
  )]

  as.data.frame(density_data)
}

#' Read a BED overlay file (e.g. reference gaps / centromeres)
#'
#' @param bed_file Path to BED file (3-column: contig, start, end)
#' @param chrom_levels Character vector of chromosome levels to keep
#' @return data.frame with chromosome, start, end columns, or NULL
read_bed_overlay <- function(bed_file, chrom_levels) {
  if (is.null(bed_file) || !nzchar(bed_file)) return(NULL)
  if (!file.exists(bed_file)) {
    stop("BED overlay file not found: ", bed_file)
  }

  cat("Reading BED file from:", bed_file, "\n")
  bed_data <- fread(bed_file, header = FALSE, sep = "\t",
                    col.names = c("contig", "start", "end"),
                    colClasses = c("character", "integer", "integer"))

  # Extract chromosome from contig (handles PREFIX#HAPLOTYPE#chr format)
  bed_data$chromosome <- sub(".*#", "", bed_data$contig)

  # Filter to our chromosome set
  bed_data <- bed_data[bed_data$chromosome %in% chrom_levels, ]
  bed_data$chromosome <- factor(bed_data$chromosome, levels = chrom_levels)
  bed_data <- bed_data[!is.na(bed_data$chromosome), ]

  cat("Read", nrow(bed_data), "BED intervals\n")
  as.data.frame(bed_data)
}

#' Build a chromosome density ideogram plot
#'
#' @param density_data data.frame from compute_density_bins()
#' @param chrom_lengths data.frame with chromosome, length columns
#' @param bed_data data.frame from read_bed_overlay() or NULL
#' @param plot_title Plot title string
#' @param scale_type Fill scale transform: "log1p", "sqrt", "log", "identity"
#' @return ggplot object
plot_ideogram <- function(density_data, chrom_lengths, bed_data, plot_title, scale_type = "log1p") {
  chrom_levels <- levels(chrom_lengths$chromosome)

  # Layout parameters
  ideogram_height  <- 0.6
  ideogram_spacing <- 1.0
  y_positions <- seq(length(chrom_levels), 1, by = -1) * ideogram_spacing

  chrom_y_map <- data.frame(
    chromosome = factor(chrom_levels, levels = chrom_levels),
    y_pos = y_positions
  )

  # Merge y positions
  density_data  <- merge(density_data, chrom_y_map, by = "chromosome")
  chrom_lengths <- merge(chrom_lengths, chrom_y_map, by = "chromosome")
  if (!is.null(bed_data) && nrow(bed_data) > 0) {
    bed_data <- merge(bed_data, chrom_y_map, by = "chromosome")
  }

  p <- ggplot() +
    # Chromosome backgrounds
    geom_rect(data = chrom_lengths,
              aes(xmin = 0, xmax = length,
                  ymin = y_pos - ideogram_height / 2,
                  ymax = y_pos + ideogram_height / 2),
              fill = "white", color = "black", linewidth = 0.5) +
    # Density tiles
    geom_rect(data = density_data,
              aes(xmin = bin_start, xmax = bin_end,
                  ymin = y_pos - ideogram_height / 2,
                  ymax = y_pos + ideogram_height / 2,
                  fill = count)) +
    # BED overlay
    { if (!is.null(bed_data) && nrow(bed_data) > 0)
        geom_rect(data = bed_data,
                  aes(xmin = start, xmax = end,
                      ymin = y_pos - ideogram_height / 2,
                      ymax = y_pos + ideogram_height / 2,
                      alpha = "Reference gaps"),
                  fill = "black")
    } +
    scale_fill_gradientn(
      colors = c("#FFFACD", "yellow", "orange", "#FF4500", "#FF0000"),
      trans  = scale_type,
      name   = paste0("Variant count\n(", scale_type, " scale)"),
      labels = comma,
      na.value = "white",
      breaks = c(0, 10, 100, 1000, 10000, 100000),
      guide  = guide_colorbar(barheight = unit(4, "in"),
                               barwidth  = unit(0.5, "cm"),
                               label.position = "right")
    ) +
    { if (!is.null(bed_data) && nrow(bed_data) > 0)
        scale_alpha_manual(values = c("Reference gaps" = 0.7),
                           name = NULL,
                           guide = guide_legend(order = 2))
    } +
    scale_y_continuous(breaks = y_positions,
                       labels = chrom_levels,
                       expand = c(0.02, 0)) +
    scale_x_continuous(labels = unit_format(unit = "Mb", scale = 1e-6),
                       expand = c(0.01, 0)) +
    labs(title = plot_title, x = "Position", y = "Chromosome") +
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
      plot.background  = element_rect(fill = "white", color = NA)
    )

  p
}

#' Save plot and print summary statistics
#'
#' @param p ggplot object
#' @param output_file Path to write PNG
#' @param data data.frame with chromosome and ref_length columns for summary
save_and_summarize <- function(p, output_file, data) {
  cat("Saving plot to:", output_file, "\n")
  ggsave(output_file, p, width = 12, height = 10, dpi = 300, bg = "white",
         device = grDevices::png, type = "cairo")
  cat("Done!\n")

  cat("\nSummary by chromosome:\n")
  summary_stats <- data %>%
    group_by(chromosome) %>%
    summarise(
      n_variants = n(),
      total_bp   = sum(ref_length),
      .groups    = "drop"
    ) %>%
    arrange(chromosome)

  print(summary_stats)
}
