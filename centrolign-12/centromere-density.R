#!/usr/bin/env Rscript
# centromere-density.R — Density plot of off-reference segments along the
# centrolign reference path (zoomed view, not a full-chromosome ideogram).
#
# Usage: Rscript centromere-density.R <augref-segs.tsv> <output.png> [title] [min_length] [bin_kb]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript centromere-density.R <input.tsv> <output.png> [title] [min_length] [bin_kb]\n")
  quit(status = 1)
}

input_file <- args[1]
output_file <- args[2]
plot_title <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "Off-Reference Segment Density"
min_length <- if (length(args) >= 4 && nzchar(args[4])) as.numeric(args[4]) else 50
bin_kb <- if (length(args) >= 5 && nzchar(args[5])) as.numeric(args[5]) else 10

# Read and compute segment lengths
dt <- fread(input_file, header = FALSE)
setnames(dt, c("source", "src_start", "src_end", "augref_path",
               "ref_contig", "ref_start", "ref_end"))
dt[, seg_len := src_end - src_start]
dt <- dt[seg_len >= min_length]
cat("Segments >= ", min_length, " bp: ", nrow(dt), "\n")

# Reference span
ref_max <- max(dt$ref_end)
ref_min <- min(dt$ref_start)
cat("Reference span: ", ref_min, " - ", ref_max, " (",
    round((ref_max - ref_min) / 1e6, 2), " Mb)\n")

# Bin by ref_start position
bin_size <- bin_kb * 1000
dt[, bin := floor(ref_start / bin_size) * bin_size]

# Count segments and total bp per bin
bin_counts <- dt[, .(n_segments = .N,
                     total_bp = sum(seg_len)),
                 by = bin]

# Fill in empty bins
all_bins <- data.table(bin = seq(floor(ref_min / bin_size) * bin_size,
                                  floor(ref_max / bin_size) * bin_size,
                                  by = bin_size))
bin_counts <- merge(all_bins, bin_counts, by = "bin", all.x = TRUE)
bin_counts[is.na(n_segments), c("n_segments", "total_bp") := 0]

# Position in kb for axis labels
bin_counts[, pos_kb := bin / 1000]

# Two-panel plot: segment count and total off-ref bp (log1p scale)
p1 <- ggplot(bin_counts, aes(x = pos_kb, y = n_segments)) +
  geom_col(fill = "#E74C3C", width = bin_kb * 0.9) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Segment count (log)",
       title = plot_title,
       subtitle = paste0("Off-reference segments in ", bin_kb, " kb bins")) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

p2 <- ggplot(bin_counts, aes(x = pos_kb, y = total_bp / 1000)) +
  geom_col(fill = "#3498DB", width = bin_kb * 0.9) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000, 100000),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Reference position (kb)", y = "Total off-ref bp, kb (log)") +
  theme_bw(base_size = 14)

p <- gridExtra::arrangeGrob(p1, p2, ncol = 1, heights = c(1, 0.8))

ggsave(output_file, p, width = 10, height = 6, dpi = 300)
cat("Saved:", output_file, "\n")
