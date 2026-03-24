#!/usr/bin/env Rscript

# contig-depth-plot.R — Scatter plot of read depth vs contig length for augref contigs
#
# Usage: Rscript scripts/contig-depth-plot.R
#          --depth <contig-depth.tsv>   (from vg depth -b 1000000000)
#          --segs <augref-segs.tsv>
#          --sample <sample_name>
#          --output <output.png>

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
depth_path  <- NULL
segs_path   <- NULL
sample_name <- NULL
output      <- NULL

i <- 1
while (i <= length(args)) {
  if (args[i] == "--depth" && i + 1 <= length(args)) {
    depth_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--segs" && i + 1 <= length(args)) {
    segs_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--sample" && i + 1 <= length(args)) {
    sample_name <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--output" && i + 1 <= length(args)) {
    output <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(depth_path) || is.null(segs_path) || is.null(output)) {
  cat("Usage: Rscript contig-depth-plot.R --depth <tsv> --segs <tsv> --sample <name> --output <png>\n")
  quit(status = 1)
}
if (is.null(sample_name)) sample_name <- "Unknown"

save_png <- function(plot, file, width = 10, height = 7) {
  tryCatch({
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(file, width = width, height = height, units = "in", res = 300)
      print(plot)
      dev.off()
    } else {
      ggsave(file, plot = plot, width = width, height = height, dpi = 300,
             device = grDevices::png, type = "cairo")
    }
  }, error = function(e) {
    grDevices::png(file, width = width * 300, height = height * 300, res = 300, type = "cairo")
    print(plot)
    dev.off()
  })
  cat("Saved:", file, "\n")
}

# ---------------------------------------------------------------------------
# Read data
# ---------------------------------------------------------------------------
depth <- fread(depth_path, header = FALSE,
               col.names = c("path", "bin_start", "bin_end", "mean_depth", "stddev"))
cat("Depth records:", nrow(depth), "\n")

segs <- fread(segs_path, select = c(1, 2, 3, 4),
              col.names = c("source", "src_start", "src_end", "augref_path"))
segs[, length := src_end - src_start]
segs[, source_sample := sub("#.*", "", source)]
segs <- unique(segs, by = "augref_path")
cat("Augref contigs:", nrow(segs), "\n")

# Get backbone (on-ref) depth for reference line
backbone_depth <- depth[grepl("^augref_", path) & !grepl("_[0-9]+_alt$", path)]
ref_depth <- if (nrow(backbone_depth) > 0) {
  # Weighted mean by bin size
  backbone_depth[, weighted.mean(mean_depth, bin_end - bin_start)]
} else NA_real_

cat("Backbone reference depth:", round(ref_depth, 1), "x\n")

# Filter depth to off-ref contigs only
depth_offref <- depth[grepl("_[0-9]+_alt$", path)]

# Left join: all off-ref contigs get depth (0 if no reads mapped)
merged <- merge(segs, depth_offref, by.x = "augref_path", by.y = "path", all.x = TRUE)
merged[is.na(mean_depth), c("mean_depth", "stddev") := 0]

n_covered <- sum(merged$mean_depth > 0)
n_total   <- nrow(merged)
cat("Contigs with depth > 0:", n_covered, "of", n_total, "\n")

if (n_total == 0) {
  cat("No augref contigs found; creating empty plot.\n")
  file.create(output)
  quit(status = 0)
}

# ---------------------------------------------------------------------------
# Plot: length vs depth scatter
# ---------------------------------------------------------------------------
merged[, has_depth := mean_depth > 0]

p <- ggplot(merged, aes(x = length, y = mean_depth, color = source_sample)) +
  geom_point(aes(shape = has_depth), alpha = 0.7, size = 2) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 4),
                     labels = c("TRUE" = "Covered", "FALSE" = "No reads"),
                     name = NULL) +
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_discrete(name = "Source Haplotype") +
  labs(title = paste("Augref Contig Read Depth:", sample_name),
       subtitle = paste0(n_covered, " of ", n_total, " off-ref contigs covered",
                         if (!is.na(ref_depth)) paste0("  |  backbone depth: ", round(ref_depth, 1), "x") else ""),
       x = "Contig Length (bp)",
       y = "Mean Read Depth") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

# Add backbone depth reference line if available
if (!is.na(ref_depth) && ref_depth > 0) {
  p <- p + geom_hline(yintercept = ref_depth, linetype = "dashed", color = "grey40")
}

save_png(p, output)
cat("Done.\n")
