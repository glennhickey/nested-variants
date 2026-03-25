#!/usr/bin/env Rscript

# mapq-dist-plot.R — MAPQ distribution: on-ref vs off-ref for GAM and BAM
#
# Usage: Rscript scripts/mapq-dist-plot.R
#          --gam-mapq <file1,file2,...>   (per-sample GAM MAPQ TSVs)
#          --bam-mapq <file1,file2,...>   (per-sample BAM MAPQ TSVs)
#          --output <output.png>
#          [--title TITLE]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
gam_paths <- NULL
bam_paths <- NULL
output    <- NULL
title     <- "Mapping Quality Distribution"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--gam-mapq" && i + 1 <= length(args)) {
    gam_paths <- strsplit(args[i + 1], ",")[[1]]; i <- i + 2
  } else if (args[i] == "--bam-mapq" && i + 1 <= length(args)) {
    bam_paths <- strsplit(args[i + 1], ",")[[1]]; i <- i + 2
  } else if (args[i] == "--output" && i + 1 <= length(args)) {
    output <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(output)) {
  cat("Usage: Rscript mapq-dist-plot.R --gam-mapq <f1,f2> --bam-mapq <f1,f2> --output <png>\n")
  quit(status = 1)
}

# ---------------------------------------------------------------------------
# Read and merge MAPQ distributions
# ---------------------------------------------------------------------------
read_mapq <- function(paths, source_label) {
  dt <- rbindlist(lapply(paths, fread))
  # Sum across samples
  dt <- dt[, .(count = sum(count)), by = .(mapq, ref_context)]
  dt[, source := source_label]
  dt
}

panels <- list()
if (!is.null(gam_paths) && length(gam_paths) > 0) {
  panels[["gam"]] <- read_mapq(gam_paths, "Graph (GAM)")
  cat("GAM:", sum(panels[["gam"]]$count), "reads\n")
}
if (!is.null(bam_paths) && length(bam_paths) > 0) {
  panels[["bam"]] <- read_mapq(bam_paths, "Linear (BAM)")
  cat("BAM:", sum(panels[["bam"]]$count), "reads\n")
}

if (length(panels) == 0) {
  cat("No MAPQ data provided.\n")
  file.create(output)
  quit(status = 0)
}

plot_dt <- rbindlist(panels)

# ---------------------------------------------------------------------------
# Plot: faceted histogram
# ---------------------------------------------------------------------------
base_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 11)
  )

ctx_colors <- c("On-reference" = "steelblue", "Off-reference" = "coral")

p <- ggplot(plot_dt, aes(x = mapq, y = count, fill = ref_context)) +
  geom_col(position = "dodge", width = 2) +
  scale_fill_manual(values = ctx_colors, name = NULL) +
  scale_y_log10(labels = scales::comma) +
  facet_wrap(~ source, scales = "free_y") +
  labs(title = title,
       subtitle = "Summed across all samples",
       x = "Mapping Quality (MAPQ)",
       y = "Number of Reads (log)") +
  base_theme

tryCatch({
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(output, width = 12, height = 5, units = "in", res = 300)
  } else {
    grDevices::png(output, width = 12 * 300, height = 5 * 300, res = 300, type = "cairo")
  }
  print(p)
  dev.off()
}, error = function(e) {
  grDevices::png(output, width = 12 * 300, height = 5 * 300, res = 300, type = "cairo")
  print(p)
  dev.off()
})
cat("Saved:", output, "\n")
cat("Done.\n")
