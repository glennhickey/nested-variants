#!/usr/bin/env Rscript

# annotation-plots.R — Annotation overlap plots for augref alt segments
#
# Usage: Rscript scripts/annotation-plots.R <per_segment.tsv> <output_prefix> [--title TITLE]
#
# Outputs:
#   {prefix}.annot-summary.png  — fraction of alt bp overlapping each annotation (source vs ref)
#   {prefix}.annot-scatter.png  — per-segment source_len vs source_overlap_frac, faceted by annotation
#   {prefix}.annot-repeats.png  — repeat class breakdown (only when repeat class data present)
#   {prefix}.annot-stats.tsv    — tabular summary

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript annotation-plots.R <per_segment.tsv> <output_prefix> [--title TITLE]\n")
  quit(status = 1)
}

input_file <- args[1]
prefix     <- args[2]
title      <- NULL

i <- 3
while (i <= length(args)) {
  if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}
if (is.null(title)) title <- "Annotation Overlap"

# ---------------------------------------------------------------------------
# Read data
# ---------------------------------------------------------------------------
cat("Reading:", input_file, "\n")
dt <- fread(input_file)

if (nrow(dt) == 0) {
  cat("No data found. Exiting.\n")
  quit(status = 0)
}

# ---------------------------------------------------------------------------
# Helper: save PNG with cairo/ragg fallback
# ---------------------------------------------------------------------------
save_png <- function(plot, file, width = 8, height = 6) {
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
# Plot 1: Summary bar chart — fraction of alt bp overlapping each annotation
# ---------------------------------------------------------------------------

# Aggregate per annotation: total overlap bp / total segment length
summary_source <- dt[, .(overlap_bp = as.numeric(sum(source_overlap_bp)),
                         total_bp = as.numeric(sum(source_len))),
                     by = .(annotation)]
summary_source[, frac := overlap_bp / total_bp]
summary_source[, coord_type := "Source"]

summary_ref <- dt[, .(overlap_bp = as.numeric(sum(ref_overlap_bp)),
                      total_bp = as.numeric(sum(ref_len))),
                  by = .(annotation)]
summary_ref[, frac := overlap_bp / total_bp]
summary_ref[, coord_type := "Reference"]

# For grouped annotations (repeats), avoid double-counting segments:
# the total_bp is summed across all classes per segment, so we need unique segments
for (ann in unique(dt$annotation)) {
  sub <- dt[annotation == ann]
  has_classes <- any(sub$annotation_class != sub$annotation)
  if (has_classes) {
    # Unique segments for this annotation
    seg_dt <- unique(sub[, .(augref_path, source_len, ref_len)])
    summary_source[annotation == ann, total_bp := sum(seg_dt$source_len)]
    summary_ref[annotation == ann, total_bp := sum(seg_dt$ref_len)]
    summary_source[annotation == ann, frac := overlap_bp / total_bp]
    summary_ref[annotation == ann, frac := overlap_bp / total_bp]
  }
}

bar_dt <- rbind(summary_source[, .(annotation, frac, coord_type)],
                summary_ref[, .(annotation, frac, coord_type)])

p1 <- ggplot(bar_dt, aes(x = annotation, y = frac, fill = coord_type)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Source" = "coral", "Reference" = "steelblue"),
                    name = NULL) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.1))) +
  labs(title = title, subtitle = "Fraction of Alt Segment bp Overlapping Annotations",
       x = "Annotation", y = "Overlap Fraction") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

save_png(p1, paste0(prefix, ".annot-summary.png"))

# ---------------------------------------------------------------------------
# Plot 2: Length vs overlap scatter, faceted by annotation
# ---------------------------------------------------------------------------

# For grouped annotations, aggregate all classes per segment before plotting
scatter_dt <- dt[, .(source_overlap_bp = sum(source_overlap_bp),
                     source_len = source_len[1],
                     ref_len = ref_len[1]),
                 by = .(augref_path, annotation)]
scatter_dt[, source_overlap_frac := pmin(source_overlap_bp / source_len, 1.0)]
scatter_dt[source_len == 0, source_overlap_frac := 0]

p2 <- ggplot(scatter_dt[source_len > 0],
             aes(x = source_len, y = source_overlap_frac)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "loess", se = FALSE, color = "coral", linewidth = 0.8) +
  facet_wrap(~ annotation, scales = "free_y") +
  scale_x_log10(labels = comma) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(title = title, subtitle = "Source Length vs Overlap Fraction",
       x = "Source Segment Length (bp, log scale)",
       y = "Source Overlap Fraction") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

save_png(p2, paste0(prefix, ".annot-scatter.png"))

# ---------------------------------------------------------------------------
# Plot 3: Repeat class breakdown (only when repeat class data present)
# ---------------------------------------------------------------------------

has_repeat_classes <- any(dt$annotation_class != dt$annotation)

if (has_repeat_classes) {
  repeat_dt <- dt[annotation_class != annotation]

  # Top 8 classes by total source + ref overlap bp
  class_totals <- repeat_dt[, .(total = as.numeric(sum(source_overlap_bp)) + as.numeric(sum(ref_overlap_bp))),
                            by = .(annotation_class)]
  setorder(class_totals, -total)
  top_classes <- head(class_totals$annotation_class, 8)
  repeat_dt[, display_class := ifelse(annotation_class %in% top_classes,
                                      annotation_class, "Other")]

  # Aggregate for stacked bar
  src_bar <- repeat_dt[, .(overlap_bp = as.numeric(sum(source_overlap_bp))),
                       by = .(display_class)]
  src_bar[, coord_type := "Source"]

  ref_bar <- repeat_dt[, .(overlap_bp = as.numeric(sum(ref_overlap_bp))),
                       by = .(display_class)]
  ref_bar[, coord_type := "Reference"]

  stack_dt <- rbind(src_bar, ref_bar)
  # Order classes: top classes first, Other last
  class_order <- c(top_classes[top_classes %in% stack_dt$display_class], "Other")
  stack_dt[, display_class := factor(display_class, levels = rev(class_order))]

  n_classes <- length(unique(stack_dt$display_class))
  if (n_classes <= 8) {
    fill_pal <- scale_fill_brewer(palette = "Set2", name = "Repeat Class")
  } else {
    fill_pal <- scale_fill_brewer(palette = "Set2", name = "Repeat Class")
  }

  p3 <- ggplot(stack_dt, aes(x = coord_type, y = overlap_bp, fill = display_class)) +
    geom_col(width = 0.6) +
    fill_pal +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
    labs(title = title, subtitle = "Repeat Class Breakdown",
         x = NULL, y = "Overlap (bp)") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )

  save_png(p3, paste0(prefix, ".annot-repeats.png"))
} else {
  cat("No repeat class data; skipping repeat breakdown plot.\n")
}

# ---------------------------------------------------------------------------
# Output TSV: per-annotation summary
# ---------------------------------------------------------------------------

# Unique segments per annotation (avoid double-counting from classes)
stats_list <- list()
for (ann in unique(dt$annotation)) {
  sub <- dt[annotation == ann]
  seg_dt <- unique(sub[, .(augref_path, source_len, ref_len)])

  total_source_bp <- sum(seg_dt$source_len)
  total_ref_bp <- sum(seg_dt$ref_len)

  src_overlap <- sum(sub$source_overlap_bp)
  ref_overlap <- sum(sub$ref_overlap_bp)

  seg_with_src <- length(unique(sub[source_overlap_bp > 0]$augref_path))
  seg_with_ref <- length(unique(sub[ref_overlap_bp > 0]$augref_path))
  total_segs <- nrow(seg_dt)

  stats_list[[length(stats_list) + 1]] <- data.table(
    annotation = ann,
    coord_type = "source",
    total_bp = total_source_bp,
    overlap_bp = src_overlap,
    overlap_frac = if (total_source_bp > 0) src_overlap / total_source_bp else 0,
    segments_with_overlap = seg_with_src,
    total_segments = total_segs
  )
  stats_list[[length(stats_list) + 1]] <- data.table(
    annotation = ann,
    coord_type = "reference",
    total_bp = total_ref_bp,
    overlap_bp = ref_overlap,
    overlap_frac = if (total_ref_bp > 0) ref_overlap / total_ref_bp else 0,
    segments_with_overlap = seg_with_ref,
    total_segments = total_segs
  )
}

stats_dt <- rbindlist(stats_list)
stats_file <- paste0(prefix, ".annot-stats.tsv")
fwrite(stats_dt, stats_file, sep = "\t")
cat("Wrote summary:", stats_file, "\n")

cat("Done.\n")
