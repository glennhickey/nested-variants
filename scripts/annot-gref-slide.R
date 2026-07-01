#!/usr/bin/env Rscript

# annot-gref-slide.R — one-slide summary of GRef annotation content and
# neighborhood breakdown. Combines two views into a single figure:
#   Left  : bar chart of "what fraction of GRef bp overlaps each annotation"
#           (Mb totals labelled). Off-reference only — the on-ref/genome-wide
#           comparator that the standalone annot-summary plot shows is dropped.
#   Right : 2x2 small-multiples of "for each GRef X category, what fraction of
#           those segments' on-ref attachment points overlap each on-ref
#           annotation" (the annot-cooccur-slide view).
#
# Inputs:
#   <prefix>.annot-stats.tsv        (annotation, coord_type, total_bp,
#                                    overlap_bp, overlap_frac, ...)
#   <prefix>.annot-per-segment.tsv  (augref_path, annotation, annotation_class,
#                                    source_overlap_bp, source_overlap_frac,
#                                    ref_overlap_bp, ref_overlap_frac, ...)
#
# Usage:
#   Rscript scripts/annot-gref-slide.R \
#       <annot-stats.tsv> <annot-per-segment.tsv> <out-prefix> [--pdf]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

.script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
if (is.na(.script_dir) || !nzchar(.script_dir)) .script_dir <- "scripts"
source(file.path(.script_dir, "plot-helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript annot-gref-slide.R",
      "<annot-stats.tsv> <annot-per-segment.tsv> <out-prefix> [--pdf]\n")
  quit(status = 1)
}
stats_in <- args[1]
segs_in  <- args[2]
prefix   <- args[3]
cli_pdf  <- "--pdf" %in% args
emit_pdf <- pdf_enabled(cli_flag = cli_pdf)

KEEP <- c("censat", "genes", "repeats", "segdups")
on_ref_cols <- c(censat  = "#a02b2b",
                 genes   = "#4daf4a",
                 repeats = "#ff7f00",
                 segdups = "#377eb8")

slide_theme <- theme_minimal(base_size = 15) +
  theme(
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.title         = element_text(face = "bold"),
    strip.text         = element_text(face = "bold", size = 14),
    strip.background   = element_rect(fill = "white", color = NA),
    legend.position    = "none"
  )

# ---------------------------------------------------------------------------
# Left panel — GRef overall totals from annot-stats.tsv
# ---------------------------------------------------------------------------

stats <- fread(stats_in)
summary_dt <- stats[coord_type == "off-reference" & annotation %in% KEEP,
                    .(annotation, frac = overlap_frac,
                      overlap_bp = as.numeric(overlap_bp),
                      segments_with_overlap)]

fmt_bp <- function(x) {
  ifelse(x >= 1e9, sprintf("%.1f Gb", x / 1e9),
  ifelse(x >= 1e6, sprintf("%.0f Mb", x / 1e6),
  ifelse(x >= 1e3, sprintf("%.0f kb", x / 1e3),
                    sprintf("%.0f bp", x))))
}
summary_dt[, bp_label := fmt_bp(overlap_bp)]
summary_dt[, annotation := factor(annotation, levels = KEEP)]

# Facet the summary the same way as the breakdown (one facet per annotation,
# stacked vertically) so gridExtra aligns rows across the two grobs. Title
# shows the capitalized annotation name plus the # of segments contributing.
cap <- function(s) paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
summary_dt[, facet_title := paste0(sapply(as.character(annotation), cap),
                                   "  (",
                                   format(segments_with_overlap, big.mark = ","),
                                   " segs)")]
facet_levels <- summary_dt[order(match(annotation, KEEP))]$facet_title
summary_dt[, facet_title := factor(facet_title, levels = facet_levels)]

p_summary <- ggplot(summary_dt,
                    aes(x = frac, y = annotation, fill = annotation)) +
  # Faint per-annotation panel background so the four rows read as visually
  # linked bands with their right-side breakdown counterparts.
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                fill = annotation),
            alpha = 0.12, color = NA, inherit.aes = FALSE,
            data = summary_dt) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.0f%%  (%s)", 100 * frac, bp_label)),
            hjust = -0.05, size = 4.4, fontface = "bold") +
  scale_fill_manual(values = on_ref_cols) +
  scale_x_continuous(labels = scales::label_percent(),
                     limits = c(0, 1.35), expand = c(0, 0)) +
  scale_y_discrete(labels = NULL, breaks = NULL) +
  facet_wrap(~ facet_title, ncol = 1, scales = "free_y") +
  labs(x = "Fraction of GRef bp overlapping annotation",
       y = NULL) +
  slide_theme +
  theme(plot.title = element_text(face = "bold", size = 16,
                                  hjust = 0.5, margin = margin(b = 10)))
p_summary <- p_summary + ggtitle("Off-Ref")

# ---------------------------------------------------------------------------
# Right panel — per-GRef-category on-ref neighborhood breakdown
# ---------------------------------------------------------------------------

cat("Reading:", segs_in, "\n")
segs <- fread(segs_in,
              select = c("augref_path", "annotation",
                         "source_overlap_frac", "ref_overlap_frac"))
segs <- segs[annotation %in% KEEP]

seg_ann <- segs[, .(source_overlap_frac = sum(source_overlap_frac),
                    ref_overlap_frac    = sum(ref_overlap_frac)),
                by = .(augref_path, annotation)]

off_segs <- lapply(KEEP, function(a)
  seg_ann[annotation == a & source_overlap_frac > 0]$augref_path)
names(off_segs) <- KEEP
on_segs  <- lapply(KEEP, function(a)
  seg_ann[annotation == a & ref_overlap_frac    > 0]$augref_path)
names(on_segs) <- KEEP

rows <- list()
for (x in KEEP) {
  denom <- length(off_segs[[x]])
  for (y in KEEP) {
    n_both <- length(intersect(off_segs[[x]], on_segs[[y]]))
    rows[[length(rows) + 1]] <- data.table(
      off_ref = x, on_ref = y,
      n_off = denom, frac = if (denom > 0) n_both / denom else 0
    )
  }
}
tab <- rbindlist(rows)

denom_by <- tab[, .(n_off = n_off[1]), by = off_ref]
# Titles here must match the summary facet titles (`paste0("GRef ", ann)`) so
# gridExtra stacks the rows in visual register. Sample counts appear in the
# body of each facet instead.
denom_by[, facet_title := paste0("GRef ", off_ref)]
tab <- merge(tab, denom_by[, .(off_ref, facet_title)], by = "off_ref")
tab[, off_ref := factor(off_ref, levels = KEEP)]
tab[, on_ref  := factor(on_ref,  levels = KEEP)]
tab[, facet_title := factor(facet_title,
                            levels = denom_by[order(match(off_ref, KEEP))]$facet_title)]

# Per-facet n-count text, placed in-panel so removing it from strip titles
# doesn't lose the info.
seg_counts <- unique(tab[, .(facet_title, on_ref = KEEP[1], n_off)])
seg_counts[, seg_lab := paste0("(", format(n_off, big.mark = ","), " segs)")]

bg_right <- unique(tab[, .(facet_title, off_ref)])

p_breakdown <- ggplot(tab, aes(x = frac, y = on_ref, fill = on_ref)) +
  # Match the left-plot per-row panel background so each row visually links
  # left ↔ right. Fill maps through the same scale as the bars because off_ref
  # and on_ref share the same category set.
  geom_rect(data = bg_right,
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                fill = off_ref),
            alpha = 0.12, color = NA, inherit.aes = FALSE) +
  geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * frac)),
            hjust = -0.15, size = 3.8, fontface = "bold") +
  scale_x_continuous(labels = scales::label_percent(),
                     limits = c(0, 1.15), expand = c(0, 0)) +
  # Pad top/bottom of each facet so the 4 right bars collectively span the
  # same vertical extent as the single left bar (which has width 0.65 with
  # ggplot's default 0.6 y-expand → ~54% of the row). Compensate with add ≈
  # (n_bars - 1) / (2 * left_frac) - default_expand.
  scale_y_discrete(limits = rev(KEEP),
                   expand = expansion(add = c(1.8, 1.8))) +
  scale_fill_manual(values = on_ref_cols) +
  facet_wrap(~ facet_title, ncol = 1) +
  labs(x = "Fraction overlapping on-ref annotation", y = NULL) +
  slide_theme +
  theme(plot.title  = element_text(face = "bold", size = 16,
                                   hjust = 0.5, margin = margin(b = 10)),
        strip.text  = element_blank(),
        panel.spacing.y = unit(0.4, "lines"))
p_breakdown <- p_breakdown + ggtitle("On-Ref Neighbourhood")

# ---------------------------------------------------------------------------
# Combine (left 1 col : right 2 cols).
# ---------------------------------------------------------------------------

combined <- gridExtra::arrangeGrob(p_summary, p_breakdown,
                                   ncol = 2, widths = c(1, 1.2))

save_plot(combined, paste0(prefix, ".annot-gref-slide.png"),
          width = 13, height = 9, pdf = emit_pdf)
