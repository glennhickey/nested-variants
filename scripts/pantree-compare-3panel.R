#!/usr/bin/env Rscript

# pantree-compare-3panel.R — three-panel GRef-vs-Pantree variant catalog plot
#
# Reads the two `*.sites.vcf-stats.tsv` summary tables that vcf-stats.R already
# emits (one per catalog), collapses to three categories — SNP / Indel+MNP / SV
# — and renders one ggplot facet per category so each gets its own y-axis. Bars
# are stacked by ref_context; Ts/Tv ratio is annotated on top of each SNP bar.
#
# Usage:
#   Rscript scripts/pantree-compare-3panel.R
#     --ours    <ours.sites.vcf-stats.tsv>
#     --pantree <pantree.vcf-stats.tsv>
#     --prefix  <out_prefix>
#     [--ours-label    NAME]   # default "GRef"
#     [--pantree-label NAME]   # default "Pantree"
#     [--title TITLE]
#
# Outputs:
#   {prefix}.pantree-compare-3panel.png
#   {prefix}.pantree-compare-3panel.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
ours_path     <- NULL
pantree_path  <- NULL
prefix        <- NULL
title         <- "GRef vs Pantree"
ours_label    <- "GRef"
pantree_label <- "Pantree"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--ours" && i + 1 <= length(args)) {
    ours_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--pantree" && i + 1 <= length(args)) {
    pantree_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--prefix" && i + 1 <= length(args)) {
    prefix <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--ours-label" && i + 1 <= length(args)) {
    ours_label <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--pantree-label" && i + 1 <= length(args)) {
    pantree_label <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}
if (is.null(ours_path) || is.null(pantree_path) || is.null(prefix)) {
  cat("Usage: Rscript pantree-compare-3panel.R --ours <tsv> --pantree <tsv> --prefix <path>\n")
  quit(status = 1)
}

# Map fine variant_type from vcf-stats.tsv to the three coarse categories.
to_category <- function(vt) {
  fcase(
    vt == "SNP",                                                "SNP",
    vt %in% c("Insertion", "Deletion", "MNP"),                  "Indel/MNP",
    vt %in% c("SV Insertion", "SV Deletion"),                   "SV",
    default = NA_character_
  )
}

read_stats <- function(path, source_label) {
  dt <- fread(path)
  dt[, source := source_label]
  dt[, category := to_category(variant_type)]
  dt <- dt[!is.na(category)]
  # ensure ts/tv numeric (file uses blanks for non-SNP rows)
  for (col in c("ts", "tv")) {
    if (!col %in% names(dt)) dt[, (col) := NA_real_]
    dt[, (col) := suppressWarnings(as.numeric(get(col)))]
  }
  dt
}

dt_o <- read_stats(ours_path,    ours_label)
dt_p <- read_stats(pantree_path, pantree_label)
all_dt <- rbindlist(list(dt_o, dt_p), use.names = TRUE, fill = TRUE)

# Aggregate counts and Ts/Tv per (source, category, ref_context).
agg <- all_dt[, .(count = sum(count, na.rm = TRUE),
                  ts    = sum(ts,    na.rm = TRUE),
                  tv    = sum(tv,    na.rm = TRUE)),
              by = .(source, category, ref_context)]

# Bar totals (over ref_context) — used for the legibility threshold on
# per-segment Ts/Tv labels.
bar_totals <- agg[, .(bar_total = sum(count)),
                  by = .(source, category)]

# Per-segment Ts/Tv (one label per (source × ref_context) inside each SNP bar)
# placed at the segment midpoint via position_stack(vjust = 0.5).  Suppress
# labels whose segment is < 12 % of the bar total — too small to read.
seg_tstv <- agg[, .(tstv = fifelse(tv > 0, ts / tv, NA_real_)),
                by = .(source, category, ref_context, count)]
seg_tstv <- merge(seg_tstv, bar_totals, by = c("source", "category"))
seg_tstv[, keep := !is.na(tstv) & bar_total > 0 &
                   count >= 0.12 * bar_total]
seg_tstv[, label := sprintf("Ts/Tv = %.2f", tstv)]

# Factor levels — fixed source order so colours stay consistent across runs.
source_levels <- c(ours_label, pantree_label)
agg[, source := factor(source, levels = source_levels)]
bar_totals[, source := factor(source, levels = source_levels)]

cat_levels <- c("SNP", "Indel/MNP", "SV")
agg[, category := factor(category, levels = cat_levels)]
bar_totals[, category := factor(category, levels = cat_levels)]
seg_tstv[, source := factor(source, levels = source_levels)]
seg_tstv[, category := factor(category, levels = cat_levels)]

ref_levels <- c("Off-reference", "On-reference")
agg[, ref_context := factor(ref_context, levels = ref_levels)]
seg_tstv[, ref_context := factor(ref_context, levels = ref_levels)]

# Persist the aggregated table for reproducibility.
fwrite(agg[order(category, source, ref_context),
           .(source, category, ref_context, count, ts, tv)],
       paste0(prefix, ".pantree-compare-3panel.tsv"), sep = "\t")
cat("Wrote:", paste0(prefix, ".pantree-compare-3panel.tsv"), "\n")

ref_colors <- c("Off-reference" = "coral", "On-reference" = "steelblue")

# Per-segment Ts/Tv labels: only on the SNP facet, one per (source × ref_context)
# placed at the segment midpoint inside the stacked bar.
snp_seg_labels <- seg_tstv[category == "SNP" & keep == TRUE]

p <- ggplot(agg, aes(x = source, y = count, fill = ref_context)) +
  geom_col(width = 0.65, alpha = 0.95) +
  geom_text(data = snp_seg_labels,
            aes(x = source, y = count, label = label, group = ref_context),
            position = position_stack(vjust = 0.5),
            size = 3.2, color = "white") +
  facet_wrap(vars(category), scales = "free_y", nrow = 1) +
  scale_fill_manual(values = ref_colors, name = NULL, drop = FALSE) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0.02, 0.12))) +
  labs(title = title,
       subtitle = paste0(ours_label, " vs ", pantree_label,
                         " — variant counts by category (stacks = ref-context",
                         "; Ts/Tv per ref-context inside SNP bars)"),
       x = NULL, y = "Count") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

save_png <- function(plot, path, w = 11, h = 5) {
  tryCatch({
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(path, width = w, height = h, units = "in", res = 300)
      print(plot); dev.off()
    } else {
      ggsave(path, plot = plot, width = w, height = h, dpi = 300,
             device = grDevices::png, type = "cairo")
    }
  }, error = function(e) {
    grDevices::png(path, width = w * 300, height = h * 300, res = 300, type = "cairo")
    print(plot); dev.off()
  })
}
save_png(p, paste0(prefix, ".pantree-compare-3panel.png"))
cat("Saved:", paste0(prefix, ".pantree-compare-3panel.png"), "\n")
cat("Done.\n")
