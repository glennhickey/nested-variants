#!/usr/bin/env Rscript

# augref-depth-slide.R — slide-ready two-panel figure from an augref-depth TSV.
# Left:  GRef segment length histogram (log x)
# Right: GRef segment depth histogram  (linear x, log y)
#
# Usage:
#   Rscript scripts/augref-depth-slide.R <augref-depth.tsv> <out-prefix> [--pdf]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

.script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
if (is.na(.script_dir) || !nzchar(.script_dir)) .script_dir <- "scripts"
source(file.path(.script_dir, "plot-helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript augref-depth-slide.R <augref-depth.tsv> <out-prefix> [--pdf]\n")
  quit(status = 1)
}
tsv_in   <- args[1]
prefix   <- args[2]
cli_pdf  <- "--pdf" %in% args
emit_pdf <- pdf_enabled(cli_flag = cli_pdf)

dt <- fread(tsv_in, header = FALSE,
            col.names = c("path", "bin_start", "bin_end", "mean_depth", "stddev"))
dt <- dt[grepl("_alt$", path)]
dt[, length     := bin_end - bin_start + 1]
dt[, depth_int  := round(mean_depth)]
cat("GRef segments:", nrow(dt), "\n")
cat("Length: min=", min(dt$length), " median=", median(dt$length),
    " max=", max(dt$length), "\n", sep = "")
cat("Depth : min=", min(dt$depth_int), " median=", median(dt$depth_int),
    " max=", max(dt$depth_int), "\n", sep = "")

slide_theme <- theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    axis.title       = element_text(face = "bold"),
    legend.position  = "none"
  )

p_len <- ggplot(dt, aes(x = length)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", linewidth = 0.2) +
  scale_x_log10(labels = scales::label_comma(),
                breaks = c(100, 1e3, 1e4, 1e5)) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "Length (bp)", y = "GRef Segments") +
  slide_theme

step <- 5
hist_dt <- dt[, .(N = .N),
              by = .(bin_lo = ((depth_int - 1) %/% step) * step + 1)]
setorder(hist_dt, bin_lo)
hist_dt[, bin_mid := bin_lo + (step - 1) / 2]

first_N  <- hist_dt[bin_lo == 1, N]
# Cap the y-axis so the first bar visually terminates at `cap`; the tick at
# `cap` is relabelled with the true singleton-bin count, and a broken-axis
# marker just below tells the reader we've compressed the top.
cap       <- 40000
top_label <- format(round(first_N / 1000) * 1000, big.mark = ",")
break_y   <- cap * 0.90
break_h   <- cap * 0.06
bar_xmin  <- 1 - 0.5
bar_xmax  <- step + 0.5 - step * 0.05  # matches geom_col width = step * 0.95

p_dep <- ggplot(hist_dt, aes(x = bin_mid, y = N)) +
  geom_col(fill = "coral", width = step * 0.95) +
  # White band across the first bar to clear space for the break marker.
  annotate("rect",
           xmin = bar_xmin, xmax = bar_xmax,
           ymin = break_y, ymax = break_y + break_h,
           fill = "white", color = NA) +
  # Two short diagonal lines: classic broken-axis symbol.
  annotate("segment",
           x = bar_xmin, xend = bar_xmax,
           y = break_y + break_h * 0.15, yend = break_y + break_h * 0.55,
           linewidth = 0.7, color = "grey30") +
  annotate("segment",
           x = bar_xmin, xend = bar_xmax,
           y = break_y + break_h * 0.45, yend = break_y + break_h * 0.85,
           linewidth = 0.7, color = "grey30") +
  scale_x_continuous(labels = scales::label_comma()) +
  scale_y_continuous(
    breaks = c(0, 10000, 20000, 30000, cap),
    labels = c("0", "10,000", "20,000", "30,000", top_label)) +
  coord_cartesian(ylim = c(0, cap), expand = FALSE) +
  labs(x = "Depth (# haplotypes)", y = NULL) +
  slide_theme

p <- gridExtra::arrangeGrob(p_len, p_dep, ncol = 2)

save_plot(p, paste0(prefix, ".augref-depth-slide.png"),
          width = 12, height = 5, pdf = emit_pdf)
