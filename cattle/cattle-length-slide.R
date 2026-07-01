#!/usr/bin/env Rscript

# cattle-length-slide.R — slide-ready GRef segment length histogram for the
# cattle pangenome.  Reads the gref-segs TSV (one row per gref segment) and
# emits a single-panel length histogram matching the styling of the HPRC
# augref-depth-slide plots.
#
# Input columns (TSV, no header):
#   source_path, source_start, source_end, gref_path, ref_path, ref_start, ref_end
# Segment length = source_end - source_start (the source-haplotype view of the
# non-reference insertion).
#
# Usage:
#   Rscript scripts/cattle/cattle-length-slide.R <gref-segs.tsv[.gz]> <out-prefix> [--pdf]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Locate plot-helpers.R in the top-level scripts/ dir (repo-relative).
.script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
if (is.na(.script_dir) || !nzchar(.script_dir)) .script_dir <- "cattle"
source(file.path(.script_dir, "..", "scripts", "plot-helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript cattle-length-slide.R <gref-segs.tsv[.gz]> <out-prefix>",
      "[--onref-bp N] [--pdf]\n")
  quit(status = 1)
}
tsv_in   <- args[1]
prefix   <- args[2]
cli_pdf  <- "--pdf" %in% args
emit_pdf <- pdf_enabled(cli_flag = cli_pdf)

# --onref-bp: total on-reference length in bp, for the totals inset.
onref_bp <- NA_real_
onref_idx <- which(args == "--onref-bp")
if (length(onref_idx) == 1 && onref_idx + 1 <= length(args)) {
  onref_bp <- as.numeric(args[onref_idx + 1])
}

read_cmd <- if (grepl("\\.gz$", tsv_in)) paste("zcat", shQuote(tsv_in)) else tsv_in
dt <- fread(cmd = read_cmd, header = FALSE,
            col.names = c("source_path", "source_start", "source_end",
                          "gref_path", "ref_path", "ref_start", "ref_end"))
dt[, length := source_end - source_start]
dt <- dt[length > 0]
cat("GRef segments:", nrow(dt), "\n")
cat("Length: min=", min(dt$length), " median=", median(dt$length),
    " max=", max(dt$length), "\n", sep = "")

slide_theme <- theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    axis.title       = element_text(face = "bold"),
    legend.position  = "none"
  )

offref_bp <- sum(as.numeric(dt$length))

fmt_bp <- function(bp) {
  if (is.na(bp)) return("—")
  if (bp >= 1e9) sprintf("%.2f Gb", bp / 1e9)
  else if (bp >= 1e6) sprintf("%.1f Mb", bp / 1e6)
  else if (bp >= 1e3) sprintf("%.0f kb", bp / 1e3)
  else sprintf("%d bp", as.integer(bp))
}

inset_lines <- c("Total length",
                 paste0("On-ref:  ", fmt_bp(onref_bp)),
                 paste0("Off-ref: ", fmt_bp(offref_bp)))
inset_label <- paste(inset_lines, collapse = "\n")
cat(inset_label, "\n")

p <- ggplot(dt, aes(x = length)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", linewidth = 0.2) +
  annotate("label",
           x = Inf, y = Inf, hjust = 1.05, vjust = 1.15,
           label = inset_label,
           family = "mono", size = 4.5, fontface = "plain",
           label.size = 0.3, fill = "white",
           label.padding = unit(0.5, "lines")) +
  scale_x_log10(labels = scales::label_comma(),
                breaks = c(100, 1e3, 1e4, 1e5, 1e6)) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "Length (bp)", y = "GRef Segments") +
  slide_theme

save_plot(p, paste0(prefix, ".cattle-length-slide.png"),
          width = 7, height = 5.5, pdf = emit_pdf)
