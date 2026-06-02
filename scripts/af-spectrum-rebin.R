#!/usr/bin/env Rscript

# af-spectrum-rebin.R — re-bin the non-reference allele-frequency spectrum
# from a deconstruct records TSV at an arbitrary af_step, and emit
# {prefix}.af-spectrum-fine.{tsv,png,pdf}.
#
# The default `deconstruct_sites_stats` rule bins at af_step = 0.05 (21 points
# across 0..1), which collapses ~24 discrete AF levels into the 0.5 bin on a
# 464-haplotype panel and exaggerates the spike there.  This script reads
# nonref_af straight from records.tsv and re-bins at a finer step so the
# spectrum reflects the underlying distribution instead of the binning.
#
# Usage:
#   Rscript scripts/af-spectrum-rebin.R <records.tsv> <out_prefix> \
#           [--af-step 0.01] [--title TITLE] [--pdf]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Shared title-suppression + save_plot helper (handles --pdf and FIGURE_TITLES)
.script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
if (is.na(.script_dir) || !nzchar(.script_dir)) .script_dir <- "scripts"
source(file.path(.script_dir, "plot-helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript af-spectrum-rebin.R <records.tsv> <out_prefix> [--af-step S] [--title T] [--pdf]\n")
  quit(status = 1)
}
records_tsv <- args[1]
out_prefix  <- args[2]
af_step     <- 0.01
title       <- NULL
cli_pdf     <- FALSE

i <- 3
while (i <= length(args)) {
  if (args[i] == "--af-step" && i + 1 <= length(args)) {
    af_step <- as.numeric(args[i + 1]); i <- i + 2
  } else if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--pdf") {
    cli_pdf <- TRUE; i <- i + 1
  } else {
    i <- i + 1
  }
}
if (is.null(title)) title <- basename(out_prefix)
emit_pdf <- pdf_enabled(cli_flag = cli_pdf)

cat("Reading records:", records_tsv, " (this is a multi-GB stream)...\n")
dt <- fread(records_tsv, select = c("ref_context", "nonref_af"),
            colClasses = c(ref_context = "character", nonref_af = "numeric"))
dt <- dt[!is.na(nonref_af)]
cat("Records with AF:", nrow(dt), "\n")

dt[, af_plot := round(nonref_af / af_step) * af_step]
dt[, af_plot := pmin(pmax(af_plot, 0), 1)]
af_counts <- dt[, .(count = .N), by = .(af_plot, ref_context)]
setorder(af_counts, ref_context, af_plot)

tsv_out <- paste0(out_prefix, ".af-spectrum-fine.tsv")
fwrite(af_counts, tsv_out, sep = "\t")
cat("Wrote:", tsv_out, "\n")

# Match the original af-spectrum panel styling.
p <- ggplot(af_counts, aes(x = af_plot, y = count, color = ref_context)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.0) +
  scale_color_manual(values = c("On-reference" = "steelblue",
                                "Off-reference" = "coral"),
                     name = NULL) +
  scale_y_log10(labels = scales::comma) +
  scale_x_continuous(limits = c(-0.02, 1.02)) +
  labs(x = "Allele Frequency", y = "Sites (log scale)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.position  = "right"
  )

save_plot(p, paste0(out_prefix, ".af-spectrum-fine.png"),
          width = 8, height = 6, pdf = emit_pdf)

cat("Done.\n")
