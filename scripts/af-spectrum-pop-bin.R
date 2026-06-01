#!/usr/bin/env Rscript

# af-spectrum-pop-bin.R — re-bin the population-AF spectrum from a 2-column
# TSV (CHROM, AF) produced by af-spectrum-pop.sh, and emit
# {prefix}.af-spectrum-pop.{tsv,png,pdf}.
#
# ref_context is derived from CHROM: any contig ending in _<N>_alt is
# off-reference, otherwise on-reference. Matches vcf-stats.R's convention.
#
# Usage:
#   Rscript scripts/af-spectrum-pop-bin.R <ac.tsv> <out_prefix> [--af-step S] [--pdf]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

.script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
if (is.na(.script_dir) || !nzchar(.script_dir)) .script_dir <- "scripts"
source(file.path(.script_dir, "plot-helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript af-spectrum-pop-bin.R <ac.tsv> <out_prefix> [--af-step S] [--pdf]\n")
  quit(status = 1)
}
in_tsv     <- args[1]
out_prefix <- args[2]
af_step    <- 0.01
cli_pdf    <- FALSE
i <- 3
while (i <= length(args)) {
  if (args[i] == "--af-step" && i + 1 <= length(args)) {
    af_step <- as.numeric(args[i + 1]); i <- i + 2
  } else if (args[i] == "--pdf") {
    cli_pdf <- TRUE; i <- i + 1
  } else {
    i <- i + 1
  }
}
emit_pdf <- pdf_enabled(cli_flag = cli_pdf)

cat("Reading:", in_tsv, "\n")
dt <- fread(in_tsv, header = FALSE, col.names = c("CHROM", "AF_str"),
            colClasses = c(CHROM = "character", AF_str = "character"))
cat("Records read:", nrow(dt), "\n")

# AF may be ".", a single value, or comma-separated (multi-allelic). We sum
# multi-allelic AFs as the non-ref total (clamped at 1).
is_multi <- grepl(",", dt$AF_str, fixed = TRUE)
is_dot   <- dt$AF_str == "."
dt[, nonref_af := NA_real_]
dt[!is_multi & !is_dot, nonref_af := pmin(as.numeric(AF_str), 1.0)]
if (any(is_multi)) {
  dt[is_multi, nonref_af := vapply(AF_str, function(x) {
    vals <- suppressWarnings(as.numeric(unlist(strsplit(x, ","))))
    min(sum(vals, na.rm = TRUE), 1.0)
  }, numeric(1))]
}
dt <- dt[!is.na(nonref_af)]
cat("Records with AF:", nrow(dt), "\n")

dt[, ref_context := fifelse(grepl("_[0-9]+_alt$", CHROM),
                            "Off-reference", "On-reference")]
dt[, af_plot := round(nonref_af / af_step) * af_step]
dt[, af_plot := pmin(pmax(af_plot, 0), 1)]
af_counts <- dt[, .(count = .N), by = .(af_plot, ref_context)]
setorder(af_counts, ref_context, af_plot)

tsv_out <- paste0(out_prefix, ".af-spectrum-pop.tsv")
fwrite(af_counts, tsv_out, sep = "\t")
cat("Wrote:", tsv_out, "\n")

p <- ggplot(af_counts, aes(x = af_plot, y = count, color = ref_context)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.0) +
  scale_color_manual(values = c("On-reference" = "steelblue",
                                "Off-reference" = "coral"),
                     name = NULL) +
  scale_y_log10(labels = scales::comma) +
  scale_x_continuous(limits = c(-0.02, 1.02)) +
  labs(x = "Non-Reference Frequency", y = "Sites (log scale)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.position  = "right"
  )

save_plot(p, paste0(out_prefix, ".af-spectrum-pop.png"),
          width = 8, height = 6, pdf = emit_pdf)
cat("Done.\n")
