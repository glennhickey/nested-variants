#!/usr/bin/env Rscript

# vcfeval-snp-slide.R — slide version of the SNP row of vcfeval-detailed.
# Reads the vcfeval-detailed TSV (sample, ref_context, variant_type,
# category, giab, count) and emits a two-panel figure:
# Off-reference SNP | On-reference SNP, each panel three stacked bars
# (Shared, Call only, FreeBayes only) coloured by GIAB region.
#
# Bars are mean per-sample counts; each panel has its own y-scale because
# on-ref counts are ~100× off-ref.
#
# Usage:
#   Rscript scripts/vcfeval-snp-slide.R <vcfeval-detailed.tsv> <out-prefix> [--pdf]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

.script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
if (is.na(.script_dir) || !nzchar(.script_dir)) .script_dir <- "scripts"
source(file.path(.script_dir, "plot-helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript vcfeval-snp-slide.R <vcfeval-detailed.tsv> <out-prefix> [--pdf]\n")
  quit(status = 1)
}
tsv_in   <- args[1]
prefix   <- args[2]
cli_pdf  <- "--pdf" %in% args
emit_pdf <- pdf_enabled(cli_flag = cli_pdf)

dt <- fread(tsv_in)
dt <- dt[variant_type == "SNP"]
n_samples <- uniqueN(dt$sample)

# Per-sample panel is defined by (ref_context, category, giab). Sum across
# samples then divide → mean per sample.
agg <- dt[, .(count = sum(count) / n_samples),
          by = .(ref_context, category, giab)]

agg[, ref_context := factor(ref_context,
                            levels = c("On-reference", "Off-reference"))]
agg[, category    := factor(category,
                            levels = c("Shared", "Call only", "FreeBayes only"))]
agg[, giab        := factor(giab, levels = c("Easy", "Segdup", "Hard"))]

slide_theme <- theme_minimal(base_size = 16) +
  theme(
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title         = element_text(face = "bold"),
    axis.text.x        = element_text(angle = 20, hjust = 1),
    strip.text         = element_text(face = "bold", size = 16),
    strip.background   = element_rect(fill = "white", color = NA),
    legend.position    = "right",
    legend.title       = element_text(face = "bold")
  )

giab_cols <- c(Easy   = "#4daf4a",
               Segdup = "#a02b2b",
               Hard   = "#ff7f00")

# Per-(ref_context, category) totals, for text labels above every bar so the
# tiny bars (e.g. on-ref FreeBayes only) still show their value.
bar_totals <- agg[, .(count = sum(count)), by = .(ref_context, category)]

p <- ggplot(agg, aes(x = category, y = count, fill = giab)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.4) +
  geom_text(data = bar_totals,
            aes(x = category, y = count,
                label = scales::label_number(accuracy = 0.1,
                                             scale_cut = scales::cut_short_scale())(count)),
            vjust = -0.4, size = 4.4, fontface = "bold",
            inherit.aes = FALSE) +
  scale_fill_manual(values = giab_cols, name = "GIAB region") +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~ ref_context, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "SNPs per Sample (mean)") +
  slide_theme

cat("Samples:", n_samples, "\n")
print(dcast(agg, ref_context + category ~ giab, value.var = "count"))

save_plot(p, paste0(prefix, ".vcfeval-snp-slide.png"),
          width = 12, height = 5.5, pdf = emit_pdf)
