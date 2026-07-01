#!/usr/bin/env Rscript

# per-sample-types-slide.R — slide-ready side-by-side per-sample variant-count
# violins.  Reads the per-sample-types TSV emitted by vcf-stats.R and produces
# a two-panel figure: on-reference | off-reference, each panel three violins
# (SNP, Indel/MNP, SV) with per-sample counts.
#
# Y-scales are free per panel because on-ref counts are ~10-30× higher than
# off-ref.
#
# Usage:
#   Rscript scripts/per-sample-types-slide.R <per-sample-types.tsv> <out-prefix> [--pdf]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

.script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
if (is.na(.script_dir) || !nzchar(.script_dir)) .script_dir <- "scripts"
source(file.path(.script_dir, "plot-helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript per-sample-types-slide.R <per-sample-types.tsv> <out-prefix> [--pdf]\n")
  quit(status = 1)
}
tsv_in   <- args[1]
prefix   <- args[2]
cli_pdf  <- "--pdf" %in% args
emit_pdf <- pdf_enabled(cli_flag = cli_pdf)

dt <- fread(tsv_in)
sv_types <- c("SV Insertion", "SV Deletion")
dt <- dt[variant_type != "Other"]
# GRCh38 is embedded in the v2.1 CHM13 graph as a haploid path, so its
# per-sample variant counts sit ~30-40% below the diploid modal (visible as
# the low outlier). Dropping it for the slide.
dt <- dt[sample != "GRCh38"]
dt[, type_class := fcase(
  variant_type == "SNP",         "SNP",
  variant_type %in% sv_types,    "SV",
  default =                      "Indel/MNP")]
dt[, type_class := factor(type_class, levels = c("SNP", "Indel/MNP", "SV"))]
dt[, ref_context := factor(ref_context,
                           levels = c("On-reference", "Off-reference"))]

# Collapse (sample, ref_context, type_class) → summed count so SVs (Ins+Del)
# and Indels/MNPs (Ins+Del+MNP) each land as one violin per sample.
ps <- dt[, .(count = sum(count)),
         by = .(sample, ref_context, type_class)]
n_samples <- uniqueN(ps$sample)

slide_theme <- theme_minimal(base_size = 16) +
  theme(
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title         = element_text(face = "bold"),
    strip.text         = element_text(face = "bold", size = 16),
    legend.position    = "none"
  )

p <- ggplot(ps, aes(x = type_class, y = count, fill = ref_context)) +
  geom_violin(width = 0.85, alpha = 0.55, scale = "width",
              color = NA) +
  geom_jitter(aes(color = ref_context),
              width = 0.12, size = 0.9, alpha = 0.55) +
  scale_fill_manual(values = c("On-reference"  = "steelblue",
                               "Off-reference" = "coral")) +
  scale_color_manual(values = c("On-reference"  = "steelblue4",
                                "Off-reference" = "coral4")) +
  scale_y_continuous(labels = scales::label_comma()) +
  facet_wrap(~ ref_context, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "Variants per Sample") +
  slide_theme

cat("Samples:", n_samples, "\n")
cat("Per-sample medians (rounded):\n")
print(ps[, .(median = format(round(median(count)), big.mark = ",")),
         by = .(ref_context, type_class)][order(ref_context, type_class)])

save_plot(p, paste0(prefix, ".per-sample-types-slide.png"),
          width = 12, height = 5.5, pdf = emit_pdf)
