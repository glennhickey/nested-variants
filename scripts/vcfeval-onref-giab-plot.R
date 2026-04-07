#!/usr/bin/env Rscript
# vcfeval-onref-giab-plot.R — Bar chart of on-ref TP/FP/FN by GIAB region
#
# Usage: Rscript vcfeval-onref-giab-plot.R <input.tsv> <output.png> [--title TITLE]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
input_tsv <- args[1]
output_png <- args[2]
title <- "Call vs DeepVariant On-Reference"

i <- 3
while (i <= length(args)) {
  if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

dt <- fread(input_tsv)
if (nrow(dt) == 0) {
  cat("No data — creating empty plot.\n")
  file.create(output_png)
  quit(status = 0)
}

snp_colors <- c("TP" = "#4DAF4A", "FP" = "#E41A1C", "FN" = "#377EB8")

for (vt in c("SNP", "Indel")) {
  sub <- dt[variant_type == vt]
  if (nrow(sub) == 0) next

  p <- ggplot(sub, aes(x = giab_region, y = count, fill = category)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = snp_colors, name = NULL) +
    scale_y_continuous(labels = comma) +
    labs(title = title,
         subtitle = paste0(vt, " TP/FP/FN by GIAB Region (on-reference)"),
         x = "GIAB Region", y = paste(vt, "Count")) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )

  suffix <- if (vt == "SNP") "" else paste0("-", tolower(vt))
  out <- sub("\\.png$", paste0(suffix, ".png"), output_png)
  ggsave(out, p, width = 8, height = 6, dpi = 300, bg = "white",
         device = grDevices::png, type = "cairo")
  cat("Saved:", out, "\n")
}
