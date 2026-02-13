#!/usr/bin/env Rscript

# vcf-stats.R — VCF variant statistics with on-ref vs off-ref breakdown
#
# Usage: Rscript scripts/vcf-stats.R <input.vcf.gz> <output_prefix> [--title TITLE]
#
# Outputs:
#   {prefix}.vcf-stats.tsv      — summary table
#   {prefix}.variant-types.png  — grouped bar chart of variant types
#   {prefix}.size-dist.png      — indel/SV size distribution
#   {prefix}.af-spectrum.png    — allele frequency histogram (only when AF present)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript vcf-stats.R <input.vcf.gz> <output_prefix> [--title TITLE]\n")
  quit(status = 1)
}

vcf    <- args[1]
prefix <- args[2]
title  <- NULL

i <- 3
while (i <= length(args)) {
  if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}
if (is.null(title)) title <- basename(vcf)

# ---------------------------------------------------------------------------
# Read VCF via bcftools (normalize multi-allelic, extract fields)
# ---------------------------------------------------------------------------
cat("Reading VCF:", vcf, "\n")

# Try with AF first (filter out all-homref sites from vg call -A)
cmd_af <- sprintf(
  "bcftools view -c1 '%s' 2>/dev/null | bcftools norm -m- 2>/dev/null | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%INFO/AF\\n' 2>/dev/null",
  vcf
)
dt <- tryCatch(
  fread(cmd = cmd_af, col.names = c("CHROM", "POS", "REF", "ALT", "AF")),
  error = function(e) NULL,
  warning = function(w) NULL
)

has_af <- !is.null(dt) && nrow(dt) > 0 && !all(is.na(dt$AF) | dt$AF == ".")

if (is.null(dt) || nrow(dt) == 0) {
  # Fallback: read without AF
  cmd_no_af <- sprintf(
    "bcftools view -c1 '%s' 2>/dev/null | bcftools norm -m- 2>/dev/null | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' 2>/dev/null",
    vcf
  )
  dt <- fread(cmd = cmd_no_af, col.names = c("CHROM", "POS", "REF", "ALT"))
  has_af <- FALSE
}

if (!has_af && "AF" %in% names(dt)) {
  dt[, AF := NULL]
}

if (has_af) {
  dt[, AF := as.numeric(AF)]
  dt <- dt[!is.na(AF)]
}

cat("Read", nrow(dt), "variant records\n")
cat("AF available:", has_af, "\n")

if (nrow(dt) == 0) {
  cat("No variants found. Exiting.\n")
  quit(status = 0)
}

# ---------------------------------------------------------------------------
# Classify variants
# ---------------------------------------------------------------------------
dt[, ref_len := nchar(REF)]
dt[, alt_len := nchar(ALT)]
dt[, size := abs(alt_len - ref_len)]

dt[, variant_type := fifelse(
  ref_len == 1L & alt_len == 1L, "SNP",
  fifelse(size < 50L, "Indel", "SV")
)]

# On-ref vs off-ref
dt[, ref_context := fifelse(
  grepl("_[0-9]+_alt$", CHROM), "Off-reference", "On-reference"
)]

# Ts/Tv classification (SNPs only)
transitions <- c("AG", "GA", "CT", "TC")
dt[variant_type == "SNP", tstv := fifelse(
  paste0(REF, ALT) %in% transitions, "Ts", "Tv"
)]

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
summary_dt <- dt[, .(count = .N), by = .(ref_context, variant_type)]

# Add Ts/Tv counts for SNPs
tstv_dt <- dt[variant_type == "SNP",
              .(ts = sum(tstv == "Ts"), tv = sum(tstv == "Tv")),
              by = .(ref_context)]
tstv_dt[, tstv_ratio := round(ts / tv, 2)]

summary_dt <- merge(summary_dt, tstv_dt, by = "ref_context", all.x = TRUE)
summary_dt[variant_type != "SNP", c("ts", "tv", "tstv_ratio") := .(NA, NA, NA)]

# Add AF summaries if available
if (has_af) {
  af_dt <- dt[, .(mean_af = round(mean(AF), 4), median_af = round(median(AF), 4)),
              by = .(ref_context, variant_type)]
  summary_dt <- merge(summary_dt, af_dt, by = c("ref_context", "variant_type"), all.x = TRUE)
}

# Sort for readability
setorder(summary_dt, ref_context, variant_type)

tsv_file <- paste0(prefix, ".vcf-stats.tsv")
fwrite(summary_dt, tsv_file, sep = "\t")
cat("Wrote summary:", tsv_file, "\n")

# Print summary to stdout
cat("\n")
print(summary_dt)
cat("\n")

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
# Plot 1: Variant type bar chart
# ---------------------------------------------------------------------------
plot_dt <- dt[, .(count = .N), by = .(ref_context, variant_type)]

# Get Ts/Tv labels for SNP bars
tstv_labels <- tstv_dt[, .(ref_context, label = paste0("Ts/Tv=", tstv_ratio))]
plot_dt <- merge(plot_dt, tstv_labels, by = "ref_context", all.x = TRUE)
plot_dt[variant_type != "SNP", label := NA_character_]

# Order variant types
plot_dt[, variant_type := factor(variant_type, levels = c("SNP", "Indel", "SV"))]

p1 <- ggplot(plot_dt, aes(x = variant_type, y = count, fill = ref_context)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = label),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                    name = NULL) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  labs(title = title, subtitle = "Variant Type Counts",
       x = "Variant Type", y = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

save_png(p1, paste0(prefix, ".variant-types.png"))

# ---------------------------------------------------------------------------
# Plot 2: Size distribution (indels + SVs only)
# ---------------------------------------------------------------------------
size_dt <- dt[variant_type %in% c("Indel", "SV") & size > 0]

if (nrow(size_dt) > 0) {
  p2 <- ggplot(size_dt, aes(x = size, fill = ref_context)) +
    geom_histogram(bins = 50, position = "identity", alpha = 0.5) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "grey40") +
    annotate("text", x = 50, y = Inf, label = "50 bp", vjust = 1.5, hjust = -0.1,
             size = 3, color = "grey40") +
    scale_x_log10(labels = scales::comma) +
    scale_fill_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                      name = NULL) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
    labs(title = title, subtitle = "Indel / SV Size Distribution",
         x = "Size (bp, log scale)", y = "Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )

  save_png(p2, paste0(prefix, ".size-dist.png"))
} else {
  cat("No indels/SVs with size > 0; skipping size distribution plot.\n")
  # Create empty file so Snakemake sees the output
  file.create(paste0(prefix, ".size-dist.png"))
}

# ---------------------------------------------------------------------------
# Plot 3: AF spectrum (only when AF is available)
# ---------------------------------------------------------------------------
if (has_af) {
  dt[, variant_type_f := factor(variant_type, levels = c("SNP", "Indel", "SV"))]

  p3 <- ggplot(dt, aes(x = AF, color = ref_context)) +
    geom_freqpoly(bins = 50, linewidth = 0.8) +
    facet_wrap(~ variant_type_f, ncol = 1, scales = "free_y") +
    scale_color_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                       name = NULL) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
    labs(title = title, subtitle = "Allele Frequency Spectrum",
         x = "Allele Frequency", y = "Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      strip.text = element_text(face = "bold"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )

  save_png(p3, paste0(prefix, ".af-spectrum.png"), height = 8)
  cat("AF spectrum plot generated.\n")
} else {
  cat("No AF field; skipping allele frequency spectrum plot.\n")
}

cat("Done.\n")
