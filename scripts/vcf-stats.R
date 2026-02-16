#!/usr/bin/env Rscript

# vcf-stats.R — VCF variant statistics with on-ref vs off-ref breakdown
#
# Usage: Rscript scripts/vcf-stats.R <input.vcf.gz> <output_prefix>
#          [--title TITLE] [--mode sites|variants]
#
# Modes:
#   sites    — (default) read VCF as-is; multi-allelic sites classified by largest allele
#   variants — pipe through `bcftools norm -m-` first to split multi-allelic records
#
# Outputs:
#   {prefix}.vcf-stats.tsv      — summary table
#   {prefix}.variant-types.png  — grouped bar chart of variant types
#   {prefix}.size-dist.png      — indel/SV size distribution (two-panel: indels + SVs)
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
  cat("Usage: Rscript vcf-stats.R <input.vcf.gz> <output_prefix> [--title TITLE] [--mode sites|variants]\n")
  quit(status = 1)
}

vcf    <- args[1]
prefix <- args[2]
title  <- NULL
mode   <- "sites"

i <- 3
while (i <= length(args)) {
  if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--mode" && i + 1 <= length(args)) {
    mode <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}
if (!mode %in% c("sites", "variants")) {
  cat("Error: --mode must be 'sites' or 'variants', got '", mode, "'\n", sep = "")
  quit(status = 1)
}
if (is.null(title)) title <- basename(vcf)

mode_label <- if (mode == "sites") "(per site)" else "(per variant)"

# ---------------------------------------------------------------------------
# Read VCF
# ---------------------------------------------------------------------------
cat("Reading VCF:", vcf, " (mode:", mode, ")\n")

# Build pipeline prefix: variants mode pipes through bcftools norm -m- first
# +fill-tags computes AF from genotypes when INFO/AF is absent (e.g., single-sample vg call)
if (mode == "variants") {
  pipe_prefix <- sprintf("bcftools norm -m- '%s' 2>/dev/null | bcftools view -c1 2>/dev/null | bcftools +fill-tags - -- -t AF 2>/dev/null", vcf)
} else {
  pipe_prefix <- sprintf("bcftools view -c1 '%s' 2>/dev/null | bcftools +fill-tags - -- -t AF 2>/dev/null", vcf)
}

# Try with AF first (filter out all-homref sites from vg call -A)
cmd_af <- sprintf(
  "%s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%INFO/AF\\n' 2>/dev/null",
  pipe_prefix
)
dt <- tryCatch(
  fread(cmd = cmd_af, col.names = c("CHROM", "POS", "REF", "ALT", "AF_str")),
  error = function(e) NULL,
  warning = function(w) NULL
)

has_af <- !is.null(dt) && nrow(dt) > 0 && !all(is.na(dt$AF_str) | dt$AF_str == ".")

if (is.null(dt) || nrow(dt) == 0) {
  # Fallback: read without AF
  cmd_no_af <- sprintf(
    "%s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' 2>/dev/null",
    pipe_prefix
  )
  dt <- fread(cmd = cmd_no_af, col.names = c("CHROM", "POS", "REF", "ALT"))
  has_af <- FALSE
}

if (!has_af && "AF_str" %in% names(dt)) {
  dt[, AF_str := NULL]
}

cat("Read", nrow(dt), "variant records\n")
cat("AF available:", has_af, "\n")

if (nrow(dt) == 0) {
  cat("No variants found. Exiting.\n")
  quit(status = 0)
}

# ---------------------------------------------------------------------------
# Classify variants (per-site, no multi-allelic splitting)
# For multi-allelic sites: SV > Indel > SNP (use largest allele)
# ---------------------------------------------------------------------------
dt[, ref_len := nchar(REF)]

# Compute max size across comma-separated ALT alleles
dt[, size := sapply(seq_len(.N), function(i) {
  alts <- unlist(strsplit(ALT[i], ","))
  max(abs(nchar(alts) - ref_len[i]))
})]

dt[, variant_type := fifelse(
  size == 0L, "SNP",
  fifelse(size < 50L, "Indel", "SV")
)]

# On-ref vs off-ref
dt[, ref_context := fifelse(
  grepl("_[0-9]+_alt$", CHROM), "Off-reference", "On-reference"
)]

# Ts/Tv classification (SNPs only — biallelic SNPs where ALT has no comma)
dt[variant_type == "SNP" & !grepl(",", ALT), tstv := {
  transitions <- c("AG", "GA", "CT", "TC")
  fifelse(paste0(REF, ALT) %in% transitions, "Ts", "Tv")
}]

# ---------------------------------------------------------------------------
# AF: compute per-site non-reference frequency (sum of alt AFs)
# ---------------------------------------------------------------------------
if (has_af) {
  # Ensure AF_str is character (fread may auto-detect as numeric after norm -m-)
  if (!is.character(dt$AF_str)) dt[, AF_str := as.character(AF_str)]
  dt[, nonref_af := sapply(AF_str, function(x) {
    vals <- as.numeric(unlist(strsplit(x, ",")))
    min(sum(vals, na.rm = TRUE), 1.0)
  })]
  dt[, AF_str := NULL]
}

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
summary_dt <- dt[, .(count = .N), by = .(ref_context, variant_type)]

# Add Ts/Tv counts for SNPs
tstv_dt <- dt[variant_type == "SNP" & !is.na(tstv),
              .(ts = sum(tstv == "Ts"), tv = sum(tstv == "Tv")),
              by = .(ref_context)]
tstv_dt[, tstv_ratio := round(ts / tv, 2)]

summary_dt <- merge(summary_dt, tstv_dt, by = "ref_context", all.x = TRUE)
summary_dt[variant_type != "SNP", c("ts", "tv", "tstv_ratio") := .(NA, NA, NA)]

# Add AF summaries if available
if (has_af) {
  af_dt <- dt[, .(mean_af = round(mean(nonref_af), 4),
                   median_af = round(median(nonref_af), 4)),
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
if (nrow(tstv_dt) > 0) {
  tstv_labels <- tstv_dt[, .(ref_context, label = paste0("Ts/Tv=", tstv_ratio))]
  plot_dt <- merge(plot_dt, tstv_labels, by = "ref_context", all.x = TRUE)
  plot_dt[variant_type != "SNP", label := NA_character_]
} else {
  plot_dt[, label := NA_character_]
}

# Order variant types
plot_dt[, variant_type := factor(variant_type, levels = c("SNP", "Indel", "SV"))]

p1 <- ggplot(plot_dt, aes(x = variant_type, y = count, fill = ref_context)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = label),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                    name = NULL) +
  scale_y_log10(labels = scales::comma) +
  labs(title = title, subtitle = paste("Variant Type Counts", mode_label),
       x = "Variant Type", y = "Count (log scale)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

save_png(p1, paste0(prefix, ".variant-types.png"))

# ---------------------------------------------------------------------------
# Plot 2: Size distribution — two-panel (Indels 1-49 bp / SVs 50-1000 bp)
# ---------------------------------------------------------------------------
size_dt <- dt[variant_type %in% c("Indel", "SV") & size > 0]

if (nrow(size_dt) > 0) {
  size_dt[, panel := fifelse(size < 50L, "Indels (1-49 bp)", "SVs (50-1000 bp)")]
  size_dt[, panel := factor(panel, levels = c("Indels (1-49 bp)", "SVs (50-1000 bp)"))]
  size_counts <- size_dt[size <= 1000, .(count = .N), by = .(size, ref_context, panel)]

  p2 <- ggplot(size_counts, aes(x = size, y = count, color = ref_context)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~panel, scales = "free") +
    scale_color_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                       name = NULL) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
    labs(title = title, subtitle = paste("Indel / SV Size Distribution", mode_label),
         x = "Size (bp)", y = "Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold")
    )

  save_png(p2, paste0(prefix, ".size-dist.png"), width = 12)
} else {
  cat("No indels/SVs with size > 0; skipping size distribution plot.\n")
  # Create empty file so Snakemake sees the output
  file.create(paste0(prefix, ".size-dist.png"))
}

# ---------------------------------------------------------------------------
# Plot 3: AF spectrum — per-site non-reference frequency (log y-axis)
# ---------------------------------------------------------------------------
if (has_af) {
  cat("AF spectrum: ", nrow(dt), " sites with non-ref AF\n")

  # Adaptive plot: bar chart for discrete AF (few samples), freqpoly for continuous
  dt[, af_rounded := round(nonref_af, 2)]
  n_unique <- length(unique(dt$af_rounded))
  cat("Unique AF values (rounded):", n_unique, "\n")

  y_label <- if (mode == "sites") "Sites (log scale)" else "Variants (log scale)"
  common_theme <- theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )

  if (n_unique <= 30) {
    # Discrete AF: count exact rounded values, use bar chart
    af_counts <- dt[, .(count = .N), by = .(af_rounded, ref_context)]

    # Bar width: 80% of minimum spacing (or 0.02 if only one value)
    af_unique <- sort(unique(af_counts$af_rounded))
    bar_width <- if (length(af_unique) > 1) 0.8 * min(diff(af_unique)) else 0.02

    p3 <- ggplot(af_counts, aes(x = af_rounded, y = count, fill = ref_context)) +
      geom_col(position = position_dodge(width = bar_width), width = bar_width) +
      scale_fill_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                        name = NULL) +
      scale_y_log10(labels = scales::comma) +
      scale_x_continuous(breaks = af_unique, limits = c(-0.02, 1.02)) +
      labs(title = title, subtitle = paste("Non-Reference Allele Frequency Spectrum", mode_label),
           x = "Non-Reference Frequency", y = y_label) +
      common_theme
  } else {
    # Continuous AF: use frequency polygon (many samples)
    p3 <- ggplot(dt, aes(x = nonref_af, color = ref_context)) +
      geom_freqpoly(bins = 50, linewidth = 0.8) +
      scale_color_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                         name = NULL) +
      scale_y_log10(labels = scales::comma) +
      labs(title = title, subtitle = paste("Non-Reference Allele Frequency Spectrum", mode_label),
           x = "Non-Reference Frequency (sum of alt AFs per site)", y = y_label) +
      common_theme
  }

  dt[, af_rounded := NULL]
  save_png(p3, paste0(prefix, ".af-spectrum.png"), height = 6)
  cat("AF spectrum plot generated.\n")
} else {
  cat("No AF field; skipping allele frequency spectrum plot.\n")
  # Create empty file so Snakemake sees the output
  file.create(paste0(prefix, ".af-spectrum.png"))
}

cat("Done.\n")
