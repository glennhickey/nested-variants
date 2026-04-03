#!/usr/bin/env Rscript
# vcf-stats-by-chrom.R — Per-chromosome breakdown plots
#
# Usage:
#   Rscript vcf-stats-by-chrom.R <output_prefix> \
#     --stats chrom1:stats.tsv chrom2:stats.tsv ... \
#     --per-sample chrom1:per-sample.tsv chrom2:per-sample.tsv ... \
#     --records chrom1:records.tsv chrom2:records.tsv ...
#
# Produces:
#   <prefix>.snp-counts.png        — SNP counts per chromosome
#   <prefix>.indel-counts.png      — MNP + Indel counts per chromosome
#   <prefix>.sv-counts.png         — SV counts per chromosome
#   <prefix>.af-spectrum.png       — AF spectrum lines per chromosome
#   <prefix>.per-sample-snp.png    — Per-sample SNP box plots per chromosome
#   <prefix>.per-sample-indel.png  — Per-sample MNP/Indel box plots per chromosome
#   <prefix>.per-sample-sv.png     — Per-sample SV box plots per chromosome

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
output_prefix <- args[1]
args <- args[-1]

# Parse --stats, --per-sample, --records sections
stats_pairs <- c()
ps_pairs <- c()
rec_pairs <- c()
current <- NULL
for (a in args) {
  if (a == "--stats") { current <- "stats"; next }
  if (a == "--per-sample") { current <- "ps"; next }
  if (a == "--records") { current <- "rec"; next }
  if (is.null(current)) next
  if (current == "stats") stats_pairs <- c(stats_pairs, a)
  else if (current == "ps") ps_pairs <- c(ps_pairs, a)
  else if (current == "rec") rec_pairs <- c(rec_pairs, a)
}

# Helper: parse chrom:path pairs
read_pairs <- function(pairs) {
  rbindlist(lapply(pairs, function(pair) {
    parts <- strsplit(pair, ":", fixed = TRUE)[[1]]
    chrom <- parts[1]
    path <- paste(parts[-1], collapse = ":")
    dt <- fread(path)
    dt[, chrom := chrom]
    dt
  }))
}

# Natural sort for chromosome factor levels
nat_sort <- function(v) {
  u <- unique(v)
  u[order(nchar(u), u)]
}

# Variant type groups
type_groups <- list(
  snp   = "SNP",
  indel = c("MNP", "Insertion", "Deletion"),
  sv    = c("SV Insertion", "SV Deletion")
)
group_labels <- c(snp = "SNP", indel = "MNP / Indel", sv = "SV")

# --- Variant type counts per chromosome (3 plots) ---
if (length(stats_pairs) > 0) {
  stats <- read_pairs(stats_pairs)
  stats[, chrom := factor(chrom, levels = nat_sort(chrom))]

  for (grp in names(type_groups)) {
    types <- type_groups[[grp]]
    grp_data <- stats[variant_type %in% types]
    # Sum across types within group
    grp_agg <- grp_data[, .(count = sum(count)), by = .(chrom, ref_context)]

    p <- ggplot(grp_agg, aes(x = chrom, y = count, fill = ref_context)) +
      geom_col(position = "dodge") +
      scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
      scale_fill_manual(values = c("Off-reference" = "#E74C3C",
                                   "On-reference" = "#3498DB")) +
      labs(title = paste(group_labels[grp], "Counts by Chromosome"),
           x = NULL, y = "Count", fill = NULL) +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(face = "bold"))

    # Add Ts/Tv labels above each SNP bar (like 2A)
    if (grp == "snp") {
      tstv <- grp_data[!is.na(tstv_ratio)]
      if (nrow(tstv) > 0) {
        p <- p + geom_text(data = tstv,
                           aes(x = chrom, y = count, fill = ref_context,
                               label = paste0("Ts/Tv=", tstv_ratio)),
                           position = position_dodge(0.9),
                           vjust = -0.3, show.legend = FALSE, size = 3)
      }
    }

    fname <- paste0(output_prefix, ".", grp, "-counts.png")
    ggsave(fname, p, width = 8, height = 5, dpi = 300,
           device = grDevices::png, type = "cairo")
    cat("Saved:", fname, "\n")
  }
}

# --- AF spectrum per chromosome ---
if (length(rec_pairs) > 0) {
  recs <- read_pairs(rec_pairs)
  recs <- recs[variant_type == "SNP" & !is.na(nonref_af)]
  recs[, chrom := factor(chrom, levels = nat_sort(chrom))]

  af_step <- 0.05
  recs[, af_bin := floor(nonref_af / af_step) * af_step + af_step / 2]
  recs[af_bin > 1, af_bin := 1]

  af_counts <- recs[, .(sites = .N), by = .(chrom, ref_context, af_bin)]

  p_c <- ggplot(af_counts, aes(x = af_bin, y = sites,
                                color = chrom, linetype = ref_context)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    scale_y_log10(labels = comma) +
    scale_linetype_manual(values = c("Off-reference" = "solid",
                                     "On-reference" = "dashed")) +
    labs(title = "SNP Allele Frequency Spectrum by Chromosome",
         x = "Non-Reference Frequency", y = "Sites (log scale)",
         color = "Chromosome", linetype = NULL) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(paste0(output_prefix, ".af-spectrum.png"), p_c,
         width = 8, height = 5, dpi = 300,
         device = grDevices::png, type = "cairo")
  cat("Saved:", paste0(output_prefix, ".af-spectrum.png"), "\n")
}

# --- Per-sample box plots per chromosome (3 plots) ---
if (length(ps_pairs) > 0) {
  ps <- read_pairs(ps_pairs)
  ps[, chrom := factor(chrom, levels = nat_sort(chrom))]

  for (grp in names(type_groups)) {
    types <- type_groups[[grp]]
    ps_grp <- ps[variant_type %in% types & ref_context == "Off-reference"]
    # Sum across types within group per sample
    ps_agg <- ps_grp[, .(count = sum(count)), by = .(sample, chrom)]

    n_samples <- length(unique(ps_agg$sample))

    p <- ggplot(ps_agg, aes(x = chrom, y = count, fill = chrom)) +
      geom_boxplot(alpha = 0.3, outlier.shape = NA) +
      geom_jitter(aes(color = chrom), width = 0.2, size = 0.8, alpha = 0.4,
                  show.legend = FALSE) +
      scale_y_continuous(labels = comma) +
      labs(title = paste(group_labels[grp], "(off-ref, n=", n_samples, "samples)"),
           x = NULL, y = "Count per Sample",
           fill = "Chromosome") +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(face = "bold"))

    fname <- paste0(output_prefix, ".per-sample-", grp, ".png")
    ggsave(fname, p, width = 8, height = 5, dpi = 300,
           device = grDevices::png, type = "cairo")
    cat("Saved:", fname, "\n")
  }
}
