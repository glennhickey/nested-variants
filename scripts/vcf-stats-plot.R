#!/usr/bin/env Rscript

# vcf-stats-plot.R — render VCF-stats figures from pre-computed TSVs
#
# Usage: Rscript scripts/vcf-stats-plot.R <prefix>
#          [--title TITLE] [--mode sites|variants] [--filter all|pass]
#          [--pdf]
#
# Inputs (read on demand from <prefix>.<name>.tsv; missing TSVs are silently
# skipped and the corresponding figure is left empty / not written):
#   <prefix>.vcf-stats.tsv               (always required)
#   <prefix>.size-dist.tsv               (optional — needed for size-dist plots)
#   <prefix>.af-spectrum.tsv             (optional — AF spectrum plot)
#   <prefix>.tr-counts.tsv               (optional — TR overlay on variant-types)
#   <prefix>.vcf-stats-by-annot.tsv      (optional — annotation breakdown)
#   <prefix>.giab-strat.tsv              (optional — GIAB stratification)
#   <prefix>.giab-strat-tstv.tsv         (optional — Ts/Tv overlay on GIAB)
#   <prefix>.per-sample-types.tsv        (optional — per-sample violins + SV)
#   <prefix>.per-sample-types-by-pop.tsv (optional — per-sample by super-pop)
#   <prefix>.per-sample-giab-strat.tsv   (optional — per-sample GIAB)
#
# Outputs (PNG always; PDF when --pdf):
#   <prefix>.variant-types.png(.pdf)
#   <prefix>.size-dist.png(.pdf)
#   <prefix>.size-dist-log.png(.pdf)
#   <prefix>.af-spectrum.png(.pdf)
#   <prefix>.variant-types-by-annot.png(.pdf)
#   <prefix>.giab-strat.png(.pdf)
#   <prefix>.per-sample-types.png(.pdf)
#   <prefix>.per-sample-types-by-pop.png(.pdf)
#   <prefix>.per-sample-sv-types.png(.pdf)
#   <prefix>.per-sample-giab-strat.png(.pdf)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript vcf-stats-plot.R <prefix> [--title T] [--mode M] [--filter F] [--pdf]\n")
  quit(status = 1)
}
prefix    <- args[1]
title     <- NULL
mode      <- "sites"
filter    <- "all"
emit_pdf  <- FALSE

i <- 2
while (i <= length(args)) {
  if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--mode" && i + 1 <= length(args)) {
    mode <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--filter" && i + 1 <= length(args)) {
    filter <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--pdf") {
    emit_pdf <- TRUE; i <- i + 1
  } else {
    i <- i + 1
  }
}
if (is.null(title)) title <- basename(prefix)

mode_label   <- if (mode == "sites") "(per site)" else "(per variant)"
filter_label <- if (filter == "pass") ", PASS only" else ""

# Shared save_plot helper
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
if (is.na(script_dir) || !nzchar(script_dir)) script_dir <- "scripts"
source(file.path(script_dir, "plot-helpers.R"))

# Paper-figure title suppression — clear the title up front so every
# `title = title` inside labs() below resolves to NULL when the user has
# disabled titles. Subtitles are wrapped with if_titles() at each call site.
if (!titles_on()) title <- NULL
emit_pdf <- pdf_enabled(cli_flag = emit_pdf)
save_png <- function(plot, file, width = 8, height = 6) {
  save_plot(plot, file, width = width, height = height, pdf = emit_pdf)
}

read_tsv <- function(suffix) {
  path <- paste0(prefix, ".", suffix)
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  tryCatch(fread(path), error = function(e) NULL)
}

theme_common <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

# ---------------------------------------------------------------------------
# 1. Variant-types bar chart
# ---------------------------------------------------------------------------
summary_dt <- read_tsv("vcf-stats.tsv")
if (is.null(summary_dt) || nrow(summary_dt) == 0) {
  cat("No vcf-stats.tsv at", prefix, "— skipping all plots.\n")
  quit(status = 0)
}

plot_dt <- summary_dt[, .(ref_context, variant_type, count)]
# Drop minor types <0.5% (matches the compute-side threshold)
total_variants <- sum(plot_dt$count)
type_totals <- plot_dt[, .(total = sum(count)), by = variant_type]
minor_types <- type_totals[total / total_variants < 0.005
                           & !variant_type %in% c("SV Insertion", "SV Deletion"), variant_type]
if (length(minor_types) > 0) {
  cat("Dropping minor variant types from plots (<0.5%):",
      paste(minor_types, collapse = ", "), "\n")
  plot_dt <- plot_dt[!variant_type %in% minor_types]
}
type_levels <- intersect(
  c("SNP", "MNP", "Insertion", "Deletion", "SV Insertion", "SV Deletion"),
  unique(plot_dt$variant_type))
plot_dt[, variant_type := factor(variant_type, levels = type_levels)]

# Ts/Tv labels for SNP bars
snp_row <- summary_dt[variant_type == "SNP" & !is.na(tstv_ratio)]
if (nrow(snp_row) > 0) {
  tstv_labels <- snp_row[, .(ref_context, label = paste0("Ts/Tv=", tstv_ratio))]
  plot_dt <- merge(plot_dt, tstv_labels, by = "ref_context", all.x = TRUE)
  plot_dt[!variant_type %in% "SNP", label := NA_character_]
} else {
  plot_dt[, label := NA_character_]
}

# Tandem-repeat overlay (faded portion atop indel bars)
tr_counts <- read_tsv("tr-counts.tsv")
has_tr <- !is.null(tr_counts) && nrow(tr_counts) > 0

if (has_tr) {
  plot_dt <- merge(plot_dt, tr_counts, by = c("ref_context", "variant_type"), all.x = TRUE)
  plot_dt[is.na(tr_count), tr_count := 0L]
  indel_types <- c("Insertion", "Deletion", "SV Insertion", "SV Deletion")
  tr_dt <- plot_dt[variant_type %in% indel_types & tr_count > 0,
                   .(ref_context, variant_type, count = tr_count)]
}

p1 <- ggplot(plot_dt, aes(x = variant_type, y = count, fill = ref_context)) +
  geom_col(position = "dodge", width = 0.7) +
  { if (has_tr) geom_col(data = tr_dt, aes(group = ref_context),
                         position = "dodge", width = 0.7, alpha = 0.35,
                         fill = "white", show.legend = FALSE) } +
  geom_text(aes(label = label),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                    name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = title,
       subtitle = if_titles(paste0("Variant Type Counts ", mode_label, filter_label,
                         if (has_tr) " (faded = tandem repeat)" else "")),
       x = "Variant Type", y = "Count") +
  theme_common
save_png(p1, paste0(prefix, ".variant-types.png"))

# ---------------------------------------------------------------------------
# 2. Size distribution (linear + log)
# ---------------------------------------------------------------------------
size_counts <- read_tsv("size-dist.tsv")
if (!is.null(size_counts) && nrow(size_counts) > 0) {
  size_counts[, panel := factor(panel, levels = c("Small (1-49 bp)", "Structural (50-1000 bp)"))]
  p2 <- ggplot(size_counts, aes(x = size, y = count, color = direction, linetype = ref_context)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~panel, scales = "free") +
    scale_color_manual(values = c("Insertion" = "coral", "Deletion" = "steelblue"),
                       name = "Direction") +
    scale_linetype_manual(values = c("On-reference" = "solid", "Off-reference" = "dashed"),
                          name = "Context") +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
    labs(title = title,
         subtitle = if_titles(paste0("Indel / SV Size Distribution ", mode_label, filter_label)),
         x = "Size (bp)", y = "Count") +
    theme_common +
    theme(axis.line = element_line(color = "black", linewidth = 0.5),
          strip.text = element_text(face = "bold"))
  save_png(p2, paste0(prefix, ".size-dist.png"), width = 12)

  p2log <- p2 + scale_y_log10(labels = scales::comma) + labs(y = "Count (log scale)")
  save_png(p2log, paste0(prefix, ".size-dist-log.png"), width = 12)
} else {
  file.create(paste0(prefix, ".size-dist.png"))
  file.create(paste0(prefix, ".size-dist-log.png"))
}

# ---------------------------------------------------------------------------
# 3. AF spectrum
# ---------------------------------------------------------------------------
af_counts <- read_tsv("af-spectrum.tsv")
if (!is.null(af_counts) && nrow(af_counts) > 0) {
  y_label <- if (mode == "sites") "Sites (log scale)" else "Variants (log scale)"
  p3 <- ggplot(af_counts, aes(x = af_plot, y = count, color = ref_context)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.2) +
    scale_color_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                       name = NULL) +
    scale_y_log10(labels = scales::comma) +
    scale_x_continuous(limits = c(-0.02, 1.02)) +
    labs(title = title,
         subtitle = if_titles(paste0("Non-Reference Allele Frequency Spectrum ",
                           mode_label, filter_label)),
         x = "Non-Reference Frequency", y = y_label) +
    theme_common
  save_png(p3, paste0(prefix, ".af-spectrum.png"), height = 6)
} else {
  file.create(paste0(prefix, ".af-spectrum.png"))
}

# ---------------------------------------------------------------------------
# 4. Variant types by annotation
# ---------------------------------------------------------------------------
ann_counts <- read_tsv("vcf-stats-by-annot.tsv")
if (!is.null(ann_counts) && nrow(ann_counts) > 0) {
  ann_counts <- ann_counts[variant_type %in% type_levels]
  ann_counts[, variant_type := factor(variant_type, levels = type_levels)]
  p_annot <- ggplot(ann_counts, aes(x = variant_type, y = count, fill = ref_context)) +
    geom_col(position = "dodge", width = 0.7) +
    facet_wrap(~annotation, scales = "free_y") +
    scale_fill_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                      name = NULL) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = title, subtitle = if_titles("Variant Types by Annotation Region")) +
    theme_common +
    coord_flip()
  save_png(p_annot, paste0(prefix, ".variant-types-by-annot.png"), width = 12, height = 8)
}

# ---------------------------------------------------------------------------
# 5. GIAB stratification
# ---------------------------------------------------------------------------
giab_counts <- read_tsv("giab-strat.tsv")
if (!is.null(giab_counts) && nrow(giab_counts) > 0) {
  region_levels <- c("Easy", "Segdup", "Hard")
  region_colors <- c("Easy" = "forestgreen", "Segdup" = "firebrick", "Hard" = "darkorange")
  giab_plot_dt <- giab_counts[variant_type %in% type_levels & !is.na(giab_region)]
  giab_plot_dt[, variant_type := factor(variant_type, levels = type_levels)]
  giab_plot_dt[, giab_region := factor(giab_region, levels = region_levels)]

  if (nrow(giab_plot_dt) > 0) {
    sv_type_names <- c("SV Insertion", "SV Deletion")
    giab_plot_dt[, type_class := fcase(
      variant_type == "SNP",                     "SNP",
      variant_type %in% sv_type_names,           "SV",
      default =                                  "Indel/MNP")]
    giab_plot_dt[, type_class := factor(type_class, levels = c("SNP", "Indel/MNP", "SV"))]

    giab_tstv <- read_tsv("giab-strat-tstv.tsv")
    if (!is.null(giab_tstv) && nrow(giab_tstv) > 0) {
      giab_tstv[, giab_region := factor(giab_region, levels = region_levels)]
      giab_tstv[, variant_type := factor("SNP", levels = type_levels)]
      giab_tstv[, type_class := factor("SNP", levels = c("SNP", "Indel/MNP", "SV"))]
      giab_tstv <- merge(
        giab_tstv,
        giab_plot_dt[variant_type == "SNP", .(giab_region, ref_context, count)],
        by = c("giab_region", "ref_context"), all.x = TRUE)
      giab_tstv[, label := paste0("Ts/Tv=", tstv_ratio)]
    }

    p_giab <- ggplot(giab_plot_dt,
                     aes(x = variant_type, y = count, fill = giab_region)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.7) +
      scale_fill_manual(values = region_colors, name = "GIAB Region") +
      scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0.02, 0.15))) +
      facet_wrap(vars(ref_context, type_class), scales = "free", ncol = 3) +
      labs(title = title,
           subtitle = if_titles(paste0("GIAB Genome Stratification ", mode_label, filter_label)),
           x = "Variant Type", y = "Count") +
      theme_common
    if (!is.null(giab_tstv) && nrow(giab_tstv) > 0) {
      p_giab <- p_giab +
        geom_text(data = giab_tstv,
                  aes(x = variant_type, y = count, label = label, group = giab_region),
                  position = position_dodge(width = 0.7),
                  vjust = -0.3, size = 2.8)
    }
    save_png(p_giab, paste0(prefix, ".giab-strat.png"), width = 12, height = 6)
  } else {
    file.create(paste0(prefix, ".giab-strat.png"))
  }
}

# ---------------------------------------------------------------------------
# 6. Per-sample types (violin + jitter); 6a SV-only sub-plot.
# ---------------------------------------------------------------------------
ps_counts <- read_tsv("per-sample-types.tsv")
if (!is.null(ps_counts) && nrow(ps_counts) > 0) {
  n_samples <- uniqueN(ps_counts$sample)
  ps_counts <- ps_counts[variant_type %in% type_levels & variant_type != "Other"]
  ps_counts[, variant_type := factor(variant_type, levels = type_levels)]
  sv_type_names <- c("SV Insertion", "SV Deletion")
  ps_counts[, type_class := fcase(
    variant_type == "SNP",                     "SNP",
    variant_type %in% sv_type_names,           "SV",
    default =                                  "Indel/MNP")]
  ps_counts[, type_class := factor(type_class, levels = c("SNP", "Indel/MNP", "SV"))]
  ps_counts[, size_class := fifelse(variant_type %in% sv_type_names,
                                    "Structural Variants", "Small Variants")]
  ps_counts[, size_class := factor(size_class,
                                   levels = c("Small Variants", "Structural Variants"))]
  ps_summary <- ps_counts[, .(mean_count = mean(count),
                              min_count = min(count),
                              max_count = max(count)),
                          by = .(variant_type, ref_context, size_class)]

  p_ps <- ggplot(ps_counts,
                 aes(x = variant_type, y = count, fill = ref_context)) +
    geom_violin(width = 0.7, alpha = 0.6, scale = "width",
                position = position_dodge(width = 0.7)) +
    geom_jitter(aes(color = ref_context),
                position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
                size = 1.2, alpha = 0.7) +
    scale_fill_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                      name = NULL) +
    scale_color_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                       name = NULL) +
    scale_y_continuous(labels = scales::comma) +
    facet_wrap(vars(ref_context, type_class), scales = "free", ncol = 3) +
    labs(title = title,
         subtitle = if_titles(paste0("Per-Sample Variant Counts ", mode_label, filter_label,
                           " (N=", n_samples, " samples)")),
         x = "Variant Type", y = "Count") +
    theme_common
  save_png(p_ps, paste0(prefix, ".per-sample-types.png"), width = 12)

  # 6a. Per-sample SV-only
  ps_sv <- ps_counts[size_class == "Structural Variants"]
  if (nrow(ps_sv) > 0 && sum(ps_sv$count) > 0) {
    ps_sv_summary <- ps_summary[size_class == "Structural Variants"]
    if (n_samples <= 20) {
      p_sv <- ggplot() +
        geom_col(data = ps_sv_summary,
                 aes(x = variant_type, y = mean_count, fill = ref_context),
                 width = 0.7, alpha = 0.6) +
        geom_errorbar(data = ps_sv_summary,
                      aes(x = variant_type, ymin = min_count, ymax = max_count),
                      width = 0.3) +
        geom_point(data = ps_sv,
                   aes(x = variant_type, y = count, color = ref_context),
                   position = position_jitter(width = 0.15),
                   size = 1.5, alpha = 0.8) +
        scale_fill_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                          name = NULL) +
        scale_color_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                           name = NULL) +
        scale_y_continuous(labels = scales::comma) +
        facet_wrap(~ ref_context, scales = "free_y") +
        labs(title = title,
             subtitle = if_titles(paste0("Per-Sample SV Counts ", mode_label, filter_label,
                               " (N=", n_samples, " samples, bars=mean)")),
             x = "Variant Type", y = "Count") +
        theme_common
    } else {
      p_sv <- ggplot(ps_sv,
                     aes(x = variant_type, y = count, fill = ref_context)) +
        geom_boxplot(width = 0.6, outlier.size = 1) +
        scale_fill_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                          name = NULL) +
        scale_y_continuous(labels = scales::comma) +
        facet_wrap(~ ref_context, scales = "free_y") +
        labs(title = title,
             subtitle = if_titles(paste0("Per-Sample SV Counts ", mode_label, filter_label,
                               " (N=", n_samples, " samples)")),
             x = "Variant Type", y = "Count") +
        theme_common
    }
    save_png(p_sv, paste0(prefix, ".per-sample-sv-types.png"), width = 10)
  } else {
    file.create(paste0(prefix, ".per-sample-sv-types.png"))
  }
}

# ---------------------------------------------------------------------------
# 7. Per-sample types by super-population
# ---------------------------------------------------------------------------
ps_pop <- read_tsv("per-sample-types-by-pop.tsv")
if (!is.null(ps_pop) && nrow(ps_pop) > 0) {
  ps_pop <- ps_pop[variant_type %in% type_levels & variant_type != "Other"]
  ps_pop[, variant_type := factor(variant_type, levels = type_levels)]
  sv_type_names <- c("SV Insertion", "SV Deletion")
  ps_pop[, type_class := fcase(
    variant_type == "SNP",                     "SNP",
    variant_type %in% sv_type_names,           "SV",
    default =                                  "Indel/MNP")]
  ps_pop[, type_class := factor(type_class, levels = c("SNP", "Indel/MNP", "SV"))]
  pop_levels <- c("AFR", "AMR", "EAS", "EUR", "SAS", "Unknown")
  pop_colors <- c("AFR"="#e31a1c", "AMR"="#ff7f00", "EAS"="#33a02c",
                  "EUR"="#1f78b4", "SAS"="#6a3d9a", "Unknown"="#000000")
  present <- intersect(pop_levels, unique(ps_pop$super_pop))
  extras  <- setdiff(unique(ps_pop$super_pop), pop_levels)
  ps_pop[, super_pop := factor(super_pop, levels = c(present, extras))]
  n_samples_pop <- uniqueN(ps_pop$sample)

  p_ps_pop <- ggplot(ps_pop,
                     aes(x = variant_type, y = count, fill = ref_context)) +
    geom_violin(width = 0.7, alpha = 0.6, scale = "width") +
    geom_jitter(aes(color = super_pop),
                position = position_jitter(width = 0.15, height = 0),
                size = 1.2, alpha = 0.85) +
    scale_fill_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                      name = "Ref Context") +
    scale_color_manual(values = pop_colors, name = "Super-Population",
                       drop = FALSE, na.value = "#000000") +
    scale_y_continuous(labels = scales::comma) +
    facet_wrap(vars(ref_context, type_class), scales = "free", ncol = 3) +
    labs(title = title,
         subtitle = if_titles(paste0("Per-Sample Variant Counts by Super-Population ",
                           mode_label, filter_label,
                           " (N=", n_samples_pop, " samples)")),
         x = "Variant Type", y = "Count") +
    theme_common
  save_png(p_ps_pop, paste0(prefix, ".per-sample-types-by-pop.png"), width = 12)

  # 7a. Per-sample violins coloured by AFR vs non-AFR, with variant types
  # collapsed to SNP / Indel/MNP / SV.  Used as figure 2G in the manuscript.
  # Drops 'Unknown' / reference-tagged samples (already filtered above) and
  # bins everything outside AFR (AMR, EAS, EUR, SAS) into "non-AFR".
  ps_afr <- copy(ps_pop)
  ps_afr[, type_class := factor(type_class, levels = c("SNP", "Indel/MNP", "SV"))]
  ps_afr <- ps_afr[, .(count = sum(count)),
                   by = .(sample, ref_context, type_class, super_pop)]
  ps_afr[, ancestry := fifelse(as.character(super_pop) == "AFR", "AFR", "non-AFR")]
  ps_afr[, ancestry := factor(ancestry, levels = c("AFR", "non-AFR"))]
  n_per_ancestry <- ps_afr[, .(n = uniqueN(sample)), by = ancestry]
  cat("Samples per ancestry bin:\n"); print(n_per_ancestry)
  ancestry_colors <- c("AFR" = "#e31a1c", "non-AFR" = "#1f78b4")

  p_ps_afr <- ggplot(ps_afr,
                     aes(x = ancestry, y = count, fill = ancestry)) +
    geom_violin(width = 0.7, alpha = 0.55, scale = "width") +
    geom_jitter(aes(color = ancestry),
                position = position_jitter(width = 0.15, height = 0),
                size = 1.0, alpha = 0.8) +
    scale_fill_manual(values = ancestry_colors, name = NULL, drop = FALSE) +
    scale_color_manual(values = ancestry_colors, name = NULL, drop = FALSE) +
    scale_y_continuous(labels = scales::comma) +
    facet_wrap(vars(ref_context, type_class), scales = "free", ncol = 3) +
    labs(x = NULL, y = "Count") +
    theme_common
  save_png(p_ps_afr, paste0(prefix, ".per-sample-types-by-afr.png"), width = 12)

  # 7b. Off-reference only, single y-axis: x = {SNP, Indel/MNP, SV},
  # fill = AFR vs non-AFR. Side-by-side violins per variant type.
  ps_afr_off <- ps_afr[ref_context == "Off-reference"]
  if (nrow(ps_afr_off) > 0) {
    p_ps_afr_off <- ggplot(ps_afr_off,
                           aes(x = type_class, y = count, fill = ancestry)) +
      geom_violin(width = 0.7, alpha = 0.55, scale = "width",
                  position = position_dodge(width = 0.8)) +
      geom_jitter(aes(color = ancestry),
                  position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                  size = 1.0, alpha = 0.85) +
      scale_fill_manual(values = ancestry_colors, name = NULL, drop = FALSE) +
      scale_color_manual(values = ancestry_colors, name = NULL, drop = FALSE) +
      scale_y_continuous(labels = scales::comma) +
      labs(x = NULL, y = "GRef Sites per Sample") +
      theme_common
    save_png(p_ps_afr_off,
             paste0(prefix, ".per-sample-types-by-afr-offref.png"),
             width = 8)
  }
}

# ---------------------------------------------------------------------------
# 8. Per-sample GIAB stratification
# ---------------------------------------------------------------------------
ps_giab <- read_tsv("per-sample-giab-strat.tsv")
if (!is.null(ps_giab) && nrow(ps_giab) > 0) {
  ps_giab <- ps_giab[variant_type %in% type_levels]
  region_levels_ps <- c("Easy", "Segdup", "Hard")
  ps_giab <- ps_giab[ps_giab_region %in% region_levels_ps]
  if (nrow(ps_giab) > 0) {
    ps_giab[, variant_type := factor(variant_type, levels = type_levels)]
    ps_giab[, ps_giab_region := factor(ps_giab_region, levels = region_levels_ps)]
    sv_type_names <- c("SV Insertion", "SV Deletion")
    ps_giab[, type_class := fcase(
      variant_type == "SNP",                     "SNP",
      variant_type %in% sv_type_names,           "SV",
      default =                                  "Indel/MNP")]
    ps_giab[, type_class := factor(type_class, levels = c("SNP", "Indel/MNP", "SV"))]
    region_colors_ps <- c("Easy" = "forestgreen", "Segdup" = "firebrick",
                          "Hard" = "darkorange")
    n_samples_giab <- uniqueN(ps_giab$sample)

    p_ps_giab <- ggplot(ps_giab,
                        aes(x = variant_type, y = count, fill = ps_giab_region)) +
      geom_violin(width = 0.7, alpha = 0.6, scale = "width",
                  position = position_dodge(width = 0.7)) +
      geom_jitter(aes(color = ps_giab_region),
                  position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
                  size = 1.2, alpha = 0.7) +
      scale_fill_manual(values = region_colors_ps, name = "GIAB Region") +
      scale_color_manual(values = region_colors_ps, name = "GIAB Region") +
      scale_y_continuous(labels = scales::comma) +
      facet_wrap(vars(ref_context, type_class), scales = "free", ncol = 3) +
      labs(title = title,
           subtitle = if_titles(paste0("Per-Sample GIAB Stratification ", mode_label, filter_label,
                             " (N=", n_samples_giab, " samples)")),
           x = "Variant Type", y = "Count") +
      theme_common
    save_png(p_ps_giab, paste0(prefix, ".per-sample-giab-strat.png"), width = 12)
  } else {
    file.create(paste0(prefix, ".per-sample-giab-strat.png"))
  }
}

cat("Done.\n")
