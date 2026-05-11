#!/usr/bin/env Rscript

# annotation-plots.R — Annotation overlap plots for augref alt segments
#
# Usage: Rscript scripts/annotation-plots.R <per_segment.tsv> <output_prefix> [--title TITLE]
#          [--min-overlap FRAC]  (minimum overlap fraction to count; default 0)
#        Rscript scripts/annotation-plots.R <per_segment.tsv> <output_prefix>
#          --vcf <biallelic_snps.vcf.gz> --augref-prefix <prefix> --filter <all|pass>
#          --title TITLE [--min-overlap FRAC]
#
# Base mode outputs:
#   {prefix}.annot-summary.png       — fraction of alt bp overlapping each annotation (source vs ref)
#   {prefix}.annot-scatter.png       — per-segment source_len vs source_overlap_frac, faceted by annotation
#   {prefix}.annot-repeats.png       — repeat class breakdown (only when repeat class data present)
#   {prefix}.annot-cooccur.png       — heatmap of annotation co-occurrence across coord types
#   {prefix}.annot-ancestry.png      — PCLAI ancestry co-occurrence heatmap (when pclai data present)
#   {prefix}.annot-pclai-summary.png — PCLAI ancestry breakdown bar chart (when pclai data present)
#   {prefix}.annot-stats.tsv         — tabular summary (excludes pclai)
#
# VCF mode outputs (when --vcf is provided):
#   {prefix}-counts.{filt}.png            — SNP count heatmap by annotation co-occurrence
#   {prefix}-tstv.{filt}.png              — Ts/Tv ratio heatmap by annotation co-occurrence
#   {pclai_prefix}-counts.{filt}.png      — ancestry SNP count heatmap (when pclai data present)
#   {pclai_prefix}-tstv.{filt}.png        — ancestry SNP Ts/Tv heatmap (when pclai data present)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# Title suppression for paper figures (FIGURE_TITLES=0 → labs() titles → NULL).
.script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
if (is.na(.script_dir) || !nzchar(.script_dir)) .script_dir <- "scripts"
source(file.path(.script_dir, "plot-helpers.R"))
.titles_on <- titles_on

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript annotation-plots.R <per_segment.tsv> <output_prefix> [--title TITLE]\n")
  quit(status = 1)
}

input_file    <- args[1]
prefix        <- args[2]
title         <- NULL
vcf_file      <- NULL
augref_prefix <- NULL
filt          <- NULL
min_overlap   <- 0
genome_coverage_file <- NULL

i <- 3
while (i <= length(args)) {
  if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--vcf" && i + 1 <= length(args)) {
    vcf_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--augref-prefix" && i + 1 <= length(args)) {
    augref_prefix <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--filter" && i + 1 <= length(args)) {
    filt <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--min-overlap" && i + 1 <= length(args)) {
    min_overlap <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--genome-coverage" && i + 1 <= length(args)) {
    genome_coverage_file <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}
if (is.null(title)) title <- "Annotation Overlap"
# Paper-figure title suppression — null the title up front so every
# `title = title` in labs() below resolves to NULL.  Subtitles are wrapped
# with if_titles() at each call site.
if (!.titles_on()) title <- NULL
vcf_mode <- !is.null(vcf_file)

# ---------------------------------------------------------------------------
# Read data
# ---------------------------------------------------------------------------
cat("Reading:", input_file, "\n")
cat("Min overlap fraction:", min_overlap, "\n")
dt <- fread(input_file)

if (nrow(dt) == 0) {
  cat("No data found. Exiting.\n")
  quit(status = 0)
}

# Read genome-wide coverage if provided
genome_cov <- NULL
if (!is.null(genome_coverage_file) && file.exists(genome_coverage_file)) {
  genome_cov <- fread(genome_coverage_file)
  cat("Read genome-wide coverage from:", genome_coverage_file,
      "(", nrow(genome_cov), "rows )\n")
}

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
# Aggregated view: _total rows for grouped annotations, regular for ungrouped.
# Used by summary, scatter, co-occurrence, stats, and VCF-mode plots.
# ---------------------------------------------------------------------------
total_anns <- unique(dt[annotation_class == "_total"]$annotation)
dt_agg <- rbind(
  dt[annotation %in% total_anns & annotation_class == "_total"],
  dt[!annotation %in% total_anns]
)
dt_shared <- dt_agg[annotation != "pclai"]

if (vcf_mode) {
  # =========================================================================
  # VCF mode: SNP count and Ts/Tv heatmaps by annotation co-occurrence
  # =========================================================================

  seg_ann <- dt_shared[, .(source_overlap_bp = sum(source_overlap_bp),
                           ref_overlap_bp = sum(ref_overlap_bp),
                           source_overlap_frac = sum(source_overlap_frac),
                           ref_overlap_frac = sum(ref_overlap_frac)),
                       by = .(augref_path, annotation)]

  # Read biallelic SNP VCF
  cmd <- sprintf("bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' '%s'", vcf_file)
  snps <- tryCatch(fread(cmd = cmd, col.names = c("CHROM", "POS", "REF", "ALT")),
                   error = function(e) data.table(CHROM = character(), POS = integer(),
                                                  REF = character(), ALT = character()))

  # Filter to off-reference SNPs
  snps <- snps[grepl("_[0-9]+_alt$", CHROM)]

  if (nrow(snps) == 0) {
    cat("No off-reference biallelic SNPs found. Creating empty placeholder files.\n")
    for (suf in c("counts", "tstv")) {
      out_file <- paste0(prefix, "-", suf, ".", filt, ".png")
      grDevices::png(out_file, width = 100, height = 100)
      par(mar = c(0, 0, 0, 0))
      plot.new()
      dev.off()
      cat("Saved placeholder:", out_file, "\n")
    }
  } else {
    # Map CHROM → augref_path
    known_paths <- unique(dt$augref_path)
    if (!any(snps$CHROM %in% known_paths)) {
      snps[, augref_path := paste0(augref_prefix, CHROM)]
    } else {
      snps[, augref_path := CHROM]
    }

    # Classify Ts/Tv (exclude multi-allelic sites)
    snps <- snps[!grepl(",", ALT)]
    transitions <- c("AG", "GA", "CT", "TC")
    snps[, tstv := fifelse(paste0(REF, ALT) %in% transitions, "Ts", "Tv")]

    # Build SNP count + Ts/Tv by annotation co-occurrence
    annotations <- sort(unique(seg_ann$annotation))
    results <- list()
    for (a_off in annotations) {
      segs_off <- seg_ann[annotation == a_off & source_overlap_frac > min_overlap]$augref_path
      for (a_on in annotations) {
        segs_on <- seg_ann[annotation == a_on & ref_overlap_frac > min_overlap]$augref_path
        segs_both <- intersect(segs_off, segs_on)
        sub <- snps[augref_path %in% segs_both]
        results[[length(results) + 1]] <- data.table(
          off_ref = a_off, on_ref = a_on,
          n_snps = nrow(sub),
          ts = sum(sub$tstv == "Ts"), tv = sum(sub$tstv == "Tv")
        )
      }
    }
    result_dt <- rbindlist(results)
    result_dt[, tstv_ratio := fifelse(tv > 0, round(ts / tv, 2), NA_real_)]

    # Factor levels
    result_dt[, off_ref := factor(off_ref, levels = annotations)]
    result_dt[, on_ref := factor(on_ref, levels = annotations)]

    filt_label <- if (!is.null(filt)) toupper(filt) else "ALL"

    # SNP count heatmap
    result_dt[, count_label := comma(n_snps)]
    p_counts <- ggplot(result_dt, aes(x = on_ref, y = off_ref, fill = n_snps)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = count_label), size = 3.2) +
      scale_fill_gradient(low = "white", high = "steelblue",
                          labels = comma, name = "SNP Count") +
      labs(title = title,
           subtitle = paste0("Off-Reference Biallelic SNP Count (", filt_label, ")"),
           x = "On-reference Annotation",
           y = "Off-reference Annotation") +
      coord_fixed() +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    save_png(p_counts, paste0(prefix, "-counts.", filt, ".png"))

    # Ts/Tv ratio heatmap
    result_dt[, tstv_label := fifelse(is.na(tstv_ratio), "NA",
                paste0(sprintf("%.2f", tstv_ratio), "\n(", comma(ts), "/", comma(tv), ")"))]
    p_tstv <- ggplot(result_dt, aes(x = on_ref, y = off_ref, fill = tstv_ratio)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = tstv_label), size = 3.0) +
      scale_fill_gradient(low = "coral", high = "steelblue",
                          name = "Ts/Tv", na.value = "grey90") +
      labs(title = title,
           subtitle = paste0("Off-Reference Biallelic SNP Ts/Tv (", filt_label, ")"),
           x = "On-reference Annotation",
           y = "Off-reference Annotation") +
      coord_fixed() +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    save_png(p_tstv, paste0(prefix, "-tstv.", filt, ".png"))
  }

  # -------------------------------------------------------------------------
  # PCLAI ancestry SNP heatmaps (dominant ancestry per segment)
  # -------------------------------------------------------------------------
  pclai_prefix <- sub("\\.annot-snp$", ".annot-pclai-snp", prefix)
  pclai_dt_vcf <- dt[annotation == "pclai" & annotation_class != "_total"]

  if (nrow(pclai_dt_vcf) == 0 || nrow(snps) == 0) {
    cat("No PCLAI or SNP data; creating empty PCLAI-SNP placeholders.\n")
    for (suf in c("counts", "tstv")) {
      out_file <- paste0(pclai_prefix, "-", suf, ".", filt, ".png")
      grDevices::png(out_file, width = 100, height = 100)
      par(mar = c(0, 0, 0, 0))
      plot.new()
      dev.off()
      cat("Saved placeholder:", out_file, "\n")
    }
  } else {
    # For each segment, find dominant ancestry off-ref and on-ref
    seg_src_anc <- pclai_dt_vcf[source_overlap_frac > min_overlap,
      .(annotation_class = annotation_class[which.max(source_overlap_bp)]),
      by = augref_path]
    setnames(seg_src_anc, "annotation_class", "off_ref_anc")

    seg_ref_anc <- pclai_dt_vcf[ref_overlap_frac > min_overlap,
      .(annotation_class = annotation_class[which.max(ref_overlap_bp)]),
      by = augref_path]
    setnames(seg_ref_anc, "annotation_class", "on_ref_anc")

    anc_pairs <- merge(seg_src_anc, seg_ref_anc, by = "augref_path")

    if (nrow(anc_pairs) == 0) {
      cat("No segments with PCLAI ancestry in both coords; creating empty placeholders.\n")
      for (suf in c("counts", "tstv")) {
        out_file <- paste0(pclai_prefix, "-", suf, ".", filt, ".png")
        grDevices::png(out_file, width = 100, height = 100)
        plot.new()
        dev.off()
        cat("Saved placeholder:", out_file, "\n")
      }
    } else {
      # Join with SNPs, group by (off_ref_anc, on_ref_anc)
      snp_anc <- merge(snps, anc_pairs, by = "augref_path")
      superpops <- sort(unique(c(anc_pairs$off_ref_anc, anc_pairs$on_ref_anc)))

      anc_results <- list()
      for (a_off in superpops) {
        for (a_on in superpops) {
          sub <- snp_anc[off_ref_anc == a_off & on_ref_anc == a_on]
          anc_results[[length(anc_results) + 1]] <- data.table(
            off_ref = a_off, on_ref = a_on,
            n_snps = nrow(sub),
            ts = sum(sub$tstv == "Ts"), tv = sum(sub$tstv == "Tv")
          )
        }
      }
      anc_result_dt <- rbindlist(anc_results)
      anc_result_dt[, tstv_ratio := fifelse(tv > 0, round(ts / tv, 2), NA_real_)]

      anc_result_dt[, off_ref := factor(off_ref, levels = superpops)]
      anc_result_dt[, on_ref := factor(on_ref, levels = superpops)]

      filt_label <- if (!is.null(filt)) toupper(filt) else "ALL"

      # Ancestry SNP count heatmap
      anc_result_dt[, count_label := comma(n_snps)]
      p_anc_counts <- ggplot(anc_result_dt, aes(x = on_ref, y = off_ref, fill = n_snps)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = count_label), size = 3.2) +
        scale_fill_gradient(low = "white", high = "steelblue",
                            labels = comma, name = "SNP Count") +
        labs(title = title,
             subtitle = paste0("Ancestry SNP Count (", filt_label, ")"),
             x = "On-reference Ancestry",
             y = "Off-reference Ancestry") +
        coord_fixed() +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )

      save_png(p_anc_counts, paste0(pclai_prefix, "-counts.", filt, ".png"))

      # Ancestry Ts/Tv ratio heatmap
      anc_result_dt[, tstv_label := fifelse(is.na(tstv_ratio), "NA",
          paste0(sprintf("%.2f", tstv_ratio), "\n(", comma(ts), "/", comma(tv), ")"))]
      p_anc_tstv <- ggplot(anc_result_dt, aes(x = on_ref, y = off_ref, fill = tstv_ratio)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = tstv_label), size = 3.0) +
        scale_fill_gradient(low = "coral", high = "steelblue",
                            name = "Ts/Tv", na.value = "grey90") +
        labs(title = title,
             subtitle = paste0("Ancestry SNP Ts/Tv (", filt_label, ")"),
             x = "On-reference Ancestry",
             y = "Off-reference Ancestry") +
        coord_fixed() +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )

      save_png(p_anc_tstv, paste0(pclai_prefix, "-tstv.", filt, ".png"))
    }
  }

  cat("Done (VCF mode).\n")
  quit(status = 0)
}

# =========================================================================
# Base mode: standard annotation overlap plots
# =========================================================================

# ---------------------------------------------------------------------------
# Plot 1: Summary bar chart — fraction of alt bp overlapping each annotation
# ---------------------------------------------------------------------------

summary_source <- dt_shared[, .(overlap_bp = as.numeric(sum(source_overlap_bp)),
                                total_bp = as.numeric(sum(source_len))),
                            by = .(annotation)]
summary_source[, frac := overlap_bp / total_bp]
summary_source[, coord_type := "Off-reference"]

if (!is.null(genome_cov)) {
  # Use genome-wide annotation coverage for reference bars
  gw <- genome_cov[annotation_class == "_total" | annotation_class == annotation]
  # For grouped annotations use _total; for ungrouped use self-named row
  gw_total <- gw[annotation_class == "_total"]
  gw_ungrouped <- gw[annotation_class == annotation & !annotation %in% gw_total$annotation]
  summary_ref <- rbind(gw_total, gw_ungrouped)[, .(
    overlap_bp = as.numeric(overlap_bp),
    frac = overlap_frac), by = .(annotation)]
  summary_ref[, coord_type := "Genome-wide"]
} else {
  summary_ref <- dt_shared[, .(overlap_bp = as.numeric(sum(ref_overlap_bp)),
                                total_bp = as.numeric(sum(ref_len))),
                            by = .(annotation)]
  summary_ref[, frac := overlap_bp / total_bp]
  summary_ref[, coord_type := "On-reference"]
}
ref_label <- summary_ref$coord_type[1]

bar_dt <- rbind(summary_source[, .(annotation, frac, overlap_bp, coord_type)],
                summary_ref[, .(annotation, frac, overlap_bp, coord_type)])

# Human-readable bp label (e.g., "474 Mb", "21 kb")
format_bp <- function(x) {
  ifelse(x >= 1e9, paste0(round(x / 1e9, 1), " Gb"),
  ifelse(x >= 1e6, paste0(round(x / 1e6, 1), " Mb"),
  ifelse(x >= 1e3, paste0(round(x / 1e3, 1), " kb"),
  paste0(x, " bp"))))
}
bar_dt[, bp_label := format_bp(overlap_bp)]

p1 <- ggplot(bar_dt, aes(x = annotation, y = frac, fill = coord_type)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = bp_label),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = setNames(c("coral", "steelblue"),
                                      c("Off-reference", ref_label)),
                    name = NULL) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.15))) +
  labs(title = title, subtitle = "Fraction of bp Overlapping Annotations",
       x = "Annotation", y = "Overlap Fraction") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

save_png(p1, paste0(prefix, ".annot-summary.png"))

# ---------------------------------------------------------------------------
# Plot 2: Length vs overlap scatter, faceted by annotation
# ---------------------------------------------------------------------------

# Use _total rows for grouped annotations (accurate class-agnostic overlap)
scatter_dt <- dt_shared[, .(source_overlap_bp = sum(source_overlap_bp),
                            source_len = source_len[1],
                            ref_len = ref_len[1]),
                        by = .(augref_path, annotation)]
scatter_dt[, source_overlap_frac := pmin(source_overlap_bp / source_len, 1.0)]
scatter_dt[source_len == 0, source_overlap_frac := 0]

p2 <- ggplot(scatter_dt[source_len > 0],
             aes(x = source_len, y = source_overlap_frac)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "loess", se = FALSE, color = "coral", linewidth = 0.8) +
  facet_wrap(~ annotation, scales = "free_y") +
  scale_x_log10(labels = comma) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(title = title, subtitle = "Source Length vs Overlap Fraction",
       x = "Source Segment Length (bp, log scale)",
       y = "Source Overlap Fraction") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

save_png(p2, paste0(prefix, ".annot-scatter.png"))

# ---------------------------------------------------------------------------
# Plot 3: Repeat class breakdown (only when repeat class data present)
# ---------------------------------------------------------------------------

has_repeat_classes <- any(dt$annotation == "repeats" & !dt$annotation_class %in% c("repeats", "_total"))

if (has_repeat_classes) {
  repeat_dt <- dt[annotation == "repeats" & !annotation_class %in% c("repeats", "_total")]

  # Collapse RepeatMasker class/family (e.g. "Simple_repeat/unknown" → "Simple_repeat")
  repeat_dt[, annotation_class := sub("/.*", "", annotation_class)]

  # Re-aggregate overlap bp after collapsing subclasses
  repeat_dt <- repeat_dt[, .(source_overlap_bp = sum(source_overlap_bp),
                             ref_overlap_bp = sum(ref_overlap_bp),
                             source_len = source_len[1],
                             ref_len = ref_len[1]),
                         by = .(augref_path, annotation, annotation_class)]

  # Total bp per coord type (unique segments to avoid double-counting across classes)
  seg_dt <- unique(repeat_dt[, .(augref_path, source_len, ref_len)])
  total_source_bp <- as.numeric(sum(seg_dt$source_len))
  total_ref_bp <- as.numeric(sum(seg_dt$ref_len))

  # Top 8 classes by combined overlap fraction
  class_totals <- repeat_dt[, .(total_frac =
    as.numeric(sum(source_overlap_bp)) / total_source_bp +
    as.numeric(sum(ref_overlap_bp)) / total_ref_bp),
    by = .(annotation_class)]
  setorder(class_totals, -total_frac)
  top_classes <- head(class_totals$annotation_class, 8)
  repeat_dt[, display_class := ifelse(annotation_class %in% top_classes,
                                      annotation_class, "Other")]

  # Aggregate overlap fractions per class
  src_bar <- repeat_dt[, .(overlap_frac = as.numeric(sum(source_overlap_bp)) / total_source_bp),
                       by = .(display_class)]
  src_bar[, coord_type := "Off-reference"]

  if (!is.null(genome_cov)) {
    # Use genome-wide per-class coverage for reference bars
    gw_rep <- genome_cov[annotation == "repeats" & !annotation_class %in% c("_total", "repeats")]
    gw_rep[, annotation_class := sub("/.*", "", annotation_class)]
    gw_rep <- gw_rep[, .(overlap_bp = sum(overlap_bp), genome_bp = genome_bp[1]),
                      by = .(annotation_class)]
    gw_rep[, display_class := ifelse(annotation_class %in% top_classes,
                                      annotation_class, "Other")]
    ref_bar <- gw_rep[, .(overlap_frac = sum(overlap_bp) / genome_bp[1]),
                       by = .(display_class)]
    ref_bar[, coord_type := "Genome-wide"]
  } else {
    ref_bar <- repeat_dt[, .(overlap_frac = as.numeric(sum(ref_overlap_bp)) / total_ref_bp),
                         by = .(display_class)]
    ref_bar[, coord_type := "On-reference"]
  }
  ref_label_repeat <- ref_bar$coord_type[1]

  repeat_bar_dt <- rbind(src_bar, ref_bar)
  # Order classes by total fraction (off + on), Other last
  class_order <- c(top_classes[top_classes %in% repeat_bar_dt$display_class], "Other")
  repeat_bar_dt[, display_class := factor(display_class, levels = rev(class_order))]

  p3 <- ggplot(repeat_bar_dt, aes(x = display_class, y = overlap_frac, fill = coord_type)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = setNames(c("coral", "steelblue"),
                                        c("Off-reference", ref_label_repeat)),
                      name = NULL) +
    scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.1))) +
    labs(title = title, subtitle = "Repeat Class Breakdown (fraction of segment bp)",
         x = NULL, y = "Overlap Fraction") +
    coord_flip() +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )

  save_png(p3, paste0(prefix, ".annot-repeats.png"))
} else {
  cat("No repeat class data; skipping repeat breakdown plot.\n")
}

# ---------------------------------------------------------------------------
# Plot 4: Annotation co-occurrence heatmap across coord types
# ---------------------------------------------------------------------------

# For each segment, determine which annotations it overlaps in off-ref and on-ref
# Use _total rows for grouped annotations to avoid cross-class double-counting
seg_ann <- dt_shared[, .(source_overlap_bp = sum(source_overlap_bp),
                         ref_overlap_bp = sum(ref_overlap_bp),
                         source_overlap_frac = sum(source_overlap_frac),
                         ref_overlap_frac = sum(ref_overlap_frac)),
                     by = .(augref_path, annotation)]

annotations <- sort(unique(seg_ann$annotation))
all_segs <- unique(seg_ann$augref_path)
n_segs <- length(all_segs)

# Build co-occurrence: for each pair (off-ref ann A, on-ref ann B),
# count segments with overlap fraction above threshold in both
cooccur_list <- list()
for (a_off in annotations) {
  segs_off <- seg_ann[annotation == a_off & source_overlap_frac > min_overlap]$augref_path
  for (a_on in annotations) {
    segs_on <- seg_ann[annotation == a_on & ref_overlap_frac > min_overlap]$augref_path
    n_both <- length(intersect(segs_off, segs_on))
    cooccur_list[[length(cooccur_list) + 1]] <- data.table(
      off_ref = a_off, on_ref = a_on,
      n_segments = n_both,
      frac_segments = n_both / n_segs
    )
  }
}
cooccur_dt <- rbindlist(cooccur_list)

# Factor levels for consistent ordering
cooccur_dt[, off_ref := factor(off_ref, levels = annotations)]
cooccur_dt[, on_ref := factor(on_ref, levels = annotations)]

# Format segment count labels
cooccur_dt[, label := paste0(n_segments, "\n(", sprintf("%.1f%%", 100 * frac_segments), ")")]

p4 <- ggplot(cooccur_dt, aes(x = on_ref, y = off_ref, fill = frac_segments)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3.2) +
  scale_fill_gradient(low = "white", high = "steelblue",
                      labels = percent, name = "Fraction of\nSegments") +
  labs(title = title,
       subtitle = paste0("Annotation Co-occurrence (", n_segs, " segments)"),
       x = "On-reference Annotation",
       y = "Off-reference Annotation") +
  coord_fixed() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save_png(p4, paste0(prefix, ".annot-cooccur.png"))

# ---------------------------------------------------------------------------
# Plot 5: PCLAI ancestry co-occurrence heatmap
# ---------------------------------------------------------------------------

pclai_dt <- dt[annotation == "pclai" & annotation_class != "_total"]
if (nrow(pclai_dt) > 0) {
  # Each row has annotation_class = super-population, with source/ref overlap bp
  # For each segment, determine dominant ancestry in each coord type
  # (the class with the most overlap bp)
  seg_src_anc <- pclai_dt[source_overlap_frac > min_overlap,
    .(annotation_class = annotation_class[which.max(source_overlap_bp)]),
    by = augref_path]
  setnames(seg_src_anc, "annotation_class", "off_ref_anc")

  seg_ref_anc <- pclai_dt[ref_overlap_frac > min_overlap,
    .(annotation_class = annotation_class[which.max(ref_overlap_bp)]),
    by = augref_path]
  setnames(seg_ref_anc, "annotation_class", "on_ref_anc")

  anc_pairs <- merge(seg_src_anc, seg_ref_anc, by = "augref_path")
  n_anc_segs <- nrow(anc_pairs)

  if (n_anc_segs > 0) {
    # Count segments per (off_ref, on_ref) ancestry pair
    anc_cooccur <- anc_pairs[, .N, by = .(off_ref_anc, on_ref_anc)]
    setnames(anc_cooccur, "N", "n_segments")
    anc_cooccur[, frac_segments := n_segments / n_anc_segs]

    # Ensure all pairs present (fill missing with 0)
    superpops <- sort(unique(c(anc_cooccur$off_ref_anc, anc_cooccur$on_ref_anc)))
    all_pairs <- CJ(off_ref_anc = superpops, on_ref_anc = superpops)
    anc_cooccur <- merge(all_pairs, anc_cooccur, by = c("off_ref_anc", "on_ref_anc"), all.x = TRUE)
    anc_cooccur[is.na(n_segments), c("n_segments", "frac_segments") := .(0, 0)]

    anc_cooccur[, off_ref_anc := factor(off_ref_anc, levels = superpops)]
    anc_cooccur[, on_ref_anc := factor(on_ref_anc, levels = superpops)]
    anc_cooccur[, label := paste0(n_segments, "\n(", sprintf("%.1f%%", 100 * frac_segments), ")")]

    p5 <- ggplot(anc_cooccur, aes(x = on_ref_anc, y = off_ref_anc, fill = frac_segments)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = label), size = 3.2) +
      scale_fill_gradient(low = "white", high = "steelblue",
                          labels = percent, name = "Fraction of\nSegments") +
      labs(title = title,
           subtitle = paste0("Local Ancestry Co-occurrence (", n_anc_segs, " segments)"),
           x = "On-reference Ancestry",
           y = "Off-reference Ancestry") +
      coord_fixed() +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    save_png(p5, paste0(prefix, ".annot-ancestry.png"))
  } else {
    cat("No segments with PCLAI ancestry in both coords; skipping ancestry heatmap.\n")
  }
} else {
  cat("No PCLAI data; skipping ancestry heatmap.\n")
}

# ---------------------------------------------------------------------------
# Plot 6: PCLAI ancestry summary bar chart
# ---------------------------------------------------------------------------

pclai_anc_dt <- dt[annotation == "pclai" & annotation_class != "_total"]
if (nrow(pclai_anc_dt) > 0) {
  # Total bp per coord type (unique segments to avoid double-counting across ancestries)
  seg_dt_pclai <- unique(pclai_anc_dt[, .(augref_path, source_len, ref_len)])
  total_source_bp_pclai <- as.numeric(sum(seg_dt_pclai$source_len))
  total_ref_bp_pclai <- as.numeric(sum(seg_dt_pclai$ref_len))

  src_anc_bar <- pclai_anc_dt[, .(overlap_frac = as.numeric(sum(source_overlap_bp)) / total_source_bp_pclai),
                               by = .(annotation_class)]
  src_anc_bar[, coord_type := "Off-reference"]

  ref_anc_bar <- pclai_anc_dt[, .(overlap_frac = as.numeric(sum(ref_overlap_bp)) / total_ref_bp_pclai),
                               by = .(annotation_class)]
  ref_anc_bar[, coord_type := "On-reference"]

  anc_bar_dt <- rbind(src_anc_bar, ref_anc_bar)
  # Order by total fraction (off + on)
  anc_order <- anc_bar_dt[, .(total = sum(overlap_frac)), by = annotation_class]
  setorder(anc_order, -total)
  anc_bar_dt[, annotation_class := factor(annotation_class, levels = rev(anc_order$annotation_class))]

  p6 <- ggplot(anc_bar_dt, aes(x = annotation_class, y = overlap_frac, fill = coord_type)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("Off-reference" = "coral", "On-reference" = "steelblue"),
                      name = NULL) +
    scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.1))) +
    labs(title = title, subtitle = "Local Ancestry Breakdown (fraction of segment bp)",
         x = NULL, y = "Overlap Fraction") +
    coord_flip() +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )

  save_png(p6, paste0(prefix, ".annot-pclai-summary.png"))
} else {
  cat("No PCLAI data; skipping ancestry summary bar chart.\n")
}

# ---------------------------------------------------------------------------
# Output TSV: per-annotation summary
# ---------------------------------------------------------------------------

# Use _total rows for grouped annotations to avoid cross-class double-counting
stats_list <- list()
for (ann in unique(dt_shared$annotation)) {
  sub <- dt_shared[annotation == ann]
  seg_dt <- unique(sub[, .(augref_path, source_len, ref_len)])

  total_source_bp <- sum(seg_dt$source_len)
  total_ref_bp <- sum(seg_dt$ref_len)

  src_overlap <- sum(sub$source_overlap_bp)
  ref_overlap <- sum(sub$ref_overlap_bp)

  seg_with_src <- length(unique(sub[source_overlap_bp > 0]$augref_path))
  seg_with_ref <- length(unique(sub[ref_overlap_bp > 0]$augref_path))
  total_segs <- nrow(seg_dt)

  stats_list[[length(stats_list) + 1]] <- data.table(
    annotation = ann,
    coord_type = "off-reference",
    total_bp = total_source_bp,
    overlap_bp = src_overlap,
    overlap_frac = if (total_source_bp > 0) src_overlap / total_source_bp else 0,
    segments_with_overlap = seg_with_src,
    total_segments = total_segs
  )
  stats_list[[length(stats_list) + 1]] <- data.table(
    annotation = ann,
    coord_type = "on-reference",
    total_bp = total_ref_bp,
    overlap_bp = ref_overlap,
    overlap_frac = if (total_ref_bp > 0) ref_overlap / total_ref_bp else 0,
    segments_with_overlap = seg_with_ref,
    total_segments = total_segs
  )
}

stats_dt <- rbindlist(stats_list)
stats_file <- paste0(prefix, ".annot-stats.tsv")
fwrite(stats_dt, stats_file, sep = "\t")
cat("Wrote summary:", stats_file, "\n")

cat("Done.\n")
