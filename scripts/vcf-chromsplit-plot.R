#!/usr/bin/env Rscript

# vcf-chromsplit-plot.R — Per-contig concordance/discordance plots
#
# Usage: Rscript scripts/vcf-chromsplit-plot.R <chromsplit.tsv> <output_prefix>
#          [--title TITLE] [--strip-prefix PREFIX]
#          [--annot annot-per-segment.tsv] [--segs augref-segs.tsv]
#          [--giab-beds bed1,bed2,...] [--giab-names name1,name2,...]
#
# Reads the merged long-format chromsplit TSV (from vcfeval_chromsplit_merge)
# and produces:
#   {prefix}.chromsplit.png           — FP vs FN scatter per contig
#   {prefix}.chromsplit-top.png       — Top 30 most discordant off-ref contigs
#   {prefix}.chromsplit-concordant.png — Top 30 most concordant off-ref contigs
#   {prefix}.chromsplit-annot.png     — SNP TP/FP/FN by annotation (off-ref, when --annot)
#   {prefix}.chromsplit-giab.png      — SNP TP/FP/FN by GIAB region (off-ref, when --giab-beds)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript vcf-chromsplit-plot.R <chromsplit.tsv> <output_prefix> [--title TITLE] [--strip-prefix PREFIX]\n")
  quit(status = 1)
}

input_tsv    <- args[1]
prefix       <- args[2]
title        <- NULL
strip_prefix <- NULL
annot_file   <- NULL
segs_file    <- NULL
giab_beds_arg  <- NULL
giab_names_arg <- NULL

i <- 3
while (i <= length(args)) {
  if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--strip-prefix" && i + 1 <= length(args)) {
    strip_prefix <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--annot" && i + 1 <= length(args)) {
    annot_file <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--segs" && i + 1 <= length(args)) {
    segs_file <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--giab-beds" && i + 1 <= length(args)) {
    giab_beds_arg <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--giab-names" && i + 1 <= length(args)) {
    giab_names_arg <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(title)) title <- "Per-Contig FP/FN Breakdown"

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
# Read data
# ---------------------------------------------------------------------------
cat("Reading:", input_tsv, "\n")
dt <- data.table::fread(input_tsv)

if (nrow(dt) == 0) {
  cat("No data — creating empty outputs.\n")
  file.create(paste0(prefix, ".chromsplit.png"))
  file.create(paste0(prefix, ".chromsplit-top.png"))
  file.create(paste0(prefix, ".chromsplit-concordant.png"))
  if (!is.null(annot_file)) file.create(paste0(prefix, ".chromsplit-annot.png"))
  if (!is.null(giab_beds_arg)) file.create(paste0(prefix, ".chromsplit-giab.png"))
  quit(status = 0)
}

n_contigs <- length(unique(dt$contig))
n_samples <- length(unique(dt$sample))
cat("Contigs:", n_contigs, "  Samples:", n_samples, "\n")

# Strip prefix from contig names for readability
if (!is.null(strip_prefix) && nzchar(strip_prefix)) {
  dt[, contig := sub(paste0("^", strip_prefix), "", contig, fixed = FALSE)]
}

# Combine SV columns into SNP/Indel if present
has_sv <- "SV_SNP_FP" %in% names(dt)
if (has_sv) {
  dt[, SNP_FP := SNP_FP + SV_SNP_FP]
  dt[, SNP_FN := SNP_FN + SV_SNP_FN]
  dt[, INDEL_FP := INDEL_FP + SV_INDEL_FP]
  dt[, INDEL_FN := INDEL_FN + SV_INDEL_FN]
}

# TP columns (may not exist in older TSVs)
has_tp <- "SNP_TP" %in% names(dt)
if (has_tp && has_sv) {
  dt[, SNP_TP := SNP_TP + SV_SNP_TP]
  dt[, INDEL_TP := INDEL_TP + SV_INDEL_TP]
}
if (!has_tp) {
  dt[, SNP_TP := 0L]
  dt[, INDEL_TP := 0L]
}
dt[, total_tp := SNP_TP + INDEL_TP]

base_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Color/label palettes
all_colors <- c("SNP_TP" = "#4DAF4A", "INDEL_TP" = "#66C2A5",
                "SNP_FP" = "#E41A1C", "SNP_FN" = "#377EB8",
                "INDEL_FP" = "#FF7F00", "INDEL_FN" = "#984EA3")
all_labels <- c("SNP_TP" = "SNP TP", "INDEL_TP" = "Indel TP",
                "SNP_FP" = "SNP FP", "SNP_FN" = "SNP FN",
                "INDEL_FP" = "Indel FP", "INDEL_FN" = "Indel FN")
snp_colors <- c("SNP_TP" = "#4DAF4A", "SNP_FP" = "#E41A1C", "SNP_FN" = "#377EB8")
snp_labels <- c("SNP_TP" = "SNP TP", "SNP_FP" = "SNP FP", "SNP_FN" = "SNP FN")

# ---------------------------------------------------------------------------
# Plot 1: FP vs FN scatter, faceted by variant type
# ---------------------------------------------------------------------------
snp_dt <- dt[, .(sample, contig, FP = SNP_FP, FN = SNP_FN)]
snp_dt[, type := "SNP"]
indel_dt <- dt[, .(sample, contig, FP = INDEL_FP, FN = INDEL_FN)]
indel_dt[, type := "Indel"]
scatter_dt <- rbind(snp_dt, indel_dt)

p_scatter <- ggplot(scatter_dt, aes(x = FP, y = FN)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", alpha = 0.5) +
  facet_wrap(~ type, scales = "free") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(title = title,
       subtitle = paste0("Each point = one contig (", n_contigs, " contigs",
                         if (n_samples > 1) paste0(", ", n_samples, " samples") else "",
                         ")"),
       x = "False Positives", y = "False Negatives") +
  base_theme

if (n_samples > 1) {
  p_scatter <- p_scatter +
    geom_point(aes(color = sample), alpha = 0.5, size = 1.8) +
    scale_color_brewer(palette = "Set1", name = "Sample")
} else {
  p_scatter <- p_scatter +
    geom_point(alpha = 0.5, size = 1.8, color = "steelblue")
}

save_png(p_scatter, paste0(prefix, ".chromsplit.png"), width = 10, height = 5)
rm(snp_dt, indel_dt, scatter_dt)

# ---------------------------------------------------------------------------
# Identify off-ref contigs (requires --segs)
# ---------------------------------------------------------------------------
offref_contigs <- NULL
if (!is.null(segs_file)) {
  cat("Reading augref segments:", segs_file, "\n")
  segs <- fread(segs_file, select = c(4, 5, 6, 7),
                col.names = c("augref_path", "ref_path", "ref_start", "ref_end"))
  segs <- unique(segs, by = "augref_path")

  if (!is.null(strip_prefix) && nzchar(strip_prefix)) {
    segs[, contig := sub(paste0("^", strip_prefix), "", augref_path)]
  } else {
    segs[, contig := augref_path]
  }

  offref_contigs <- unique(segs$contig)
  cat("Off-ref contigs in chromsplit:", length(intersect(unique(dt$contig), offref_contigs)),
      "of", n_contigs, "\n")
}

# ---------------------------------------------------------------------------
# Plot 2 & 3: Top discordant / concordant off-ref contigs (stacked bar)
# ---------------------------------------------------------------------------
# Use off-ref subset if segs available, otherwise all contigs
dt_bar <- if (!is.null(offref_contigs)) dt[contig %in% offref_contigs] else dt
offref_label <- if (!is.null(offref_contigs)) " off-ref" else ""

agg <- dt_bar[, .(SNP_FP = sum(SNP_FP), SNP_FN = sum(SNP_FN),
                   INDEL_FP = sum(INDEL_FP), INDEL_FN = sum(INDEL_FN),
                   SNP_TP = sum(SNP_TP), INDEL_TP = sum(INDEL_TP),
                   total_errors = sum(total_errors), total_tp = sum(total_tp)),
              by = contig]

bar_measures <- c("SNP_TP", "INDEL_TP", "SNP_FP", "SNP_FN", "INDEL_FP", "INDEL_FN")

# --- Discordant ---
setorder(agg, -total_errors)
n_show <- min(30, nrow(agg))
top_disc <- agg[1:n_show]

bar_disc <- melt(top_disc, id.vars = "contig",
                 measure.vars = bar_measures,
                 variable.name = "error_type", value.name = "count")
bar_disc[, contig := factor(contig, levels = rev(top_disc$contig))]

p_disc <- ggplot(bar_disc, aes(x = contig, y = count, fill = error_type)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = all_colors, labels = all_labels, name = NULL) +
  scale_y_continuous(labels = comma) +
  labs(title = title,
       subtitle = paste0("Top ", n_show, " most discordant", offref_label, " contigs",
                         if (n_samples > 1) " (summed across samples)" else ""),
       x = NULL, y = "Variant Count") +
  base_theme

bar_height <- max(6, n_show * 0.25)
save_png(p_disc, paste0(prefix, ".chromsplit-top.png"), width = 10, height = bar_height)

# --- Concordant ---
setorder(agg, -total_tp)
top_conc <- agg[total_tp > 0][1:min(30, sum(agg$total_tp > 0))]

if (nrow(top_conc) > 0) {
  n_show_c <- nrow(top_conc)
  bar_conc <- melt(top_conc, id.vars = "contig",
                   measure.vars = bar_measures,
                   variable.name = "error_type", value.name = "count")
  bar_conc[, contig := factor(contig, levels = rev(top_conc$contig))]

  p_conc <- ggplot(bar_conc, aes(x = contig, y = count, fill = error_type)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = all_colors, labels = all_labels, name = NULL) +
    scale_y_continuous(labels = comma) +
    labs(title = title,
         subtitle = paste0("Top ", n_show_c, " most concordant", offref_label, " contigs",
                           if (n_samples > 1) " (summed across samples)" else ""),
         x = NULL, y = "Variant Count") +
    base_theme

  bar_height_c <- max(6, n_show_c * 0.25)
  save_png(p_conc, paste0(prefix, ".chromsplit-concordant.png"), width = 10, height = bar_height_c)
} else {
  cat("No TP data — creating empty concordant plot.\n")
  file.create(paste0(prefix, ".chromsplit-concordant.png"))
}
rm(agg, dt_bar)

# ---------------------------------------------------------------------------
# Plot 4 & 5: Top discordant / concordant on-ref contigs (stacked bar)
# ---------------------------------------------------------------------------
if (!is.null(offref_contigs)) {
  dt_onref <- dt[!contig %in% offref_contigs]
  cat("On-ref contigs:", uniqueN(dt_onref$contig), "\n")

  agg_on <- dt_onref[, .(SNP_FP = sum(SNP_FP), SNP_FN = sum(SNP_FN),
                          INDEL_FP = sum(INDEL_FP), INDEL_FN = sum(INDEL_FN),
                          SNP_TP = sum(SNP_TP), INDEL_TP = sum(INDEL_TP),
                          total_errors = sum(total_errors), total_tp = sum(total_tp)),
                     by = contig]

  # --- On-ref Discordant ---
  setorder(agg_on, -total_errors)
  n_show_on <- min(30, nrow(agg_on))
  top_disc_on <- agg_on[1:n_show_on]

  bar_disc_on <- melt(top_disc_on, id.vars = "contig",
                       measure.vars = bar_measures,
                       variable.name = "error_type", value.name = "count")
  bar_disc_on[, contig := factor(contig, levels = rev(top_disc_on$contig))]

  p_disc_on <- ggplot(bar_disc_on, aes(x = contig, y = count, fill = error_type)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = all_colors, labels = all_labels, name = NULL) +
    scale_y_continuous(labels = comma) +
    labs(title = title,
         subtitle = paste0("Top ", n_show_on, " most discordant on-ref contigs",
                           if (n_samples > 1) " (summed across samples)" else ""),
         x = NULL, y = "Variant Count") +
    base_theme

  bar_height_on <- max(6, n_show_on * 0.25)
  save_png(p_disc_on, paste0(prefix, ".chromsplit-top-onref.png"), width = 10, height = bar_height_on)

  # --- On-ref Concordant ---
  setorder(agg_on, -total_tp)
  top_conc_on <- agg_on[total_tp > 0][1:min(30, sum(agg_on$total_tp > 0))]

  if (nrow(top_conc_on) > 0) {
    n_show_on_c <- nrow(top_conc_on)
    bar_conc_on <- melt(top_conc_on, id.vars = "contig",
                         measure.vars = bar_measures,
                         variable.name = "error_type", value.name = "count")
    bar_conc_on[, contig := factor(contig, levels = rev(top_conc_on$contig))]

    p_conc_on <- ggplot(bar_conc_on, aes(x = contig, y = count, fill = error_type)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = all_colors, labels = all_labels, name = NULL) +
      scale_y_continuous(labels = comma) +
      labs(title = title,
           subtitle = paste0("Top ", n_show_on_c, " most concordant on-ref contigs",
                             if (n_samples > 1) " (summed across samples)" else ""),
           x = NULL, y = "Variant Count") +
      base_theme

    bar_height_on_c <- max(6, n_show_on_c * 0.25)
    save_png(p_conc_on, paste0(prefix, ".chromsplit-concordant-onref.png"), width = 10, height = bar_height_on_c)
  } else {
    file.create(paste0(prefix, ".chromsplit-concordant-onref.png"))
  }
  rm(agg_on, dt_onref)
} else {
  file.create(paste0(prefix, ".chromsplit-top-onref.png"))
  file.create(paste0(prefix, ".chromsplit-concordant-onref.png"))
}

# ---------------------------------------------------------------------------
# Annotation / GIAB stratification (off-ref contigs only, SNPs only)
# ---------------------------------------------------------------------------
if (!is.null(segs_file)) {
  dt_offref <- dt[contig %in% offref_contigs]

  # --- Annotation-stratified TP/FP/FN bars ---
  if (!is.null(annot_file)) {
    cat("Reading annotations:", annot_file, "\n")
    annot <- fread(annot_file)

    # Use _total aggregated rows for grouped annotations (e.g., repeats)
    annot_agg <- annot[annotation_class == "_total" |
                       !annotation %in% annot[annotation_class == "_total"]$annotation]

    # Filter out pclai
    annot_agg <- annot_agg[annotation != "pclai"]

    # Classify by reference-position overlap
    contig_annots <- annot_agg[ref_overlap_frac > 0, .(augref_path, annotation)]
    contig_annots <- unique(contig_annots, by = c("augref_path", "annotation"))

    # Strip prefix to match chromsplit contig names
    if (!is.null(strip_prefix) && nzchar(strip_prefix)) {
      contig_annots[, contig := sub(paste0("^", strip_prefix), "", augref_path)]
    } else {
      contig_annots[, contig := augref_path]
    }
    contig_annots[, augref_path := NULL]

    dt_annot <- merge(dt_offref, contig_annots, by = "contig", allow.cartesian = TRUE)
    rm(annot, annot_agg, contig_annots)

    if (nrow(dt_annot) > 0) {
      annot_sum <- dt_annot[, .(SNP_TP = sum(SNP_TP), SNP_FP = sum(SNP_FP),
                                SNP_FN = sum(SNP_FN)),
                            by = annotation]
      rm(dt_annot)
      annot_bar <- melt(annot_sum, id.vars = "annotation",
                        measure.vars = c("SNP_TP", "SNP_FP", "SNP_FN"),
                        variable.name = "error_type", value.name = "count")

      p_annot <- ggplot(annot_bar, aes(x = annotation, y = count, fill = error_type)) +
        geom_col(position = "dodge", width = 0.7) +
        scale_fill_manual(values = snp_colors, labels = snp_labels, name = NULL) +
        scale_y_continuous(labels = comma) +
        labs(title = title,
             subtitle = paste0("SNP TP/FP/FN by Annotation (off-ref contigs",
                               if (n_samples > 1) paste0(", ", n_samples, " samples") else "",
                               ")"),
             x = "Annotation", y = "SNP Count") +
        base_theme

      save_png(p_annot, paste0(prefix, ".chromsplit-annot.png"), width = 8, height = 6)
    } else {
      cat("No annotation matches — creating empty annot plot.\n")
      file.create(paste0(prefix, ".chromsplit-annot.png"))
    }
  }

  # --- GIAB-stratified TP/FP/FN bars ---
  if (!is.null(giab_beds_arg)) {
    giab_bed_files <- strsplit(giab_beds_arg, ",")[[1]]
    giab_names <- if (!is.null(giab_names_arg)) {
      gsub("_", " ", strsplit(giab_names_arg, ",")[[1]])
    } else {
      paste0("Region", seq_along(giab_bed_files))
    }

    cat("GIAB stratification from", length(giab_bed_files), "BEDs\n")

    # Write reference-space BED for contigs in chromsplit data only
    segs_used <- segs[contig %in% unique(dt$contig)]
    cat("GIAB BED contigs:", nrow(segs_used), "of", nrow(segs), "segs\n")
    tmp_bed <- tempfile(fileext = ".bed")
    fwrite(segs_used[, .(ref_path, ref_start, ref_end, contig)],
           tmp_bed, sep = "\t", col.names = FALSE)
    rm(segs_used)
    tmp_sorted <- tempfile(fileext = ".sorted.bed")
    system(sprintf("LC_ALL=C sort -k1,1 -k2,2n '%s' > '%s'", tmp_bed, tmp_sorted))
    unlink(tmp_bed)

    # bedtools intersect each GIAB BED
    contig_giab <- data.table()
    for (k in seq_along(giab_bed_files)) {
      cmd <- sprintf("bedtools intersect -a '%s' -b '%s' -wa -u",
                     tmp_sorted, giab_bed_files[k])
      hits <- tryCatch(
        fread(cmd = cmd, header = FALSE, select = 4, col.names = "contig"),
        error = function(e) data.table(contig = character(0)),
        warning = function(w) data.table(contig = character(0))
      )
      if (nrow(hits) > 0) {
        hits <- unique(hits)
        hits[, giab_region := giab_names[k]]
        contig_giab <- rbind(contig_giab, hits)
      }
    }
    unlink(tmp_sorted)

    if (nrow(contig_giab) > 0) {
      dt_giab <- merge(dt_offref, contig_giab, by = "contig", allow.cartesian = TRUE)

      giab_sum <- dt_giab[, .(SNP_TP = sum(SNP_TP), SNP_FP = sum(SNP_FP),
                              SNP_FN = sum(SNP_FN)),
                          by = giab_region]
      giab_bar <- melt(giab_sum, id.vars = "giab_region",
                       measure.vars = c("SNP_TP", "SNP_FP", "SNP_FN"),
                       variable.name = "error_type", value.name = "count")
      giab_bar[, giab_region := factor(giab_region, levels = giab_names)]

      p_giab <- ggplot(giab_bar, aes(x = giab_region, y = count, fill = error_type)) +
        geom_col(position = "dodge", width = 0.7) +
        scale_fill_manual(values = snp_colors, labels = snp_labels, name = NULL) +
        scale_y_continuous(labels = comma) +
        labs(title = title,
             subtitle = paste0("SNP TP/FP/FN by GIAB Region (off-ref contigs",
                               if (n_samples > 1) paste0(", ", n_samples, " samples") else "",
                               ")"),
             x = "GIAB Region", y = "SNP Count") +
        base_theme

      save_png(p_giab, paste0(prefix, ".chromsplit-giab.png"), width = 8, height = 6)
    } else {
      cat("No GIAB hits — creating empty GIAB plot.\n")
      file.create(paste0(prefix, ".chromsplit-giab.png"))
    }
  }
}

cat("Done.\n")
