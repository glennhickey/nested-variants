#!/usr/bin/env Rscript

# call-summary-panel.R — Compact per-sample call summary with annotation stacking
#
# Produces a single figure with 3 horizontal stacked bars (mean per sample):
#   Off-ref SNPs | On-ref SNPs | On-ref SVs
# Each bar is stacked by exclusive annotation (genes, repeats, segdups, censat, other).
# SNP bar segments are labelled with the Ts/Tv ratio for that annotation region.
#
# Usage: Rscript scripts/call-summary-panel.R
#          --per-sample <per-sample-types.tsv>
#          --annot <annot-exclusive.tsv>
#          --output <output.png>
#          [--title TITLE]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
ps_path     <- NULL
annot_path  <- NULL
vcf_path    <- NULL
min_sv_size <- 50L
output      <- NULL
title       <- "vg call Summary"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--per-sample" && i + 1 <= length(args)) {
    ps_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--annot" && i + 1 <= length(args)) {
    annot_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--vcf" && i + 1 <= length(args)) {
    vcf_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--min-sv-size" && i + 1 <= length(args)) {
    min_sv_size <- as.integer(args[i + 1]); i <- i + 2
  } else if (args[i] == "--output" && i + 1 <= length(args)) {
    output <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(ps_path) || is.null(output)) {
  cat("Usage: Rscript call-summary-panel.R --per-sample <tsv> --annot <tsv> --output <png>\n")
  quit(status = 1)
}

save_png <- function(plot, file, width = 10, height = 5) {
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
ps <- fread(ps_path)
n_samples <- uniqueN(ps$sample)
cat("Per-sample types:", nrow(ps), "rows,", n_samples, "samples\n")

# Compute mean per-sample counts for three metrics
mean_counts <- ps[, .(mean_count = mean(count)), by = .(ref_context, variant_type)]

# Define bar metrics
offref_snp   <- mean_counts[ref_context == "Off-reference" & variant_type == "SNP", mean_count]
offref_mnp   <- mean_counts[ref_context == "Off-reference" & variant_type == "MNP", mean_count]
offref_indel <- sum(mean_counts[ref_context == "Off-reference" &
                                variant_type %in% c("Insertion", "Deletion"), mean_count])

if (length(offref_snp) == 0)   offref_snp   <- 0
if (length(offref_mnp) == 0)   offref_mnp   <- 0
if (length(offref_indel) == 0) offref_indel <- 0

# On-ref SV count: use VCF with size threshold if --vcf given, else per-sample-types.tsv
if (!is.null(vcf_path)) {
  cat("Computing on-ref SVs from VCF (size >=", min_sv_size, "bp):", vcf_path, "\n")
  # Pre-filter to indels only (skip SNPs/MNPs) and PASS, then read CHROM/REF/ALT/GTs
  sv_cmd <- sprintf(
    "bcftools view -f PASS -v indels '%s' 2>/dev/null | bcftools query -f '%%CHROM\\t%%REF\\t%%ALT[\\t%%GT]\\n' 2>/dev/null",
    vcf_path)
  sv_raw <- fread(cmd = sv_cmd, header = FALSE)
  if (nrow(sv_raw) == 0) {
    onref_sv_ins <- 0; onref_sv_del <- 0
    cat("No indels found in VCF\n")
  } else {
  sv_chrom <- sv_raw[[1]]
  sv_ref   <- sv_raw[[2]]
  sv_alt   <- sv_raw[[3]]
  # On-ref: CHROM does NOT end in _NNN_alt
  is_onref <- !grepl("_[0-9]+_alt$", sv_chrom)
  # Fast path: biallelic (no comma) — compute signed size directly
  is_multi <- grepl(",", sv_alt, fixed = TRUE)
  sv_size_signed <- integer(length(sv_alt))
  sv_size_signed[!is_multi] <- nchar(sv_alt[!is_multi]) - nchar(sv_ref[!is_multi])
  # Slow path: multi-allelic — use largest allele
  multi_idx <- which(is_multi)
  if (length(multi_idx) > 0) {
    sv_size_signed[multi_idx] <- vapply(multi_idx, function(i) {
      alts <- unlist(strsplit(sv_alt[i], ","))
      alts <- alts[!is.na(alts) & nzchar(alts) & alts != "*" & alts != "."]
      if (length(alts) == 0) return(0L)
      diffs <- nchar(alts) - nchar(sv_ref[i])
      diffs[which.max(abs(diffs))]
    }, integer(1))
  }
  sv_size <- abs(sv_size_signed)
  is_sv <- is_onref & sv_size >= min_sv_size
  is_sv_ins <- is_sv & sv_size_signed > 0
  is_sv_del <- is_sv & sv_size_signed < 0
  cat("On-ref SVs >=", min_sv_size, "bp:", sum(is_sv), "(", sum(is_sv_ins), "ins,", sum(is_sv_del), "del )\n")
  gt_cols <- 4:ncol(sv_raw)
  if (sum(is_sv_ins) > 0) {
    onref_sv_ins <- mean(vapply(gt_cols, function(col) sum(grepl("[1-9]", sv_raw[[col]][is_sv_ins])), numeric(1)))
  } else { onref_sv_ins <- 0 }
  if (sum(is_sv_del) > 0) {
    onref_sv_del <- mean(vapply(gt_cols, function(col) sum(grepl("[1-9]", sv_raw[[col]][is_sv_del])), numeric(1)))
  } else { onref_sv_del <- 0 }
  cat("Mean per sample — On-ref SV Ins:", round(onref_sv_ins), "  On-ref SV Del:", round(onref_sv_del), "\n")
  } # end else (sv_raw has rows)
} else {
  onref_sv_ins <- sum(mean_counts[ref_context == "On-reference" & variant_type == "SV Insertion", mean_count])
  onref_sv_del <- sum(mean_counts[ref_context == "On-reference" & variant_type == "SV Deletion", mean_count])
}

cat("Mean per sample — Off-ref SNPs:", round(offref_snp),
    "  Off-ref MNPs:", round(offref_mnp),
    "  Off-ref Indels:", round(offref_indel),
    "  On-ref SV Ins:", round(onref_sv_ins),
    "  On-ref SV Del:", round(onref_sv_del), "\n")

# ---------------------------------------------------------------------------
# Build bar data with annotation stacking
# ---------------------------------------------------------------------------
annot_colors <- c("genes" = "forestgreen", "repeats" = "orange",
                  "segdups" = "firebrick", "censat" = "mediumpurple",
                  "Other" = "grey70")
annot_order  <- c("genes", "repeats", "segdups", "censat", "Other")

# Define bar categories
bar_cats <- c("Off-ref SNPs", "Off-ref MNPs", "Off-ref Indels", "On-ref SV Ins", "On-ref SV Del")

if (!is.null(annot_path) && file.exists(annot_path)) {
  annot <- fread(annot_path)
  cat("Annotation-exclusive:", nrow(annot), "rows\n")

  # Map variant_type to our bar categories
  annot[, bar_cat := fifelse(
    ref_context == "Off-reference" & variant_type == "SNP", "Off-ref SNPs",
    fifelse(ref_context == "Off-reference" & variant_type == "MNP", "Off-ref MNPs",
    fifelse(ref_context == "Off-reference" & variant_type %in% c("Insertion", "Deletion"), "Off-ref Indels",
    fifelse(ref_context == "On-reference" & variant_type == "SV Insertion", "On-ref SV Ins",
    fifelse(ref_context == "On-reference" & variant_type == "SV Deletion", "On-ref SV Del",
    NA_character_)))))]

  annot <- annot[!is.na(bar_cat)]

  # Aggregate by bar_cat x annotation
  annot_agg <- annot[, .(count = sum(count),
                         ts = sum(ts, na.rm = TRUE),
                         tv = sum(tv, na.rm = TRUE)),
                     by = .(bar_cat, annot_exclusive)]

  # Compute proportion within each bar_cat
  annot_agg[, total := sum(count), by = bar_cat]
  annot_agg[, prop := count / total]

  # Mean per-sample count for each bar_cat
  bar_totals <- data.table(
    bar_cat = bar_cats,
    mean_total = c(offref_snp, offref_mnp, offref_indel, onref_sv_ins, onref_sv_del))

  annot_agg <- merge(annot_agg, bar_totals, by = "bar_cat")
  annot_agg[, bar_value := prop * mean_total]

  # Ts/Tv ratio per segment
  annot_agg[tv > 0, tstv_ratio := round(ts / tv, 2)]

  # Ensure all annotations present
  annot_agg[, annot_exclusive := factor(annot_exclusive,
    levels = intersect(annot_order, unique(annot_agg$annot_exclusive)))]

  plot_dt <- annot_agg[, .(bar_cat, annotation = annot_exclusive,
                            bar_value, tstv_ratio)]
} else {
  # No annotation data — single "Other" segment per bar
  plot_dt <- data.table(
    bar_cat = bar_cats,
    annotation = factor("Other", levels = "Other"),
    bar_value = c(offref_snp, offref_mnp, offref_indel, onref_sv_ins, onref_sv_del),
    tstv_ratio = NA_real_)
}

# Filter out bars with 0 count
bar_keep <- plot_dt[, .(total = sum(bar_value)), by = bar_cat][total > 0, bar_cat]
plot_dt <- plot_dt[bar_cat %in% bar_keep]

plot_dt[, bar_cat := factor(bar_cat, levels = rev(intersect(bar_cats, bar_keep)))]

# ---------------------------------------------------------------------------
# Build the plot
# ---------------------------------------------------------------------------

# Label for each segment: Ts/Tv only for SNP bars
plot_dt[, label := fifelse(
  !is.na(tstv_ratio) & bar_value > 0 & grepl("SNP", bar_cat),
  paste0("Ts/Tv=", tstv_ratio),
  "")]

# Suppress labels for very small segments (< 5% of bar total)
bar_totals_dt <- plot_dt[, .(total = sum(bar_value)), by = bar_cat]
plot_dt <- merge(plot_dt, bar_totals_dt, by = "bar_cat", suffixes = c("", "_total"))
plot_dt[bar_value / total < 0.05, label := ""]

p <- ggplot(plot_dt, aes(x = bar_cat, y = bar_value, fill = annotation)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 2.8, lineheight = 0.85, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = annot_colors,
                    name = "Annotation") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  labs(title = title,
       subtitle = paste0("Mean per sample (N=", n_samples, ", PASS)"),
       x = NULL, y = "Mean Variant Count per Sample") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "bottom"
  )

save_png(p, output, width = 10, height = 4)
cat("Done.\n")
