#!/usr/bin/env Rscript

# pantree-compare.R — Side-by-side comparison of two variant catalogs
#
# Usage: Rscript scripts/pantree-compare.R
#          --ours <records.tsv> --pantree <records.tsv> --prefix <path>
#          [--title TITLE] [--ours-label LABEL] [--pantree-label LABEL]
#
# Both TSVs must have columns:
#   CHROM POS ref_context variant_type size size_signed nonref_af is_repeat
#
# Outputs:
#   {prefix}.pantree-types.png      — variant type counts side-by-side
#   {prefix}.pantree-types-pct.png  — same but percentages
#   {prefix}.pantree-size-dist.png  — indel size distribution overlay
#   {prefix}.pantree-af.png         — AF spectrum overlay
#   {prefix}.pantree-compare.tsv    — summary comparison table

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
ours_path    <- NULL
pantree_path <- NULL
prefix       <- NULL
title        <- "Variant Catalog Comparison"
ours_label   <- "Deconstruct"
pantree_label <- "Pantree"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--ours" && i + 1 <= length(args)) {
    ours_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--pantree" && i + 1 <= length(args)) {
    pantree_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--prefix" && i + 1 <= length(args)) {
    prefix <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--ours-label" && i + 1 <= length(args)) {
    ours_label <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--pantree-label" && i + 1 <= length(args)) {
    pantree_label <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(ours_path) || is.null(pantree_path) || is.null(prefix)) {
  cat("Usage: Rscript pantree-compare.R --ours <tsv> --pantree <tsv> --prefix <path>\n")
  quit(status = 1)
}

save_png <- function(plot, file, width = 10, height = 7) {
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
cat("Reading:", ours_path, "\n")
dt_ours <- fread(ours_path)
dt_ours[, source := ours_label]

cat("Reading:", pantree_path, "\n")
dt_pt <- fread(pantree_path)
dt_pt[, source := pantree_label]

# Ensure is_repeat is logical
for (d in list(dt_ours, dt_pt)) {
  if ("is_repeat" %in% names(d) && is.character(d$is_repeat)) {
    d[, is_repeat := (is_repeat == "TRUE")]
  }
}

# Normalize CHROM to base chromosome for fair comparison.
# Our deconstruct records have augref path names like "augref_CHM13#0#chr1" or
# "augref_CHM13#0#chr10_10001_alt"; pantree has plain "chr1", "chr10".
# Strip augref prefix and _NNN_alt suffix to get the base chromosome.
normalize_chrom <- function(x) {
  x <- sub("^augref_[^#]+#[^#]+#", "", x)   # augref_REF#HAP#chrom → chrom
  x <- sub("_[0-9]+_alt$", "", x)            # chr10_10001_alt → chr10
  x
}

dt_ours[, base_chrom := normalize_chrom(CHROM)]
dt_pt[, base_chrom := normalize_chrom(CHROM)]

# Restrict to shared base chromosomes so the comparison is fair
# (e.g. pantree may exclude chrY)
shared_chroms <- intersect(unique(dt_ours$base_chrom), unique(dt_pt$base_chrom))
n_ours_chroms <- uniqueN(dt_ours$base_chrom)
n_pt_chroms   <- uniqueN(dt_pt$base_chrom)

if (length(shared_chroms) < n_ours_chroms) {
  dropped <- setdiff(unique(dt_ours$base_chrom), shared_chroms)
  cat("Dropping", length(dropped), "base chromosomes from", ours_label,
      "not in", pantree_label, ":", paste(dropped, collapse = ", "), "\n")
  dt_ours <- dt_ours[base_chrom %in% shared_chroms]
}
if (length(shared_chroms) < n_pt_chroms) {
  dropped <- setdiff(unique(dt_pt$base_chrom), shared_chroms)
  cat("Dropping", length(dropped), "base chromosomes from", pantree_label,
      "not in", ours_label, ":", paste(dropped, collapse = ", "), "\n")
  dt_pt <- dt_pt[base_chrom %in% shared_chroms]
}

cat(ours_label, ":", nrow(dt_ours), "records\n")
cat(pantree_label, ":", nrow(dt_pt), "records\n")

# Variant type levels
type_levels <- c("SNP", "MNP", "Insertion", "Deletion", "SV Insertion", "SV Deletion", "Other")

# ---------------------------------------------------------------------------
# Plot 1: Variant type counts (side-by-side)
# ---------------------------------------------------------------------------
counts_ours <- dt_ours[, .(count = .N), by = .(source, ref_context, variant_type)]
counts_pt   <- dt_pt[, .(count = .N), by = .(source, ref_context, variant_type)]
counts_all  <- rbind(counts_ours, counts_pt)
counts_all[, variant_type := factor(variant_type, levels = intersect(type_levels, unique(variant_type)))]

fill_breaks <- c(paste0(ours_label, ".On-reference"), paste0(ours_label, ".Off-reference"),
                 paste0(pantree_label, ".On-reference"), paste0(pantree_label, ".Off-reference"))
fill_values <- setNames(
  c("steelblue", "coral", "dodgerblue3", "tomato3"), fill_breaks)
fill_labels <- setNames(
  c(paste(ours_label, "On-ref"), paste(ours_label, "Off-ref"),
    paste(pantree_label, "On-ref"), paste(pantree_label, "Off-ref")), fill_breaks)

p1 <- ggplot(counts_all, aes(x = variant_type, y = count, fill = interaction(source, ref_context))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    breaks = fill_breaks,
    values = fill_values,
    labels = fill_labels,
    name = NULL
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = title, subtitle = "Variant Type Counts",
       x = "Variant Type", y = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
save_png(p1, paste0(prefix, ".pantree-types.png"))

# ---------------------------------------------------------------------------
# Plot 2: Variant type percentages
# ---------------------------------------------------------------------------
counts_all[, total := sum(count), by = source]
counts_all[, pct := 100 * count / total]

p2 <- ggplot(counts_all, aes(x = variant_type, y = pct, fill = interaction(source, ref_context))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    breaks = fill_breaks,
    values = fill_values,
    labels = fill_labels,
    name = NULL
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = title, subtitle = "Variant Type Proportions",
       x = "Variant Type", y = "Percent of Total") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
save_png(p2, paste0(prefix, ".pantree-types-pct.png"))

# ---------------------------------------------------------------------------
# Plot 3: Size distribution overlay
# ---------------------------------------------------------------------------
indel_types <- c("Insertion", "Deletion", "SV Insertion", "SV Deletion")
dt_indels <- rbind(
  dt_ours[variant_type %in% indel_types, .(source, ref_context, size, size_signed, variant_type)],
  dt_pt[variant_type %in% indel_types, .(source, ref_context, size, size_signed, variant_type)]
)

if (nrow(dt_indels) > 0) {
  dt_indels[, direction := fifelse(size_signed > 0, "Insertion", "Deletion")]

  # Small indels (1-49bp)
  dt_small <- dt_indels[size > 0 & size < 50,
                        .(count = .N), by = .(source, ref_context, size, direction)]
  # SVs (50-1000bp)
  dt_sv <- dt_indels[size >= 50 & size <= 1000,
                     .(count = .N), by = .(source, ref_context, size, direction)]

  dt_small[, size_panel := "Small Indels (1-49 bp)"]
  dt_sv[, size_panel := "Structural Variants (50-1000 bp)"]
  dt_size <- rbind(dt_small, dt_sv)

  if (nrow(dt_size) > 0) {
    p3 <- ggplot(dt_size, aes(x = size, y = count, color = source, linetype = direction)) +
      geom_line(linewidth = 0.7) +
      scale_color_manual(values = setNames(c("steelblue", "tomato3"), c(ours_label, pantree_label))) +
      scale_y_continuous(labels = scales::comma) +
      facet_grid(ref_context ~ size_panel, scales = "free") +
      labs(title = title, subtitle = "Indel Size Distribution",
           x = "Size (bp)", y = "Count") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA)
      )
    save_png(p3, paste0(prefix, ".pantree-size-dist.png"), width = 12, height = 8)
  } else {
    file.create(paste0(prefix, ".pantree-size-dist.png"))
  }
} else {
  file.create(paste0(prefix, ".pantree-size-dist.png"))
}

# ---------------------------------------------------------------------------
# Plot 4: AF spectrum overlay
# ---------------------------------------------------------------------------
has_af_ours <- "nonref_af" %in% names(dt_ours) && !all(is.na(dt_ours$nonref_af))
has_af_pt   <- "nonref_af" %in% names(dt_pt)   && !all(is.na(dt_pt$nonref_af))

if (has_af_ours && has_af_pt) {
  af_step <- 0.05
  dt_af <- rbind(
    dt_ours[!is.na(nonref_af), .(source, nonref_af, ref_context)],
    dt_pt[!is.na(nonref_af), .(source, nonref_af, ref_context)]
  )
  dt_af[, af_bin := round(nonref_af / af_step) * af_step]
  dt_af[, af_bin := pmin(af_bin, 1.0)]
  af_counts <- dt_af[, .(count = .N), by = .(source, ref_context, af_bin)]

  p4 <- ggplot(af_counts, aes(x = af_bin, y = count, color = source, linetype = ref_context)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.5) +
    scale_color_manual(values = setNames(c("steelblue", "tomato3"), c(ours_label, pantree_label))) +
    scale_y_log10(labels = scales::comma) +
    labs(title = title, subtitle = "Allele Frequency Spectrum",
         x = "Non-reference Allele Frequency", y = "Count (log scale)") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
  save_png(p4, paste0(prefix, ".pantree-af.png"))
} else {
  file.create(paste0(prefix, ".pantree-af.png"))
  cat("Skipping AF spectrum: AF data not available in both sources.\n")
}

# ---------------------------------------------------------------------------
# Summary TSV
# ---------------------------------------------------------------------------
summary_ours <- dt_ours[, .(count = .N), by = .(ref_context, variant_type)]
summary_ours[, source := ours_label]
summary_ours[, pct := 100 * count / sum(count)]

summary_pt <- dt_pt[, .(count = .N), by = .(ref_context, variant_type)]
summary_pt[, source := pantree_label]
summary_pt[, pct := 100 * count / sum(count)]

summary_all <- rbind(summary_ours, summary_pt)
setcolorder(summary_all, c("source", "ref_context", "variant_type", "count", "pct"))

fwrite(summary_all, paste0(prefix, ".pantree-compare.tsv"), sep = "\t")
cat("Wrote:", paste0(prefix, ".pantree-compare.tsv"), "\n")
cat("Done.\n")
