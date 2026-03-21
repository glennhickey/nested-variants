#!/usr/bin/env Rscript

# vcf-compare.R — Compare two multisample VCFs per-sample (shared vs unique)
#
# Usage: Rscript scripts/vcf-compare.R <call.vcf.gz> <dv.vcf.gz> <output_prefix>
#          [--mode sites|variants] [--filter all|pass]
#          [--label-a LABEL] [--label-b LABEL] [--title TITLE]
#
# Outputs:
#   {prefix}.compare.png  — grouped bar chart of shared / A-only / B-only
#   {prefix}.compare.tsv  — per-sample comparison counts

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript vcf-compare.R <vcf_a.vcf.gz> <vcf_b.vcf.gz> <output_prefix> [options]\n")
  quit(status = 1)
}

vcf_a   <- args[1]
vcf_b   <- args[2]
prefix  <- args[3]
title        <- NULL
mode         <- "sites"
filter       <- "all"
label_a      <- "Call"
label_b      <- "DeepVariant"
strip_prefix <- NULL
no_sv        <- FALSE

i <- 4
while (i <= length(args)) {
  if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--mode" && i + 1 <= length(args)) {
    mode <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--filter" && i + 1 <= length(args)) {
    filter <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--label-a" && i + 1 <= length(args)) {
    label_a <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--label-b" && i + 1 <= length(args)) {
    label_b <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--strip-prefix" && i + 1 <= length(args)) {
    strip_prefix <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--no-sv") {
    no_sv <- TRUE; i <- i + 1
  } else {
    i <- i + 1
  }
}
if (!mode %in% c("sites", "variants")) {
  cat("Error: --mode must be 'sites' or 'variants'\n"); quit(status = 1)
}
if (!filter %in% c("all", "pass")) {
  cat("Error: --filter must be 'all' or 'pass'\n"); quit(status = 1)
}
if (is.null(title)) title <- paste(label_a, "vs", label_b)

mode_label <- if (mode == "sites") "(per site)" else "(per variant)"
filter_label <- if (filter == "pass") ", PASS only" else ""

# ---------------------------------------------------------------------------
# Helper: build bcftools pipeline prefix (same pattern as vcf-stats.R)
# ---------------------------------------------------------------------------
build_pipe_prefix <- function(vcf_path, mode, filter) {
  filter_cmd <- if (filter == "pass") "bcftools view -f PASS 2>/dev/null |" else ""
  if (mode == "variants") {
    sprintf("bcftools norm -m- '%s' 2>/dev/null | bcftools view -c1 2>/dev/null | %s bcftools +fill-tags - -- -t AF 2>/dev/null",
            vcf_path, filter_cmd)
  } else {
    sprintf("bcftools view -c1 '%s' 2>/dev/null | %s bcftools +fill-tags - -- -t AF 2>/dev/null",
            vcf_path, filter_cmd)
  }
}

# ---------------------------------------------------------------------------
# Helper: read VCF with per-sample GTs
# ---------------------------------------------------------------------------
read_vcf_with_gt <- function(vcf_path, mode, filter) {
  # Get sample names
  sample_names <- system(sprintf("bcftools query -l '%s' 2>/dev/null", vcf_path), intern = TRUE)
  n_samples <- length(sample_names)
  if (n_samples == 0) stop("No samples found in ", vcf_path)

  pipe_prefix <- build_pipe_prefix(vcf_path, mode, filter)
  cmd <- sprintf("%s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT]\\n' 2>/dev/null",
                 pipe_prefix)
  col_names <- c("CHROM", "POS", "REF", "ALT", paste0("GT_", seq_len(n_samples)))
  dt <- fread(cmd = cmd, col.names = col_names)
  list(dt = dt, samples = sample_names, n = n_samples)
}

# ---------------------------------------------------------------------------
# Helper: classify variant types (same logic as vcf-stats.R)
# ---------------------------------------------------------------------------
classify_variants <- function(dt) {
  dt[, ref_len := nchar(REF)]
  is_multi <- grepl(",", dt$ALT, fixed = TRUE)
  dt[, size_signed := 0L]
  dt[(!is_multi), size_signed := nchar(ALT) - ref_len]
  dt[, size := abs(size_signed)]
  multi_idx <- which(is_multi)
  if (length(multi_idx) > 0) {
    dt[multi_idx, c("size", "size_signed") := {
      res <- vapply(seq_len(.N), function(i) {
        alts <- unlist(strsplit(ALT[i], ","))
        alts <- alts[alts != "*" & alts != "."]
        if (length(alts) == 0) return(c(NA_real_, NA_real_))
        diffs <- nchar(alts) - ref_len[i]
        idx <- which.max(abs(diffs))
        c(abs(diffs[idx]), diffs[idx])
      }, numeric(2))
      if (is.null(dim(res))) res <- matrix(res, nrow = 2)
      list(res[1,], res[2,])
    }]
  }
  type_levels <- c("SNP", "MNP", "Insertion", "Deletion", "SV Insertion", "SV Deletion", "Other")
  dt[, variant_type := fifelse(
    is.na(size), "Other",
    fifelse(size == 0L & ref_len == 1L, "SNP",
    fifelse(size == 0L, "MNP",
    fifelse(size < 50L & size_signed > 0L, "Insertion",
    fifelse(size < 50L, "Deletion",
    fifelse(size_signed > 0L, "SV Insertion", "SV Deletion")))))
  )]
  dt[, ref_context := fifelse(
    grepl("_[0-9]+_alt$", CHROM), "Off-reference", "On-reference"
  )]
  dt
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
# Read both VCFs
# ---------------------------------------------------------------------------
cat("Reading VCF A:", vcf_a, "\n")
vcf_a_data <- read_vcf_with_gt(vcf_a, mode, filter)
cat("  Samples:", vcf_a_data$n, " Records:", nrow(vcf_a_data$dt), "\n")

cat("Reading VCF B:", vcf_b, "\n")
vcf_b_data <- read_vcf_with_gt(vcf_b, mode, filter)
cat("  Samples:", vcf_b_data$n, " Records:", nrow(vcf_b_data$dt), "\n")

# Strip augref prefix from CHROM values so both VCFs use the same names
# (vg call -S strips the prefix but DeepVariant keeps the FASTA contig names)
if (!is.null(strip_prefix)) {
  cat("Stripping prefix '", strip_prefix, "' from CHROM values\n", sep = "")
  vcf_a_data$dt[, CHROM := sub(strip_prefix, "", CHROM, fixed = TRUE)]
  vcf_b_data$dt[, CHROM := sub(strip_prefix, "", CHROM, fixed = TRUE)]
}

# Find overlapping samples
shared_samples <- intersect(vcf_a_data$samples, vcf_b_data$samples)
n_shared <- length(shared_samples)
cat("Shared samples:", n_shared, "\n")

if (n_shared == 0) {
  cat("No overlapping samples between the two VCFs. Exiting.\n")
  file.create(paste0(prefix, ".compare.png"))
  fwrite(data.table(sample = character(), variant_type = character(),
                    ref_context = character(), category = character(), count = integer()),
         paste0(prefix, ".compare.tsv"), sep = "\t")
  quit(status = 0)
}

only_a <- setdiff(vcf_a_data$samples, vcf_b_data$samples)
only_b <- setdiff(vcf_b_data$samples, vcf_a_data$samples)
if (length(only_a) > 0) cat("  Warning: samples only in A:", paste(only_a, collapse = ", "), "\n")
if (length(only_b) > 0) cat("  Warning: samples only in B:", paste(only_b, collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# Restrict to shared contigs so caller-specific contig sets don't inflate X-only
# ---------------------------------------------------------------------------
chroms_a <- unique(vcf_a_data$dt$CHROM)
chroms_b <- unique(vcf_b_data$dt$CHROM)
shared_chroms <- intersect(chroms_a, chroms_b)
only_a_chroms <- setdiff(chroms_a, shared_chroms)
only_b_chroms <- setdiff(chroms_b, shared_chroms)
if (length(only_a_chroms) > 0) {
  n_drop_a <- nrow(vcf_a_data$dt[CHROM %in% only_a_chroms])
  cat("Dropping", length(only_a_chroms), "contigs (", n_drop_a, "records) from",
      label_a, "not in", label_b, "\n")
  vcf_a_data$dt <- vcf_a_data$dt[CHROM %in% shared_chroms]
}
if (length(only_b_chroms) > 0) {
  n_drop_b <- nrow(vcf_b_data$dt[CHROM %in% only_b_chroms])
  cat("Dropping", length(only_b_chroms), "contigs (", n_drop_b, "records) from",
      label_b, "not in", label_a, "\n")
  vcf_b_data$dt <- vcf_b_data$dt[CHROM %in% shared_chroms]
}
cat("Shared contigs:", length(shared_chroms), "\n")

# ---------------------------------------------------------------------------
# Classify variants in both
# ---------------------------------------------------------------------------
dt_a <- classify_variants(vcf_a_data$dt)
dt_b <- classify_variants(vcf_b_data$dt)

# Drop SV categories if requested (callers like DeepVariant don't call SVs)
if (no_sv) {
  dt_a <- dt_a[!variant_type %in% c("SV Insertion", "SV Deletion")]
  dt_b <- dt_b[!variant_type %in% c("SV Insertion", "SV Deletion")]
}

# ---------------------------------------------------------------------------
# For each shared sample: determine shared / A-only / B-only
# ---------------------------------------------------------------------------
sv_types <- if (no_sv) character(0) else c("SV Insertion", "SV Deletion")
type_levels <- c("SNP", "MNP", "Insertion", "Deletion", sv_types, "Other")

results_list <- vector("list", n_shared)
for (si in seq_along(shared_samples)) {
  sname <- shared_samples[si]

  # Find GT column index in each VCF
  idx_a <- which(vcf_a_data$samples == sname)
  idx_b <- which(vcf_b_data$samples == sname)
  gt_col_a <- paste0("GT_", idx_a)
  gt_col_b <- paste0("GT_", idx_b)

  # Variant keys for carriers in A
  carriers_a <- dt_a[grepl("[1-9]", get(gt_col_a)),
                     .(key = paste(CHROM, POS, REF, ALT, sep = ":"),
                       variant_type, ref_context)]
  # Variant keys for carriers in B
  carriers_b <- dt_b[grepl("[1-9]", get(gt_col_b)),
                     .(key = paste(CHROM, POS, REF, ALT, sep = ":"),
                       variant_type, ref_context)]

  # Shared keys
  shared_keys <- intersect(carriers_a$key, carriers_b$key)

  # Classify
  carriers_a[, category := fifelse(key %in% shared_keys, "Shared",
                                   paste0(label_a, " only"))]
  carriers_b[!key %in% shared_keys, category := paste0(label_b, " only")]
  carriers_b <- carriers_b[!key %in% shared_keys]

  combined <- rbind(carriers_a[, .(variant_type, ref_context, category)],
                    carriers_b[, .(variant_type, ref_context, category)])
  counts <- combined[, .(count = .N), by = .(variant_type, ref_context, category)]
  counts[, sample := sname]
  results_list[[si]] <- counts
}

compare_dt <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

# Ensure all combinations exist
cat_levels <- c("Shared", paste0(label_a, " only"), paste0(label_b, " only"))
all_combos <- CJ(sample = shared_samples,
                  variant_type = unique(compare_dt$variant_type),
                  ref_context = unique(compare_dt$ref_context),
                  category = cat_levels)
compare_dt <- merge(all_combos, compare_dt,
                     by = c("sample", "variant_type", "ref_context", "category"),
                     all.x = TRUE)
compare_dt[is.na(count), count := 0L]

# ---------------------------------------------------------------------------
# Write TSV
# ---------------------------------------------------------------------------
setorder(compare_dt, sample, ref_context, variant_type, category)
fwrite(compare_dt, paste0(prefix, ".compare.tsv"), sep = "\t")
cat("Wrote comparison:", paste0(prefix, ".compare.tsv"), "\n")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
# Summary across samples
compare_summary <- compare_dt[, .(mean_count = mean(count),
                                   min_count = min(count),
                                   max_count = max(count)),
                               by = .(variant_type, ref_context, category)]

compare_dt[, variant_type := factor(variant_type,
  levels = intersect(type_levels, unique(variant_type)))]
compare_summary[, variant_type := factor(variant_type,
  levels = intersect(type_levels, unique(variant_type)))]
compare_dt[, category := factor(category, levels = cat_levels)]
compare_summary[, category := factor(category, levels = cat_levels)]

cat_colors <- c("Shared" = "forestgreen",
                setNames("steelblue", paste0(label_a, " only")),
                setNames("coral", paste0(label_b, " only")))

if (n_shared <= 20) {
  p_compare <- ggplot() +
    geom_col(data = compare_summary,
             aes(x = variant_type, y = mean_count, fill = category),
             position = position_dodge(width = 0.7), width = 0.7, alpha = 0.6) +
    geom_errorbar(data = compare_summary,
                  aes(x = variant_type, ymin = min_count, ymax = max_count,
                      group = category),
                  position = position_dodge(width = 0.7), width = 0.3) +
    geom_point(data = compare_dt,
               aes(x = variant_type, y = count, color = category),
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
               size = 1.5, alpha = 0.8) +
    facet_wrap(~ref_context, scales = "free_y") +
    scale_fill_manual(values = cat_colors, name = NULL) +
    scale_color_manual(values = cat_colors, name = NULL) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = title,
         subtitle = paste0(label_a, " vs ", label_b, " ", mode_label, filter_label,
                           " (N=", n_shared, " samples, bars=mean)"),
         x = "Variant Type", y = "Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
} else {
  p_compare <- ggplot(compare_dt,
                      aes(x = variant_type, y = count, fill = category)) +
    geom_boxplot(position = position_dodge(width = 0.7), width = 0.6,
                 outlier.size = 1) +
    facet_wrap(~ref_context, scales = "free_y") +
    scale_fill_manual(values = cat_colors, name = NULL) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = title,
         subtitle = paste0(label_a, " vs ", label_b, " ", mode_label, filter_label,
                           " (N=", n_shared, " samples)"),
         x = "Variant Type", y = "Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
}

save_png(p_compare, paste0(prefix, ".compare.png"), width = 12)

cat("Done.\n")
