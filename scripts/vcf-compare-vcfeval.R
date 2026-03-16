#!/usr/bin/env Rscript

# vcf-compare-vcfeval.R — Plot vcfeval TP/FP/FN counts by variant type
#
# Usage: Rscript scripts/vcf-compare-vcfeval.R <output_prefix>
#          --vcfeval-dirs dir1,dir2,...   (per-sample vcfeval output directories)
#          --samples sample1,sample2,... (corresponding sample names)
#          [--label-a LABEL] [--label-b LABEL] [--title TITLE]
#          [--filter all|pass] [--no-sv]
#
# Each vcfeval directory must contain tp.vcf.gz, fp.vcf.gz, fn.vcf.gz.
# TP = shared variants, FN = base-only (label-a only), FP = call-only (label-b only)
#
# Outputs:
#   {prefix}.vcfeval-compare.png  — grouped bar chart of TP / base-only / call-only
#   {prefix}.vcfeval-compare.tsv  — per-sample comparison counts

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript vcf-compare-vcfeval.R <output_prefix> [options]\n")
  quit(status = 1)
}

prefix       <- args[1]
title        <- NULL
label_a      <- "Call"
label_b      <- "DeepVariant"
vcfeval_dirs <- NULL
sample_names <- NULL
no_sv        <- FALSE
filter       <- "all"

i <- 2
while (i <= length(args)) {
  if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--label-a" && i + 1 <= length(args)) {
    label_a <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--label-b" && i + 1 <= length(args)) {
    label_b <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--vcfeval-dirs" && i + 1 <= length(args)) {
    vcfeval_dirs <- unlist(strsplit(args[i + 1], ",", fixed = TRUE)); i <- i + 2
  } else if (args[i] == "--samples" && i + 1 <= length(args)) {
    sample_names <- unlist(strsplit(args[i + 1], ",", fixed = TRUE)); i <- i + 2
  } else if (args[i] == "--filter" && i + 1 <= length(args)) {
    filter <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--no-sv") {
    no_sv <- TRUE; i <- i + 1
  } else {
    i <- i + 1
  }
}

if (is.null(vcfeval_dirs) || is.null(sample_names)) {
  cat("Error: --vcfeval-dirs and --samples are required\n")
  quit(status = 1)
}
if (length(vcfeval_dirs) != length(sample_names)) {
  cat("Error: --vcfeval-dirs and --samples must have the same number of entries\n")
  quit(status = 1)
}
if (!filter %in% c("all", "pass")) {
  cat("Error: --filter must be 'all' or 'pass'\n"); quit(status = 1)
}
filter_label <- if (filter == "pass") ", PASS only" else ""
if (is.null(title)) title <- paste(label_a, "vs", label_b, "(vcfeval)")

n_samples <- length(sample_names)
cat("Samples:", n_samples, "\n")

# ---------------------------------------------------------------------------
# Helper: classify variant types (same logic as vcf-compare.R / vcf-stats.R)
# ---------------------------------------------------------------------------
classify_variants <- function(dt) {
  if (nrow(dt) == 0) return(dt[, .(variant_type = character(), ref_context = character())])
  dt[, ref_len := nchar(REF)]
  is_multi <- grepl(",", dt$ALT, fixed = TRUE)
  dt[, size_signed := 0L]
  dt[(!is_multi), size_signed := nchar(ALT) - ref_len]
  dt[, size := abs(size_signed)]
  multi_idx <- which(is_multi)
  if (length(multi_idx) > 0) {
    dt[multi_idx, c("size", "size_signed") := {
      res <- sapply(seq_len(.N), function(i) {
        alts <- unlist(strsplit(ALT[i], ","))
        diffs <- nchar(alts) - ref_len[i]
        idx <- which.max(abs(diffs))
        c(abs(diffs[idx]), diffs[idx])
      })
      list(res[1,], res[2,])
    }]
  }
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
# Helper: read a VCF via bcftools query (no sample GT needed — just records)
# ---------------------------------------------------------------------------
read_vcf_records <- function(vcf_path, filter = "all") {
  if (!file.exists(vcf_path)) {
    cat("  Warning: missing", vcf_path, "\n")
    return(data.table(CHROM = character(), POS = integer(), REF = character(), ALT = character()))
  }
  if (filter == "pass") {
    cmd <- sprintf("bcftools view -f PASS '%s' 2>/dev/null | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' 2>/dev/null", vcf_path)
  } else {
    cmd <- sprintf("bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' '%s' 2>/dev/null", vcf_path)
  }
  dt <- tryCatch(
    fread(cmd = cmd, col.names = c("CHROM", "POS", "REF", "ALT")),
    error = function(e) data.table(CHROM = character(), POS = integer(), REF = character(), ALT = character())
  )
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
# Read vcfeval outputs for each sample
# ---------------------------------------------------------------------------
cat_levels <- c("Shared", paste0(label_a, " only"), paste0(label_b, " only"))

results_list <- vector("list", n_samples)
for (si in seq_len(n_samples)) {
  sname <- sample_names[si]
  ve_dir <- vcfeval_dirs[si]
  cat("Reading vcfeval results for", sname, "from", ve_dir, "\n")

  # Use tp-baseline.vcf.gz for TP counts: it preserves the baseline (label_a)
  # representation so MNPs are not decomposed into SNPs.  tp.vcf.gz contains
  # the call-side view where --decompose splits MNPs into individual SNPs.
  tp_dt <- read_vcf_records(file.path(ve_dir, "tp-baseline.vcf.gz"), filter)
  fn_dt <- read_vcf_records(file.path(ve_dir, "fn.vcf.gz"), filter)
  fp_dt <- read_vcf_records(file.path(ve_dir, "fp.vcf.gz"), filter)

  tp_dt <- classify_variants(tp_dt)
  fn_dt <- classify_variants(fn_dt)
  fp_dt <- classify_variants(fp_dt)

  # Drop SV categories if requested
  if (no_sv) {
    tp_dt <- tp_dt[!variant_type %in% c("SV Insertion", "SV Deletion")]
    fn_dt <- fn_dt[!variant_type %in% c("SV Insertion", "SV Deletion")]
    fp_dt <- fp_dt[!variant_type %in% c("SV Insertion", "SV Deletion")]
  }

  tp_counts <- tp_dt[, .(count = .N), by = .(variant_type, ref_context)][, category := "Shared"]
  fn_counts <- fn_dt[, .(count = .N), by = .(variant_type, ref_context)][, category := paste0(label_a, " only")]
  fp_counts <- fp_dt[, .(count = .N), by = .(variant_type, ref_context)][, category := paste0(label_b, " only")]

  combined <- rbind(tp_counts, fn_counts, fp_counts)
  combined[, sample := sname]
  results_list[[si]] <- combined

  cat("  TP:", nrow(tp_dt), " FN:", nrow(fn_dt), " FP:", nrow(fp_dt), "\n")
}

compare_dt <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

if (nrow(compare_dt) == 0) {
  cat("No variants found in vcfeval outputs. Creating empty outputs.\n")
  file.create(paste0(prefix, ".vcfeval-compare.png"))
  fwrite(data.table(sample = character(), variant_type = character(),
                    ref_context = character(), category = character(), count = integer()),
         paste0(prefix, ".vcfeval-compare.tsv"), sep = "\t")
  quit(status = 0)
}

# Ensure all combinations exist
all_combos <- CJ(sample = sample_names,
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
fwrite(compare_dt, paste0(prefix, ".vcfeval-compare.tsv"), sep = "\t")
cat("Wrote comparison:", paste0(prefix, ".vcfeval-compare.tsv"), "\n")

# ---------------------------------------------------------------------------
# Plot (same style as vcf-compare.R)
# ---------------------------------------------------------------------------
sv_types <- if (no_sv) character(0) else c("SV Insertion", "SV Deletion")
type_levels <- c("SNP", "MNP", "Insertion", "Deletion", sv_types, "Other")

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

if (n_samples <= 20) {
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
         subtitle = paste0(label_a, " vs ", label_b, " via vcfeval",
                           filter_label, " (N=", n_samples, " samples, bars=mean)"),
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
         subtitle = paste0(label_a, " vs ", label_b, " via vcfeval",
                           filter_label, " (N=", n_samples, " samples)"),
         x = "Variant Type", y = "Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
}

save_png(p_compare, paste0(prefix, ".vcfeval-compare.png"), width = 12)

cat("Done.\n")
