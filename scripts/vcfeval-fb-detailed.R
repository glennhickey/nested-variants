#!/usr/bin/env Rscript

# vcfeval-fb-detailed.R
#
# Extended vcfeval comparison panel for the FreeBayes-vs-vg-call short-read
# figure (4b). Parallel to scripts/vcf-compare-vcfeval.R but:
#
#   1. Splits FB-only (FP) calls into two buckets:
#        - fb_only_in_graph    — FP (CHROM, POS) appears in that sample's call
#                                VCF (vg call -a -A emits every graph site, so
#                                "in sample's call VCF" = "in graph")
#        - fb_only_not_in_graph — FP (CHROM, POS) absent from the sample's call
#                                VCF entirely
#   2. Stacks each bar by GIAB region (Easy / Segdup / Other Difficult / Other)
#      using the pipeline's augref-space GIAB BEDs.
#   3. Adds a Ts/Tv annotation above each SNP bar section, averaged across
#      samples within that (category × ref_context × giab) cell.
#
# Usage:
#   Rscript scripts/vcfeval-fb-detailed.R <output_prefix>
#     --vcfeval-dirs d1,d2,...
#     --samples s1,s2,...
#     --giab-beds  easy.bed,segdup.bed,otherdif.bed
#     --giab-names Easy,Segdup,Other_Difficult
#     [--filter-label STR] [--title TITLE]
#
# The "in-graph" lookup uses the baseline VCF that vcfeval was given, which
# lives in each vcfeval dir as `truth.vcf.gz` (prefix-renamed to match the
# augref-space CHROMs of the other vcfeval outputs).
#
# Note: vcfeval strips the FILTER field on its output VCFs (all records become
# FILTER='.'), so we do NOT re-apply PASS filtering here. The caller picks
# filter-ness by choosing the vcfeval-fb/{all,pass}/ dir; --filter-label is
# only for the plot subtitle.
#
# Outputs:
#   {prefix}.vcfeval-detailed.png
#   {prefix}.vcfeval-detailed.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript vcfeval-fb-detailed.R <output_prefix> [options]\n")
  quit(status = 1)
}

prefix        <- args[1]
title         <- NULL
vcfeval_dirs  <- NULL
sample_names  <- NULL
giab_beds     <- NULL
giab_names    <- NULL
filter_label_arg <- ""

i <- 2
while (i <= length(args)) {
  if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--vcfeval-dirs" && i + 1 <= length(args)) {
    vcfeval_dirs <- unlist(strsplit(args[i + 1], ",", fixed = TRUE)); i <- i + 2
  } else if (args[i] == "--samples" && i + 1 <= length(args)) {
    sample_names <- unlist(strsplit(args[i + 1], ",", fixed = TRUE)); i <- i + 2
  } else if (args[i] == "--giab-beds" && i + 1 <= length(args)) {
    giab_beds <- unlist(strsplit(args[i + 1], ",", fixed = TRUE)); i <- i + 2
  } else if (args[i] == "--giab-names" && i + 1 <= length(args)) {
    giab_names <- unlist(strsplit(args[i + 1], ",", fixed = TRUE)); i <- i + 2
  } else if (args[i] == "--filter-label" && i + 1 <= length(args)) {
    filter_label_arg <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

for (req in c("vcfeval_dirs", "sample_names", "giab_beds", "giab_names")) {
  if (is.null(get(req))) { cat("Error: --", req, " is required\n", sep = ""); quit(status = 1) }
}
if (length(vcfeval_dirs) != length(sample_names)) {
  cat("Error: --vcfeval-dirs and --samples must be the same length\n")
  quit(status = 1)
}
if (length(giab_beds) != length(giab_names)) {
  cat("Error: --giab-beds and --giab-names must be the same length\n")
  quit(status = 1)
}
n_samples <- length(sample_names)
filter_label <- if (nzchar(filter_label_arg)) paste0(", ", filter_label_arg) else ""
if (is.null(title)) title <- "Call vs FreeBayes — detailed (vcfeval)"

# Canonicalise the GIAB names to the display values used throughout the
# codebase (Easy / Segdup / Other Difficult), matching vcf-stats.R's palette.
normalize_giab_name <- function(nm) {
  nm <- gsub("_", " ", nm, fixed = TRUE)
  nm
}
giab_display <- normalize_giab_name(giab_names)
giab_order   <- c("Easy", "Segdup", "Other Difficult", "Other")
giab_colors  <- c("Easy" = "forestgreen",
                  "Segdup" = "firebrick",
                  "Other Difficult" = "darkorange",
                  "Other" = "grey60")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
read_vcf_records <- function(vcf) {
  if (!file.exists(vcf)) {
    return(data.table(CHROM = character(), POS = integer(),
                      REF = character(), ALT = character()))
  }
  cmd <- sprintf("bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' '%s' 2>/dev/null", vcf)
  tryCatch(
    fread(cmd = cmd, col.names = c("CHROM", "POS", "REF", "ALT"),
          colClasses = c(CHROM = "character", POS = "integer",
                         REF = "character", ALT = "character")),
    error = function(e) data.table(CHROM = character(), POS = integer(),
                                   REF = character(), ALT = character())
  )
}

in_graph_positions <- function(fp_dt, truth_vcf) {
  # For each FP record, return a logical vector indicating whether its
  # (CHROM, POS) is present in truth_vcf.
  #
  # Naive approach: fread the full 50 M-row truth VCF into a keyed table and
  # join. On HPRC that's minutes per sample. Instead, narrow the truth read
  # to *only* the ~10k FP positions via a BED + `bcftools view -R`. tabix
  # handles the region-tree lookup in a single indexed pass.
  if (nrow(fp_dt) == 0 || !file.exists(truth_vcf)) {
    return(rep(FALSE, nrow(fp_dt)))
  }
  bed <- tempfile(fileext = ".bed")
  on.exit(unlink(bed), add = TRUE)
  # BED is half-open 0-based; one row per unique FP position.
  uniq_pos <- unique(fp_dt[, .(CHROM, POS)])
  setorder(uniq_pos, CHROM, POS)
  fwrite(uniq_pos[, .(CHROM, POS - 1L, POS)], bed, sep = "\t", col.names = FALSE)
  cmd <- sprintf(
    "bcftools view -R '%s' '%s' 2>/dev/null | bcftools query -f '%%CHROM\\t%%POS\\n' 2>/dev/null",
    bed, truth_vcf)
  hits <- tryCatch(
    fread(cmd = cmd, col.names = c("CHROM", "POS"),
          colClasses = c(CHROM = "character", POS = "integer")),
    error = function(e) data.table(CHROM = character(), POS = integer())
  )
  if (nrow(hits) == 0) return(rep(FALSE, nrow(fp_dt)))
  hits <- unique(hits, by = c("CHROM", "POS"))
  setkey(hits, CHROM, POS)
  # data.table by-key membership lookup: returns row indices; NA where missing.
  idx <- hits[fp_dt[, .(CHROM, POS)], on = c("CHROM", "POS"),
              which = TRUE, mult = "first"]
  !is.na(idx)
}

classify_variants <- function(dt) {
  if (nrow(dt) == 0) {
    return(dt[, `:=`(variant_type = character(), ref_context = character(),
                     tstv = character())])
  }
  dt[, ref_len := nchar(REF)]
  is_multi <- grepl(",", dt$ALT, fixed = TRUE)
  dt[, size_signed := 0L]
  dt[(!is_multi), size_signed := nchar(ALT) - ref_len]
  dt[, size := abs(size_signed)]
  mi <- which(is_multi)
  if (length(mi) > 0) {
    dt[mi, c("size", "size_signed") := {
      res <- sapply(seq_len(.N), function(i) {
        alts <- unlist(strsplit(ALT[i], ","))
        diffs <- nchar(alts) - ref_len[i]
        idx <- which.max(abs(diffs))
        c(abs(diffs[idx]), diffs[idx])
      })
      list(res[1, ], res[2, ])
    }]
  }
  dt[, variant_type := fifelse(
    is.na(size), "Other",
    fifelse(size == 0L & ref_len == 1L, "SNP",
    fifelse(size == 0L,                 "MNP",
    fifelse(size < 50L & size_signed > 0L, "Insertion",
    fifelse(size < 50L,                 "Deletion",
    fifelse(size_signed > 0L, "SV Insertion", "SV Deletion")))))
  )]
  dt[, ref_context := fifelse(
    grepl("_[0-9]+_alt$", CHROM), "Off-reference", "On-reference"
  )]
  # Ts/Tv classification for bi-allelic SNPs (single-base REF & ALT)
  dt[, tstv := NA_character_]
  snp_mask <- dt$variant_type == "SNP" & nchar(dt$REF) == 1L & nchar(dt$ALT) == 1L
  dt[snp_mask, tstv := fifelse(
    paste0(REF, ALT) %in% c("AG", "GA", "CT", "TC"), "Ts", "Tv"
  )]
  dt
}

# Spatial-join a (CHROM, POS) data.table against a BED file using bedtools
# intersect. Returns a CHROM/POS-keyed data.table marking which rows hit.
bed_membership <- function(dt, bed_file) {
  if (nrow(dt) == 0 || !file.exists(bed_file)) return(rep(FALSE, nrow(dt)))
  tmp_in  <- tempfile(fileext = ".bed")
  tmp_sorted <- tempfile(fileext = ".bed")
  tmp_out <- tempfile(fileext = ".bed")
  on.exit(unlink(c(tmp_in, tmp_sorted, tmp_out)), add = TRUE)
  fwrite(dt[, .(CHROM, POS - 1L, POS)], tmp_in, sep = "\t", col.names = FALSE)
  system(sprintf("LC_ALL=C sort -k1,1 -k2,2n '%s' > '%s'", tmp_in, tmp_sorted))
  # -u: unique hits only; -wa: print A record (ensures 1 row per input)
  system(sprintf("bedtools intersect -u -a '%s' -b '%s' > '%s' 2>/dev/null",
                 tmp_sorted, bed_file, tmp_out))
  hits <- tryCatch(
    fread(tmp_out, header = FALSE, col.names = c("CHROM", "bed_start", "POS")),
    error = function(e) data.table(CHROM = character(), bed_start = integer(),
                                    POS = integer())
  )
  if (nrow(hits) == 0) return(rep(FALSE, nrow(dt)))
  hits[, hit := TRUE]
  merged <- merge(dt[, .(CHROM, POS, .I)], hits[, .(CHROM, POS, hit)],
                  by = c("CHROM", "POS"), all.x = TRUE, sort = FALSE)
  setorder(merged, I)
  !is.na(merged$hit)
}

assign_giab_region <- function(dt, giab_beds, giab_display_names) {
  # One pass per BED; first hit wins in the priority order the user provided,
  # defaulting to "Other" when nothing matches.
  if (nrow(dt) == 0) return(character(0))
  out <- rep("Other", nrow(dt))
  assigned <- rep(FALSE, nrow(dt))
  for (k in seq_along(giab_beds)) {
    nm <- giab_display_names[k]
    hits <- bed_membership(dt[!assigned], giab_beds[k])
    idx <- which(!assigned)[hits]
    if (length(idx) > 0) {
      out[idx] <- nm
      assigned[idx] <- TRUE
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Per-sample: build a long table (sample × category × ref_context × variant_type
# × giab_region) → count + Ts/Tv tallies.
# ---------------------------------------------------------------------------
rows <- vector("list", n_samples)
for (si in seq_len(n_samples)) {
  sname  <- sample_names[si]
  ve_dir <- vcfeval_dirs[si]
  cat("[", sname, "] reading vcfeval from ", ve_dir, "\n", sep = "")

  tp <- classify_variants(read_vcf_records(file.path(ve_dir, "tp-baseline.vcf.gz")))
  fn <- classify_variants(read_vcf_records(file.path(ve_dir, "fn.vcf.gz")))
  fp <- classify_variants(read_vcf_records(file.path(ve_dir, "fp.vcf.gz")))

  # --- Classify FP as in-graph vs not-in-graph by querying the baseline
  # truth VCF only at the ~few-thousand FP positions (via bcftools view -R).
  # Avoids a full fread of the 50 M-row truth.vcf.gz on HPRC.
  if (nrow(fp) > 0) {
    fp[, in_graph := in_graph_positions(fp, file.path(ve_dir, "truth.vcf.gz"))]
  } else {
    fp[, in_graph := logical(0)]
  }
  cat("  TP:", nrow(tp), " FN:", nrow(fn),
      " FP:", nrow(fp),
      " (in_graph:", sum(fp$in_graph),
      " / not_in_graph:", sum(!fp$in_graph), ")\n")

  # --- GIAB region assignment for all three subsets ----------------------
  tp[, giab := assign_giab_region(tp, giab_beds, giab_display)]
  fn[, giab := assign_giab_region(fn, giab_beds, giab_display)]
  fp[, giab := assign_giab_region(fp, giab_beds, giab_display)]

  # --- Tag each subset with its comparison category ---------------------
  tp[, category := "Shared"]
  fn[, category := "Call only"]
  if (nrow(fp) > 0) {
    fp[, category := fifelse(in_graph, "FB only (in graph)", "FB only (not in graph)")]
  } else {
    fp[, category := character(0)]
  }

  cols <- c("variant_type", "ref_context", "giab", "category", "tstv")
  combined <- rbindlist(list(tp[, ..cols], fn[, ..cols], fp[, ..cols]),
                        use.names = TRUE, fill = TRUE)
  combined[, sample := sname]
  rows[[si]] <- combined
}

all_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
if (nrow(all_dt) == 0) {
  cat("No records found. Writing empty outputs.\n")
  file.create(paste0(prefix, ".vcfeval-detailed.png"))
  fwrite(data.table(), paste0(prefix, ".vcfeval-detailed.tsv"), sep = "\t")
  quit(status = 0)
}

# ---------------------------------------------------------------------------
# Aggregate counts: (sample × ref_context × variant_type × category × giab)
# ---------------------------------------------------------------------------
counts_dt <- all_dt[, .(count = .N,
                        ts = sum(tstv == "Ts", na.rm = TRUE),
                        tv = sum(tstv == "Tv", na.rm = TRUE)),
                    by = .(sample, ref_context, variant_type, category, giab)]

# Ensure all combinations exist with zeros
cat_levels  <- c("Shared", "Call only", "FB only (in graph)", "FB only (not in graph)")
type_levels <- c("SNP", "MNP", "Insertion", "Deletion", "SV Insertion", "SV Deletion", "Other")
present_giab <- intersect(giab_order, unique(counts_dt$giab))
if (length(present_giab) == 0) present_giab <- giab_order

all_combos <- CJ(sample = sample_names,
                 ref_context = unique(counts_dt$ref_context),
                 variant_type = intersect(type_levels, unique(counts_dt$variant_type)),
                 category = cat_levels,
                 giab = present_giab)
counts_dt <- merge(all_combos, counts_dt,
                   by = c("sample", "ref_context", "variant_type", "category", "giab"),
                   all.x = TRUE)
counts_dt[is.na(count), `:=`(count = 0L, ts = 0L, tv = 0L)]

# Ordering
counts_dt[, variant_type := factor(variant_type,
  levels = intersect(type_levels, unique(variant_type)))]
counts_dt[, category := factor(category, levels = cat_levels)]
counts_dt[, giab := factor(giab, levels = giab_order)]
setorder(counts_dt, sample, ref_context, variant_type, category, giab)

fwrite(counts_dt, paste0(prefix, ".vcfeval-detailed.tsv"), sep = "\t")
cat("Wrote:", paste0(prefix, ".vcfeval-detailed.tsv"), "\n")

# ---------------------------------------------------------------------------
# Ts/Tv mean per (ref_context × variant_type × category × giab) across samples.
# Compute each sample's per-stratum Ts/Tv, then average those ratios.
# ---------------------------------------------------------------------------
snp_sample <- counts_dt[variant_type == "SNP" & (ts + tv) > 0,
                        .(ratio = ts / tv),
                        by = .(sample, ref_context, category, giab)]
snp_mean <- snp_sample[, .(tstv = mean(ratio)),
                       by = .(ref_context, category, giab)]

# ---------------------------------------------------------------------------
# Plot: stacked bars (fill = GIAB region) per (variant_type × category),
#       faceted by ref_context. Height = mean across samples.
# ---------------------------------------------------------------------------
bar_summary <- counts_dt[, .(mean_count = mean(count)),
                         by = .(ref_context, variant_type, category, giab)]

# Compute stacked y-positions for Ts/Tv labels on SNP bars.
#   x-axis position of each bar = category (within a given variant_type panel
#   of the facet). We'll place text at the top of each giab segment.
snp_bar <- bar_summary[variant_type == "SNP"]
setorder(snp_bar, ref_context, category, giab)
snp_bar[, y_top := cumsum(mean_count), by = .(ref_context, category)]
snp_bar[, y_center := y_top - mean_count / 2]
# Join the Ts/Tv ratio for the label
snp_bar <- merge(snp_bar, snp_mean,
                 by = c("ref_context", "category", "giab"),
                 all.x = TRUE)
# Only label bars that have meaningful SNP counts (> ~1% of panel total)
snp_bar[, panel_total := sum(mean_count), by = .(ref_context, category)]
snp_bar[, keep := !is.na(tstv) & mean_count > 0.05 * panel_total]

cat_colors_shape <- c("Shared" = "forestgreen",
                      "Call only" = "steelblue",
                      "FB only (in graph)" = "coral",
                      "FB only (not in graph)" = "firebrick4")

p <- ggplot(bar_summary,
            aes(x = category, y = mean_count, fill = giab)) +
  geom_col(position = position_stack(), width = 0.7, alpha = 0.9) +
  facet_grid(ref_context ~ variant_type, scales = "free_y", switch = "y") +
  scale_fill_manual(values = giab_colors, name = "GIAB region", drop = FALSE) +
  scale_y_continuous(labels = scales::comma) +
  geom_text(data = snp_bar[keep == TRUE],
            aes(x = category, y = y_center,
                label = sprintf("%.2f", tstv)),
            inherit.aes = FALSE, size = 2.8, color = "white") +
  labs(title = title,
       subtitle = paste0("Call vs FreeBayes via vcfeval", filter_label,
                         " — FB-only split by graph membership",
                         ", stacks = GIAB region, mean across N=",
                         n_samples, " samples",
                         "; text on SNP bars = Ts/Tv"),
       x = NULL, y = "Mean count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    strip.placement = "outside"
  )

save_png <- function(plot, path, w = 14, h = 8) {
  tryCatch({
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(path, width = w, height = h, units = "in", res = 300)
      print(plot); dev.off()
    } else {
      ggsave(path, plot = plot, width = w, height = h, dpi = 300,
             device = grDevices::png, type = "cairo")
    }
  }, error = function(e) {
    grDevices::png(path, width = w * 300, height = h * 300, res = 300, type = "cairo")
    print(plot); dev.off()
  })
}
save_png(p, paste0(prefix, ".vcfeval-detailed.png"))
cat("Saved:", paste0(prefix, ".vcfeval-detailed.png"), "\n")
cat("Done.\n")
