#!/usr/bin/env Rscript

# vcf-stats.R — VCF variant statistics with on-ref vs off-ref breakdown
#
# Usage: Rscript scripts/vcf-stats.R <input.vcf.gz> <output_prefix>
#          [--title TITLE] [--mode sites|variants] [--filter all|pass]
#          [--af-step SIZE]       (AF rounding step for spectrum plot; default 0.1)
#          [--annot-beds F1,F2]   (comma-separated augref-space annotation BED paths)
#          [--annot-names N1,N2]  (comma-separated display names, same order)
#
# Modes:
#   sites    — (default) read VCF as-is; multi-allelic sites classified by largest allele
#   variants — pipe through `bcftools norm -m-` first to split multi-allelic records
#
# Outputs:
#   {prefix}.vcf-stats.tsv      — summary table
#   {prefix}.variant-types.png  — grouped bar chart of variant types
#   {prefix}.size-dist.png      — indel/SV size distribution (two-panel)
#   {prefix}.size-dist-log.png  — same as above with log y-axis
#   {prefix}.af-spectrum.png    — allele frequency histogram (only when AF present)
#   {prefix}.variant-types-by-annot.png — variant types faceted by annotation (when --annot-beds)
#   {prefix}.vcf-stats-by-annot.tsv     — annotation × variant type counts (when --annot-beds)

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
title   <- NULL
mode    <- "sites"
filter  <- "all"
af_step <- 0.1
annot_beds     <- NULL
annot_names_arg <- NULL
giab_strat_beds_arg  <- NULL
giab_strat_names_arg <- NULL
per_sample <- FALSE
ref_sample <- NULL
no_sv      <- FALSE
tsv_input  <- FALSE
dump_records <- FALSE
records_only <- FALSE

i <- 3
while (i <= length(args)) {
  if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--mode" && i + 1 <= length(args)) {
    mode <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--filter" && i + 1 <= length(args)) {
    filter <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--af-step" && i + 1 <= length(args)) {
    af_step <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--annot-beds" && i + 1 <= length(args)) {
    annot_beds <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--annot-names" && i + 1 <= length(args)) {
    annot_names_arg <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--giab-strat-beds" && i + 1 <= length(args)) {
    giab_strat_beds_arg <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--giab-strat-names" && i + 1 <= length(args)) {
    giab_strat_names_arg <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--per-sample") {
    per_sample <- TRUE
    i <- i + 1
  } else if (args[i] == "--ref-sample" && i + 1 <= length(args)) {
    ref_sample <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--no-sv") {
    no_sv <- TRUE
    i <- i + 1
  } else if (args[i] == "--tsv") {
    tsv_input <- TRUE
    i <- i + 1
  } else if (args[i] == "--dump-records") {
    dump_records <- TRUE
    i <- i + 1
  } else if (args[i] == "--records-only") {
    records_only <- TRUE
    i <- i + 1
  } else {
    i <- i + 1
  }
}
if (!mode %in% c("sites", "variants")) {
  cat("Error: --mode must be 'sites' or 'variants', got '", mode, "'\n", sep = "")
  quit(status = 1)
}
if (!filter %in% c("all", "pass")) {
  cat("Error: --filter must be 'all' or 'pass', got '", filter, "'\n", sep = "")
  quit(status = 1)
}
if (is.null(title)) title <- basename(vcf)

mode_label <- if (mode == "sites") "(per site)" else "(per variant)"
filter_label <- if (filter == "pass") ", PASS only" else ""

# ---------------------------------------------------------------------------
# Read VCF (or pre-extracted TSV)
# ---------------------------------------------------------------------------
if (tsv_input) {
  cat("Reading pre-extracted TSV:", vcf, "\n")
  dt <- fread(vcf)
  has_af <- "nonref_af" %in% names(dt) && !all(is.na(dt$nonref_af))
  has_tr <- "is_repeat" %in% names(dt) && any(dt$is_repeat == TRUE | dt$is_repeat == "TRUE")
  # Convert is_repeat to logical if character
  if ("is_repeat" %in% names(dt) && is.character(dt$is_repeat)) {
    dt[, is_repeat := (is_repeat == "TRUE")]
  }
  # Ts/Tv classification for on-reference SNPs with single-base REF and ALT
  if ("REF" %in% names(dt) && "ALT" %in% names(dt)) {
    dt[variant_type == "SNP" & nchar(REF) == 1 & nchar(ALT) == 1, tstv := {
      transitions <- c("AG", "GA", "CT", "TC")
      fifelse(paste0(REF, ALT) %in% transitions, "Ts", "Tv")
    }]
  }
  cat("Read", nrow(dt), "variant records\n")
  cat("AF available:", has_af, "\n")
  if (has_tr) cat("TR annotation detected:", sum(dt$is_repeat), "tandem repeat indels\n")
  if (nrow(dt) == 0) { cat("No variants found. Exiting.\n"); quit(status = 0) }
} else {
cat("Reading VCF:", vcf, " (mode:", mode, ", filter:", filter, ")\n")

# Build pipeline prefix: variants mode pipes through bcftools norm -m- first
# +fill-tags computes AF from genotypes when INFO/AF is absent (e.g., single-sample vg call)
# When filter=="pass", insert bcftools view -f PASS to keep only PASS variants
filter_cmd <- if (filter == "pass") "bcftools view -f PASS 2>/dev/null |" else ""

if (mode == "variants") {
  pipe_prefix <- sprintf("bcftools norm -m- '%s' 2>/dev/null | bcftools view -c1 2>/dev/null | %s bcftools +fill-tags - -- -t AF 2>/dev/null", vcf, filter_cmd)
} else {
  pipe_prefix <- sprintf("bcftools view -c1 '%s' 2>/dev/null | %s bcftools +fill-tags - -- -t AF 2>/dev/null", vcf, filter_cmd)
}

# Try with AF + TR_MOTIF first (filter out all-homref sites from vg call -A)
cmd_af <- sprintf(
  "%s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%INFO/AF\\t%%INFO/TR_MOTIF\\n' 2>/dev/null",
  pipe_prefix
)
dt <- tryCatch(
  fread(cmd = cmd_af, col.names = c("CHROM", "POS", "REF", "ALT", "AF_str", "TR_MOTIF")),
  error = function(e) NULL,
  warning = function(w) NULL
)

has_af <- !is.null(dt) && nrow(dt) > 0 && !all(is.na(dt$AF_str) | dt$AF_str == ".")

if (is.null(dt) || nrow(dt) == 0) {
  # Fallback: read without AF (still include TR_MOTIF)
  cmd_no_af <- sprintf(
    "%s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%INFO/TR_MOTIF\\n' 2>/dev/null",
    pipe_prefix
  )
  dt <- tryCatch(
    fread(cmd = cmd_no_af, col.names = c("CHROM", "POS", "REF", "ALT", "TR_MOTIF")),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if (is.null(dt) || nrow(dt) == 0) {
    # Final fallback: no AF, no TR_MOTIF
    cmd_bare <- sprintf(
      "%s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' 2>/dev/null",
      pipe_prefix
    )
    dt <- fread(cmd = cmd_bare, col.names = c("CHROM", "POS", "REF", "ALT"))
  }
  has_af <- FALSE
}

if (!has_af && "AF_str" %in% names(dt)) {
  dt[, AF_str := NULL]
}

# Detect TR_MOTIF availability
has_tr <- "TR_MOTIF" %in% names(dt) && !all(is.na(dt$TR_MOTIF) | dt$TR_MOTIF == ".")
if ("TR_MOTIF" %in% names(dt) && !has_tr) {
  dt[, TR_MOTIF := NULL]
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

# Compute max size and direction across comma-separated ALT alleles.
# Fast path: biallelic (no comma in ALT) — the vast majority of rows.
is_multi <- grepl(",", dt$ALT, fixed = TRUE)
dt[, size_signed := 0L]
dt[(!is_multi), size_signed := nchar(ALT) - ref_len]
dt[, size := abs(size_signed)]

# Slow path: multi-allelic only (typically <5% of rows)
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

# NA size can arise from unusual ALT fields (e.g., '*' spanning deletions)
dt[, variant_type := fifelse(
  is.na(size), "Other",
  fifelse(size == 0L & ref_len == 1L, "SNP",
  fifelse(size == 0L, "MNP",
  fifelse(size < 50L & size_signed > 0L, "Insertion",
  fifelse(size < 50L, "Deletion",
  fifelse(size_signed > 0L, "SV Insertion", "SV Deletion")))))
)]

# On-ref vs off-ref
dt[, ref_context := fifelse(
  grepl("_[0-9]+_alt$", CHROM), "Off-reference", "On-reference"
)]

# Drop SV categories if requested
if (no_sv) {
  dt <- dt[!variant_type %in% c("SV Insertion", "SV Deletion")]
}

# Ts/Tv classification (SNPs only — biallelic SNPs where ALT has no comma)
dt[variant_type == "SNP" & !grepl(",", ALT), tstv := {
  transitions <- c("AG", "GA", "CT", "TC")
  fifelse(paste0(REF, ALT) %in% transitions, "Ts", "Tv")
}]

# Tandem repeat flag for indels (from TR_MOTIF annotation)
if (has_tr) {
  cat("TR_MOTIF field detected — flagging tandem repeat indels\n")
  # For multi-allelic (sites mode): any non-"." value among comma-separated motifs
  dt[, is_repeat := FALSE]
  dt[variant_type %in% c("Insertion", "Deletion", "SV Insertion", "SV Deletion"),
     is_repeat := grepl("[^.,]", TR_MOTIF)]
  dt[, TR_MOTIF := NULL]
  cat("Tandem repeat indels:", sum(dt$is_repeat), "of",
      sum(dt$variant_type %in% c("Insertion", "Deletion", "SV Insertion", "SV Deletion")),
      "indel records\n")
}

# ---------------------------------------------------------------------------
# AF: compute per-site non-reference frequency (sum of alt AFs)
# ---------------------------------------------------------------------------
if (has_af) {
  # Ensure AF_str is character (fread may auto-detect as numeric after norm -m-)
  if (!is.character(dt$AF_str)) dt[, AF_str := as.character(AF_str)]
  # Fast path: single AF value (biallelic, no comma) — vast majority
  # Pre-filter "." to avoid millions of as.numeric() warnings
  is_multi_af <- grepl(",", dt$AF_str, fixed = TRUE)
  is_valid_af <- !is_multi_af & dt$AF_str != "."
  dt[, nonref_af := NA_real_]
  dt[(is_valid_af), nonref_af := pmin(as.numeric(AF_str), 1.0)]
  # Slow path: multi-allelic AF (sum comma-separated values)
  if (any(is_multi_af)) {
    dt[is_multi_af, nonref_af := sapply(AF_str, function(x) {
      vals <- as.numeric(unlist(strsplit(x, ",")))
      min(sum(vals, na.rm = TRUE), 1.0)
    })]
  }
  dt[, AF_str := NULL]
}

}  # end of else (VCF reading path)

# ---------------------------------------------------------------------------
# Dump per-record TSV if requested
# ---------------------------------------------------------------------------
if (dump_records) {
  records_path <- paste0(prefix, ".records.tsv")
  cols <- intersect(c("CHROM", "POS", "ref_context", "variant_type",
                       "size", "size_signed", "nonref_af", "is_repeat"), names(dt))
  fwrite(dt[, ..cols], records_path, sep = "\t", nThread = 1)
  cat("Wrote per-record TSV:", records_path, "\n")
  if (records_only) {
    cat("--records-only: skipping plots.\n")
    quit(status = 0)
  }
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
# Plot 1: Variant type bar chart (with tandem repeat overlay when available)
# ---------------------------------------------------------------------------
plot_dt <- dt[, .(count = .N), by = .(ref_context, variant_type)]

# Get Ts/Tv labels for SNP bars
if (nrow(tstv_dt) > 0) {
  tstv_labels <- tstv_dt[, .(ref_context, label = paste0("Ts/Tv=", tstv_ratio))]
  plot_dt <- merge(plot_dt, tstv_labels, by = "ref_context", all.x = TRUE)
  plot_dt[!variant_type %in% "SNP", label := NA_character_]
} else {
  plot_dt[, label := NA_character_]
}

# Order variant types (drop empty levels)
sv_types <- if (no_sv) character(0) else c("SV Insertion", "SV Deletion")
type_levels <- c("SNP", "MNP", "Insertion", "Deletion", sv_types, "Other")
plot_dt[, variant_type := factor(variant_type, levels = intersect(type_levels, unique(variant_type)))]

if (has_tr) {
  # Tandem repeat overlay: solid full bar + faded TR portion on top
  indel_types <- c("Insertion", "Deletion", "SV Insertion", "SV Deletion")
  tr_counts <- dt[variant_type %in% indel_types,
                  .(tr_count = sum(is_repeat)), by = .(ref_context, variant_type)]
  plot_dt <- merge(plot_dt, tr_counts, by = c("ref_context", "variant_type"), all.x = TRUE)
  plot_dt[is.na(tr_count), tr_count := 0L]

  # TR overlay data (only indel types with tr_count > 0)
  tr_dt <- plot_dt[variant_type %in% indel_types & tr_count > 0,
                   .(ref_context, variant_type, count = tr_count)]

  p1 <- ggplot(plot_dt, aes(x = variant_type, y = count, fill = ref_context)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_col(data = tr_dt, position = "dodge", width = 0.7, alpha = 0.35,
             fill = "white") +
    geom_text(aes(label = label),
              position = position_dodge(width = 0.7),
              vjust = -0.5, size = 3, na.rm = TRUE) +
    scale_fill_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                      name = NULL) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = title, subtitle = paste0("Variant Type Counts ", mode_label, filter_label,
                                          " (faded = tandem repeat)"),
         x = "Variant Type", y = "Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
} else {
  p1 <- ggplot(plot_dt, aes(x = variant_type, y = count, fill = ref_context)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_text(aes(label = label),
              position = position_dodge(width = 0.7),
              vjust = -0.5, size = 3, na.rm = TRUE) +
    scale_fill_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                      name = NULL) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = title, subtitle = paste0("Variant Type Counts ", mode_label, filter_label),
         x = "Variant Type", y = "Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
}

save_png(p1, paste0(prefix, ".variant-types.png"))

# ---------------------------------------------------------------------------
# Plot 2: Size distribution — two-panel (Small 1-49 bp / Structural 50-1000 bp)
# ---------------------------------------------------------------------------
size_dt <- dt[variant_type %in% c("Insertion", "Deletion", "SV Insertion", "SV Deletion") & size > 0]

if (nrow(size_dt) > 0) {
  size_dt[, direction := fifelse(size_signed > 0L, "Insertion", "Deletion")]
  size_dt[, panel := fifelse(size < 50L, "Small (1-49 bp)", "Structural (50-1000 bp)")]
  size_dt[, panel := factor(panel, levels = c("Small (1-49 bp)", "Structural (50-1000 bp)"))]
  size_counts <- size_dt[size <= 1000, .(count = .N), by = .(size, direction, ref_context, panel)]

  p2 <- ggplot(size_counts, aes(x = size, y = count, color = direction, linetype = ref_context)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~panel, scales = "free") +
    scale_color_manual(values = c("Insertion" = "coral", "Deletion" = "steelblue"),
                       name = "Direction") +
    scale_linetype_manual(values = c("On-reference" = "solid", "Off-reference" = "dashed"),
                          name = "Context") +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
    labs(title = title, subtitle = paste0("Indel / SV Size Distribution ", mode_label, filter_label),
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

  p2log <- p2 +
    scale_y_log10(labels = scales::comma) +
    labs(y = "Count (log scale)")
  save_png(p2log, paste0(prefix, ".size-dist-log.png"), width = 12)
} else {
  cat("No indels/SVs with size > 0; skipping size distribution plot.\n")
  # Create empty files so Snakemake sees the outputs
  file.create(paste0(prefix, ".size-dist.png"))
  file.create(paste0(prefix, ".size-dist-log.png"))
}

# ---------------------------------------------------------------------------
# Plot 3: AF spectrum — per-site non-reference frequency (log y-axis)
#
# AF values are inherently discrete (multiples of 1/(2N) for N diploid samples),
# so we always count exact values rather than binning.  When there are many
# unique values (large N), we round to ~100 evenly-spaced bins first.
# ---------------------------------------------------------------------------
if (has_af) {
  cat("AF spectrum: ", nrow(dt), " sites with non-ref AF\n")

  n_raw <- length(unique(dt$nonref_af))
  cat("Unique raw AF values:", n_raw, "\n")

  # Round AF to nearest step (default 0.1 = 10%)
  dt[, af_plot := round(nonref_af / af_step) * af_step]
  n_unique <- length(unique(dt$af_plot))
  cat("Unique AF values for plot:", n_unique, "\n")

  af_counts <- dt[, .(count = .N), by = .(af_plot, ref_context)]
  af_unique <- sort(unique(af_counts$af_plot))

  y_label <- if (mode == "sites") "Sites (log scale)" else "Variants (log scale)"
  common_theme <- theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )

  p3 <- ggplot(af_counts, aes(x = af_plot, y = count, color = ref_context)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.2) +
    scale_color_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                       name = NULL) +
    scale_y_log10(labels = scales::comma) +
    scale_x_continuous(limits = c(-0.02, 1.02)) +
    labs(title = title, subtitle = paste0("Non-Reference Allele Frequency Spectrum ", mode_label, filter_label),
         x = "Non-Reference Frequency", y = y_label) +
    common_theme

  dt[, af_plot := NULL]
  save_png(p3, paste0(prefix, ".af-spectrum.png"), height = 6)
  cat("AF spectrum plot generated.\n")
} else {
  cat("No AF field; skipping allele frequency spectrum plot.\n")
  # Create empty file so Snakemake sees the output
  file.create(paste0(prefix, ".af-spectrum.png"))
}

# ---------------------------------------------------------------------------
# Plot 4: Annotation-stratified variant type counts (when --annot provided)
# ---------------------------------------------------------------------------
if (!is.null(annot_beds)) {
  bed_files <- strsplit(annot_beds, ",")[[1]]
  anames   <- strsplit(annot_names_arg, ",")[[1]]

  cat("Annotating variants from", length(bed_files), "augref BED files (bedtools)\n")

  # Check bedtools is available
  if (system("bedtools --version >/dev/null 2>&1") != 0) {
    stop("bedtools is required for --annot-beds but is not found in PATH")
  }

  # Write variant positions as sorted BED for bedtools intersect -sorted.
  # This avoids loading annotation BEDs into R memory entirely.
  tmp_unsorted <- tempfile(fileext = ".bed")
  fwrite(dt[, .(CHROM, POS - 1L, POS - 1L + ref_len)],
         tmp_unsorted, sep = "\t", col.names = FALSE)
  tmp_bed <- tempfile(fileext = ".sorted.bed")
  system(sprintf("LC_ALL=C sort -k1,1 -k2,2n '%s' > '%s'", tmp_unsorted, tmp_bed))
  unlink(tmp_unsorted)

  # For each annotation: bedtools intersect → hit positions → flag in dt
  dt[, bed_start := POS - 1L]
  setkey(dt, CHROM, bed_start)
  for (k in seq_along(bed_files)) {
    cat("  bedtools intersect:", anames[k], "\n")
    hit_col <- paste0(anames[k], "_hit")
    cmd <- sprintf(
      "bedtools intersect -a '%s' -b '%s' -u -sorted 2>/dev/null",
      tmp_bed, bed_files[k]
    )
    hits <- tryCatch(
      fread(cmd = cmd, select = 1:2, col.names = c("CHROM", "bed_start"), header = FALSE),
      error = function(e) data.table(CHROM = character(0), bed_start = integer(0)),
      warning = function(w) data.table(CHROM = character(0), bed_start = integer(0))
    )
    set(dt, j = hit_col, value = FALSE)
    if (nrow(hits) > 0) {
      hits <- unique(hits)
      dt[hits, (hit_col) := TRUE, on = .(CHROM, bed_start)]
    }
  }
  dt[, bed_start := NULL]
  unlink(tmp_bed)

  # Accumulate annotation counts per annotation (avoids melt on full table)
  ann_counts_list <- vector("list", length(anames))
  for (k in seq_along(anames)) {
    col <- paste0(anames[k], "_hit")
    sub_counts <- dt[get(col) == TRUE, .(count = .N), by = .(variant_type, ref_context)]
    if (nrow(sub_counts) > 0) sub_counts[, annotation := anames[k]]
    ann_counts_list[[k]] <- sub_counts
  }
  ann_counts <- rbindlist(ann_counts_list, use.names = TRUE, fill = TRUE)

  if (nrow(ann_counts) > 0) {
    ann_counts[, variant_type := factor(variant_type, levels = type_levels)]

    p_annot <- ggplot(ann_counts, aes(x = variant_type, y = count, fill = ref_context)) +
      geom_col(position = "dodge", width = 0.7) +
      facet_wrap(~annotation, scales = "free_y") +
      scale_fill_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                        name = NULL) +
      scale_y_continuous(labels = scales::comma) +
      labs(title = title, subtitle = "Variant Types by Annotation Region") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA)
      ) +
      coord_flip()

    save_png(p_annot, paste0(prefix, ".variant-types-by-annot.png"), width = 12, height = 8)

    # Write annotation-stratified summary TSV
    setorder(ann_counts, annotation, ref_context, variant_type)
    fwrite(ann_counts, paste0(prefix, ".vcf-stats-by-annot.tsv"), sep = "\t")
    cat("Wrote annotation-stratified stats:", paste0(prefix, ".vcf-stats-by-annot.tsv"), "\n")
  } else {
    cat("No variants in annotated regions; creating empty annotation plot.\n")
    file.create(paste0(prefix, ".variant-types-by-annot.png"))
    fwrite(data.table(annotation = character(), ref_context = character(),
                      variant_type = character(), count = integer()),
           paste0(prefix, ".vcf-stats-by-annot.tsv"), sep = "\t")
  }
}

# ---------------------------------------------------------------------------
# Plot 5: GIAB genome stratification (standalone, when --giab-strat-beds)
# ---------------------------------------------------------------------------
if (!is.null(giab_strat_beds_arg)) {
  strat_files  <- strsplit(giab_strat_beds_arg, ",")[[1]]
  strat_names  <- gsub("_", " ", strsplit(giab_strat_names_arg, ",")[[1]])

  cat("GIAB stratification from", length(strat_files), "BEDs\n")

  # Write variant BED
  tmp_bed <- tempfile(fileext = ".sorted.bed")
  tmp_unsorted_giab <- tempfile(fileext = ".bed")
  fwrite(dt[, .(CHROM, POS - 1L, POS - 1L + ref_len)],
         tmp_unsorted_giab, sep = "\t", col.names = FALSE)
  system(sprintf("LC_ALL=C sort -k1,1 -k2,2n '%s' > '%s'", tmp_unsorted_giab, tmp_bed))
  unlink(tmp_unsorted_giab)

  # bedtools intersect each partition
  dt[, giab_bed_start := POS - 1L]
  setkey(dt, CHROM, giab_bed_start)
  for (k in seq_along(strat_files)) {
    col <- paste0("giab_", k)
    cmd <- sprintf("bedtools intersect -a '%s' -b '%s' -u -sorted 2>/dev/null",
                   tmp_bed, strat_files[k])
    hits <- tryCatch(
      fread(cmd = cmd, select = 1:2, col.names = c("CHROM", "giab_bed_start"), header = FALSE),
      error = function(e) data.table(CHROM = character(0), giab_bed_start = integer(0)),
      warning = function(w) data.table(CHROM = character(0), giab_bed_start = integer(0))
    )
    set(dt, j = col, value = FALSE)
    if (nrow(hits) > 0) {
      hits <- unique(hits)
      dt[hits, (col) := TRUE, on = .(CHROM, giab_bed_start)]
    }
  }
  dt[, giab_bed_start := NULL]
  unlink(tmp_bed)

  # Classify each variant into exactly one category (priority: first hit wins)
  # Off-ref variants get "Off-reference" since GIAB BEDs are reference-only
  dt[, giab_region := "Unclassified"]
  dt[ref_context == "Off-reference", giab_region := "Off-reference"]
  for (k in rev(seq_along(strat_names))) {
    col <- paste0("giab_", k)
    dt[get(col) == TRUE, giab_region := strat_names[k]]
  }

  # Summary table
  giab_counts <- dt[, .(count = .N), by = .(variant_type, giab_region)]
  giab_counts[, variant_type := factor(variant_type, levels = type_levels)]
  region_levels <- c(strat_names, "Off-reference")
  giab_counts[, giab_region := factor(giab_region, levels = region_levels)]

  # Standalone grouped bar chart
  region_colors <- c("Easy" = "forestgreen", "Segdup" = "firebrick",
                     "Other Difficult" = "darkorange", "Off-reference" = "steelblue")
  p_giab <- ggplot(giab_counts[giab_region != "Unclassified"],
                   aes(x = variant_type, y = count, fill = giab_region)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = region_colors, name = "GIAB Region") +
    scale_y_continuous(labels = scales::comma) +
    labs(title = title,
         subtitle = paste0("GIAB Genome Stratification ", mode_label, filter_label),
         x = "Variant Type", y = "Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )

  save_png(p_giab, paste0(prefix, ".giab-strat.png"), width = 10, height = 6)

  # Write TSV
  setorder(giab_counts, giab_region, variant_type)
  fwrite(giab_counts, paste0(prefix, ".giab-strat.tsv"), sep = "\t")
  cat("Wrote GIAB strat:", paste0(prefix, ".giab-strat.tsv"), "\n")

  # Clean up columns
  for (k in seq_along(strat_names)) set(dt, j = paste0("giab_", k), value = NULL)
  dt[, giab_region := NULL]
}

# ---------------------------------------------------------------------------
# Per-sample variant stats (when --per-sample)
# ---------------------------------------------------------------------------
if (per_sample) {
  cat("Per-sample mode enabled\n")

  # 1. Read sample names
  sample_cmd <- sprintf("bcftools query -l '%s' 2>/dev/null", vcf)
  all_sample_names <- system(sample_cmd, intern = TRUE)
  n_all_samples <- length(all_sample_names)
  # Exclude reference sample from per-sample analysis (it always has 0 variants)
  if (!is.null(ref_sample) && ref_sample %in% all_sample_names) {
    sample_names <- setdiff(all_sample_names, ref_sample)
    cat("Found", n_all_samples, "samples, excluding reference sample:", ref_sample, "\n")
  } else {
    sample_names <- all_sample_names
  }
  n_samples <- length(sample_names)
  cat("Per-sample analysis on", n_samples, "samples:", paste(sample_names, collapse = ", "), "\n")

  if (n_samples == 0) {
    cat("No samples found in VCF; skipping per-sample stats.\n")
    file.create(paste0(prefix, ".per-sample-types.png"))
    fwrite(data.table(sample = character(), variant_type = character(),
                      ref_context = character(), count = integer()),
           paste0(prefix, ".per-sample-types.tsv"), sep = "\t")
    if (!is.null(giab_strat_beds_arg)) {
      file.create(paste0(prefix, ".per-sample-giab-strat.png"))
      fwrite(data.table(sample = character(), variant_type = character(),
                        ps_giab_region = character(), count = integer()),
             paste0(prefix, ".per-sample-giab-strat.tsv"), sep = "\t")
    }
  } else {
    # 2. Read per-sample GTs through same pipeline as site-level
    #    (bcftools outputs all samples; we filter to non-ref samples after melting)
    gt_cmd <- sprintf(
      "%s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT]\\n' 2>/dev/null",
      pipe_prefix
    )
    gt_cols <- c("CHROM", "POS", "REF", "ALT", paste0("GT_", seq_len(n_all_samples)))
    gt_dt <- fread(cmd = gt_cmd, col.names = gt_cols)
    cat("Read", nrow(gt_dt), "variant records with GTs\n")

    if (nrow(gt_dt) > 0) {
      # 3. Classify variants (same logic as site-level)
      gt_dt[, ref_len := nchar(REF)]
      is_multi_gt <- grepl(",", gt_dt$ALT, fixed = TRUE)
      gt_dt[, size_signed := 0L]
      gt_dt[(!is_multi_gt), size_signed := nchar(ALT) - ref_len]
      gt_dt[, size := abs(size_signed)]
      multi_idx_gt <- which(is_multi_gt)
      if (length(multi_idx_gt) > 0) {
        gt_dt[multi_idx_gt, c("size", "size_signed") := {
          res <- sapply(seq_len(.N), function(i) {
            alts <- unlist(strsplit(ALT[i], ","))
            diffs <- nchar(alts) - ref_len[i]
            idx <- which.max(abs(diffs))
            c(abs(diffs[idx]), diffs[idx])
          })
          list(res[1,], res[2,])
        }]
      }
      gt_dt[, variant_type := fifelse(
        is.na(size), "Other",
        fifelse(size == 0L & ref_len == 1L, "SNP",
        fifelse(size == 0L, "MNP",
        fifelse(size < 50L & size_signed > 0L, "Insertion",
        fifelse(size < 50L, "Deletion",
        fifelse(size_signed > 0L, "SV Insertion", "SV Deletion")))))
      )]
      gt_dt[, ref_context := fifelse(
        grepl("_[0-9]+_alt$", CHROM), "Off-reference", "On-reference"
      )]

      # Drop SV categories if requested
      if (no_sv) {
        gt_dt <- gt_dt[!variant_type %in% c("SV Insertion", "SV Deletion")]
      }

      # 4. Count carriers per sample without melting (avoids exceeding R's 2^31
      #    vector limit when n_variants × n_samples is very large).
      gt_sample_cols <- paste0("GT_", seq_len(n_all_samples))
      ps_list <- vector("list", n_samples)
      total_carriers <- 0L
      for (si in seq_along(sample_names)) {
        sname <- sample_names[si]
        col_idx <- match(sname, all_sample_names)
        gt_col <- gt_sample_cols[col_idx]
        carriers <- grepl("[1-9]", gt_dt[[gt_col]])
        total_carriers <- total_carriers + sum(carriers)
        ps_list[[si]] <- gt_dt[carriers, .(count = .N), by = .(variant_type, ref_context)
                               ][, sample := sname]
      }
      ps_counts <- rbindlist(ps_list)
      cat("Carrier genotype rows:", total_carriers, "\n")

      # Ensure all sample × type × context combinations exist (fill with 0)
      all_combos <- CJ(sample = sample_names,
                        variant_type = unique(ps_counts$variant_type),
                        ref_context = unique(ps_counts$ref_context))
      ps_counts <- merge(all_combos, ps_counts,
                          by = c("sample", "variant_type", "ref_context"), all.x = TRUE)
      ps_counts[is.na(count), count := 0L]

      # Summary stats across samples
      ps_summary <- ps_counts[, .(mean_count = mean(count),
                                   min_count = min(count),
                                   max_count = max(count),
                                   sd_count = sd(count)),
                               by = .(variant_type, ref_context)]

      # Write per-sample TSV
      setorder(ps_counts, sample, ref_context, variant_type)
      fwrite(ps_counts, paste0(prefix, ".per-sample-types.tsv"), sep = "\t")
      cat("Wrote per-sample stats:", paste0(prefix, ".per-sample-types.tsv"), "\n")

      # 6. Plot
      ps_counts[, variant_type := factor(variant_type,
        levels = intersect(type_levels, unique(variant_type)))]
      ps_summary[, variant_type := factor(variant_type,
        levels = intersect(type_levels, unique(variant_type)))]

      if (n_samples <= 20) {
        # Bar at mean + jittered dots + min/max error bars
        p_ps <- ggplot() +
          geom_col(data = ps_summary,
                   aes(x = variant_type, y = mean_count, fill = ref_context),
                   position = position_dodge(width = 0.7), width = 0.7, alpha = 0.6) +
          geom_errorbar(data = ps_summary,
                        aes(x = variant_type, ymin = min_count, ymax = max_count,
                            group = ref_context),
                        position = position_dodge(width = 0.7), width = 0.3) +
          geom_point(data = ps_counts,
                     aes(x = variant_type, y = count, color = ref_context),
                     position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
                     size = 1.5, alpha = 0.8) +
          scale_fill_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                            name = NULL) +
          scale_color_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                             name = NULL) +
          scale_y_continuous(labels = scales::comma) +
          labs(title = title,
               subtitle = paste0("Per-Sample Variant Counts ", mode_label, filter_label,
                                 " (N=", n_samples, " samples, bars=mean)"),
               x = "Variant Type", y = "Count") +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background  = element_rect(fill = "white", color = NA)
          )
      } else {
        # Boxplot for large sample counts
        p_ps <- ggplot(ps_counts,
                       aes(x = variant_type, y = count, fill = ref_context)) +
          geom_boxplot(position = position_dodge(width = 0.7), width = 0.6,
                       outlier.size = 1) +
          scale_fill_manual(values = c("On-reference" = "steelblue", "Off-reference" = "coral"),
                            name = NULL) +
          scale_y_continuous(labels = scales::comma) +
          labs(title = title,
               subtitle = paste0("Per-Sample Variant Counts ", mode_label, filter_label,
                                 " (N=", n_samples, " samples)"),
               x = "Variant Type", y = "Count") +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background  = element_rect(fill = "white", color = NA)
          )
      }

      save_png(p_ps, paste0(prefix, ".per-sample-types.png"))

      # -----------------------------------------------------------------------
      # Per-sample GIAB stratification (when --giab-strat-beds + --per-sample)
      # -----------------------------------------------------------------------
      if (!is.null(giab_strat_beds_arg)) {
        strat_files_ps  <- strsplit(giab_strat_beds_arg, ",")[[1]]
        strat_names_ps  <- gsub("_", " ", strsplit(giab_strat_names_arg, ",")[[1]])

        cat("Per-sample GIAB stratification from", length(strat_files_ps), "BEDs\n")

        # Write variant BED from gt_dt
        tmp_bed_ps <- tempfile(fileext = ".sorted.bed")
        tmp_unsorted_ps <- tempfile(fileext = ".bed")
        fwrite(gt_dt[, .(CHROM, POS - 1L, POS - 1L + ref_len)],
               tmp_unsorted_ps, sep = "\t", col.names = FALSE)
        system(sprintf("LC_ALL=C sort -k1,1 -k2,2n '%s' > '%s'", tmp_unsorted_ps, tmp_bed_ps))
        unlink(tmp_unsorted_ps)

        # bedtools intersect each partition
        gt_dt[, ps_bed_start := POS - 1L]
        setkey(gt_dt, CHROM, ps_bed_start)
        for (k in seq_along(strat_files_ps)) {
          col <- paste0("ps_giab_", k)
          cmd <- sprintf("bedtools intersect -a '%s' -b '%s' -u -sorted 2>/dev/null",
                         tmp_bed_ps, strat_files_ps[k])
          hits <- tryCatch(
            fread(cmd = cmd, select = 1:2, col.names = c("CHROM", "ps_bed_start"), header = FALSE),
            error = function(e) data.table(CHROM = character(0), ps_bed_start = integer(0)),
            warning = function(w) data.table(CHROM = character(0), ps_bed_start = integer(0))
          )
          set(gt_dt, j = col, value = FALSE)
          if (nrow(hits) > 0) {
            hits <- unique(hits)
            gt_dt[hits, (col) := TRUE, on = .(CHROM, ps_bed_start)]
          }
        }
        gt_dt[, ps_bed_start := NULL]
        unlink(tmp_bed_ps)

        # Classify each variant into one GIAB region
        gt_dt[, ps_giab_region := "Unclassified"]
        gt_dt[ref_context == "Off-reference", ps_giab_region := "Off-reference"]
        for (k in rev(seq_along(strat_names_ps))) {
          col <- paste0("ps_giab_", k)
          gt_dt[get(col) == TRUE, ps_giab_region := strat_names_ps[k]]
        }

        # Count per sample × variant_type × giab_region (loop to avoid melt)
        ps_giab_list <- vector("list", n_samples)
        for (si in seq_along(sample_names)) {
          sname <- sample_names[si]
          col_idx <- match(sname, all_sample_names)
          gt_col <- gt_sample_cols[col_idx]
          carriers <- grepl("[1-9]", gt_dt[[gt_col]])
          ps_giab_list[[si]] <- gt_dt[carriers, .(count = .N),
                                       by = .(variant_type, ps_giab_region)
                                       ][, sample := sname]
        }
        ps_giab_counts <- rbindlist(ps_giab_list)

        # Ensure all combos exist
        region_levels_ps <- c(strat_names_ps, "Off-reference")
        all_giab_combos <- CJ(sample = sample_names,
                               variant_type = unique(ps_giab_counts$variant_type),
                               ps_giab_region = region_levels_ps)
        ps_giab_counts <- merge(all_giab_combos, ps_giab_counts,
                                 by = c("sample", "variant_type", "ps_giab_region"), all.x = TRUE)
        ps_giab_counts[is.na(count), count := 0L]

        # Write TSV
        setorder(ps_giab_counts, sample, ps_giab_region, variant_type)
        fwrite(ps_giab_counts, paste0(prefix, ".per-sample-giab-strat.tsv"), sep = "\t")
        cat("Wrote per-sample GIAB strat:", paste0(prefix, ".per-sample-giab-strat.tsv"), "\n")

        # Summary for plotting
        ps_giab_summary <- ps_giab_counts[, .(mean_count = mean(count),
                                               min_count = min(count),
                                               max_count = max(count)),
                                           by = .(variant_type, ps_giab_region)]

        ps_giab_counts[, variant_type := factor(variant_type,
          levels = intersect(type_levels, unique(variant_type)))]
        ps_giab_summary[, variant_type := factor(variant_type,
          levels = intersect(type_levels, unique(variant_type)))]
        ps_giab_counts[, ps_giab_region := factor(ps_giab_region, levels = region_levels_ps)]
        ps_giab_summary[, ps_giab_region := factor(ps_giab_region, levels = region_levels_ps)]

        # Filter out Unclassified
        ps_giab_counts_plot <- ps_giab_counts[ps_giab_region != "Unclassified"]
        ps_giab_summary_plot <- ps_giab_summary[ps_giab_region != "Unclassified"]

        region_colors_ps <- c("Easy" = "forestgreen", "Segdup" = "firebrick",
                               "Other Difficult" = "darkorange", "Off-reference" = "steelblue")

        if (n_samples <= 20) {
          p_ps_giab <- ggplot() +
            geom_col(data = ps_giab_summary_plot,
                     aes(x = variant_type, y = mean_count, fill = ps_giab_region),
                     position = position_dodge(width = 0.7), width = 0.7, alpha = 0.6) +
            geom_errorbar(data = ps_giab_summary_plot,
                          aes(x = variant_type, ymin = min_count, ymax = max_count,
                              group = ps_giab_region),
                          position = position_dodge(width = 0.7), width = 0.3) +
            geom_point(data = ps_giab_counts_plot,
                       aes(x = variant_type, y = count, color = ps_giab_region),
                       position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
                       size = 1.5, alpha = 0.8) +
            scale_fill_manual(values = region_colors_ps, name = "GIAB Region") +
            scale_color_manual(values = region_colors_ps, name = "GIAB Region") +
            scale_y_continuous(labels = scales::comma) +
            labs(title = title,
                 subtitle = paste0("Per-Sample GIAB Stratification ", mode_label, filter_label,
                                   " (N=", n_samples, " samples, bars=mean)"),
                 x = "Variant Type", y = "Count") +
            theme_minimal() +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5),
              panel.background = element_rect(fill = "white", color = NA),
              plot.background  = element_rect(fill = "white", color = NA)
            )
        } else {
          p_ps_giab <- ggplot(ps_giab_counts_plot,
                              aes(x = variant_type, y = count, fill = ps_giab_region)) +
            geom_boxplot(position = position_dodge(width = 0.7), width = 0.6,
                         outlier.size = 1) +
            scale_fill_manual(values = region_colors_ps, name = "GIAB Region") +
            scale_y_continuous(labels = scales::comma) +
            labs(title = title,
                 subtitle = paste0("Per-Sample GIAB Stratification ", mode_label, filter_label,
                                   " (N=", n_samples, " samples)"),
                 x = "Variant Type", y = "Count") +
            theme_minimal() +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5),
              panel.background = element_rect(fill = "white", color = NA),
              plot.background  = element_rect(fill = "white", color = NA)
            )
        }

        save_png(p_ps_giab, paste0(prefix, ".per-sample-giab-strat.png"), width = 10)

        # Clean up
        for (k in seq_along(strat_names_ps)) set(gt_dt, j = paste0("ps_giab_", k), value = NULL)
        gt_dt[, ps_giab_region := NULL]
      }
    } else {
      cat("No variant records; creating empty per-sample outputs.\n")
      file.create(paste0(prefix, ".per-sample-types.png"))
      fwrite(data.table(sample = character(), variant_type = character(),
                        ref_context = character(), count = integer()),
             paste0(prefix, ".per-sample-types.tsv"), sep = "\t")
      if (!is.null(giab_strat_beds_arg)) {
        file.create(paste0(prefix, ".per-sample-giab-strat.png"))
        fwrite(data.table(sample = character(), variant_type = character(),
                          ps_giab_region = character(), count = integer()),
               paste0(prefix, ".per-sample-giab-strat.tsv"), sep = "\t")
      }
    }
  }
}

cat("Done.\n")
