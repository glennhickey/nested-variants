#!/usr/bin/env Rscript
# vcf-stats-combine.R — Combine per-chromosome vcf-stats TSVs into
# aggregate and per-chromosome summary plots.
#
# Usage:
#   Rscript vcf-stats-combine.R <output_prefix> \
#     [--per-sample chrom1:per-sample.tsv chrom2:per-sample.tsv ...] \
#     <chrom1:stats.tsv> [chrom2:stats.tsv] ...
#
# Produces:
#   <prefix>.variant-types.png          — aggregate variant type counts
#   <prefix>.by-chrom.variant-types.png — variant types grouped by chromosome
#   <prefix>.per-sample-types.png       — per-sample box plots (if --per-sample)
#   <prefix>.vcf-stats.tsv              — aggregate stats TSV

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# When EMIT_SVG=1, every ggsave(...png) also produces a sibling .svg
local({
  if (!nzchar(Sys.getenv("EMIT_SVG"))) return(invisible(NULL))
  .ggsave <- ggplot2::ggsave
  assign("ggsave", function(filename, plot = ggplot2::last_plot(), ...) {
    .ggsave(filename, plot = plot, ...)
    if (grepl("\\.png$", filename)) {
      pdf_file <- sub("\\.png$", ".pdf", filename)
      args <- list(...); args$device <- grDevices::cairo_pdf
      args$type <- NULL; args$dpi <- NULL
      tryCatch(do.call(.ggsave, c(list(pdf_file, plot), args)),
               error = function(e) {
                 cat("PDF emit failed for ", pdf_file, ": ",
                     conditionMessage(e), "\n", sep = "")
                 if (file.exists(pdf_file) && file.info(pdf_file)$size == 0)
                   file.remove(pdf_file)
               })
    }
  }, envir = .GlobalEnv)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript vcf-stats-combine.R <output_prefix> [--per-sample chrom:file ...] <chrom:file> ...\n")
  quit(status = 1)
}

output_prefix <- args[1]
args <- args[-1]

# Parse --per-sample flag: collects chrom:path pairs until next non-pair arg
per_sample_pairs <- c()
stats_pairs <- c()
i <- 1
in_per_sample <- FALSE
while (i <= length(args)) {
  if (args[i] == "--per-sample") {
    in_per_sample <- TRUE
    i <- i + 1
  } else if (in_per_sample && grepl(":", args[i])) {
    per_sample_pairs <- c(per_sample_pairs, args[i])
    i <- i + 1
  } else {
    in_per_sample <- FALSE
    stats_pairs <- c(stats_pairs, args[i])
    i <- i + 1
  }
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

# Read all per-chromosome stats (may be empty if only --per-sample)
all_stats <- if (length(stats_pairs) > 0) read_pairs(stats_pairs) else NULL

# Variant type ordering
type_order <- c("SNP", "MNP", "Insertion", "Deletion", "SV Insertion", "SV Deletion")

if (!is.null(all_stats)) {

# Natural sort chromosomes: numeric first (chr1..chr22), then X, Y, M, other
chrom_order <- unique(all_stats$chrom)
.sfx <- sub("^chr", "", chrom_order)
.key <- suppressWarnings(as.numeric(.sfx))
.key[.sfx == "X"] <- 100
.key[.sfx == "Y"] <- 101
.key[.sfx == "M"] <- 102
.key[is.na(.key)] <- 1000
chrom_order <- chrom_order[order(.key, chrom_order)]
all_stats[, chrom := factor(chrom, levels = chrom_order)]

# --- Aggregate variant types ---
agg <- all_stats[, .(count = sum(count),
                      ts = sum(ts, na.rm = TRUE),
                      tv = sum(tv, na.rm = TRUE)),
                 by = .(ref_context, variant_type)]
agg[, tstv_ratio := ifelse(tv > 0, round(ts / tv, 2), NA)]
agg[, variant_type := factor(variant_type, levels = type_order)]

p_agg <- ggplot(agg, aes(x = variant_type, y = count, fill = ref_context)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = comma(count)), position = position_dodge(0.9),
            vjust = -0.3, size = 3) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c("Off-reference" = "#E74C3C", "On-reference" = "#3498DB")) +
  labs(title = "Variant Type Counts (all chromosomes)",
       x = "Variant Type", y = "Count", fill = NULL) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1))

snp_row <- agg[variant_type == "SNP" & ref_context == "Off-reference"]
if (nrow(snp_row) > 0 && !is.na(snp_row$tstv_ratio)) {
  p_agg <- p_agg + annotate("text", x = "SNP", y = max(agg$count, na.rm = TRUE) * 0.95,
                             label = paste0("Ts/Tv=", snp_row$tstv_ratio), size = 3.5)
}

ggsave(paste0(output_prefix, ".variant-types.png"), p_agg,
       width = 8, height = 5, dpi = 300,
       device = grDevices::png, type = "cairo")
cat("Saved:", paste0(output_prefix, ".variant-types.png"), "\n")

# --- Per-chromosome variant types ---
by_chrom <- all_stats[, .(count = sum(count)), by = .(chrom, ref_context, variant_type)]
by_chrom[, variant_type := factor(variant_type, levels = type_order)]

p_chrom <- ggplot(by_chrom[ref_context == "Off-reference"],
                  aes(x = variant_type, y = count, fill = chrom)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Off-Reference Variant Types by Chromosome",
       x = "Variant Type", y = "Count", fill = "Chromosome") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(paste0(output_prefix, ".by-chrom.variant-types.png"), p_chrom,
       width = 8, height = 5, dpi = 300,
       device = grDevices::png, type = "cairo")
cat("Saved:", paste0(output_prefix, ".by-chrom.variant-types.png"), "\n")

# --- Write aggregate TSV ---
fwrite(agg, paste0(output_prefix, ".vcf-stats.tsv"), sep = "\t")
cat("Saved:", paste0(output_prefix, ".vcf-stats.tsv"), "\n")

} # end if (!is.null(all_stats))

# --- Per-sample box plots (if --per-sample provided) ---
if (length(per_sample_pairs) > 0) {
  ps <- read_pairs(per_sample_pairs)

  # Sum counts for samples appearing in multiple chromosomes
  ps_agg <- ps[, .(count = sum(count)), by = .(sample, variant_type, ref_context)]
  ps_agg[, variant_type := factor(variant_type, levels = type_order)]

  # Box + jittered per-sample points, split SNP from the rest since SNP
  # counts are ~20x larger and compress the other categories flat.
  ps_off <- ps_agg[ref_context == "Off-reference"]
  ps_off[, facet_group := ifelse(variant_type == "SNP",
                                 "SNP", "Indel / MNP / SV")]
  ps_off[, facet_group := factor(facet_group,
                                 levels = c("SNP", "Indel / MNP / SV"))]

  p_ps <- ggplot(ps_off, aes(x = variant_type, y = count)) +
    geom_boxplot(fill = "#E74C3C", alpha = 0.3, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.6, alpha = 0.35, colour = "#7B241C") +
    scale_y_continuous(labels = comma) +
    facet_wrap(~ facet_group, scales = "free") +
    labs(title = "Per-Sample Variant Counts (all chromosomes)",
         subtitle = paste0("Off-ref, ", length(unique(ps_agg$sample)), " samples"),
         x = "Variant Type", y = "Count") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(paste0(output_prefix, ".per-sample-types.png"), p_ps,
         width = 8, height = 5, dpi = 300,
         device = grDevices::png, type = "cairo")
  cat("Saved:", paste0(output_prefix, ".per-sample-types.png"), "\n")
}
