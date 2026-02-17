#!/usr/bin/env Rscript

# segment-polymorphism.R — Per-segment polymorphism table from deconstruct VCF
#
# Usage: Rscript scripts/segment-polymorphism.R \
#          --vcf <deconstruct.vcf.gz> \
#          --augref-prefix <prefix> \
#          [--annot <annot-per-segment.tsv>] \
#          --output <output.tsv>

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

vcf_file      <- NULL
augref_prefix <- NULL
annot_file    <- NULL
output_file   <- NULL

i <- 1
while (i <= length(args)) {
  if (args[i] == "--vcf" && i + 1 <= length(args)) {
    vcf_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--augref-prefix" && i + 1 <= length(args)) {
    augref_prefix <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--annot" && i + 1 <= length(args)) {
    annot_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--output" && i + 1 <= length(args)) {
    output_file <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(vcf_file) || is.null(augref_prefix) || is.null(output_file)) {
  cat("Usage: Rscript segment-polymorphism.R --vcf <vcf> --augref-prefix <prefix> --output <out.tsv>\n")
  quit(status = 1)
}

# ---------------------------------------------------------------------------
# 1. Parse segment lengths from VCF header
# ---------------------------------------------------------------------------
cat("Parsing contig lengths from VCF header...\n")
header_cmd <- sprintf("bcftools view -h '%s' | grep '^##contig'", vcf_file)
header_lines <- system(header_cmd, intern = TRUE)

contigs <- data.table(
  augref_path    = sub(".*ID=([^,>]+).*", "\\1", header_lines),
  segment_length = as.integer(sub(".*length=([0-9]+).*", "\\1", header_lines))
)
contigs <- contigs[grepl("_[0-9]+_alt$", augref_path)]
contigs[, short_name := sub(paste0("^", augref_prefix), "", augref_path)]
cat(sprintf("  %d alt segments in header\n", nrow(contigs)))

# ---------------------------------------------------------------------------
# 2. Read off-ref variants from VCF (streaming via bcftools query)
# ---------------------------------------------------------------------------
cat("Reading variants from VCF...\n")
cmd <- sprintf(
  "bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%INFO/AF\\t%%INFO/NS\\t%%INFO/AN\\n' '%s'",
  vcf_file)
dt <- fread(cmd = cmd, col.names = c("CHROM", "POS", "REF", "ALT", "AF_str", "NS", "AN"))
cat(sprintf("  %d total variants\n", nrow(dt)))

# Keep only off-ref (alt segment) variants
dt <- dt[grepl("_[0-9]+_alt$", CHROM)]
dt[, short_name := sub(paste0("^", augref_prefix), "", CHROM)]
cat(sprintf("  %d off-ref variants\n", nrow(dt)))

# ---------------------------------------------------------------------------
# 3. Classify variants
# ---------------------------------------------------------------------------
cat("Classifying variants...\n")

# Parse AF (take first value for multi-allelic)
dt[, af := as.numeric(sub(",.*", "", AF_str))]

# Polymorphic = AF >= 5% and AF <= 95% (present in some but not all haplotypes)
dt[, polymorphic := af >= 0.05 & af <= 0.95]

# Variant type
dt[, var_type := fifelse(nchar(REF) == 1 & nchar(ALT) == 1, "SNP",
                  fifelse(abs(nchar(ALT) - nchar(REF)) < 50, "Indel", "SV"))]

# Ts/Tv for SNPs
transitions <- c("AG", "GA", "CT", "TC")
dt[var_type == "SNP", tstv := fifelse(paste0(REF, ALT) %in% transitions, "Ts", "Tv")]

# ---------------------------------------------------------------------------
# 4. Aggregate per segment
# ---------------------------------------------------------------------------
cat("Aggregating per segment...\n")
seg_stats <- dt[, .(
  total_variants       = .N,
  polymorphic_variants = sum(polymorphic),
  snp_count            = sum(var_type == "SNP"),
  indel_count          = sum(var_type == "Indel"),
  sv_count             = sum(var_type == "SV"),
  mean_af              = round(mean(af, na.rm = TRUE), 4),
  mean_ns              = round(mean(NS, na.rm = TRUE), 1),
  snp_ts               = sum(tstv == "Ts", na.rm = TRUE),
  snp_tv               = sum(tstv == "Tv", na.rm = TRUE)
), by = short_name]

seg_stats[, snp_tstv := fifelse(snp_tv > 0, round(snp_ts / snp_tv, 2), NA_real_)]

# ---------------------------------------------------------------------------
# 5. Join with contig lengths, compute density
# ---------------------------------------------------------------------------
cat("Joining with contig lengths...\n")
master <- merge(contigs, seg_stats, by = "short_name", all.x = TRUE)

# Zero-fill segments with no variants
count_cols <- c("total_variants", "polymorphic_variants", "snp_count",
                "indel_count", "sv_count", "snp_ts", "snp_tv")
for (col in count_cols) master[is.na(get(col)), (col) := 0]

master[, variant_density_per_kb := round(total_variants / (segment_length / 1000), 2)]
master[, polymorphic_density_per_kb := round(polymorphic_variants / (segment_length / 1000), 2)]

# ---------------------------------------------------------------------------
# 6. Join annotation overlaps (optional)
# ---------------------------------------------------------------------------
if (!is.null(annot_file)) {
  cat("Joining annotation overlaps from:", annot_file, "\n")
  annot <- fread(annot_file)

  # For grouped annotations where annotation_class differs from annotation
  # (e.g., pclai with AFR/EUR/EAS/SAS/AMR), create per-class columns.
  # Other annotations (genes, segdups, censat) and repeats get a single column
  # with the max overlap fraction across sub-classes.
  annot[, pivot_col := fifelse(
    annotation == "pclai" & annotation_class != annotation,
    paste0(annotation, "_", annotation_class),
    annotation)]

  annot_wide <- dcast(annot, augref_path ~ pivot_col,
                      value.var = "source_overlap_frac",
                      fun.aggregate = max, fill = 0)
  master <- merge(master, annot_wide, by = "augref_path", all.x = TRUE)
}

# ---------------------------------------------------------------------------
# 7. Write TSV (sorted by polymorphic_density_per_kb descending)
# ---------------------------------------------------------------------------
setorder(master, -polymorphic_density_per_kb)
fwrite(master, output_file, sep = "\t")
cat(sprintf("Wrote %d segments to %s\n", nrow(master), output_file))
