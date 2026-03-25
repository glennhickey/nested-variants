#!/usr/bin/env Rscript

# contig-depth-summary.R — Averaged contig depth: pack vs BAM
#
# Reads multiple per-sample depth files (from vg depth and/or samtools coverage),
# averages across samples, and produces a two-panel scatter of length vs mean depth.
#
# Usage: Rscript scripts/contig-depth-summary.R
#          --pack-depths <file1,file2,...>    (from vg depth -b 1000000000)
#          --bam-depths <file1,file2,...>     (from samtools coverage)
#          --segs <augref-segs.tsv>
#          --output <output.png>
#          [--title TITLE]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
pack_paths <- NULL
bam_paths  <- NULL
segs_path  <- NULL
output     <- NULL
title      <- "Augref Contig Read Depth"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--pack-depths" && i + 1 <= length(args)) {
    pack_paths <- strsplit(args[i + 1], ",")[[1]]; i <- i + 2
  } else if (args[i] == "--bam-depths" && i + 1 <= length(args)) {
    bam_paths <- strsplit(args[i + 1], ",")[[1]]; i <- i + 2
  } else if (args[i] == "--segs" && i + 1 <= length(args)) {
    segs_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--output" && i + 1 <= length(args)) {
    output <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--title" && i + 1 <= length(args)) {
    title <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(segs_path) || is.null(output)) {
  cat("Usage: Rscript contig-depth-summary.R --pack-depths <f1,f2> --bam-depths <f1,f2> --segs <tsv> --output <png>\n")
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
# Read contig metadata
# ---------------------------------------------------------------------------
segs <- fread(segs_path, select = c(1, 2, 3, 4),
              col.names = c("source", "src_start", "src_end", "augref_path"))
segs[, length := src_end - src_start]
segs <- unique(segs, by = "augref_path")
cat("Augref contigs:", nrow(segs), "\n")

# ---------------------------------------------------------------------------
# Helper: read and average depth across samples
# ---------------------------------------------------------------------------
read_pack_depths <- function(paths) {
  all_dt <- rbindlist(lapply(paths, function(f) {
    dt <- fread(f, header = FALSE,
                col.names = c("path", "bin_start", "bin_end", "mean_depth", "stddev"))
    # Keep only augref off-ref contigs
    dt[grepl("_[0-9]+_alt$", path), .(path, mean_depth)]
  }))
  if (nrow(all_dt) == 0) return(data.table(augref_path = character(), mean_depth = numeric()))
  # Average across samples per contig
  all_dt[, .(mean_depth = mean(mean_depth)), by = .(augref_path = path)]
}

read_bam_depths <- function(paths) {
  all_dt <- rbindlist(lapply(paths, function(f) {
    dt <- fread(f, header = TRUE)
    setnames(dt, 1, "rname")
    # Keep only augref off-ref contigs
    dt[grepl("_[0-9]+_alt$", rname), .(rname, meandepth)]
  }))
  if (nrow(all_dt) == 0) return(data.table(augref_path = character(), mean_depth = numeric()))
  all_dt[, .(mean_depth = mean(meandepth)), by = .(augref_path = rname)]
}

# ---------------------------------------------------------------------------
# Read pack depths
# ---------------------------------------------------------------------------
panels <- list()

if (!is.null(pack_paths) && length(pack_paths) > 0) {
  pack_avg <- read_pack_depths(pack_paths)
  pack_merged <- merge(segs, pack_avg, by = "augref_path", all.x = TRUE)
  pack_merged[is.na(mean_depth), mean_depth := 0]

  # Get backbone depth for reference line
  backbone_depths <- rbindlist(lapply(pack_paths, function(f) {
    dt <- fread(f, header = FALSE,
                col.names = c("path", "bin_start", "bin_end", "mean_depth", "stddev"))
    dt[grepl("^augref_", path) & !grepl("_[0-9]+_alt$", path)]
  }))
  ref_depth <- if (nrow(backbone_depths) > 0) {
    backbone_depths[, weighted.mean(mean_depth, bin_end - bin_start)]
  } else NA_real_

  n_covered <- sum(pack_merged$mean_depth > 0)
  cat("Pack: ", n_covered, "of", nrow(pack_merged), "contigs covered, backbone:", round(ref_depth, 1), "x\n")

  pack_merged[, source := "Graph (vg pack)"]
  panels[["pack"]] <- list(dt = pack_merged, ref_depth = ref_depth, n_covered = n_covered)
}

# ---------------------------------------------------------------------------
# Read BAM depths
# ---------------------------------------------------------------------------
if (!is.null(bam_paths) && length(bam_paths) > 0) {
  bam_avg <- read_bam_depths(bam_paths)
  bam_merged <- merge(segs, bam_avg, by = "augref_path", all.x = TRUE)
  bam_merged[is.na(mean_depth), mean_depth := 0]

  # Backbone from BAM
  backbone_bam <- rbindlist(lapply(bam_paths, function(f) {
    dt <- fread(f, header = TRUE)
    setnames(dt, 1, "rname")
    dt[grepl("^augref_", rname) & !grepl("_[0-9]+_alt$", rname)]
  }))
  ref_depth_bam <- if (nrow(backbone_bam) > 0) {
    backbone_bam[, weighted.mean(meandepth, endpos - startpos + 1)]
  } else NA_real_

  n_covered_bam <- sum(bam_merged$mean_depth > 0)
  cat("BAM: ", n_covered_bam, "of", nrow(bam_merged), "contigs covered, backbone:", round(ref_depth_bam, 1), "x\n")

  bam_merged[, source := "Linear (surject)"]
  panels[["bam"]] <- list(dt = bam_merged, ref_depth = ref_depth_bam, n_covered = n_covered_bam)
}

if (length(panels) == 0) {
  cat("No depth data provided.\n")
  file.create(output)
  quit(status = 0)
}

# ---------------------------------------------------------------------------
# Build combined plot
# ---------------------------------------------------------------------------
n_samples_pack <- if (!is.null(pack_paths)) length(pack_paths) else 0
n_samples_bam  <- if (!is.null(bam_paths))  length(bam_paths)  else 0
n_samples <- max(n_samples_pack, n_samples_bam)

plot_dt <- rbindlist(lapply(panels, function(p) p$dt), fill = TRUE)
plot_dt[, covered := mean_depth > 0]

# Build subtitle per facet
ref_lines <- rbindlist(lapply(names(panels), function(nm) {
  p <- panels[[nm]]
  data.table(source = p$dt$source[1], ref_depth = p$ref_depth)
}))

p <- ggplot(plot_dt, aes(x = length, y = mean_depth)) +
  geom_point(aes(shape = covered), alpha = 0.6, size = 1.5, color = "steelblue") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 4),
                     labels = c("TRUE" = "Covered", "FALSE" = "No reads"),
                     name = NULL) +
  geom_hline(data = ref_lines, aes(yintercept = ref_depth),
             linetype = "dashed", color = "grey40") +
  facet_wrap(~ source, scales = "free_y") +
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = title,
       subtitle = paste0("Mean across ", n_samples, " samples  |  dashed line = backbone depth"),
       x = "Contig Length (bp)",
       y = "Mean Read Depth") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 11)
  )

save_png(p, output, width = 12, height = 5)
cat("Done.\n")
