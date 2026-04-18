#!/usr/bin/env Rscript

# augref-frequency.R
#
# Reads per-segment depth from `vg depth -P augref_<REF>` and produces two
# summary plots on the off-reference (_alt) contigs:
#   - <prefix>.augref-frequency.png         histogram of haplotype counts
#   - <prefix>.augref-frequency-by-size.png length × depth scatter
#
# Input columns (vg depth output): path, bin_start, bin_end, mean_depth, stddev
# Depth = number of graph paths traversing each node in the segment
#         (= augref path itself + any sample haplotypes that share the content)

library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript augref-frequency.R <depth.tsv> <out-prefix>\n")
  quit(status = 1)
}
input_tsv <- args[1]
prefix    <- args[2]

cat("Reading:", input_tsv, "\n")
dt <- fread(input_tsv, header = FALSE,
            col.names = c("path", "bin_start", "bin_end", "mean_depth", "stddev"))
cat("Rows read:", nrow(dt), "\n")

# Keep only off-reference (_alt) contigs; main chromosomes always saturate
# depth and aren't informative for population frequency.
dt_alt <- dt[grepl("_alt$", path)]
dt_alt[, length := bin_end - bin_start + 1]
# Round depth to nearest integer; the mean across bases can have tiny
# fractional jitter where a few nodes lie off the shared subpath.
dt_alt[, depth_int := round(mean_depth)]
cat("Off-reference (_alt) segments:", nrow(dt_alt), "\n")

if (nrow(dt_alt) == 0) {
  cat("No _alt segments found; creating empty plots.\n")
  file.create(paste0(prefix, ".augref-frequency.png"))
  file.create(paste0(prefix, ".augref-frequency-by-size.png"))
  quit(status = 0)
}

save_png <- function(plot, path, w = 8, h = 6) {
  tryCatch({
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(path, width = w, height = h, units = "in", res = 300)
      print(plot)
      dev.off()
    } else {
      ggsave(path, plot = plot, width = w, height = h, dpi = 300,
             device = grDevices::png, type = "cairo")
    }
  }, error = function(e) {
    grDevices::png(path, width = w * 300, height = h * 300, res = 300, type = "cairo")
    print(plot)
    dev.off()
  })
}

# ---- Histogram: number of haplotypes covering each segment ----
max_depth <- max(dt_alt$depth_int, na.rm = TRUE)
p_hist <- ggplot(dt_alt, aes(x = depth_int)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.8) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = min(10, max_depth + 1))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Off-Reference Segment Frequency",
       subtitle = paste0("Number of graph paths traversing each augref segment (N=",
                         format(nrow(dt_alt), big.mark = ","), " segments)"),
       x = "# paths covering segment (including augref)",
       y = "# segments") +
  theme_minimal() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )
save_png(p_hist, paste0(prefix, ".augref-frequency.png"))
cat("Saved:", paste0(prefix, ".augref-frequency.png"), "\n")

# ---- Scatter: segment length vs frequency ----
p_scatter <- ggplot(dt_alt, aes(x = length, y = depth_int)) +
  geom_point(alpha = 0.3, size = 0.8, color = "steelblue") +
  scale_x_log10(labels = scales::comma, breaks = scales::breaks_log(n = 8)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = min(10, max_depth + 1))) +
  labs(title = "Off-Reference Segment: Length vs Population Frequency",
       subtitle = paste0("N=", format(nrow(dt_alt), big.mark = ","), " off-ref segments"),
       x = "Segment length (bp, log scale)",
       y = "# paths covering segment") +
  theme_minimal() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )
save_png(p_scatter, paste0(prefix, ".augref-frequency-by-size.png"))
cat("Saved:", paste0(prefix, ".augref-frequency-by-size.png"), "\n")

# Text summary
cat("\nPopulation-frequency summary:\n")
by_depth <- dt_alt[, .(n_segments = .N, total_bp = sum(length)), by = depth_int][order(depth_int)]
print(by_depth)
