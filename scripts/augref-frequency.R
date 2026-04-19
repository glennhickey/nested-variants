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

# ---- CCDF: # segments covered by AT LEAST N paths ----
max_depth <- max(dt_alt$depth_int, na.rm = TRUE)
ccdf_dt <- data.table(
  depth = seq_len(max_depth),
  n_segments = sapply(seq_len(max_depth),
                      function(k) sum(dt_alt$depth_int >= k))
)
ccdf_dt[, total_bp := sapply(depth,
                              function(k) sum(dt_alt$length[dt_alt$depth_int >= k]))]
ccdf_dt[, frac_segments := n_segments / nrow(dt_alt)]
ccdf_dt[, frac_bp       := total_bp / sum(dt_alt$length)]

p_hist <- ggplot(ccdf_dt, aes(x = depth, y = n_segments)) +
  geom_step(direction = "hv", linewidth = 0.9, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = min(10, max_depth))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Off-Reference Segment Frequency (CCDF)",
       subtitle = paste0("# segments supported by at least N graph paths (N=",
                         format(nrow(dt_alt), big.mark = ","), " off-ref segments)"),
       x = "# paths covering segment ≥ N (including augref)",
       y = "# segments with depth ≥ N") +
  theme_minimal() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )
save_png(p_hist, paste0(prefix, ".augref-frequency.png"))
cat("Saved:", paste0(prefix, ".augref-frequency.png"), "\n")

# ---- CCDF of total bp: # bp in segments covered by AT LEAST N paths ----
p_scatter <- ggplot(ccdf_dt, aes(x = depth, y = total_bp)) +
  geom_step(direction = "hv", linewidth = 0.9, color = "firebrick") +
  geom_point(size = 2, color = "firebrick") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = min(10, max_depth))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Off-Reference Sequence Content (CCDF)",
       subtitle = paste0("Total bp in segments supported by at least N graph paths"),
       x = "# paths covering segment ≥ N (including augref)",
       y = "Total bp (segments with depth ≥ N)") +
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
