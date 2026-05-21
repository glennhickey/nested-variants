#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(RColorBrewer)

# Title suppression for paper figures (FIGURE_TITLES=0 → labs() titles → NULL).
.script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
if (is.na(.script_dir) || !nzchar(.script_dir)) .script_dir <- "scripts"
source(file.path(.script_dir, "plot-helpers.R"))

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check minimum number of arguments
if (length(args) < 2) {
  cat("Usage: Rscript offref-length-hist.R <output.png> <data1.tsv|bed> [data2.tsv|bed] [...] [min_threshold] [cumulative]\n")
  cat("Example: Rscript offref-length-hist.R result.png file1.tsv file2.bed file3.tsv 50 TRUE\n")
  cat("  output.png: Output filename for the plot\n")
  cat("  data files: One or more TSV/BED files where length is calculated as COLUMN3 - COLUMN2\n")
  cat("  min_threshold: Optional minimum value threshold (values below this will be excluded)\n")
  cat("  cumulative: Optional TRUE/FALSE for cumulative distribution (default: FALSE)\n")
  cat("  Labels are auto-generated from filenames (basename without extension)\n")
  quit(status = 1)
}

output <- args[1]

# Find where data files end and optional parameters begin
# Data files are those that exist and end with .tsv, .txt, or .bed
data_files <- c()
opt_start <- length(args) + 1

for (i in 2:length(args)) {
  arg <- args[i]
  # Check if this looks like a file (ends with .tsv, .txt, or .bed and exists)
  if ((grepl("\\.(tsv|txt|bed)$", arg, ignore.case = TRUE)) && file.exists(arg)) {
    data_files <- c(data_files, arg)
  } else {
    opt_start <- i
    break
  }
}

# Check we have at least one data file
if (length(data_files) == 0) {
  cat("Error: No valid data files found\n")
  quit(status = 1)
}

# Parse optional parameters
min_threshold <- NULL
cumulative <- FALSE

if (opt_start <= length(args)) {
  # First optional arg could be min_threshold
  first_opt <- args[opt_start]
  if (grepl("^[0-9.]+$", first_opt)) {
    min_threshold <- as.numeric(first_opt)
    opt_start <- opt_start + 1
  }

  # Second optional arg could be cumulative
  if (opt_start <= length(args)) {
    cumulative <- as.logical(toupper(args[opt_start]))
  }
}

# Read all data files and calculate lengths
all_lengths <- list()
all_labels <- c()

for (i in 1:length(data_files)) {
  file <- data_files[i]
  cat("Reading", file, "...\n")

  # Read TSV file
  data <- read.delim(file, header = FALSE)

  # Calculate lengths as COLUMN3 - COLUMN2
  lengths <- data$V3 - data$V2
  all_lengths[[i]] <- lengths

  # Generate label from filename (basename without extension)
  label <- sub("\\.[^.]*$", "", basename(file))
  all_labels <- c(all_labels, label)

  cat("  Found", length(lengths), "intervals\n")
}

# Combine all data into a single dataframe
df <- data.frame(
  value = unlist(all_lengths),
  dataset = rep(all_labels, times = sapply(all_lengths, length))
)

# Apply minimum threshold filter if specified
if (!is.null(min_threshold)) {
  original_count <- nrow(df)
  df <- df[df$value >= min_threshold, ]
  filtered_count <- original_count - nrow(df)
  cat("Filtered out", filtered_count, "values below threshold", min_threshold, "\n")
  cat("Remaining values:", nrow(df), "\n")
}

# For many chromosomes we drop per-dataset rainbow colors and prefer
# a muted-lines + summary scheme (cumulative) or a box-per-chrom plot
# (histogram).  Below 10 datasets we keep the original coloured style.
many <- length(all_labels) >= 10

# Derive clean chrom labels by stripping the ".augref-segs" suffix when
# present, then natural-sort (chr1, chr2, ..., chr22, chrX, chrY).
chrom_labels <- sub("\\.augref-segs$", "", all_labels)
nat_levels <- function(v) {
  u <- unique(v)
  suffix <- sub("^chr", "", u)
  key <- suppressWarnings(as.numeric(suffix))
  key[suffix == "X"] <- 100
  key[suffix == "Y"] <- 101
  key[suffix == "M"] <- 102
  key[is.na(key)] <- 1000
  u[order(key, u)]
}
df$chrom <- factor(sub("\\.augref-segs$", "", df$dataset),
                   levels = nat_levels(chrom_labels))

# Colours (only used in the <10-dataset path)
n_datasets <- length(all_labels)
if (n_datasets <= 2) {
  colors <- c("#0072B2", "#D55E00")[1:n_datasets]
} else if (n_datasets <= 8) {
  colors <- RColorBrewer::brewer.pal(max(3, n_datasets), "Set2")[1:n_datasets]
} else {
  colors <- rainbow(n_datasets)
}
color_map <- setNames(colors, all_labels)

# Create plot based on cumulative flag
if (cumulative) {
  df_sorted <- df[order(df$value, decreasing = TRUE), ]
  df_cumulative <- do.call(rbind, lapply(unique(df$dataset), function(ds) {
    subset_data <- df_sorted[df_sorted$dataset == ds, ]
    data.frame(
      value = subset_data$value,
      count = seq_len(nrow(subset_data)),
      dataset = ds
    )
  }))

  if (many) {
    # Combined cumulative across ALL chromosomes (one bold line)
    all_sorted <- sort(df$value, decreasing = TRUE)
    df_combined <- data.frame(value = all_sorted, count = seq_along(all_sorted))

    p <- ggplot() +
      geom_line(data = df_cumulative,
                aes(x = value, y = count, group = dataset),
                colour = "grey60", linewidth = 0.4, alpha = 0.5) +
      geom_line(data = df_combined,
                aes(x = value, y = count),
                colour = "black", linewidth = 1.2) +
      scale_x_log10(labels = scales::comma,
                    breaks = scales::breaks_log(n = 10)) +
      scale_y_log10(labels = scales::comma) +
      labs(title = if_titles("Off-Reference Segment Lengths (Cumulative Count)"),
           subtitle = if_titles(sprintf("Grey: each of %d chromosomes.  Black: all chromosomes combined.",
                              n_datasets)),
           x = "Length (log scale)",
           y = "Count >= Length (log scale)") +
      theme_minimal()
  } else {
    p <- ggplot(df_cumulative, aes(x = value, y = count, color = dataset)) +
      geom_line(linewidth = 0.8, alpha = 0.8) +
      geom_point(size = 1.5, alpha = 0.6) +
      scale_x_log10(labels = scales::comma,
                    breaks = scales::breaks_log(n = 10)) +
      scale_y_log10(labels = scales::comma) +
      scale_color_manual(values = color_map) +
      labs(title = if_titles("Off-Reference Segment Lengths (Cumulative Count)"),
           x = "Length (log scale)",
           y = "Count >= Length (log scale)",
           color = "Dataset") +
      theme_minimal()
  }
} else {
  if (many) {
    # One boxplot per chromosome on a log y-axis, vertical x-labels.
    p <- ggplot(df, aes(x = chrom, y = value)) +
      geom_boxplot(fill = "#6baed6", colour = "#08519c",
                   outlier.size = 0.4, outlier.alpha = 0.3,
                   linewidth = 0.4) +
      scale_y_log10(labels = scales::comma,
                    breaks = scales::breaks_log(n = 10)) +
      labs(title = if_titles("Off-Reference Interval Lengths by Chromosome"),
           x = NULL, y = "Length (log scale)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  } else {
    min_val <- min(df$value)
    max_val <- max(df$value)
    log_breaks <- 10^seq(log10(min_val), log10(max_val), length.out = 21)

    p <- ggplot(df, aes(x = value, fill = dataset, color = dataset)) +
      geom_histogram(alpha = 0.5, position = "identity",
                     breaks = log_breaks) +
      scale_x_log10(labels = scales::comma,
                    breaks = scales::breaks_log(n = 10)) +
      scale_y_continuous(trans = "log1p", labels = scales::comma) +
      scale_fill_manual(values = color_map) +
      scale_color_manual(values = color_map) +
      labs(title = if_titles("Off-Reference Interval Lengths (Log Scale)"),
           x = "Value (log scale)",
           y = "Frequency (log scale)",
           fill = "Dataset") +
      theme_minimal() +
      guides(color = "none")
  }
}

save_plot <- function(plot, path, w = 8, h = 6) {
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
    grDevices::png(path, width = w*300, height = h*300, res = 300, type = "cairo")
    print(plot)
    dev.off()
  })
}

save_plot(p, output)

# Log-log cumulative CCDF: for each dataset, plot the number of segments
# with length >= x on both log-scaled axes. Reveals the tail distribution.
loglog_output <- sub("\\.png$", "-loglog.png", output)
ccdf_list <- lapply(data_files, function(f) {
  d <- read.delim(f, header = FALSE)
  vals <- sort(d[[3]] - d[[2]], decreasing = TRUE)
  data.frame(length = vals, rank = seq_along(vals),
             dataset = tools::file_path_sans_ext(basename(f)))
})
ccdf_df <- do.call(rbind, ccdf_list)
p_loglog <- ggplot(ccdf_df, aes(x = length, y = rank, color = dataset)) +
  geom_step(linewidth = 0.6) +
  scale_x_log10(labels = scales::comma, breaks = scales::breaks_log(n = 10)) +
  scale_y_log10(labels = scales::comma) +
  scale_color_manual(values = color_map) +
  labs(x = "GRef Segment Length", y = "Count >= Length") +
  guides(color = "none") +
  theme_minimal()
save_plot(p_loglog, loglog_output)

# Backing TSV (figure 1B): (dataset, length, rank) — one row per segment,
# sorted by length descending; rank = count of segments with length >= row's length.
loglog_tsv <- sub("\\.png$", ".tsv", loglog_output)
tryCatch({
  utils::write.table(ccdf_df[, c("dataset", "length", "rank")],
                     file = loglog_tsv, sep = "\t",
                     quote = FALSE, row.names = FALSE)
  cat("TSV saved to", loglog_tsv, "\n")
}, error = function(e) {
  cat("TSV emit failed for ", loglog_tsv, ": ", conditionMessage(e), "\n", sep = "")
})

# Always emit a PDF companion of the log-log CCDF — used as figure 1B in
# the manuscript, so vector output matters.
loglog_pdf_output <- sub("\\.png$", ".pdf", loglog_output)
tryCatch({
  ggsave(loglog_pdf_output, plot = p_loglog, width = 8, height = 6,
         device = grDevices::cairo_pdf)
  cat("PDF saved to", loglog_pdf_output, "\n")
}, error = function(e) {
  cat("PDF emit failed for ", loglog_pdf_output, ": ", conditionMessage(e),
      "\n", sep = "")
  if (file.exists(loglog_pdf_output) && file.info(loglog_pdf_output)$size == 0) {
    file.remove(loglog_pdf_output)
  }
})

cat("Histogram saved to", output, "\n")

# Emit SVG alongside PNG when EMIT_SVG=1
if (nzchar(Sys.getenv("EMIT_SVG")) && grepl("\\.png$", output)) {
  pdf_output <- sub("\\.png$", ".pdf", output)
  tryCatch({
    ggsave(pdf_output, plot = p, width = 8, height = 6, device = grDevices::cairo_pdf)
    cat("PDF saved to", pdf_output, "\n")
  }, error = function(e) {
    cat("PDF emit failed for ", pdf_output, ": ", conditionMessage(e), "\n", sep = "")
    if (file.exists(pdf_output) && file.info(pdf_output)$size == 0) file.remove(pdf_output)
  })
}
