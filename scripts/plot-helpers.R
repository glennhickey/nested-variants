# plot-helpers.R — shared ggplot save helper used by every plotting script in
# this repo.  Callers `source()` this file and then use save_plot() instead of
# rolling their own ggsave/agg_png wrapper.
#
# save_plot(p, file, width, height, pdf = FALSE)
#   - Writes <file> as PNG (preferring ragg, falling back to cairo).
#   - When pdf = TRUE, also writes <file with .png replaced by .pdf>
#     using grDevices::cairo_pdf (vector output).
#
# The calling script decides PDF policy. Conventions:
#   1. CLI flag `--pdf` enables PDF emission for this run.
#   2. Env var EMIT_PDF=1 enables it globally (handy for ad-hoc reruns).
# Use `pdf_enabled()` to resolve both into a single boolean.

pdf_enabled <- function(cli_flag = FALSE) {
  env <- Sys.getenv("EMIT_PDF", unset = "")
  cli_flag || env %in% c("1", "true", "TRUE", "yes")
}

save_plot <- function(plot, file, width = 8, height = 6, pdf = FALSE) {
  # PNG (ragg-preferred → cairo fallback)
  tryCatch({
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(file, width = width, height = height, units = "in", res = 300)
      print(plot); dev.off()
    } else {
      ggplot2::ggsave(file, plot = plot, width = width, height = height, dpi = 300,
                      device = grDevices::png, type = "cairo")
    }
  }, error = function(e) {
    grDevices::png(file, width = width * 300, height = height * 300, res = 300, type = "cairo")
    print(plot); dev.off()
  })
  cat("Saved:", file, "\n")

  if (isTRUE(pdf)) {
    pdf_file <- sub("\\.png$", ".pdf", file)
    tryCatch({
      grDevices::cairo_pdf(pdf_file, width = width, height = height)
      print(plot); dev.off()
      cat("Saved:", pdf_file, "\n")
    }, error = function(e) {
      cat("PDF emit failed for ", pdf_file, ": ", conditionMessage(e), "\n", sep = "")
      if (file.exists(pdf_file) && file.info(pdf_file)$size == 0) {
        file.remove(pdf_file)
      }
    })
  }
}
