# =============================================================================
# theme_p3d.R  —  House ggplot theme + type system for the P3D / GEDI figure set
# -----------------------------------------------------------------------------
# Source AFTER palette.R:   source("palette.R"); source("theme_p3d.R")
#
# Environment notes (this cluster):
#   * No Cairo graphics device. Everything here is device-agnostic; the saver
#     uses plain ggsave() (no cairo_pdf / cairo_png).
#   * Custom font is OPTIONAL via showtext. If showtext/sysfonts are missing OR
#     Google Fonts can't be reached, the theme silently falls back to the
#     device default sans family, so scripts stay portable and never error.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

# --- One type scale (points) -------------------------------------------------
p3d_type <- list(
  axis_title  = 11,
  tick        = 9,
  annotation  = 8.5,
  panel_label = 12,   # a / b / c ... tags on multi-panel figures
  legend      = 10,
  strip       = 10    # family / facet labels
)

# --- Optional custom font (Source Sans 3). Safe no-op if unavailable. ---------
p3d_font         <- "sans"     # fallback family
p3d_use_showtext <- FALSE
if (requireNamespace("showtext", quietly = TRUE) &&
    requireNamespace("sysfonts",  quietly = TRUE)) {
  .ok <- tryCatch({
    sysfonts::font_add_google("Source Sans 3", "p3d_sans"); TRUE
  }, error = function(e) FALSE)
  if (isTRUE(.ok)) {
    showtext::showtext_auto()
    p3d_font         <- "p3d_sans"
    p3d_use_showtext <- TRUE
  }
}

# --- The theme ---------------------------------------------------------------
# grid = "x" (default) draws only vertical reference lines — right for a
# horizontal coefficient plot. Use "y", "xy", or "none" elsewhere.
theme_p3d <- function(base_size = 11, base_family = p3d_font, grid = "x") {
  th <- theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text                 = element_text(family = base_family, color = p3d_gray$text_hi),
      axis.title           = element_text(size = p3d_type$axis_title, color = p3d_gray$text_hi),
      axis.title.y         = element_blank(),
      axis.text            = element_text(size = p3d_type$tick, color = p3d_gray$text_lo),
      axis.ticks.x         = element_line(color = p3d_gray$grid, linewidth = 0.3),
      axis.ticks.y         = element_blank(),
      legend.position      = "top",
      legend.direction     = "horizontal",
      legend.title         = element_blank(),
      # margin on the text gives inter-key spacing portably (avoids
      # legend.key.spacing.x, which only exists in ggplot2 >= 3.5).
      legend.text          = element_text(size = p3d_type$legend, color = p3d_gray$text_lo,
                                          margin = margin(r = 14)),
      legend.margin        = margin(b = 2),
      plot.title           = element_blank(),   # caption carries the title
      plot.subtitle        = element_blank(),
      plot.caption         = element_text(size = p3d_type$annotation, hjust = 0,
                                          color = p3d_gray$text_lo, lineheight = 1.15,
                                          margin = margin(t = 10)),
      plot.caption.position = "plot",
      plot.tag             = element_text(size = p3d_type$panel_label, face = "bold",
                                          color = p3d_gray$text_hi),
      strip.text           = element_text(size = p3d_type$strip, color = p3d_gray$strip,
                                          hjust = 0),
      panel.grid.minor     = element_blank(),
      panel.grid.major     = element_blank(),
      plot.margin          = margin(6, 12, 6, 6)
    )
  if (grepl("x", grid)) th <- th + theme(
    panel.grid.major.x = element_line(color = p3d_gray$grid, linewidth = 0.3))
  if (grepl("y", grid)) th <- th + theme(
    panel.grid.major.y = element_line(color = p3d_gray$band_b, linewidth = 0.25))
  th
}

# --- Device-agnostic saver ---------------------------------------------------
# Writes PNG (and optionally PDF) with NO cairo dependency. If showtext is
# active, it sets the showtext DPI to match the export so text scales correctly.
# House rule: >= 300 dpi (600 for line-heavy figures) at a fixed target width.
save_p3d <- function(plot, basename, width_in, height_in,
                     dpi = 600, dir = ".", formats = c("png", "pdf")) {
  if (isTRUE(p3d_use_showtext)) showtext::showtext_opts(dpi = dpi)
  # Prefer ragg for PNG (cairo-free, crisp text, no X server needed); fall back
  # to the default device if ragg isn't installed.
  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else NULL
  paths <- character(0)
  for (fmt in formats) {
    f <- file.path(dir, paste0(basename, ".", fmt))
    if (fmt == "png" && !is.null(png_device)) {
      ggplot2::ggsave(filename = f, plot = plot, device = png_device,
                      width = width_in, height = height_in, units = "in", dpi = dpi)
    } else {
      ggplot2::ggsave(filename = f, plot = plot,
                      width = width_in, height = height_in, units = "in", dpi = dpi)
    }
    paths <- c(paths, f)
    message(sprintf("  wrote %s  (%.2f x %.2f in @ %d dpi)", f, width_in, height_in, dpi))
  }
  invisible(paths)
}

invisible(TRUE)
