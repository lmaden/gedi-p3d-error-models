# =====================================================================
# fig_common.R
# Shared helpers for the main-figure scripts.
#
# Sourced by every section_G_fig*.R and section_G_suppfig_*.R.
# Holds:
#   - Tracker -> manuscript site ID mapping (per story-lock v2.2 §0)
#   - Screened-site classification (CHM-screened sites 1, 2, 3;
#     DTM-screened site 10)
#   - Lab-aligned color palette and ggplot theme add-ons
#   - Output directories under plots/section_G/
#   - Cluster-safe graphics device setup
# =====================================================================

Sys.setenv(DISPLAY = "")
options(device = pdf, bitmapType = "cairo")

if (!exists("log_progress"))  source("analysis_utils.R")
if (!exists("PROJECT_ROOT"))  source("analysis_config.R")

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(data.table)
  library(patchwork); library(cowplot); library(scales)
})

# ---------------------------------------------------------------------
# Output directories
# ---------------------------------------------------------------------
SECTION_G_PLOTS <- file.path(OUT_ROOT, "plots", "section_G")
dir.create(SECTION_G_PLOTS, recursive = TRUE, showWarnings = FALSE)

manuscript_tables_dir <- file.path(PROJECT_ROOT, "manuscript_tables")

# ---------------------------------------------------------------------
# Site ID mapping (per story-lock v2.2 §0)
# Tracker 1..15 = Manuscript 1..15
# Tracker 17 = Manuscript 16 (neon_ster2022)
# Tracker 18 = Manuscript 17 (neon_cper2021)
# Tracker 19 = Manuscript 18 (ltbmu_20180710)
# Tracker 20 = Manuscript 19 (neon_sjer)
# Tracker 16 (neon_nogp) and 21 (nasa_sonoma_whole) are OMITTED.
# ---------------------------------------------------------------------
site_id_lookup_path <- file.path(manuscript_tables_dir, "site_id_lookup.csv")

build_site_lookup_fallback <- function() {
  data.table(
    tracker_id    = c(1:15, 17, 18, 19, 20),
    manuscript_id = 1:19
  )
}

load_site_lookup <- function() {
  if (file.exists(site_id_lookup_path)) {
    sl <- as.data.table(fread(site_id_lookup_path))
    nms <- tolower(names(sl))
    setnames(sl, names(sl), nms)
    # Normalize tracker column name -> tracker_id
    if ("tracker_site"    %in% nms && !"tracker_id" %in% nms) {
      setnames(sl, "tracker_site",    "tracker_id")
    } else if ("tracker_site_id" %in% nms && !"tracker_id" %in% nms) {
      setnames(sl, "tracker_site_id", "tracker_id")
    }
    # Normalize manuscript column name -> manuscript_id
    if ("manuscript_site" %in% nms && !"manuscript_id" %in% nms) {
      setnames(sl, "manuscript_site", "manuscript_id")
    }
    # Optional pretty columns; rename for downstream convenience
    if ("site_short_name" %in% names(sl)) setnames(sl, "site_short_name", "short_name")
    if ("flag_status"     %in% names(sl)) sl[, flag_status     := as.character(flag_status)]
    if ("dominant_lc"     %in% names(sl)) sl[, dominant_lc     := as.character(dominant_lc)]
    return(sl)
  }
  log_progress("  site_id_lookup.csv not found; using built-in fallback mapping")
  build_site_lookup_fallback()
}

SITE_LOOKUP <- load_site_lookup()

tracker_to_manuscript <- function(x) {
  # Accept character or numeric; return integer manuscript IDs.
  xi <- suppressWarnings(as.integer(as.character(x)))
  m  <- match(xi, SITE_LOOKUP$tracker_id)
  SITE_LOOKUP$manuscript_id[m]
}

manuscript_to_tracker <- function(x) {
  xi <- suppressWarnings(as.integer(as.character(x)))
  m  <- match(xi, SITE_LOOKUP$manuscript_id)
  SITE_LOOKUP$tracker_id[m]
}

# ---------------------------------------------------------------------
# Site frame membership flags (canonical per story-lock §0)
# Both spaces provided:
#   TRACKER_* - for brms fit data (factor labels are tracker IDs)
#   MANUSCRIPT_* - for CSV-keyed data (post-renumbering)
# Site 16-19 manuscript == tracker 17-20; sites 1-15 are identical.
# ---------------------------------------------------------------------
TRACKER_CHM_SCREENED <- c(1L, 2L, 3L)
TRACKER_DTM_SCREENED <- c(10L)
TRACKER_CHM_16       <- c(4:9, 10:15, 17, 18, 19, 20)
TRACKER_DTM_18       <- c(1:9, 11:15, 17, 18, 19, 20)
TRACKER_PHASE25_15   <- intersect(TRACKER_CHM_16, TRACKER_DTM_18)

MANUSCRIPT_CHM_SCREENED <- c(1L, 2L, 3L)
MANUSCRIPT_DTM_SCREENED <- c(10L)
MANUSCRIPT_CHM_16       <- setdiff(1:19, MANUSCRIPT_CHM_SCREENED)
MANUSCRIPT_DTM_18       <- setdiff(1:19, MANUSCRIPT_DTM_SCREENED)
MANUSCRIPT_PHASE25_15   <- intersect(MANUSCRIPT_CHM_16, MANUSCRIPT_DTM_18)

site_frame_class <- function(tracker_id) {
  ti <- suppressWarnings(as.integer(as.character(tracker_id)))
  out <- character(length(ti))
  out[ti %in% TRACKER_PHASE25_15]                              <- "both"
  out[ti %in% TRACKER_DTM_18 & !(ti %in% TRACKER_PHASE25_15)]  <- "dtm_only"
  out[ti %in% TRACKER_CHM_16 & !(ti %in% TRACKER_PHASE25_15)]  <- "chm_only"
  out[ti %in% TRACKER_CHM_SCREENED]                             <- "chm_screened"
  out[ti %in% TRACKER_DTM_SCREENED & !(ti %in% TRACKER_CHM_16)] <- "dtm_screened"
  out[out == ""] <- "unknown"
  out
}

# Manuscript-space frame classifier (preferred when working from CSVs).
manuscript_frame_class <- function(manuscript_id) {
  mi <- suppressWarnings(as.integer(as.character(manuscript_id)))
  out <- character(length(mi))
  out[mi %in% MANUSCRIPT_PHASE25_15]                              <- "both"
  out[mi %in% MANUSCRIPT_DTM_18 & !(mi %in% MANUSCRIPT_PHASE25_15)] <- "dtm_only"
  out[mi %in% MANUSCRIPT_CHM_16 & !(mi %in% MANUSCRIPT_PHASE25_15)] <- "chm_only"
  out[mi %in% MANUSCRIPT_CHM_SCREENED]                              <- "chm_screened"
  out[mi %in% MANUSCRIPT_DTM_SCREENED & !(mi %in% MANUSCRIPT_CHM_16)] <- "dtm_screened"
  out[out == ""] <- "unknown"
  out
}

# ---------------------------------------------------------------------
# Lab-aligned color palette
# Restrained: categorical for product (CHM blue, DTM orange);
# screening status uses gray-scale + open symbols.
# ---------------------------------------------------------------------
COLOR_CHM             <- "#1f78b4"    # CHM product (blue)
COLOR_DTM             <- "#ff7f00"    # DTM product (orange)
COLOR_19SITE          <- "#9e9e9e"    # 19-site sensitivity context (gray)
COLOR_CHM_SCREENED    <- "#7f7f7f"    # CHM-screened site marker (medium gray)
COLOR_DTM_SCREENED    <- "#c0c0c0"    # DTM-screened site marker (light gray)
COLOR_REF_LINE        <- "#d62728"    # red, for 1:1 and regression overlays
COLOR_ZERO_LINE       <- "#7f7f7f"    # gray for zero / median guidelines

# ---------------------------------------------------------------------
# Minus-sign helper (U+2212): ensures captions/labels use proper minus
# rather than ASCII hyphen in numeric strings.
# ---------------------------------------------------------------------
MINUS <- "\u2212"

fmt_signed <- function(x, digits = 2) {
  # Convert numeric to character with U+2212 minus and digits decimals.
  s <- formatC(abs(x), format = "f", digits = digits)
  ifelse(is.na(x), NA_character_,
         ifelse(x < 0, paste0(MINUS, s), s))
}

# ---------------------------------------------------------------------
# Lab-aligned ggplot theme (cowplot half_open base, tighter strips)
# ---------------------------------------------------------------------
theme_section_G <- function(base_size = 11) {
  cowplot::theme_cowplot(font_size = base_size) +
    theme(
      panel.grid.major = element_line(color = "gray92", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray97", color = NA),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      plot.title       = element_text(size = base_size, face = "bold"),
      plot.subtitle    = element_text(size = base_size - 1, color = "gray30"),
      legend.position  = "right",
      legend.title     = element_text(size = base_size - 1),
      legend.text      = element_text(size = base_size - 2),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1)
    )
}

# ---------------------------------------------------------------------
# Save figure as both PDF (publication) and PNG (preview).
# Args:
#   plot:      ggplot or patchwork object
#   stem:      filename stem (no extension)
#   width_in:  width in inches
#   height_in: height in inches
#   dpi:       PNG resolution
# Returns:    invisible(c(pdf_path, png_path))
#
# Device selection:
#   PDF: prefers cairo_pdf (better font handling) when capabilities()
#        report cairo support; falls back to base pdf() otherwise.
#        Set the env var SECTION_G_FORCE_BASE_PDF=1 to skip cairo.
#   PNG: prefers ragg::agg_png, then cairo PNG, then base png().
# ---------------------------------------------------------------------
.have_cairo <- function() {
  if (identical(Sys.getenv("SECTION_G_FORCE_BASE_PDF"), "1")) return(FALSE)
  ok <- tryCatch(isTRUE(capabilities("cairo")), error = function(e) FALSE)
  if (!ok) return(FALSE)
  # capabilities() can report TRUE while the actual device errors;
  # do a tiny smoke test to confirm.
  tf <- tempfile(fileext = ".pdf")
  res <- tryCatch({
    grDevices::cairo_pdf(tf, width = 1, height = 1)
    grDevices::dev.off()
    TRUE
  }, error = function(e) FALSE)
  if (file.exists(tf)) try(file.remove(tf), silent = TRUE)
  res
}

.pdf_device <- function() {
  if (.have_cairo()) {
    log_progress("  PDF device: cairo_pdf")
    return(grDevices::cairo_pdf)
  }
  log_progress("  PDF device: base pdf() (cairo unavailable on this R build)")
  function(filename, width, height, ...) {
    grDevices::pdf(file = filename, width = width, height = height, ...)
  }
}

.png_device <- function() {
  if (requireNamespace("ragg", quietly = TRUE)) {
    log_progress("  PNG device: ragg::agg_png")
    return(ragg::agg_png)
  }
  if (isTRUE(capabilities("cairo"))) {
    log_progress("  PNG device: cairo png()")
    return(function(filename, width, height, units, res, ...) {
      grDevices::png(filename, width = width, height = height,
                     units = units, res = res, type = "cairo", ...)
    })
  }
  log_progress("  PNG device: base png()")
  function(filename, width, height, units, res, ...) {
    grDevices::png(filename, width = width, height = height,
                   units = units, res = res, ...)
  }
}

save_figure <- function(plot, stem, width_in = 6.5, height_in = 6.0,
                        dpi = 300, family = NULL) {
  pdf_path <- file.path(SECTION_G_PLOTS, paste0(stem, ".pdf"))
  png_path <- file.path(SECTION_G_PLOTS, paste0(stem, ".png"))
  pdev <- .pdf_device()
  qdev <- .png_device()
  ggsave(pdf_path, plot = plot, device = pdev,
         width = width_in, height = height_in, units = "in")
  ggsave(png_path, plot = plot, device = qdev,
         width = width_in, height = height_in, units = "in", dpi = dpi)
  log_progress(sprintf("  Wrote: %s", pdf_path))
  log_progress(sprintf("  Wrote: %s", png_path))
  invisible(c(pdf_path, png_path))
}

# ---------------------------------------------------------------------
# Banner for each section_G script.
# ---------------------------------------------------------------------
fig_banner <- function(figure_label, purpose) {
  log_progress("================================================================")
  log_progress(sprintf("%s: %s", figure_label, purpose))
  log_progress(sprintf("  Outputs: %s", SECTION_G_PLOTS))
  log_progress("================================================================")
}

# Quick sanity print when sourced directly.
if (!interactive() && sys.nframe() == 0) {
  log_progress("fig_common.R loaded.")
  log_progress(sprintf("  SITE_LOOKUP rows: %d", nrow(SITE_LOOKUP)))
  log_progress(sprintf("  Output dir: %s", SECTION_G_PLOTS))
}
