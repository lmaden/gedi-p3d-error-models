# =============================================================================
# palette.R  —  House color vocabulary for the P3D / GEDI figure set (Chapter 1)
# -----------------------------------------------------------------------------
# Source at the top of every figure script:   source("palette.R")
# Colorblind-safe throughout (Okabe-Ito base; Paul Tol 'muted' for high-N cats).
# =============================================================================

# --- LOCKED model + reference colors (do NOT change: these define identity) ---
p3d_col <- list(
  chm  = "#0072B2",   # CHM  — blue        (Okabe-Ito)
  dtm  = "#D55E00",   # DTM  — vermillion  (Okabe-Ito)
  zero = "#4D4D4D"    # zero / reference lines — neutral gray
)

# Named vector for scale_color_manual() / scale_shape_manual() on the two models.
# Keys match the series labels used in the figures.
p3d_series <- c("CHM 16-site (primary)" = p3d_col$chm,
                "DTM 18-site"           = p3d_col$dtm)

# --- Neutral grays used for chrome (bands, gridlines, text) -------------------
p3d_gray <- list(
  band_a  = "#FFFFFF",  # unshaded family band
  band_b  = "#F4F5F6",  # subtle alternating shaded family band
  grid    = "#E9EBED",  # vertical reference gridlines
  strip   = "#2E2E2E",  # family label text
  text_hi = "#262626",  # primary text
  text_lo = "#595959"   # secondary text / caption
)

# --- Curated CATEGORICAL palette (for later figures) --------------------------
# Okabe-Ito base. Note the CHM blue and DTM vermillion are RESERVED for model
# identity, so the general-purpose categorical helper skips them by default to
# avoid a category accidentally reading as "CHM" or "DTM".
okabe_ito <- c(black      = "#000000", orange     = "#E69F00",
               skyblue    = "#56B4E9", green      = "#009E73",
               yellow     = "#F0E442", blue       = "#0072B2",
               vermillion = "#D55E00", purple     = "#CC79A7")

# Forest functional types (Table 2: BDF, ENF, EBF, DNF) — 4 fixed hues, none of
# which collide with the reserved CHM/DTM colors.
p3d_fft <- c("BDF" = "#009E73",   # broadleaf deciduous   — green
             "ENF" = "#E69F00",   # evergreen needleleaf  — orange
             "EBF" = "#CC79A7",   # evergreen broadleaf   — reddish purple
             "DNF" = "#56B4E9")   # deciduous needleleaf  — sky blue

# Ecoregions (~9 EPA Level-II): Okabe-Ito only yields 8 distinct usable hues, so
# use Paul Tol 'muted' (9 colors, colorblind-safe) for this higher-cardinality
# case. Assign by name in your script for a stable mapping across figures.
p3d_ecoregion <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
                   "#DDCC77", "#CC6677", "#882255", "#AA4499")

# --- Accessor -----------------------------------------------------------------
# palette_p3d("series") / "fft" / "ecoregion" / "cat"
#   "cat" = general qualitative ramp that SKIPS the reserved CHM/DTM hues.
palette_p3d <- function(kind = c("series", "fft", "ecoregion", "cat"), n = NULL) {
  kind <- match.arg(kind)
  out <- switch(
    kind,
    series    = p3d_series,
    fft       = p3d_fft,
    ecoregion = p3d_ecoregion,
    cat       = unname(okabe_ito[!names(okabe_ito) %in% c("blue", "vermillion")])
  )
  if (!is.null(n)) out <- out[seq_len(n)]
  out
}

# Convenience ggplot scales (optional).
scale_color_p3d_series <- function(...) ggplot2::scale_color_manual(values = p3d_series, ...)
scale_fill_p3d_series  <- function(...) ggplot2::scale_fill_manual(values  = p3d_series, ...)

invisible(TRUE)
