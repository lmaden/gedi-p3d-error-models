# =============================================================================
# section_H_fig3_coef_twoseries.R
#
# Manuscript Figure 3, TWO-series version: 16-site CHM and
# 18-site DTM only. The 19-site CHM sensitivity series is dropped from the
# MAIN figure and retained on the supplementary comparison figure (which
# reuses the existing word/media/image3.png -- no regeneration needed).
#
# This is a FRESH, NON-MUTATING generator. The DATA-PREP stage (Sections 1-4:
# CSV load, detect_col standardization, DTM fixef extraction, the 17 main-effect
# list, |est_16| ordering, and BOTH sign-reversal counts) and the VERIFICATION
# stage (Section 8) are reproduced VERBATIM from the committed script so the
# story-lock 2F fingerprints and the 14/17 & 13/17 counts are unchanged.
#
# ONLY the presentation diverges (Sections 5-7), per the figure-upgrade brief:
#   - source the reusable house style (palette.R + theme_p3d.R)
#   - group the 17 predictors into 4 labeled families with subtle alternating
#     background shading; order by |max effect| WITHIN each family
#   - dumbbell connectors between the CHM and DTM points so sign reversals
#     (opposite-sign pairs crossing the zero line) read at a glance
#   - locked colors, thin dashed zero line, consistent point/bar weight
#   - fixed-effect R^2 + reversal count carried in the caption (not the panel)
#   - width/height UNCHANGED (7.5 x 7.3) so the docx <wp:extent> aspect
#     (5943600 / 5781712 = 1.028) still swaps cleanly; now exported at 600 dpi
#     (dpi does not affect aspect)
#
# -----------------------------------------------------------------------------
# NOTE -- point estimate = posterior MEDIAN (matches manuscript convention)
#   The manuscript reports "posterior medians with 95% equal-tailed quantile
#   credible intervals throughout the paper and supplement," and Table 6 says
#   the same, so this figure now plots MEDIANS. The credible intervals were
#   already equal-tailed Q2.5/Q97.5 quantiles, so only the point moved.
#     - DTM: fixef(fit_dtm, robust = TRUE, ...)  -> median  (done in Section 2)
#     - CHM: est_16 read from chm_sensitivity_fixef_comparison.csv. CONFIRM that
#       section_10b_chm_sensitivity_refit.R wrote est_16 as a MEDIAN (Table 6's
#       published medians match it, so it almost certainly did). If it wrote a
#       mean, regenerate that column as a median.
#   RE-FINGERPRINT after this change (Section 8): the anchors move only at the
#   3rd decimal, but a coefficient sitting on zero (relative geoaccuracy, CHM
#   ~= -0.00) could in principle flip sign and change the 14/17 count. Confirm
#   the counts before any docx media swap; if a count moves, update story-lock
#   2F and the in-text "14 of 17" claim to match.
#   To revert to the committed posterior-mean behavior: robust = FALSE above.
#
# Sources (identical to the committed generator):
#   - manuscript_tables/chm_sensitivity_fixef_comparison.csv  (CHM 16/19)
#   - load_checkpoint("10_models_stage2")$fit_dtm_s2 -> fixef()  (DTM 18)
# =============================================================================

source("fig_common.R")     # load_checkpoint, save_figure, logging, dirs
source("palette.R")                  # p3d_col, p3d_series, p3d_gray   (place next
source("theme_p3d.R")                # theme_p3d(), save_p3d()          to common)

fig_banner("Figure 3 (two-series)",
                 "CHM 16-site + DTM 18-site coefficient forest plot (D4)")

suppressPackageStartupMessages({
  library(brms); library(data.table); library(ggplot2); library(dplyr)
})

# -----------------------------------------------------------------------------
# 1. CHM 16-site and 19-site coefficients from comparison CSV   (VERBATIM)
# -----------------------------------------------------------------------------
log_subsection("Loading CHM fixef comparison (16-site vs 19-site)")

chm_cmp_path <- file.path(manuscript_tables_dir,
                          "chm_sensitivity_fixef_comparison.csv")
if (!file.exists(chm_cmp_path)) {
  stop("Required CHM comparison CSV not found: ", chm_cmp_path,
       "\n  Run section_10b_chm_sensitivity_refit.R first.")
}
chm_cmp <- as.data.table(fread(chm_cmp_path))
log_progress(sprintf("  Loaded %d rows, %d cols",
                     nrow(chm_cmp), ncol(chm_cmp)))
cat("  Columns: ", paste(names(chm_cmp), collapse = ", "), "\n")

detect_col <- function(dt, patterns) {
  hits <- unlist(lapply(patterns, function(p) grep(p, names(dt),
                                                   value = TRUE, ignore.case = TRUE)))
  if (!length(hits)) return(NA_character_)
  hits[1]
}

col_pred <- detect_col(chm_cmp, c("^predictor$", "^term$", "^param$", "^variable$"))
col_e19  <- detect_col(chm_cmp, c("est_full$", "estimate_full", "estimate.full",
                                  "Estimate_full", "est_19", "Estimate.full"))
col_l19  <- detect_col(chm_cmp, c("q2\\.?5_full$", "lower_full", "lo_full",
                                  "Q2\\.5.full", "ci_lo_full", "ci\\.low\\.full"))
col_h19  <- detect_col(chm_cmp, c("q97\\.?5_full$", "upper_full", "hi_full",
                                  "Q97\\.5.full", "ci_hi_full", "ci\\.high\\.full"))
col_e16  <- detect_col(chm_cmp, c("est_16$", "estimate_16", "estimate.16",
                                  "Estimate_16", "Estimate.16"))
col_l16  <- detect_col(chm_cmp, c("q2\\.?5_16$", "lower_16", "lo_16",
                                  "Q2\\.5.16", "ci_lo_16", "ci\\.low\\.16"))
col_h16  <- detect_col(chm_cmp, c("q97\\.?5_16$", "upper_16", "hi_16",
                                  "Q97\\.5.16", "ci_hi_16", "ci\\.high\\.16"))
needed <- c(col_pred, col_e19, col_l19, col_h19, col_e16, col_l16, col_h16)
if (any(is.na(needed))) {
  stop("Could not locate all expected columns in chm_sensitivity_fixef_comparison.csv.\n",
       "Columns present: ", paste(names(chm_cmp), collapse = ", "),
       "\nDetected: pred=", col_pred, " e19=", col_e19, " l19=", col_l19,
       " h19=", col_h19, " e16=", col_e16, " l16=", col_l16, " h16=", col_h16)
}

chm <- data.table(
  predictor = as.character(chm_cmp[[col_pred]]),
  est_19    = as.numeric(chm_cmp[[col_e19]]),
  lo_19     = as.numeric(chm_cmp[[col_l19]]),
  hi_19     = as.numeric(chm_cmp[[col_h19]]),
  est_16    = as.numeric(chm_cmp[[col_e16]]),
  lo_16     = as.numeric(chm_cmp[[col_l16]]),
  hi_16     = as.numeric(chm_cmp[[col_h16]])
)
log_progress(sprintf("  CHM comparison standardized to %d predictor rows",
                     nrow(chm)))

# -----------------------------------------------------------------------------
# 2. DTM 18-site coefficients from fit   (VERBATIM)
# -----------------------------------------------------------------------------
log_subsection("Loading DTM 18-site fit and extracting fixef")

s2 <- load_checkpoint("10_models_stage2")
if (!"fit_dtm_s2" %in% names(s2)) {
  stop("fit_dtm_s2 not found in 10_models_stage2 checkpoint.")
}
fit_dtm <- s2$fit_dtm_s2
# robust = TRUE -> Estimate is the posterior MEDIAN (matches the manuscript
# convention; see header NOTE). Q2.5/Q97.5 are the same equal-tailed quantiles
# either way. Set robust = FALSE to revert to the posterior mean.
dtm_fx_raw <- as.data.table(fixef(fit_dtm, robust = TRUE, probs = c(0.025, 0.975)),
                            keep.rownames = "predictor")
setnames(dtm_fx_raw,
         c("Estimate", "Est.Error", "Q2.5", "Q97.5"),
         c("est_dtm", "se_dtm", "lo_dtm", "hi_dtm"),
         skip_absent = TRUE)
if (!"est_dtm" %in% names(dtm_fx_raw) && "Estimate" %in% names(dtm_fx_raw)) {
  setnames(dtm_fx_raw, "Estimate", "est_dtm")
}
dtm <- dtm_fx_raw[, .(predictor, est_dtm, lo_dtm, hi_dtm)]
log_progress(sprintf("  DTM fixef extracted: %d coefficient rows", nrow(dtm)))

# -----------------------------------------------------------------------------
# 3. The 17 main-effect continuous predictors   (VERBATIM, story-lock 2F)
#    PRETTY labels reconciled to manuscript Table 2/6 wording (display only;
#    does not touch any value or count). Changes vs committed labels:
#      Mean terrain slope           -> Terrain slope (mean)
#      Waveform structural complexity -> Waveform complexity (WSCI)
#      Leaf-on image ratio          -> Leaf-on fraction
#      Absolute geolocation         -> Absolute geoaccuracy
#      Relative geolocation         -> Relative geoaccuracy
#      Slope SD                     -> Slope variability (SD)
#      Aspect (cos) / Aspect (sin)  -> Aspect (N-S) / Aspect (E-W)
#      Stereo angle                 -> Stereo ratio
# -----------------------------------------------------------------------------
MAIN_EFFECTS_17 <- c(
  "slope_mean_z",  "rh_98_z",       "wsci_z",
  "cover_z",       "meta_leafon_z", "meta_sunel_z",
  "meta_offnad_z", "meta_relgeo_z", "view_az_cos_z",
  "meta_absgeo_z", "meta_az_conc_z","meta_fwdrev_z",
  "view_az_sin_z", "slope_sd_z",    "aspect_cos_z",
  "aspect_sin_z",  "meta_stereo_z"
)
PRETTY <- c(
  slope_mean_z  = "Terrain slope (mean)",
  rh_98_z       = "Canopy height (RH98)",
  wsci_z        = "Waveform complexity (WSCI)",
  cover_z       = "Canopy cover",
  meta_leafon_z = "Leaf-on fraction",
  meta_sunel_z  = "Sun elevation",
  meta_offnad_z = "Off-nadir angle",
  meta_relgeo_z = "Relative geoaccuracy",
  view_az_cos_z = "Viewing azimuth (cos)",
  meta_absgeo_z = "Absolute geoaccuracy",
  meta_az_conc_z= "Azimuth concentration",
  meta_fwdrev_z = "Forward/reverse ratio",
  view_az_sin_z = "Viewing azimuth (sin)",
  slope_sd_z    = "Slope variability (SD)",
  aspect_cos_z  = "Aspect (N-S)",
  aspect_sin_z  = "Aspect (E-W)",
  meta_stereo_z = "Stereo ratio"
)

chm <- chm[predictor %in% MAIN_EFFECTS_17]
dtm <- dtm[predictor %in% MAIN_EFFECTS_17]
if (nrow(chm) != 17L) {
  log_progress(sprintf("  WARNING: expected 17 CHM main-effect rows, got %d. Predictors found: %s",
                       nrow(chm), paste(chm$predictor, collapse = ",")))
}
if (nrow(dtm) != 17L) {
  log_progress(sprintf("  WARNING: expected 17 DTM main-effect rows, got %d. Predictors found: %s",
                       nrow(dtm), paste(dtm$predictor, collapse = ",")))
}

all_fx <- merge(chm, dtm, by = "predictor", all = FALSE)

# Sign-reversal counts (BOTH computed; only 16-vs-18 goes on the figure).
all_fx[, signrev_16_dtm := sign(est_16) != sign(est_dtm)]
all_fx[, signrev_19_dtm := sign(est_19) != sign(est_dtm)]
n_rev_16 <- sum(all_fx$signrev_16_dtm, na.rm = TRUE)
n_rev_19 <- sum(all_fx$signrev_19_dtm, na.rm = TRUE)
log_progress(sprintf("  Sign reversals CHM 16-site vs DTM 18-site: %d / %d",
                     n_rev_16, nrow(all_fx)))
log_progress(sprintf("  Sign reversals CHM 19-site vs DTM 18-site: %d / %d",
                     n_rev_19, nrow(all_fx)))

# -----------------------------------------------------------------------------
# 4. Order by |CHM 16-site estimate| descending   (VERBATIM)
#    Kept so all_fx stays in its canonical order for the Section 8 fingerprint.
#    The FIGURE uses a family-grouped order built on a copy in Section 5.
# -----------------------------------------------------------------------------
ord <- order(-abs(all_fx$est_16))
all_fx <- all_fx[ord]
all_fx[, pretty := factor(PRETTY[predictor], levels = rev(PRETTY[predictor]))]

# =============================================================================
# 5. FIGURE DATA -- family grouping + within-family magnitude order  (UPGRADE)
# =============================================================================
log_subsection("Composing two-series forest plot (grouped, house style)")

# Predictor -> family map (anchored to Table 2 categories; geolocation split out
# of the acquisition block to foreground the misregistration narrative).
FAMILY <- c(
  rh_98_z       = "Canopy structure", wsci_z = "Canopy structure",
  cover_z       = "Canopy structure",
  slope_mean_z  = "Terrain", slope_sd_z = "Terrain",
  aspect_cos_z  = "Terrain", aspect_sin_z = "Terrain",
  meta_absgeo_z = "Acquisition (geolocation)",
  meta_relgeo_z = "Acquisition (geolocation)",
  meta_leafon_z = "Acquisition (viewing and illumination)",
  view_az_cos_z = "Acquisition (viewing and illumination)",
  meta_offnad_z = "Acquisition (viewing and illumination)",
  meta_fwdrev_z = "Acquisition (viewing and illumination)",
  meta_sunel_z  = "Acquisition (viewing and illumination)",
  view_az_sin_z = "Acquisition (viewing and illumination)",
  meta_stereo_z = "Acquisition (viewing and illumination)",
  meta_az_conc_z= "Acquisition (viewing and illumination)"
)
# Top-to-bottom family order on the plot = Table 2 narrative order
# (Terrain -> Vegetation/canopy -> Acquisition). One-line change to reorder.
FAMILY_ORDER <- c("Terrain", "Canopy structure",
                  "Acquisition (geolocation)",
                  "Acquisition (viewing and illumination)")

# Fixed-effect R^2 for the caption (Table 3). Override or wire to a table.
fe_r2_chm <- 0.128
fe_r2_dtm <- 0.587

pf <- copy(all_fx)                                   # plotting copy; all_fx untouched
pf[, family := factor(FAMILY[predictor], levels = FAMILY_ORDER)]
pf[, mag    := pmax(abs(est_16), abs(est_dtm))]
setorder(pf, family, -mag)
pf[, yrow  := (.N + 1L) - seq_len(.N)]               # first (top) family = highest y
pf[, label := PRETTY[predictor]]

# alternating family bands
bands <- pf[, .(ymin = min(yrow) - 0.5, ymax = max(yrow) + 0.5,
                ytop = max(yrow) + 0.5), by = family]
bands[, fill := rep(c(p3d_gray$band_b, p3d_gray$band_a), length.out = .N)]

# long form. dodge = 0 -> both series sit EVEN on the row center, so the
# connector between them is purely horizontal (movement only along x).
dodge <- 0
long2 <- rbind(
  pf[, .(predictor, label, yrow, ypos = yrow + dodge,
         series = "CHM 16-site (primary)", est = est_16, lo = lo_16, hi = hi_16)],
  pf[, .(predictor, label, yrow, ypos = yrow - dodge,
         series = "DTM 18-site",           est = est_dtm, lo = lo_dtm, hi = hi_dtm)]
)
long2[, series := factor(series, levels = c("CHM 16-site (primary)", "DTM 18-site"))]

# CHM<->DTM connectors (horizontal, since both points share yrow); darker when
# the sign reverses so the zero-crossing links stand out.
conn <- pf[, .(predictor, x1 = est_16, y1 = yrow + dodge,
               x2 = est_dtm, y2 = yrow - dodge, rev = signrev_16_dtm)]
conn[, col := ifelse(rev, "#6E6E6E", "#CFCFCF")]

# -----------------------------------------------------------------------------
# 6. Build forest plot   (UPGRADE)
# -----------------------------------------------------------------------------
# Horizontal CHM<->DTM connectors ON by default. Because both series sit on the
# same row center (dodge = 0), the link moves only along x, so it reads as the
# coefficient shifting between models -- and crosses the zero line exactly for
# the reversals. Set FALSE to drop the links.
show_connectors <- TRUE

xlim <- c(-1.38, 1.22)   # widen if a future CI exceeds this; whiskers must fit

p <- ggplot() +
  geom_rect(data = bands,
            aes(xmin = xlim[1], xmax = xlim[2], ymin = ymin, ymax = ymax),
            fill = bands$fill, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, color = p3d_col$zero, linetype = "22", linewidth = 0.4) +
  (if (show_connectors)
     geom_segment(data = conn, aes(x = x1, xend = x2, y = y1, yend = y2),
                  color = conn$col, linewidth = 0.5, lineend = "round") else NULL) +
  geom_segment(data = long2, aes(x = lo, xend = hi, y = ypos, yend = ypos, color = series),
               linewidth = 0.55, lineend = "butt", alpha = 0.85) +
  geom_point(data = long2, aes(x = est, y = ypos, color = series, shape = series),
             size = 2.3) +
  geom_text(data = bands, aes(x = xlim[1] + 0.02, y = ytop - 0.06, label = family),
            hjust = 0, vjust = 1, size = 3.1, fontface = "bold",
            color = p3d_gray$strip, inherit.aes = FALSE) +
  scale_color_manual(values = p3d_series,
                     labels = c("CHM 16-site (primary)" = "CHM 16-site", "DTM 18-site" = "DTM 18-site")) +
  scale_shape_manual(values = c("CHM 16-site (primary)" = 16L, "DTM 18-site" = 17L)) +
  scale_y_continuous(breaks = pf$yrow, labels = pf$label,
                     expand = expansion(mult = c(0.01, 0.03))) +
  scale_x_continuous(limits = xlim, breaks = seq(-1, 1, 0.5),
                     expand = expansion(mult = 0)) +
  labs(x = "Standardized coefficient  (m per SD of predictor)",
       caption = sprintf(paste0(
         "%d of %d shared predictors reverse sign between the CHM and DTM models.\n",
         "Fixed-effect share of R\u00b2: CHM %d%%, DTM %d%%.   Bars: 95%% credible intervals."),
         n_rev_16, nrow(all_fx), round(100 * fe_r2_chm), round(100 * fe_r2_dtm))) +
  guides(color = guide_legend(override.aes = list(linetype = 0, size = 2.6)),
         shape = "none") +
  theme_p3d(base_size = 11, grid = "x") +
  coord_cartesian(clip = "off")

# Optional: fit stats in a panel corner instead of the caption. To use, drop
# the "Fixed-effect R^2" line from the caption above and uncomment:
# p <- p + annotate("text", x = xlim[2], y = 0.7, hjust = 1, vjust = 0, size = 3,
#                   color = p3d_gray$text_lo,
#                   label = sprintf("Fixed-effect R\u00b2\nCHM %d%%  \u00b7  DTM %d%%",
#                                   round(100*fe_r2_chm), round(100*fe_r2_dtm)))

# -----------------------------------------------------------------------------
# 7. Save (SAME dimensions as the committed three-series figure)
# -----------------------------------------------------------------------------
outdir <- if (exists("out_plots")) out_plots else "."
save_p3d(p, "fig03_coef_forest_twoseries_v117", width_in = 7.5, height_in = 7.3,
         dpi = 600, dir = outdir, formats = c("png", "pdf"))
# (Prefer your verified saver? swap the line above for:
#  save_figure(p, "fig03_coef_forest_twoseries", width_in = 7.5, height_in = 7.3))

# -----------------------------------------------------------------------------
# 8. VERIFICATION (console anchors) -- fingerprint before any swap   (VERBATIM)
# -----------------------------------------------------------------------------
log_subsection("VERIFICATION (console anchors)")
cat("\nMust match story-lock 2F headline values:\n")
cat("  slope_mean_z: 16-site -0.495, 19-site -0.439, DTM 18-site ?\n")
cat("  rh_98_z:      16-site -0.441, 19-site +0.069, DTM 18-site +0.99\n")
cat("  wsci_z:       16-site +0.317, 19-site -0.868, DTM 18-site ?\n")
cat("  cover_z:      16-site -0.083, 19-site -0.161, DTM 18-site ?\n\n")
cat("Observed (rounded to 3 dp):\n")
print(all_fx[, .(
  predictor = predictor,
  CHM_16    = round(est_16, 3),
  CHM_19    = round(est_19, 3),
  DTM_18    = round(est_dtm, 3),
  flip16    = signrev_16_dtm,
  flip19    = signrev_19_dtm
)])
cat(sprintf("\nHeadline counts (should match story-lock 2F):\n"))
cat(sprintf("  Sign reversals CHM 16 vs DTM 18: %d of %d (expected 14 of 17)\n",
            n_rev_16, nrow(all_fx)))
cat(sprintf("  Sign reversals CHM 19 vs DTM 18: %d of %d (expected 13 of 17)\n",
            n_rev_19, nrow(all_fx)))
cat("\nMAIN figure plots TWO series only (CHM 16-site, DTM 18-site);",
    "\nfigure caption shows the 14/17 line only.",
    "\nThe 19-site context remains on the supplementary three-series figure",
    "\n(reuses existing word/media/image3.png; no regeneration).\n")
log_progress("section_H_fig3_coef_twoseries.R complete.")
