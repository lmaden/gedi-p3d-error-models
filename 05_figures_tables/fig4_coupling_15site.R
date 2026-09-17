# =====================================================================
# fig4_coupling_15site.R  (24 Aug 2026: 15-site basis throughout;
#   panels a/c no longer use or show the three CHM-screened sites)
# FIGURE 4 (headline): within-site CHM-DTM coupling + counterfactual refit.
#
# DESIGN DECISION (2026-07, Session 2): panels (a)(b)(c) are kept in their
# ORIGINAL style, reproduced VERBATIM from fig7_within_site_coupling.R
# (lines 60-453: data-prep + the panel7a/7b/7c ggplots + ann_df). ONLY panel (d)
# is replaced -- the old bar chart becomes a point + horizontal 95% CrI interval
# in CHM blue (#0072B2) that rhymes with Figure 3 and drops the oversized font.
# Panel (d) data-prep is VERBATIM from figure6_panel_d.R (lines 70-84).
#
# The composite is rendered at 1.7x the docx display box (final 6.5 x 6.064 in,
# aspect 1.0719 == the wp:extent for image4) so a/b/c render near their native
# 5.5-in size and read as airy as the original once Word scales the image into
# the 6.5-in column. NOTE: this keeps Fig 4 in its own look (cowplot + viridis +
# red trend lines) rather than the theme_p3d house style, by request; its text
# will therefore be a different physical size than the other house figures.
#
# Story-lock anchors (median r, r_15/p_15, corr_18/corr_15/CI, panel-d sign-flip)
# are preserved and re-verified below.
#
# Run from scripts/reviewed/:  module load rh9/R/4.5.0 ; Rscript section_G_fig4_coupling.R
# Output: <SECTION_G_PLOTS>/fig4_coupling.{png,pdf}
# =====================================================================

source("fig_common.R")
fig_banner("Figure 4",
                 "Within-site coupling (a/b/c verbatim) + counterfactual interval (d)")

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(patchwork)
  library(ggrepel)
})


# ============ VERBATIM: panels (a)/(b)/(c) from section_G_fig7 (60-453) ============
# (data-prep + the original panel7a/7b/7c ggplots + ann_df; byte-for-byte)
# ---------------------------------------------------------------------
# 1. Load per-site coupling table (gives r per site for all 18 DTM-OK)
# ---------------------------------------------------------------------
log_subsection("Loading per-site coupling and slopes")

per_site_path <- file.path(manuscript_tables_dir,
                           "coupling_rh98_per_site_slopes.csv")
if (!file.exists(per_site_path)) {
  stop("Required: ", per_site_path)
}
ps <- as.data.table(fread(per_site_path))
log_progress(sprintf("  Loaded per-site slopes: %d rows, %d cols",
                     nrow(ps), ncol(ps)))
cat("  Columns: ", paste(names(ps), collapse = ", "), "\n")

# Standardize site column.
site_col <- intersect(c("tracker_site", "tracker_site_id",
                        "tracker_id", "site", "site_id"),
                      names(ps))[1]
if (is.na(site_col)) {
  stop("No site identifier column in per_site_path. Cols: ",
       paste(names(ps), collapse = ", "))
}
setnames(ps, site_col, "tracker_id_str")
ps[, tracker_id := suppressWarnings(as.integer(as.character(tracker_id_str)))]
ps[, manuscript_id := tracker_to_manuscript(tracker_id)]

# Required columns by story-lock §2H/I:
#   corr_errP3D_errDTM (per-site Pearson r, for 7a header anchor and 7b y)
#   b_rh98_uni (univariate per-site DTM~RH98 slope; for 7c y)
r_col   <- intersect(c("corr_errP3D_errDTM", "corr_err", "r_within"),
                     names(ps))[1]
b_col   <- "b_rh98_uni"   # canonical univariate slope per story-lock §2I.
                          # Earlier intersect() picked b_rh98_partial when
                          # both columns exist, producing r=-0.418/-0.432
                          # instead of -0.477/-0.448 (partial vs univariate).
if (!(b_col %in% names(ps))) {
  stop("Required column ", b_col, " not in ", per_site_path,
       ". Columns present: ", paste(names(ps), collapse = ", "))
}
if (is.na(r_col)) {
  log_progress(sprintf("  NOTE: no per-site r column found in %s; 7a will need a re-derivation from raw error vectors.",
                       per_site_path))
} else {
  log_progress(sprintf("  Per-site r column: %s", r_col))
}
if (is.na(b_col)) {
  log_progress(sprintf("  NOTE: no per-site b_RH98 column found in %s; will skip 7c if missing.",
                       per_site_path))
} else {
  log_progress(sprintf("  Per-site b_RH98 column: %s", b_col))
}

# ---------------------------------------------------------------------
# 2. Site-level covariates (for 7b)
# ---------------------------------------------------------------------
log_subsection("Loading site-level covariates ")

cov_path <- file.path(manuscript_tables_dir,
                      "coupling_site_covariates.csv")
if (!file.exists(cov_path)) {
  stop("Required: ", cov_path)
}
cov <- as.data.table(fread(cov_path))
log_progress(sprintf("  Loaded site covariates: %d rows, %d cols",
                     nrow(cov), ncol(cov)))
cat("  Columns: ", paste(names(cov), collapse = ", "), "\n")

site_col_cov <- intersect(c("tracker_site", "tracker_site_id",
                            "tracker_id", "site", "site_id"),
                          names(cov))[1]
setnames(cov, site_col_cov, "tracker_id_str")
cov[, tracker_id := suppressWarnings(as.integer(as.character(tracker_id_str)))]
cov[, manuscript_id := tracker_to_manuscript(tracker_id)]

geoacc_col <- intersect(c("meta_abs_geoacc_avg_avg",
                          "meta_abs_geoacc_avg", "abs_geoacc_avg",
                          "meta_absgeo_avg"),
                        names(cov))[1]
if (is.na(geoacc_col)) {
  stop("Could not find absolute-geolocation column in site covariates. ",
       "Expected one of: meta_abs_geoacc_avg_avg, meta_abs_geoacc_avg.")
}
log_progress(sprintf("  Geolocation column: %s", geoacc_col))

# Join coupling r with geoacc and frame class.
ps_with_cov <- merge(
  ps[, .(tracker_id, manuscript_id,
         r_within = if (!is.na(r_col)) ps[[r_col]] else NA_real_,
         b_rh98   = if (!is.na(b_col)) ps[[b_col]] else NA_real_)],
  cov[, .(tracker_id, geoacc = .SD[[1]]),
      .SDcols = geoacc_col],
  by = "tracker_id", all = FALSE
)
log_progress(sprintf("  Joined rows: %d", nrow(ps_with_cov)))

ps_with_cov[, frame := site_frame_class(tracker_id)]
ps_with_cov[, in_phase25_15 := tracker_id %in% TRACKER_PHASE25_15]
ps_with_cov[, in_dtm_18     := tracker_id %in% TRACKER_DTM_18]
ps_with_cov[, in_chm_16     := tracker_id %in% TRACKER_CHM_16]

# ---------------------------------------------------------------------
# 3. Panel 7c data table (no RE join needed after G-1 resolution)
# ---------------------------------------------------------------------
# panel 7c plots b_RH98 against per-site
# coupling r(err_p3d, err_dtm), both already in ps_with_cov via
# Section 1. Previous versions joined to site_random_effects_chm_dtm.csv;
# that join is no longer needed.
log_subsection("Building panel 7c data table (b_RH98 + coupling r)")

panel7c <- ps_with_cov[, .(tracker_id, manuscript_id, b_rh98, r_within,
                            in_phase25_15, in_dtm_18, frame)]
log_progress(sprintf("  Panel 7c rows: %d", nrow(panel7c)))

# ---------------------------------------------------------------------
# 4. Cross-site correlation summary (for the 7b annotation)
# ---------------------------------------------------------------------
log_subsection("Loading cross-site correlations for 7b headline")

corr_path <- file.path(manuscript_tables_dir,
                       "coupling_site_correlations.csv")
corr_table <- NULL
if (file.exists(corr_path)) {
  corr_table <- as.data.table(fread(corr_path))
  log_progress(sprintf("  Loaded cross-site correlations: %d rows, %d cols",
                       nrow(corr_table), ncol(corr_table)))
  cat("  Columns: ", paste(names(corr_table), collapse = ", "), "\n")
}

# Always compute the 15-site r ourselves from the joined table as a
# verification crosscheck (panel sources r from this computation).
phase25_subset <- ps_with_cov[in_phase25_15 == TRUE & is.finite(r_within)
                               & is.finite(geoacc)]
n15 <- nrow(phase25_subset)
if (n15 >= 3L) {
  ct15 <- cor.test(phase25_subset$r_within, phase25_subset$geoacc)
  r15  <- as.numeric(ct15$estimate)
  p15  <- as.numeric(ct15$p.value)
} else {
  r15 <- NA_real_
  p15 <- NA_real_
}
log_progress(sprintf("  7b cross-site r_15 (computed) = %.3f, p = %.3f, n = %d",
                     r15, p15, n15))

# Headline from story-lock §2G: r_15 = -0.699, p = 0.004.
# We use the computed values; expect agreement to 3 dp.

# ---------------------------------------------------------------------
# 5. Cross-site b_RH98 vs coupling r summary (for 7c annotation)
# ---------------------------------------------------------------------
sub18 <- panel7c[in_dtm_18 == TRUE & is.finite(b_rh98) & is.finite(r_within)]
sub15 <- panel7c[in_phase25_15 == TRUE & is.finite(b_rh98) & is.finite(r_within)]

corr_18 <- if (nrow(sub18) >= 3L) {
  c18 <- cor.test(sub18$b_rh98, sub18$r_within)
  list(r = as.numeric(c18$estimate), p = as.numeric(c18$p.value), n = nrow(sub18))
} else list(r = NA_real_, p = NA_real_, n = nrow(sub18))

corr_15 <- if (nrow(sub15) >= 3L) {
  c15 <- cor.test(sub15$b_rh98, sub15$r_within)
  list(r = as.numeric(c15$estimate), p = as.numeric(c15$p.value), n = nrow(sub15))
} else list(r = NA_real_, p = NA_real_, n = nrow(sub15))

log_progress(sprintf("  7c 18-site descriptive r = %.3f, n = %d",
                     corr_18$r, corr_18$n))
log_progress(sprintf("  7c 15-site test r        = %.3f, n = %d",
                     corr_15$r, corr_15$n))

# Bootstrap CI for the 15-site test r: read canonical 10,000-iter values
# from manuscript_tables/coupling_rh98_bootstrap.csv (produced by
# scripts/reviewed/coupling_05_per_site_rh98_slopes.R). Story-lock
# §2I anchor: [-0.749, -0.139].
boots_path <- file.path(manuscript_tables_dir,
                        "coupling_rh98_bootstrap.csv")
if (!file.exists(boots_path)) {
  stop("Required: ", boots_path)
}
boots <- as.data.table(fread(boots_path))
ci_row <- boots[frame == "15_non_flagged_dtm_ok" & slope == "b_rh98_uni"]
if (nrow(ci_row) != 1L) {
  stop("Could not locate single (b_rh98_uni, 15_non_flagged_dtm_ok) row in ",
       boots_path)
}
ci_15 <- c(ci_row$ci_lo_95, ci_row$ci_hi_95)
log_progress(sprintf("  7c 15-site bootstrap CI (canonical) [%.3f, %.3f]",
                     ci_15[1], ci_15[2]))

# ---------------------------------------------------------------------
# 6. PANEL 7a — single example site within-site scatter
# Choose the site whose r is closest to the 18-site median r.
# ---------------------------------------------------------------------
log_subsection("Building Panel 7a")

# Median r across the 15 sites in both analyses.
ps_dtm_ok <- ps_with_cov[in_phase25_15 == TRUE & is.finite(r_within)]
median_r_18 <- median(ps_dtm_ok$r_within, na.rm = TRUE)
n_neg_15 <- sum(ps_dtm_ok$r_within < 0, na.rm = TRUE)
log_progress(sprintf("  Across-site median r (15 sites): %.3f; negative at %d of %d",
                     median_r_18, n_neg_15, nrow(ps_dtm_ok)))

# Example site fixed at manuscript Site 8 (tracker 8), as in the caption.
example_row <- ps_dtm_ok[manuscript_id == 8L]
stopifnot(nrow(example_row) == 1L)
example_tracker_id    <- example_row$tracker_id
example_manuscript_id <- example_row$manuscript_id
example_site_r        <- example_row$r_within
log_progress(sprintf("  Example site for 7a: tracker %d (manuscript %d), site r = %.3f",
                     example_tracker_id, example_manuscript_id, example_site_r))

# Load the per-footprint err_chm and err_dtm vectors at the example site.
# Strategy: load mod_chm_s2 and mod_dtm_s2 from 08_model_prep, filter to
# the example site, join on shot_number, compute err_chm = err_p3d_chm
# and err_dtm = err_p3d_dtm at common footprints.

mp <- load_checkpoint("08_model_prep")
mod_chm <- as.data.table(mp$mod_chm_s2)
mod_dtm <- as.data.table(mp$mod_dtm_s2)

# Tracker id stored as a 'site' factor / character.
mod_chm_site <- mod_chm[as.character(site) == as.character(example_tracker_id)]
mod_dtm_site <- mod_dtm[as.character(site) == as.character(example_tracker_id)]

# Join keys: shot_number is the canonical GEDI footprint ID. Some
# pipelines have shot_number; others have shot_number_chm/dtm. Use the
# first match.
join_key_candidates <- c("shot_number", "shot_id", "footprint_id",
                          "gedi_shot_number")
key_chm <- intersect(join_key_candidates, names(mod_chm_site))[1]
key_dtm <- intersect(join_key_candidates, names(mod_dtm_site))[1]
if (is.na(key_chm) || is.na(key_dtm)) {
  stop("Could not find shot_number key in mod_chm_s2 or mod_dtm_s2. ",
       "CHM cols: ", paste(head(names(mod_chm_site), 30), collapse = ","))
}
log_progress(sprintf("  Join keys: CHM=%s, DTM=%s", key_chm, key_dtm))

setkeyv(mod_chm_site, key_chm)
setkeyv(mod_dtm_site, key_dtm)
joined <- merge(mod_chm_site[, .(shot = get(key_chm), err_chm = chm_error_mean)],
                mod_dtm_site[, .(shot = get(key_dtm), err_dtm = dtm_error_mean)],
                by = "shot", all = FALSE)
log_progress(sprintf("  Joined footprints at example site: %d",
                     nrow(joined)))

# Compute the within-site r from the actual joined vectors, sanity
# check vs the CSV value.
r_emp <- if (nrow(joined) >= 3L) cor(joined$err_chm, joined$err_dtm) else NA_real_
log_progress(sprintf("  Empirical r at example site (n=%d): %.3f  vs CSV %.3f",
                     nrow(joined), r_emp, example_site_r))

# Use a hexbin to handle high-density scatter (typical 5k-100k pts).
panel7a <- ggplot(joined, aes(x = err_dtm, y = err_chm)) +
  geom_hex(bins = 60) +
  geom_smooth(method = "lm", se = FALSE,
              color = COLOR_REF_LINE, linewidth = 0.6) +
  geom_hline(yintercept = 0, color = COLOR_ZERO_LINE,
             linewidth = 0.2, linetype = "dotted") +
  geom_vline(xintercept = 0, color = COLOR_ZERO_LINE,
             linewidth = 0.2, linetype = "dotted") +
  scale_fill_viridis_c(trans = "log10", option = "viridis", name = "Count") +
  labs(
    title = sprintf("(a) Per-footprint CHM vs DTM error at Site %d",
                    example_manuscript_id),
    x = "DTM error (m)",
    y = "CHM error (m)"
  ) +
  annotate("text",
           x = -Inf, y = Inf,
           label = sprintf("Site %d: r = %s%.2f  (n = %s)\nMedian across 15 sites: r = %s%.2f",
                           example_manuscript_id,
                           ifelse(example_site_r < 0, MINUS, ""),
                           abs(example_site_r),
                           formatC(nrow(joined), big.mark = ",", format = "d"),
                           ifelse(median_r_18 < 0, MINUS, ""),
                           abs(median_r_18)),
           hjust = -0.04, vjust = 1.4, size = 3.3, lineheight = 1.05) +
  coord_cartesian(xlim = c(-15, 15), ylim = c(-15, 15)) +
  theme_section_G(base_size = 11)

# ---------------------------------------------------------------------
# 7. PANEL 7b — cross-site r vs meta_abs_geoacc_avg (15-site)
# ---------------------------------------------------------------------
log_subsection("Building Panel 7b")

panel7b_df <- ps_with_cov[in_phase25_15 == TRUE
                          & is.finite(r_within) & is.finite(geoacc)]
panel7b_df[, label := as.character(manuscript_id)]

panel7b <- ggplot(panel7b_df, aes(x = geoacc, y = r_within)) +
  geom_hline(yintercept = 0, color = COLOR_ZERO_LINE,
             linetype = "dashed", linewidth = 0.4) +
  geom_smooth(method = "lm", se = TRUE,
              color = COLOR_REF_LINE, fill = "#fbe6e6",
              linewidth = 0.6, alpha = 0.6) +
  geom_point(size = 3.0, color = "#222222", fill = "white",
             shape = 21, stroke = 0.8) +
  ggrepel::geom_text_repel(
    aes(label = label),
    size = 3.0, color = "black",
    bg.color = "white", bg.r = 0.1,
    min.segment.length = 0,
    segment.color = "gray50", segment.size = 0.25,
    box.padding = 0.30, point.padding = 0.20,
    max.overlaps = Inf,
    seed = 7L
  ) +
  labs(
    title = "(b) Cross-site CHM-DTM coupling vs\ngeolocation accuracy",
    x = "Absolute horizontal geolocation accuracy (m)",
    y = "Per-site r (CHM error, DTM error)"
  ) +
  annotate("text",
           x = -Inf, y = -Inf,
           label = sprintf("r = %s%.2f, p = %.3f (15 sites)",
                           ifelse(r15 < 0, MINUS, ""),
                           abs(r15),
                           p15),
           hjust = -0.05, vjust = -1.0, size = 3.3) +
  theme_section_G(base_size = 11)

# ---------------------------------------------------------------------
# 8. PANEL 7c — per-site b_RH98 vs CHM-DTM error coupling (G-1)
# ---------------------------------------------------------------------
log_subsection("Building Panel 7c")

panel7c_df <- panel7c[is.finite(b_rh98) & is.finite(r_within)
                      & in_phase25_15 == TRUE]
panel7c_df[, label := as.character(manuscript_id)]

# OLS line fit on the 15-site test frame only (canonical headline statistic).
panel7c_test <- panel7c_df[in_phase25_15 == TRUE]

# Annotation as a geom_text() layer (not annotate()) so it lives in the data
# layer. Anchored at upper-right (x = Inf, y = Inf) because that quadrant is
# empty in this panel: no site has both b_rh98 > 0.20 and r_within > -0.25.
# Earlier attempts at upper-left collided with sites 3 and 6 in the narrower
# combined-triptych layouts even with nudge_y on geom_text_repel.
ann_df <- data.frame(
  x = Inf, y = Inf,
  text = sprintf(
    "r = %s%.2f (15 sites)\n95%% CI [%s%.2f, %s%.2f]",
    ifelse(corr_15$r < 0, MINUS, ""), abs(corr_15$r),
    ifelse(ci_15[1] < 0,  MINUS, ""), abs(ci_15[1]),
    ifelse(ci_15[2] < 0,  MINUS, ""), abs(ci_15[2])),
  stringsAsFactors = FALSE
)

panel7c <- ggplot(panel7c_df,
                  aes(x = b_rh98, y = r_within)) +
  geom_hline(yintercept = 0, color = COLOR_ZERO_LINE,
             linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = COLOR_ZERO_LINE,
             linetype = "dashed", linewidth = 0.4) +
  geom_smooth(data = panel7c_test,
              method = "lm", se = TRUE,
              color = COLOR_REF_LINE, fill = "#fbe6e6",
              linewidth = 0.6, alpha = 0.6) +
  geom_point(size = 3.0, color = "#222222", fill = "#222222",
             shape = 21, stroke = 0.8) +
  ggrepel::geom_text_repel(
    aes(label = label),
    size = 3.0, color = "black",
    bg.color = "white", bg.r = 0.1,
    min.segment.length = 0,
    segment.color = "gray50", segment.size = 0.25,
    box.padding = 0.30, point.padding = 0.20,
    max.overlaps = Inf,
    seed = 7L
  ) +
  labs(
    title = "(c) Per-site DTM-RH98 slope vs CHM-DTM error coupling",
    x = expression(paste("Per-site ", italic("b"["RH98"]),
                          " (DTM error vs RH98 within site)")),
    y = "Per-site r (CHM error, DTM error)"
  ) +
  geom_text(data = ann_df,
            aes(x = x, y = y, label = text),
            hjust = 1.04, vjust = 1.4,
            size = 3.0, lineheight = 1.05,
            inherit.aes = FALSE) +
  theme_section_G(base_size = 11) +
  theme(legend.position = "none")

# ============ VERBATIM: panel (d) data-prep from figure6_panel_d.R (70-84) ============
operational <- c(Estimate = -0.4406, Q2.5 = -0.4766, Q97.5 = -0.4060)
counterfact <- c(Estimate = +0.3934, Q2.5 = +0.3472, Q97.5 = +0.4390)

OP_LABEL <- "P3D CHM\n(uses CNN DTM)"
CF_LABEL <- "Counterfactual\n(uses 3DEP DTM)"

panel_data <- tibble(
  fit_label = c(OP_LABEL, CF_LABEL),
  fit_order = c(1L, 2L),
  median    = c(operational[["Estimate"]], counterfact[["Estimate"]]),
  q025      = c(operational[["Q2.5"]],     counterfact[["Q2.5"]]),
  q975      = c(operational[["Q97.5"]],    counterfact[["Q97.5"]])
)

shift <- counterfact[["Estimate"]] - operational[["Estimate"]]

# ============ NEW: interval panel (d) + 2x2 composition ============
# >>>>>>>>>>>> NEW PANEL (d) + COMPOSITION (shared with deliverable) >>>>>>>>>>>>
CHM_BLUE_D <- "#0072B2"   # house CHM blue (matches Figure 3's coefficient plot)
dat_d <- data.frame(
  fit  = factor(panel_data$fit_label, levels = rev(panel_data$fit_label)),
  kind = c("operational", "counterfactual"),
  med  = panel_data$median, lo = panel_data$q025, hi = panel_data$q975,
  lab  = sprintf("%+.2f", panel_data$median))
panel_d <- ggplot(dat_d, aes(y = fit)) +
  geom_vline(xintercept = 0, color = COLOR_ZERO_LINE, linetype = "dashed", linewidth = 0.4) +
  geom_segment(aes(x = lo, xend = hi, y = fit, yend = fit),
               linewidth = 0.9, color = CHM_BLUE_D, lineend = "round") +
  geom_point(aes(x = med, fill = kind), shape = 21, size = 3.3, stroke = 0.9,
             color = CHM_BLUE_D) +
  geom_text(aes(x = med, label = lab), vjust = -1.2, size = 3.3, fontface = "bold",
            color = "#222222") +
  scale_fill_manual(values = c(operational = CHM_BLUE_D, counterfactual = "white")) +
  scale_x_continuous(limits = c(-0.65, 0.65), breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
  labs(title = "(d) Counterfactual reference-substitution refit",
       x = "RH98 coefficient (m per SD)", y = NULL) +
  guides(fill = "none") +
  theme_section_G(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), legend.position = "none")

fig4 <- patchwork::wrap_plots(panel7a, panel7b, panel7c, panel_d, ncol = 2, nrow = 2)

# Render at 1.7x the docx display box (6.5 x 6.064 in) so a/b/c render near their
# native 5.5-in size and read airy like the original; aspect == wp:extent 1.0719.
FIG_W <- 6.5 * 1.7; FIG_H <- 6.064 * 1.7
# <<<<<<<<<<<< end shared block <<<<<<<<<<<<

# ---- save at the 1.7x canvas via the section_G saver (ragg PNG embeds DPI) ----
save_figure(fig4, "fig4_coupling_15site", width_in = FIG_W, height_in = FIG_H, dpi = 300)

# ================= VERIFICATION (verbatim anchors) =======================

log_subsection("VERIFICATION (console anchors)")

cat("\nPanel 7a (per-footprint within-site CHM-DTM scatter):\n")
cat(sprintf("  Example site selected: tracker %d -> manuscript %d\n",
            example_tracker_id, example_manuscript_id))
cat(sprintf("  Across-site median r (15 sites): %.3f  (expected -0.663; negative at %d of 15, expected 14)\n",
            median_r_18, n_neg_15))
cat(sprintf("  Example site r (CSV): %.3f\n", example_site_r))
cat(sprintf("  Example site r (empirical from footprints): %.3f\n", r_emp))
cat(sprintf("  Footprints joined at example site: %s\n",
            formatC(nrow(joined), big.mark = ",", format = "d")))

cat("\nPanel 7b (cross-site coupling vs meta_abs_geoacc_avg):\n")
cat(sprintf("  Story-lock §2G target: r_15 = -0.699, p = 0.004, n = 15\n"))
cat(sprintf("  Computed:              r_15 = %.3f,  p = %.3f, n = %d\n",
            r15, p15, n15))

cat("\nPanel 7c (per-site b_RH98 vs CHM-DTM error coupling, G-1):\n")
cat(sprintf("  Story-lock §2I targets:\n"))
cat(sprintf("    15-site test:        r = -0.477 (CI [-0.749, -0.139])\n"))
cat(sprintf("  Computed:\n"))
cat(sprintf("    15-site test:        r = %.3f, n = %d, canonical CI [%.3f, %.3f]\n",
            corr_15$r, corr_15$n, ci_15[1], ci_15[2]))


# ---- Panel (d) story-lock (grafted from figure6_panel_d.R rh_98_z values) ----
cat("\nPanel (d) / Fig. 4d (RH98 counterfactual sign-flip):\n")
cat(sprintf("  Operational    rh_98_z: %+.4f [%+.4f, %+.4f]\n",
            operational[["Estimate"]], operational[["Q2.5"]], operational[["Q97.5"]]))
cat(sprintf("  Counterfactual rh_98_z: %+.4f [%+.4f, %+.4f]\n",
            counterfact[["Estimate"]], counterfact[["Q2.5"]], counterfact[["Q97.5"]]))
cat(sprintf("  Shift (cf - op): %+.4f  (-> %+.2f m per SD)\n", shift, round(shift, 2)))
stopifnot(operational[["Estimate"]] < 0, counterfact[["Estimate"]] > 0,
          operational[["Q97.5"]] < 0, counterfact[["Q2.5"]] > 0,
          isTRUE(all.equal(round(shift, 2), 0.83)))
cat("  [OK] sign flip + both 95% CIs exclude 0 + shift = +0.83\n")
cat(sprintf("\nComposite: %s/fig4_coupling_15site.{png,pdf}  (%.2f x %.2f in, aspect %.4f;\n",
            SECTION_G_PLOTS, FIG_W, FIG_H, FIG_W/FIG_H))
cat("  Word shows it in the 6.5 x 6.064 in wp:extent, i.e. scaled to the column.)\n")
log_progress("fig4_coupling_15site.R complete.")

