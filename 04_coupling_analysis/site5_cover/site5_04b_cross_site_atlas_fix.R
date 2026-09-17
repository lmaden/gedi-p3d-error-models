#!/usr/bin/env Rscript
# =============================================================================
# site5_04b_cross_site_atlas_fix.R
#
# Patch for the step 04 atlas script. The step 04 per-site processing
# was correct, but the cross-site b_cover correlation block produced
# all-NA results because the value-column candidate lists for the three
# RE CSVs didn't match the actual column names.
#
# This patch reads the already-written per-site CSVs (no recomputation),
# resolves RE columns with KNOWN correct names, and overwrites only the
# affected outputs:
#   - groundwork_task6_phase2_cross_site.csv     (overwrite)
#   - groundwork_task6_phase2_correlations.csv   (overwrite)
#   - task6_phase2_bcover_vs_re_panels.pdf       (regenerate, now 4 panels)
#
# The decomposition output carries both RE_full and RE_alt18
# is the err_alt 15-site variant, so we grab BOTH and report each.
#
# Run: source("site5_04b_cross_site_atlas_fix.R")
# Wall-clock: ~5 seconds.
# =============================================================================

source("analysis_config.R")
source("analysis_utils.R")

log_section("Step 04 Cross-Site Patch")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

MANUSCRIPT_TB <- file.path(PROJECT_ROOT, "manuscript_tables")
PLOT_DIR      <- file.path(PROJECT_ROOT, "plots", "groundwork")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Read existing per-site CSV (no recomputation)
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Reading existing per-site outputs")

site_summary_dt <- fread(file.path(MANUSCRIPT_TB,
                                    "groundwork_task6_phase2_site_summary.csv"))
log_progress(sprintf("  Site summary: %d rows", nrow(site_summary_dt)))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Read RE CSVs with confirmed column names
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Reading RE CSVs (with corrected column names)")

# Ch1 site REs
ch1_path <- file.path(MANUSCRIPT_TB, "site_random_effects_chm_dtm.csv")
ch1_raw  <- fread(ch1_path)
log_progress(sprintf("  Ch1: %s", paste(names(ch1_raw), collapse = ", ")))
ch1_re <- data.table(
  tracker_site = as.integer(ch1_raw$tracker_site_id),
  ch1_chm_re   = as.numeric(ch1_raw$chm_re_mean)
)
ch1_re <- ch1_re[!is.na(tracker_site) & !is.na(ch1_chm_re)]
log_progress(sprintf("    -> %d Ch1 site REs loaded (key=tracker_site_id, value=chm_re_mean)",
                     nrow(ch1_re)))

# 16-site REs
trackb_path <- file.path(MANUSCRIPT_TB, "track_b_re_intercepts.csv")
trackb_raw  <- fread(trackb_path)
log_progress(sprintf("  16-site: %s", paste(names(trackb_raw), collapse = ", ")))
trackb_re <- data.table(
  tracker_site = as.integer(trackb_raw$tracker_site),
  track_b_re   = as.numeric(trackb_raw$track_b_chm_re)
)
trackb_re <- trackb_re[!is.na(tracker_site) & !is.na(track_b_re)]
log_progress(sprintf("    -> %d 16-site site REs loaded (key=tracker_site, value=track_b_chm_re)",
                     nrow(trackb_re)))

# Decomposition REs, both RE_full and RE_alt18
phase3_path <- file.path(MANUSCRIPT_TB, "coupling_random_effect_shifts.csv")
phase3_raw  <- fread(phase3_path)
log_progress(sprintf("  decomposition: %s", paste(names(phase3_raw), collapse = ", ")))
log_progress(sprintf("    -> decomposition output has %d rows (15 or 18 sites)",
                     nrow(phase3_raw)))
phase3_re <- data.table(
  tracker_site     = as.integer(phase3_raw$site),  # confirmed: site = tracker site #
  phase3_re_full   = as.numeric(phase3_raw$RE_full),
  phase3_re_alt18  = as.numeric(phase3_raw$RE_alt18),
  phase3_delta     = if ("delta" %in% names(phase3_raw)) as.numeric(phase3_raw$delta) else NA_real_
)
phase3_re <- phase3_re[!is.na(tracker_site)]
log_progress(sprintf("    -> %d decomposition site REs loaded (key=site, both RE_full and RE_alt18 grabbed)",
                     nrow(phase3_re)))

# Sanity check: print which sites the decomposition covers
log_progress(sprintf("    -> decomposition tracker_site values present: %s",
                     paste(sort(phase3_re$tracker_site), collapse = ", ")))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Build cross-site table
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Building cross-site table")

cross_site <- site_summary_dt[, .(
  manuscript_site, tracker_site, site_short, flag_status, dominant_lc,
  n_tiles_valid, tile_offset_mean, tile_offset_sd,
  b_cover, se_b_cover, t_b_cover, p_b_cover, reg_r2_adj
)]

cross_site <- merge(cross_site, ch1_re,    by = "tracker_site", all.x = TRUE)
cross_site <- merge(cross_site, trackb_re, by = "tracker_site", all.x = TRUE)
cross_site <- merge(cross_site, phase3_re, by = "tracker_site", all.x = TRUE)

setorder(cross_site, manuscript_site)

# Print the cross-site table for visual inspection
cat("\nCross-site table (RE columns merged):\n")
print(cross_site[, .(
  ms = manuscript_site, site = site_short, flag = flag_status,
  b_cover = round(b_cover, 2),
  ch1 = round(ch1_chm_re, 3),
  trk_b = round(track_b_re, 3),
  ph3_full = round(phase3_re_full, 3),
  ph3_alt18 = round(phase3_re_alt18, 3)
)])

fwrite(cross_site, file.path(MANUSCRIPT_TB, "groundwork_task6_phase2_cross_site.csv"))
log_progress(sprintf("  Wrote groundwork_task6_phase2_cross_site.csv (%d rows)",
                     nrow(cross_site)))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Cross-site correlations across frames × predictors × responses
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Running cross-site correlations")

frames <- list(
  list(name = "19_sites",
       mask = rep(TRUE, nrow(cross_site))),
  list(name = "16_non_flagged",
       mask = cross_site$flag_status != "FLAGGED"),
  list(name = "15_err_alt",
       mask = cross_site$flag_status != "FLAGGED" &
              cross_site$flag_status != "DTM_excluded"),
  list(name = "forest_only_no_flagged",
       mask = cross_site$flag_status != "FLAGGED" &
              cross_site$dominant_lc %in% c("BDF", "ENF", "DNF", "EBF", "MFT", "IWL"))
)

predictors <- c("b_cover", "tile_offset_sd", "tile_offset_mean")
responses  <- c("ch1_chm_re", "track_b_re", "phase3_re_full", "phase3_re_alt18")

corr_rows <- list()
for (fr in frames) {
  sub <- cross_site[fr$mask]
  for (p in predictors) {
    for (r in responses) {
      x <- sub[[p]]; y <- sub[[r]]
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 4L) {
        corr_rows[[length(corr_rows) + 1]] <- data.table(
          frame = fr$name, predictor = p, response = r,
          n = sum(ok),
          r_pearson = NA_real_, p_pearson = NA_real_,
          r_spearman = NA_real_, p_spearman = NA_real_
        )
        next
      }
      ct_p <- tryCatch(cor.test(x[ok], y[ok], method = "pearson"),
                        error = function(e) NULL)
      ct_s <- tryCatch(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE),
                        error = function(e) NULL)
      corr_rows[[length(corr_rows) + 1]] <- data.table(
        frame = fr$name, predictor = p, response = r,
        n = sum(ok),
        r_pearson  = if (!is.null(ct_p)) unname(ct_p$estimate) else NA_real_,
        p_pearson  = if (!is.null(ct_p)) ct_p$p.value          else NA_real_,
        r_spearman = if (!is.null(ct_s)) unname(ct_s$estimate) else NA_real_,
        p_spearman = if (!is.null(ct_s)) ct_s$p.value          else NA_real_
      )
    }
  }
}
corr_dt <- rbindlist(corr_rows)

fwrite(corr_dt, file.path(MANUSCRIPT_TB, "groundwork_task6_phase2_correlations.csv"))
log_progress(sprintf("  Wrote groundwork_task6_phase2_correlations.csv (%d rows)",
                     nrow(corr_dt)))

# Print headline: b_cover vs each RE in each frame
cat("\nb_cover correlation against each RE column (Pearson):\n")
print(corr_dt[predictor == "b_cover",
              .(frame, response, n,
                r = round(r_pearson, 3),
                p = round(p_pearson, 4))])

cat("\ntile_offset_sd correlation against each RE column (Pearson):\n")
print(corr_dt[predictor == "tile_offset_sd",
              .(frame, response, n,
                r = round(r_pearson, 3),
                p = round(p_pearson, 4))])

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Regenerate scatter panels (now 4 panels)
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Regenerating bcover_vs_re_panels figure (now 4 panels)")

mk_scatter <- function(dt, re_col, re_label, frame_label) {
  d <- dt[is.finite(b_cover) & is.finite(get(re_col))]
  if (nrow(d) < 3L) {
    return(ggplot() + theme_void() +
             labs(title = sprintf("%s\n(insufficient data, n=%d)", re_label, nrow(d))))
  }
  ct <- cor.test(d$b_cover, d[[re_col]])
  r_pe <- unname(ct$estimate); p_pe <- ct$p.value
  ggplot(d, aes(x = b_cover, y = get(re_col), color = flag_status)) +
    geom_point(size = 2.5) +
    geom_text(aes(label = manuscript_site), nudge_x = 0.5, size = 3,
              color = "grey20") +
    geom_smooth(method = "lm", se = TRUE, color = "grey30",
                linetype = "dashed", inherit.aes = FALSE,
                aes(x = b_cover, y = get(re_col))) +
    geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
    geom_vline(xintercept = 0, color = "grey50", linetype = "dotted") +
    scale_color_manual(values = c("FLAGGED" = "#d6604d", "ok" = "#4393c3",
                                   "DTM_excluded" = "#bf812d"),
                       name = "Flag status") +
    labs(x = "b_cover (m / unit cover)",
         y = re_label,
         title = sprintf("%s\n%s: r = %+.2f, p = %.3f, n = %d",
                          re_label, frame_label, r_pe, p_pe, nrow(d))) +
    theme_cowplot(10) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10))
}

panel_ch1   <- mk_scatter(cross_site,
                           "ch1_chm_re",       "Ch1 CHM site RE",
                           "19 sites")
panel_trkb  <- mk_scatter(cross_site[flag_status != "FLAGGED"],
                           "track_b_re",       "16-site RE",
                           "16 non-flagged")
panel_ph3f  <- mk_scatter(cross_site,
                           "phase3_re_full",   "Decomposition RE_full",
                           sprintf("n = %d", sum(is.finite(cross_site$phase3_re_full))))
panel_ph3a  <- mk_scatter(cross_site,
                           "phase3_re_alt18",  "Decomposition RE_alt18",
                           sprintf("n = %d", sum(is.finite(cross_site$phase3_re_alt18))))

legend_dummy <- ggplot(cross_site, aes(x = b_cover, y = b_cover, color = flag_status)) +
  geom_point() +
  scale_color_manual(values = c("FLAGGED" = "#d6604d", "ok" = "#4393c3",
                                 "DTM_excluded" = "#bf812d"),
                     name = "Flag status") +
  theme_cowplot(10) + theme(legend.position = "top")
leg <- get_legend(legend_dummy)

panels_grid <- plot_grid(panel_ch1, panel_trkb, panel_ph3f, panel_ph3a,
                          ncol = 2, align = "hv")
panels_combined <- plot_grid(leg, panels_grid, ncol = 1, rel_heights = c(0.05, 1))

ggsave(file.path(PLOT_DIR, "task6_phase2_bcover_vs_re_panels.pdf"),
       panels_combined, width = 12, height = 10, bg = "white")
log_progress("  Wrote task6_phase2_bcover_vs_re_panels.pdf (4-panel)")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Final summary
# ─────────────────────────────────────────────────────────────────────────────

log_section("Cross-site patch — Complete")

cat("\nHEADLINE TEST: r(b_cover, RE) — predicted sign NEGATIVE\n")
cat("\nThe four candidate RE responses, in the relevant frame for each:\n\n")

headline <- corr_dt[predictor == "b_cover" &
  ((frame == "19_sites"        & response == "ch1_chm_re") |
   (frame == "16_non_flagged"  & response == "track_b_re") |
   (frame == "15_err_alt"      & response %in% c("phase3_re_full", "phase3_re_alt18")))]
print(headline[, .(frame, response, n,
                   r_pearson  = round(r_pearson, 3),
                   p_pearson  = round(p_pearson, 4),
                   r_spearman = round(r_spearman, 3))])

log_progress("")
log_progress("Outputs:")
log_progress("  - groundwork_task6_phase2_cross_site.csv")
log_progress("  - groundwork_task6_phase2_correlations.csv")
log_progress("  - task6_phase2_bcover_vs_re_panels.pdf")
log_progress("  - This script's stdout/stderr log")
log_progress("")
log_progress("(All other step 04 outputs are unchanged from the original run")
log_progress(" and don't need to be re-uploaded.)")
