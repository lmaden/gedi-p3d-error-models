#!/usr/bin/env Rscript
# =============================================================================
# site5_02_bootstrap_correlations.R
#
# step 02 robustness check on Step 04 cross-site correlations.
#
# step 04 produced point-estimate correlations between per-site b_cover
# and Ch1 / 16-site / decomposition site REs. With n = 13-17 in the relevant
# frames, parametric Pearson p-values are unreliable for small samples
# and the 16-site test came back at p = 0.14 (non-significant
# under the strict variance-puzzle decision rule). This script:
#
#   1. Computes site-level bootstrap 95% CIs for every (frame x predictor
#      x response) combination from step 04 (10,000 iterations).
#   2. Runs leave-one-out influence analysis on three headline
#      correlations to identify leverage points.
#   3. Reports BCa CI where it converges (small n cases sometimes don't).
#   4. Generates a forest plot of correlations + CIs and a LOO
#      diagnostic figure.
#
# Notes:
#   - Site-level resampling (sample 19/16/13 sites with replacement, treat
#     each site's b_cover as a fixed point estimate). Hierarchical
#     bootstrap (resample footprints within sites) was considered but
#     within-site b_cover SEs are tiny relative to cross-site spread
#     (e.g. Site 1: SE 0.09 on b_cover -13.91), so site-level is correct.
#   - tile_offset_mean is included for completeness but flagged in the
#     output: its correlation with site REs is mechanistically circular
#     because the per-site RE is essentially the per-site offset.
#   - phase3_re_full is identical to ch1_chm_re (confirmed in inspection
#     of the cross-site CSV); both kept for completeness but flagged.
#
# Inputs:
#   - manuscript_tables/groundwork_task6_phase2_cross_site.csv
#   - manuscript_tables/groundwork_task6_phase2_correlations.csv
#
# Outputs (manuscript_tables/):
#   - groundwork_task6_phase2_correlations_with_ci.csv
#   - groundwork_task6_phase2_loo_influence.csv
#
# Outputs (plots/groundwork/):
#   - task6_phase1_5_bootstrap_forest.pdf      forest plot of CIs
#   - task6_phase1_5_loo_influence.pdf         LOO diagnostic 3-panel
#
# Run: source("site5_02_bootstrap_correlations.R")
# Wall-clock: ~30 seconds.
# =============================================================================

source("analysis_config.R")
source("analysis_utils.R")

log_section("Step 02: Bootstrap CIs + LOO Influence")

set.seed(2026)

suppressPackageStartupMessages({
  library(data.table)
  library(boot)
  library(ggplot2)
  library(cowplot)
})

MANUSCRIPT_TB <- file.path(PROJECT_ROOT, "manuscript_tables")
PLOT_DIR      <- file.path(PROJECT_ROOT, "plots", "groundwork")

N_BOOT <- 10000L

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Read inputs
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Reading step 04 outputs")

cross_site <- fread(file.path(MANUSCRIPT_TB, "groundwork_task6_phase2_cross_site.csv"))
log_progress(sprintf("  Cross-site table: %d rows", nrow(cross_site)))

orig_corr <- fread(file.path(MANUSCRIPT_TB, "groundwork_task6_phase2_correlations.csv"))
log_progress(sprintf("  Original correlations: %d rows", nrow(orig_corr)))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Bootstrap helpers
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Defining bootstrap helpers")

# Bootstrap statistic for a single (predictor, response) pair within a frame mask.
boot_cor_stat <- function(data, indices, x_col, y_col, method) {
  d <- data[indices, ]
  ok <- is.finite(d[[x_col]]) & is.finite(d[[y_col]])
  if (sum(ok) < 3L) return(NA_real_)
  v <- tryCatch(
    cor(d[[x_col]][ok], d[[y_col]][ok], method = method),
    error = function(e) NA_real_,
    warning = function(w) suppressWarnings(
      cor(d[[x_col]][ok], d[[y_col]][ok], method = method))
  )
  v
}

# Run bootstrap and extract CIs. Returns a one-row data.table.
boot_one_pair <- function(sub, frame_name, predictor, response, n_boot) {
  ok_orig <- is.finite(sub[[predictor]]) & is.finite(sub[[response]])
  n_orig <- sum(ok_orig)
  if (n_orig < 4L) {
    return(data.table(
      frame = frame_name, predictor = predictor, response = response,
      n = n_orig,
      r_pearson = NA_real_,    r_spearman = NA_real_,
      ci_perc_lo_pe = NA_real_, ci_perc_hi_pe = NA_real_,
      ci_bca_lo_pe  = NA_real_, ci_bca_hi_pe  = NA_real_,
      ci_perc_lo_sp = NA_real_, ci_perc_hi_sp = NA_real_,
      ci_bca_lo_sp  = NA_real_, ci_bca_hi_sp  = NA_real_,
      boot_mean_pe  = NA_real_, boot_se_pe    = NA_real_,
      boot_mean_sp  = NA_real_, boot_se_sp    = NA_real_,
      bca_converged = FALSE
    ))
  }

  sub_complete <- sub[ok_orig]

  # Pearson
  bp <- boot(data = as.data.frame(sub_complete),
             statistic = boot_cor_stat, R = n_boot,
             x_col = predictor, y_col = response, method = "pearson")
  ci_perc_pe <- tryCatch(boot.ci(bp, type = "perc")$percent[4:5],
                          error = function(e) c(NA_real_, NA_real_))
  ci_bca_pe  <- tryCatch(boot.ci(bp, type = "bca")$bca[4:5],
                          error = function(e) c(NA_real_, NA_real_))

  # Spearman
  bs <- boot(data = as.data.frame(sub_complete),
             statistic = boot_cor_stat, R = n_boot,
             x_col = predictor, y_col = response, method = "spearman")
  ci_perc_sp <- tryCatch(boot.ci(bs, type = "perc")$percent[4:5],
                          error = function(e) c(NA_real_, NA_real_))
  ci_bca_sp  <- tryCatch(boot.ci(bs, type = "bca")$bca[4:5],
                          error = function(e) c(NA_real_, NA_real_))

  data.table(
    frame = frame_name, predictor = predictor, response = response,
    n = n_orig,
    r_pearson  = bp$t0,
    r_spearman = bs$t0,
    ci_perc_lo_pe = ci_perc_pe[1], ci_perc_hi_pe = ci_perc_pe[2],
    ci_bca_lo_pe  = ci_bca_pe[1],  ci_bca_hi_pe  = ci_bca_pe[2],
    ci_perc_lo_sp = ci_perc_sp[1], ci_perc_hi_sp = ci_perc_sp[2],
    ci_bca_lo_sp  = ci_bca_sp[1],  ci_bca_hi_sp  = ci_bca_sp[2],
    boot_mean_pe  = mean(bp$t, na.rm = TRUE),
    boot_se_pe    = sd(bp$t,   na.rm = TRUE),
    boot_mean_sp  = mean(bs$t, na.rm = TRUE),
    boot_se_sp    = sd(bs$t,   na.rm = TRUE),
    bca_converged = !any(is.na(c(ci_bca_pe, ci_bca_sp)))
  )
}

# Build masks (must match the step 04 patch definitions)
make_mask <- function(dt, frame_name) {
  switch(frame_name,
    "19_sites"               = rep(TRUE, nrow(dt)),
    "16_non_flagged"         = dt$flag_status != "FLAGGED",
    "15_err_alt"             = dt$flag_status != "FLAGGED" &
                               dt$flag_status != "DTM_excluded",
    "forest_only_no_flagged" = dt$flag_status != "FLAGGED" &
                               dt$dominant_lc %in% c("BDF", "ENF", "DNF",
                                                      "EBF", "MFT", "IWL"),
    stop("Unknown frame: ", frame_name)
  )
}

# Annotation: is this combination mechanistically circular or redundant?
flag_combo <- function(predictor, response) {
  if (predictor == "tile_offset_mean") {
    return("circular_predictor_is_offset_intercept")
  }
  if (response == "phase3_re_full") {
    return("redundant_identical_to_ch1_chm_re")
  }
  ""
}

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Run bootstrap on every combination from the original CSV
# ─────────────────────────────────────────────────────────────────────────────

log_subsection(sprintf("Running bootstrap (R = %d) on %d combinations",
                        N_BOOT, nrow(orig_corr)))

boot_rows <- vector("list", nrow(orig_corr))
t0 <- Sys.time()
for (i in seq_len(nrow(orig_corr))) {
  row <- orig_corr[i]
  mask <- make_mask(cross_site, row$frame)
  sub  <- cross_site[mask]
  res  <- boot_one_pair(sub, row$frame, row$predictor, row$response, N_BOOT)
  res[, note := flag_combo(row$predictor, row$response)]
  boot_rows[[i]] <- res
  if (i %% 12 == 0L || i == nrow(orig_corr)) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    log_progress(sprintf("    %d/%d combinations done (%.1fs)",
                          i, nrow(orig_corr), elapsed))
  }
}
boot_dt <- rbindlist(boot_rows)

# Merge in the original parametric p-values for comparison
boot_dt <- merge(boot_dt,
                  orig_corr[, .(frame, predictor, response,
                                p_pearson, p_spearman)],
                  by = c("frame", "predictor", "response"),
                  all.x = TRUE,
                  sort = FALSE)

# Bootstrap-based "approximate p" via percentile method:
# fraction of bootstrap distribution beyond zero (two-sided)
# (Not strictly equivalent to test p-value, but useful for quick read.)

fwrite(boot_dt, file.path(MANUSCRIPT_TB,
                            "groundwork_task6_phase2_correlations_with_ci.csv"))
log_progress(sprintf("  Wrote correlations_with_ci.csv (%d rows)", nrow(boot_dt)))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Headline summary table
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Headline summary")

headline_combos <- data.table(
  frame     = c("19_sites",      "16_non_flagged", "16_non_flagged"),
  predictor = c("b_cover",       "b_cover",        "b_cover"),
  response  = c("ch1_chm_re",    "track_b_re",     "phase3_re_alt18")
)

headline_dt <- merge(headline_combos, boot_dt,
                     by = c("frame", "predictor", "response"),
                     sort = FALSE)

# Pretty-print
fmt_ci <- function(lo, hi) {
  if (is.na(lo) || is.na(hi)) return("NA")
  sprintf("[%+.2f, %+.2f]", lo, hi)
}

cat("\nHEADLINE: bootstrap CIs on b_cover correlations\n")
cat(strrep("-", 90), "\n", sep = "")
cat(sprintf("%-30s %-3s %-7s %-22s %-22s %s\n",
            "frame x response", "n", "r_pe",
            "95% perc CI (Pe)", "95% BCa CI (Pe)", "param p"))
cat(strrep("-", 90), "\n", sep = "")
for (i in seq_len(nrow(headline_dt))) {
  r <- headline_dt[i]
  cat(sprintf("%-30s %-3d %+.3f  %-22s %-22s %.4f\n",
              sprintf("%s x %s", r$frame, r$response),
              r$n, r$r_pearson,
              fmt_ci(r$ci_perc_lo_pe, r$ci_perc_hi_pe),
              fmt_ci(r$ci_bca_lo_pe,  r$ci_bca_hi_pe),
              r$p_pearson))
}
cat(strrep("-", 90), "\n", sep = "")

cat("\nSpearman versions:\n")
cat(strrep("-", 90), "\n", sep = "")
cat(sprintf("%-30s %-3s %-7s %-22s %-22s %s\n",
            "frame x response", "n", "r_sp",
            "95% perc CI (Sp)", "95% BCa CI (Sp)", "param p"))
cat(strrep("-", 90), "\n", sep = "")
for (i in seq_len(nrow(headline_dt))) {
  r <- headline_dt[i]
  cat(sprintf("%-30s %-3d %+.3f  %-22s %-22s %.4f\n",
              sprintf("%s x %s", r$frame, r$response),
              r$n, r$r_spearman,
              fmt_ci(r$ci_perc_lo_sp, r$ci_perc_hi_sp),
              fmt_ci(r$ci_bca_lo_sp,  r$ci_bca_hi_sp),
              r$p_spearman))
}
cat(strrep("-", 90), "\n\n", sep = "")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Leave-one-out influence on the three headline correlations
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Leave-one-out influence analysis")

loo_one <- function(sub, predictor, response, frame_name) {
  ok <- is.finite(sub[[predictor]]) & is.finite(sub[[response]])
  s  <- sub[ok]
  n  <- nrow(s)
  if (n < 5L) return(NULL)

  r_full <- cor(s[[predictor]], s[[response]], method = "pearson")

  rows <- vector("list", n)
  for (i in seq_len(n)) {
    s_loo <- s[-i]
    r_loo <- cor(s_loo[[predictor]], s_loo[[response]], method = "pearson")
    rows[[i]] <- data.table(
      frame = frame_name,
      predictor = predictor,
      response = response,
      removed_ms          = s$manuscript_site[i],
      removed_site_short  = s$site_short[i],
      removed_flag        = s$flag_status[i],
      removed_lc          = s$dominant_lc[i],
      removed_b_cover     = s$b_cover[i],
      removed_response_val = s[[response]][i],
      r_full = r_full,
      r_loo  = r_loo,
      delta_r = r_loo - r_full
    )
  }
  rbindlist(rows)
}

loo_list <- list()
for (i in seq_len(nrow(headline_combos))) {
  hc <- headline_combos[i]
  mask <- make_mask(cross_site, hc$frame)
  sub  <- cross_site[mask]
  loo_list[[i]] <- loo_one(sub, hc$predictor, hc$response, hc$frame)
}
loo_dt <- rbindlist(loo_list)
fwrite(loo_dt, file.path(MANUSCRIPT_TB,
                           "groundwork_task6_phase2_loo_influence.csv"))
log_progress(sprintf("  Wrote loo_influence.csv (%d rows)", nrow(loo_dt)))

# Sites with the largest absolute delta_r in each frame
cat("\nLeave-one-out top-3 most-influential sites per headline correlation:\n")
cat(strrep("-", 90), "\n", sep = "")
for (i in seq_len(nrow(headline_combos))) {
  hc <- headline_combos[i]
  ldt <- loo_dt[frame == hc$frame & predictor == hc$predictor &
                 response == hc$response]
  if (nrow(ldt) == 0L) next
  ldt[, abs_delta := abs(delta_r)]
  setorder(ldt, -abs_delta)
  cat(sprintf("\n%s x %s (n=%d, r_full=%+.3f):\n",
              hc$frame, hc$response, nrow(ldt), ldt$r_full[1]))
  cat(sprintf("  %-3s %-18s %-12s %-7s %-7s %-7s %s\n",
              "ms", "site", "flag", "b_cov", "resp", "r_loo", "delta_r"))
  for (j in seq_len(min(3, nrow(ldt)))) {
    rj <- ldt[j]
    cat(sprintf("  %-3d %-18s %-12s %+.2f  %+.2f  %+.3f %+.3f\n",
                rj$removed_ms,
                substr(rj$removed_site_short, 1, 18),
                rj$removed_flag,
                rj$removed_b_cover,
                rj$removed_response_val,
                rj$r_loo, rj$delta_r))
  }
}
cat(strrep("-", 90), "\n\n", sep = "")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Forest plot of all bootstrap CIs (b_cover predictor focus)
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Generating forest plot")

# Focus the forest plot on b_cover correlations (the mechanistic predictor),
# excluding the redundant phase3_re_full and circular tile_offset_mean.

forest_dt <- boot_dt[predictor == "b_cover" & response != "phase3_re_full"]
forest_dt[, response_label := dplyr::recode(response,
  "ch1_chm_re"      = "Ch1 CHM site RE",
  "track_b_re"      = "16-site RE",
  "phase3_re_alt18" = "Decomposition RE_alt18 (err_alt)")]
forest_dt[, frame_label := factor(frame,
  levels = c("19_sites", "16_non_flagged", "15_err_alt",
             "forest_only_no_flagged"),
  labels = c("19 sites", "16 non-flagged", "15 err_alt",
             "forest-only non-flagged"))]
forest_dt[, response_label := factor(response_label,
  levels = c("Ch1 CHM site RE", "16-site RE", "Decomposition RE_alt18 (err_alt)"))]

forest_dt[, label_text := sprintf("r = %+.2f, n = %d", r_pearson, n)]

p_forest <- ggplot(forest_dt,
                    aes(x = r_pearson, y = response_label, color = frame_label)) +
  geom_vline(xintercept = 0, color = "grey40", linetype = "dashed") +
  geom_errorbarh(aes(xmin = ci_perc_lo_pe, xmax = ci_perc_hi_pe),
                 height = 0.18, position = position_dodge(width = 0.7),
                 linewidth = 0.8) +
  geom_point(size = 2.5, position = position_dodge(width = 0.7)) +
  geom_text(aes(label = label_text, x = ci_perc_hi_pe + 0.05),
            position = position_dodge(width = 0.7),
            hjust = 0, size = 2.8, show.legend = FALSE) +
  scale_color_manual(values = c("19 sites" = "#1b6cb0",
                                 "16 non-flagged" = "#d6604d",
                                 "15 err_alt" = "#762a83",
                                 "forest-only non-flagged" = "#2c7c2c"),
                     name = "Frame") +
  scale_x_continuous(limits = c(-1.05, 1.05),
                     breaks = seq(-1, 1, 0.25),
                     expand = expansion(mult = c(0.02, 0.18))) +
  labs(x = "Pearson r (b_cover vs site RE) with 95% percentile bootstrap CI",
       y = NULL,
       title = "Step 02 — Bootstrap CIs on cross-site b_cover correlations",
       subtitle = sprintf(
         "%d-iter site-level bootstrap. Predicted sign: NEGATIVE (steeper b_cover -> stronger ALS underestimation -> positive site RE).",
         N_BOOT)) +
  theme_cowplot(11) +
  theme(legend.position = "top",
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3))

ggsave(file.path(PLOT_DIR, "task6_phase1_5_bootstrap_forest.pdf"),
       p_forest, width = 12, height = 6, bg = "white")
log_progress("  Wrote task6_phase1_5_bootstrap_forest.pdf")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: LOO diagnostic 3-panel figure
# ─────────────────────────────────────────────────────────────────────────────

log_subsection("Generating LOO diagnostic figure")

mk_loo_panel <- function(loo_sub, title_str) {
  if (nrow(loo_sub) == 0L) {
    return(ggplot() + theme_void() +
             labs(title = sprintf("%s (no data)", title_str)))
  }
  loo_sub <- copy(loo_sub)
  loo_sub[, lab := sprintf("%d %s", removed_ms,
                            substr(removed_site_short, 1, 12))]
  setorder(loo_sub, delta_r)
  loo_sub[, lab := factor(lab, levels = lab)]

  r_full <- loo_sub$r_full[1]

  ggplot(loo_sub, aes(x = lab, y = r_loo, fill = removed_flag)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = r_full, color = "black", linetype = "dashed",
               linewidth = 0.5) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    annotate("text", x = nrow(loo_sub), y = r_full,
             label = sprintf("r_full = %+.3f", r_full),
             hjust = 1, vjust = -0.6, size = 3, color = "black") +
    scale_fill_manual(values = c("FLAGGED" = "#d6604d", "ok" = "#4393c3",
                                  "DTM_excluded" = "#bf812d"),
                       name = "Flag status") +
    coord_flip() +
    labs(x = NULL, y = "Pearson r after removing this site",
         title = title_str) +
    theme_cowplot(10) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10),
          axis.text.y = element_text(size = 8))
}

loo_panel_1 <- mk_loo_panel(
  loo_dt[frame == "19_sites" & response == "ch1_chm_re"],
  "19 sites x Ch1 CHM site RE")
loo_panel_2 <- mk_loo_panel(
  loo_dt[frame == "16_non_flagged" & response == "track_b_re"],
  "16 non-flagged x 16-site RE")
loo_panel_3 <- mk_loo_panel(
  loo_dt[frame == "16_non_flagged" & response == "phase3_re_alt18"],
  "16 non-flagged x Decomposition RE_alt18")

# Shared legend
legend_dummy <- ggplot(loo_dt, aes(x = 1, y = 1, fill = removed_flag)) +
  geom_col() +
  scale_fill_manual(values = c("FLAGGED" = "#d6604d", "ok" = "#4393c3",
                                "DTM_excluded" = "#bf812d"),
                     name = "Flag status of removed site") +
  theme_cowplot(10) + theme(legend.position = "top")
leg <- get_legend(legend_dummy)

loo_combined <- plot_grid(loo_panel_1, loo_panel_2, loo_panel_3,
                           ncol = 3, align = "hv")
loo_full <- plot_grid(leg, loo_combined, ncol = 1, rel_heights = c(0.06, 1))

ggsave(file.path(PLOT_DIR, "task6_phase1_5_loo_influence.pdf"),
       loo_full, width = 16, height = 6, bg = "white")
log_progress("  Wrote task6_phase1_5_loo_influence.pdf")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Final summary
# ─────────────────────────────────────────────────────────────────────────────

log_section("step 02 — Complete")

log_progress("Output files:")
log_progress(sprintf("  %s/groundwork_task6_phase2_correlations_with_ci.csv", MANUSCRIPT_TB))
log_progress(sprintf("  %s/groundwork_task6_phase2_loo_influence.csv",        MANUSCRIPT_TB))
log_progress(sprintf("  %s/task6_phase1_5_bootstrap_forest.pdf",              PLOT_DIR))
log_progress(sprintf("  %s/task6_phase1_5_loo_influence.pdf",                 PLOT_DIR))
log_progress("")
log_progress("Decision rule update with bootstrap CIs:")
log_progress("  * If the headline 16_non_flagged x track_b_re 95% CI EXCLUDES 0:")
log_progress("    -> b_cover variation explains residual ~12% V_site;")
log_progress("       variance puzzle closed (modulo sample size).")
log_progress("  * If the CI INCLUDES 0 but trends correctly:")
log_progress("    -> mechanism is partial; report as 'cover-slope variation")
log_progress("       contributes to but does not fully explain the residual'.")
log_progress("  * The Decomposition RE_alt18 result is a corroborating secondary")
log_progress("    test under a different post-hoc filtering regime.")
log_progress("")
log_progress("Outputs:")
log_progress("  - The two CSVs above")
log_progress("  - The two PDFs above")
log_progress("  - This script's stdout/stderr log")
