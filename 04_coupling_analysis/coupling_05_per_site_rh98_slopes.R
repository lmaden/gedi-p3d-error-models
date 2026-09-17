# =====================================================================
# coupling_05_per_site_rh98_slopes.R
#
# Step 05: per-site CNN-DTM RH98 slope variation, and its cross-site
# correlation with per-site compensation strength corr_errP3D_errDTM
# (computed in step 02).
#
# Mechanism context. step 03 established that, at the site-mean level,
# canopy-rich sites have stronger compensation (|r_15| ~ 0.37 - 0.45 for
# cover, RH98, leaf-on, WSCI; all carrying the same negative sign).
# step 03's strongest single predictor was meta_abs_geoacc_avg
# (r_15 = -0.699, p = 0.004), interpreted as a shared-noise channel.
# Step 05 refines the canopy-side test from "more canopy" to "stronger
# DTM response to canopy": for each site, fit b_RH98 (slope of CNN-DTM
# error on rh_98) and ask whether sites with steeper b_RH98 have more
# negative corr_errP3D_errDTM. The mechanistic prediction is a negative
# cross-site correlation between b_RH98 and corr_errP3D_errDTM.
#
# Inputs (all on disk, no cluster-side raster extraction):
#   checkpoints/groundwork_task4_phase2.rds                  (step 02)
#     -> $fp_data   (608,891 forest footprints across 18 sites,
#                    with err_p3d, err_alt, err_dtm precomputed)
#   data/enriched_by_site/site_NN_enriched.csv[.gz]          (18 sites)
#     -> read just (shot_number, rh_98, cover, slope_mean) per site,
#        join onto fp_data by (tracker_site, shot_number)
#   manuscript_tables/coupling_site_summary.csv       (step 02)
#     -> per-site corr_errP3D_errDTM
#
# Outputs:
#   manuscript_tables/coupling_rh98_per_site_slopes.csv
#     - per-site b_RH98 (univariate and partial), SE, R^2, n,
#       joined corr_errP3D_errDTM, flag_status, dominant_lc
#   manuscript_tables/coupling_rh98_correlations.csv
#     - cross-site Pearson + Spearman of b_RH98 vs corr, on 18-site
#       and 15-non-flagged frames, for both univariate and partial
#       slopes
#   manuscript_tables/coupling_rh98_bootstrap.csv
#     - 10,000-iter cluster bootstrap CIs on all four cross-site r
#       values
#   manuscript_tables/coupling_rh98_loo.csv
#     - leave-one-site-out r values for the headline 15-site test
#   plots/groundwork/groundwork_task4b_per_site_fits.pdf
#     - 18-panel scatter of err_dtm vs rh_98 per site, OLS line,
#       slope and R^2 annotated
#   plots/groundwork/groundwork_task4b_cross_site_scatter.pdf
#     - cross-site scatter b_RH98 vs corr_errP3D_errDTM (univariate
#       and partial; 18 + 15 frames), with bootstrap CI annotation
#   checkpoints/groundwork_task4b.rds
#
# Runtime estimate: ~1-2 minutes. Per-site OLS on ~600k footprints +
# 10,000 cluster bootstrap iterations + figures.
# =====================================================================

# ---- 1. Environment -------------------------------------------------

if (!exists("PROJECT_ROOT") || !exists("log_progress")) {
  cfg_candidates <- c(
    file.path(getwd(), "analysis_config.R"),
    "/gpfs/data1/vclgp/lmaden/chpt1/scripts/reviewed/analysis_config.R"
  )
  util_candidates <- c(
    file.path(getwd(), "analysis_utils.R"),
    "/gpfs/data1/vclgp/lmaden/chpt1/scripts/reviewed/analysis_utils.R"
  )
  cfg  <- cfg_candidates[file.exists(cfg_candidates)][1]
  util <- util_candidates[file.exists(util_candidates)][1]
  if (is.na(cfg) || is.na(util)) {
    stop("Could not locate analysis_config.R and analysis_utils.R. ",
         "Source them manually from scripts/reviewed/ before running this.")
  }
  source(cfg)
  source(util)
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

MT_DIR     <- file.path(PROJECT_ROOT, "manuscript_tables")
ENRICH_DIR <- file.path(PROJECT_ROOT, "data", "enriched_by_site")
PLOT_DIR   <- file.path(PROJECT_ROOT, "plots", "groundwork")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

log_progress(strrep("=", 70))
log_progress("Step 05: per-site CNN-DTM RH98 slopes, cross-site test")
log_progress(strrep("=", 70))

# ---- 2. Load step 02 footprint data ---------------------------------

stopifnot(checkpoint_exists("groundwork_task4_phase2"))
phase2 <- load_checkpoint("groundwork_task4_phase2")

fp_data      <- as.data.table(phase2$fp_data)
site_summary <- as.data.table(phase2$site_summary)
lookup_ms    <- as.data.table(phase2$lookup_ms)

stopifnot("err_dtm"   %in% names(fp_data))
stopifnot("shot_number" %in% names(fp_data))
stopifnot("tracker_site" %in% names(fp_data))
stopifnot("manuscript_site" %in% names(fp_data))
stopifnot("flag_status" %in% names(fp_data))
stopifnot("corr_errP3D_errDTM" %in% names(site_summary))

log_progress(sprintf(
  "step 02 fp_data: %s rows across %d sites (range manuscript %s)",
  format(nrow(fp_data), big.mark = ","),
  uniqueN(fp_data$tracker_site),
  paste(range(fp_data$manuscript_site), collapse = "-")
))

# Make sure shot_number is character on both sides for joining.
fp_data[, shot_number := as.character(shot_number)]

# ---- 3. Attach rh_98, cover, slope_mean from enriched CSVs ----------

trackers <- sort(unique(fp_data$tracker_site))

read_site_covs <- function(tid) {
  patterns <- c(
    file.path(ENRICH_DIR, sprintf("site_%02d_enriched.csv.gz", tid)),
    file.path(ENRICH_DIR, sprintf("site_%d_enriched.csv.gz",   tid)),
    file.path(ENRICH_DIR, sprintf("site_%02d_enriched.csv",    tid)),
    file.path(ENRICH_DIR, sprintf("site_%d_enriched.csv",      tid))
  )
  path <- patterns[file.exists(patterns)][1]
  if (is.na(path)) {
    log_progress(sprintf("  WARN: no enriched CSV for tracker %d", tid))
    return(NULL)
  }
  hdr   <- names(fread(path, nrows = 0L))
  want  <- c("shot_number", "rh_98", "cover", "slope_mean")
  have  <- intersect(want, hdr)
  if (!("shot_number" %in% have) || !("rh_98" %in% have)) {
    log_progress(sprintf("  WARN: tracker %d missing shot_number or rh_98",
                         tid))
    return(NULL)
  }
  DT <- fread(path, select = have,
              colClasses = list(character = "shot_number"),
              showProgress = FALSE)
  DT[, tracker_site := tid]
  DT[]
}

log_progress("Reading rh_98, cover, slope_mean per site...")
covs_list <- lapply(trackers, read_site_covs)
covs_all  <- rbindlist(covs_list, fill = TRUE, use.names = TRUE)
log_progress(sprintf("  Read %s footprint-level rows of covariates",
                     format(nrow(covs_all), big.mark = ",")))

# Merge onto fp_data. Keep only fp_data rows that match (inner join
# semantics) so the per-site fits use exactly the step 02 forest
# subset.
n_before <- nrow(fp_data)
fp_data  <- merge(
  fp_data, covs_all,
  by = c("tracker_site", "shot_number"),
  all.x = FALSE, all.y = FALSE,
  sort = FALSE
)
log_progress(sprintf(
  "  Joined: %s -> %s rows (drop = %.2f%%) after rh_98 attach",
  format(n_before, big.mark = ","),
  format(nrow(fp_data), big.mark = ","),
  100 * (1 - nrow(fp_data) / n_before)
))

# Drop rows with non-finite rh_98 (downstream OLS would drop them
# anyway; surface the count here so the per-site n is honest).
n_pre_rh <- nrow(fp_data)
fp_data  <- fp_data[is.finite(rh_98)]
log_progress(sprintf("  Dropped %s rows with non-finite rh_98",
                     format(n_pre_rh - nrow(fp_data), big.mark = ",")))

# ---- 4. Per-site OLS: err_dtm ~ rh_98 (univariate) ------------------

per_site_fits <- function(d_site, ms_id, tid, short, flag, dlc) {
  e <- d_site$err_dtm
  x <- d_site$rh_98
  ok <- is.finite(e) & is.finite(x)
  n  <- sum(ok)
  if (n < 50L) {
    return(data.table(
      manuscript_site = ms_id, tracker_site = tid,
      site_short_name = short, flag_status = flag, dominant_lc = dlc,
      n = n,
      b_rh98_uni = NA_real_, se_rh98_uni = NA_real_, r2_uni = NA_real_,
      b_rh98_partial = NA_real_, se_rh98_partial = NA_real_,
      r2_partial = NA_real_,
      n_partial = NA_integer_
    ))
  }
  # Univariate fit
  fit_u <- lm(err_dtm ~ rh_98, data = d_site, subset = ok)
  cu    <- summary(fit_u)$coefficients
  b_u   <- cu["rh_98", "Estimate"]
  se_u  <- cu["rh_98", "Std. Error"]
  r2_u  <- summary(fit_u)$r.squared

  # Partial fit: rh_98 controlling for cover and slope_mean. Drop the
  # site if cover or slope_mean has no variance or is missing for too
  # many rows.
  has_cover <- "cover" %in% names(d_site) && any(is.finite(d_site$cover))
  has_slope <- "slope_mean" %in% names(d_site) &&
                any(is.finite(d_site$slope_mean))
  b_p   <- NA_real_
  se_p  <- NA_real_
  r2_p  <- NA_real_
  n_p   <- NA_integer_
  if (has_cover && has_slope) {
    ok_p <- ok & is.finite(d_site$cover) & is.finite(d_site$slope_mean)
    n_p  <- sum(ok_p)
    if (n_p >= 50L &&
        sd(d_site$cover[ok_p])      > 0 &&
        sd(d_site$slope_mean[ok_p]) > 0) {
      fit_p <- lm(err_dtm ~ rh_98 + cover + slope_mean,
                  data = d_site, subset = ok_p)
      cp    <- summary(fit_p)$coefficients
      if ("rh_98" %in% rownames(cp)) {
        b_p  <- cp["rh_98", "Estimate"]
        se_p <- cp["rh_98", "Std. Error"]
      }
      r2_p <- summary(fit_p)$r.squared
    }
  }

  data.table(
    manuscript_site = ms_id, tracker_site = tid,
    site_short_name = short, flag_status = flag, dominant_lc = dlc,
    n = n,
    b_rh98_uni = b_u, se_rh98_uni = se_u, r2_uni = r2_u,
    b_rh98_partial = b_p, se_rh98_partial = se_p,
    r2_partial = r2_p,
    n_partial = n_p
  )
}

log_progress("Fitting per-site OLS (univariate and partial)...")

site_keys <- unique(fp_data[, .(manuscript_site, tracker_site,
                                site_short_name, flag_status, dominant_lc)])
site_keys <- site_keys[order(manuscript_site)]

per_site <- rbindlist(lapply(seq_len(nrow(site_keys)), function(i) {
  k <- site_keys[i]
  ds <- fp_data[tracker_site == k$tracker_site]
  per_site_fits(ds, k$manuscript_site, k$tracker_site,
                k$site_short_name, k$flag_status, k$dominant_lc)
}))

# Join step 02's compensation strength.
per_site <- merge(
  per_site,
  site_summary[, .(tracker_site, corr_errP3D_errDTM,
                   var_reduction_frac, p3d_rmse, alt_rmse,
                   dtm_rmse)],
  by = "tracker_site", all.x = TRUE, sort = FALSE
)

setorder(per_site, manuscript_site)

out_per_site <- file.path(MT_DIR, "coupling_rh98_per_site_slopes.csv")
fwrite(per_site, out_per_site)
log_progress(sprintf("Wrote %s", out_per_site))

# Console echo: per-site b_RH98 (univariate) sorted by magnitude.
log_progress("Per-site b_RH98 (univariate), sorted by magnitude:")
echo <- per_site[order(-abs(b_rh98_uni)),
                 .(manuscript_site, site_short_name, flag_status,
                   n, b_rh98_uni, se_rh98_uni, r2_uni,
                   corr_errP3D_errDTM)]
print(echo, digits = 3, row.names = FALSE)

# ---- 5. Cross-site Pearson and Spearman -----------------------------
#
# Two slope columns x two frames = four headline correlations.

corr_cell <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4L) {
    return(list(n = sum(ok), pearson_r = NA_real_, pearson_p = NA_real_,
                spearman_rho = NA_real_, spearman_p = NA_real_))
  }
  pp <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
  ss <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman",
                                  exact = FALSE))
  list(n = sum(ok),
       pearson_r = unname(pp$estimate), pearson_p = pp$p.value,
       spearman_rho = unname(ss$estimate), spearman_p = ss$p.value)
}

frames <- list(
  "18_sites_dtm_ok"      = per_site,
  "15_non_flagged_dtm_ok" = per_site[flag_status != "FLAGGED"]
)

slopes <- list(
  "b_rh98_uni"     = "b_rh98_uni",
  "b_rh98_partial" = "b_rh98_partial"
)

cor_rows <- list()
for (fr_name in names(frames)) {
  fr <- frames[[fr_name]]
  for (sl_name in names(slopes)) {
    sl_col <- slopes[[sl_name]]
    cc <- corr_cell(fr[[sl_col]], fr$corr_errP3D_errDTM)
    cor_rows[[length(cor_rows) + 1L]] <- data.table(
      frame = fr_name,
      slope = sl_name,
      n = cc$n,
      pearson_r = cc$pearson_r, pearson_p = cc$pearson_p,
      spearman_rho = cc$spearman_rho, spearman_p = cc$spearman_p
    )
  }
}
corr_tbl <- rbindlist(cor_rows)

out_corr <- file.path(MT_DIR, "coupling_rh98_correlations.csv")
fwrite(corr_tbl, out_corr)
log_progress(sprintf("Wrote %s", out_corr))

log_progress("Cross-site correlations:")
print(corr_tbl, digits = 3, row.names = FALSE)

# ---- 6. Cluster bootstrap CIs (10,000 iterations) -------------------
#
# Resample sites with replacement at the site level, recompute the
# Pearson r between b_RH98 and corr_errP3D_errDTM. Report 95%
# percentile CI on each of the four cells.

set.seed(2025)
n_boot <- 10000L

bootstrap_r <- function(x, y, n_boot) {
  ok <- is.finite(x) & is.finite(y)
  x  <- x[ok]; y <- y[ok]
  n  <- length(x)
  if (n < 4L) return(rep(NA_real_, n_boot))
  rs <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    xb  <- x[idx]; yb <- y[idx]
    if (sd(xb) == 0 || sd(yb) == 0) {
      rs[b] <- NA_real_
    } else {
      rs[b] <- cor(xb, yb)
    }
  }
  rs
}

log_progress(sprintf("Cluster bootstrap, %s iterations...",
                     format(n_boot, big.mark = ",")))

boot_rows <- list()
for (fr_name in names(frames)) {
  fr <- frames[[fr_name]]
  for (sl_name in names(slopes)) {
    sl_col <- slopes[[sl_name]]
    rs <- bootstrap_r(fr[[sl_col]], fr$corr_errP3D_errDTM, n_boot)
    rs_ok <- rs[is.finite(rs)]
    boot_rows[[length(boot_rows) + 1L]] <- data.table(
      frame = fr_name,
      slope = sl_name,
      n_sites = sum(is.finite(fr[[sl_col]]) & is.finite(fr$corr_errP3D_errDTM)),
      n_iter_ok = length(rs_ok),
      r_mean = mean(rs_ok),
      r_median = median(rs_ok),
      ci_lo_95 = quantile(rs_ok, 0.025, names = FALSE),
      ci_hi_95 = quantile(rs_ok, 0.975, names = FALSE),
      ci_excludes_zero = (quantile(rs_ok, 0.025, names = FALSE) > 0) |
                        (quantile(rs_ok, 0.975, names = FALSE) < 0)
    )
  }
}
boot_tbl <- rbindlist(boot_rows)

out_boot <- file.path(MT_DIR, "coupling_rh98_bootstrap.csv")
fwrite(boot_tbl, out_boot)
log_progress(sprintf("Wrote %s", out_boot))

log_progress("Bootstrap 95% CIs:")
print(boot_tbl, digits = 3, row.names = FALSE)

# ---- 7. Leave-one-site-out leverage ---------------------------------
#
# Headline test = 15 non-flagged x univariate b_RH98. Recompute Pearson
# r dropping each site in turn; flag sites whose removal moves r by
# more than 0.10.

headline_frame <- per_site[flag_status != "FLAGGED"]
headline_x     <- headline_frame$b_rh98_uni
headline_y     <- headline_frame$corr_errP3D_errDTM
ok_hl          <- is.finite(headline_x) & is.finite(headline_y)
r_full         <- cor(headline_x[ok_hl], headline_y[ok_hl])

loo_rows <- list()
for (i in which(ok_hl)) {
  keep_idx <- which(ok_hl)
  keep_idx <- keep_idx[keep_idx != i]
  r_loo <- cor(headline_x[keep_idx], headline_y[keep_idx])
  loo_rows[[length(loo_rows) + 1L]] <- data.table(
    manuscript_site = headline_frame$manuscript_site[i],
    site_short_name = headline_frame$site_short_name[i],
    b_rh98_uni      = headline_x[i],
    corr_errP3D_errDTM = headline_y[i],
    r_full          = r_full,
    r_loo           = r_loo,
    delta_r         = r_loo - r_full
  )
}
loo_tbl <- rbindlist(loo_rows)
loo_tbl <- loo_tbl[order(-abs(delta_r))]

out_loo <- file.path(MT_DIR, "coupling_rh98_loo.csv")
fwrite(loo_tbl, out_loo)
log_progress(sprintf("Wrote %s", out_loo))

log_progress("LOO leverage (headline 15-site x univariate b_RH98):")
print(loo_tbl, digits = 3, row.names = FALSE)

# ---- 8. Figure 1: per-site OLS panels -------------------------------

log_progress("Figure 1: per-site err_dtm vs rh_98 panels...")

panel_df <- copy(fp_data)
panel_df[, site_label := sprintf("Site %02d %s%s",
                                 manuscript_site, site_short_name,
                                 ifelse(flag_status == "FLAGGED",
                                        " *", ""))]
panel_df[, site_label := factor(
  site_label,
  levels = unique(panel_df[order(manuscript_site)]$site_label)
)]

# Subsample to keep file size manageable.
set.seed(2025)
panel_sub <- panel_df[, .SD[sample(.N, min(.N, 12000L))],
                      by = manuscript_site]

annot_df <- per_site[, .(
  manuscript_site, site_short_name, flag_status,
  site_label = sprintf("Site %02d %s%s",
                       manuscript_site, site_short_name,
                       ifelse(flag_status == "FLAGGED", " *", "")),
  txt = sprintf("b = %+.3f\nR^2 = %.2f", b_rh98_uni, r2_uni)
)]
annot_df[, site_label := factor(site_label,
                                levels = levels(panel_df$site_label))]

p_fig1 <- ggplot(panel_sub, aes(x = rh_98, y = err_dtm)) +
  geom_hex(bins = 40) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey40") +
  geom_smooth(method = "lm", se = FALSE,
              colour = "red", linewidth = 0.4, formula = y ~ x) +
  geom_text(data = annot_df, aes(x = -Inf, y = Inf, label = txt),
            hjust = -0.05, vjust = 1.2, size = 2.6,
            inherit.aes = FALSE) +
  facet_wrap(~ site_label, ncol = 5, scales = "free") +
  scale_fill_viridis_c(trans = "log10", name = "count") +
  labs(
    x = "GEDI rh_98 (m)",
    y = "CNN-DTM error (err_dtm = P3D DTM - 3DEP DTM, m)",
    title = "Step 05 Fig 1: per-site OLS of CNN-DTM error on RH98",
    subtitle = "Red line = univariate OLS fit. Slope b = m of DTM error per m of canopy height. Flagged sites marked with *."
  ) +
  theme_cowplot(font_size = 10) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(size = 8, lineheight = 0.9),
    panel.spacing = unit(0.3, "lines")
  )

out_fig1 <- file.path(PLOT_DIR, "groundwork_task4b_per_site_fits.pdf")
ggsave(out_fig1, p_fig1, width = 14, height = 10)
log_progress(sprintf("Wrote %s", out_fig1))

# ---- 9. Figure 2: cross-site scatter --------------------------------

log_progress("Figure 2: cross-site b_RH98 vs corr_errP3D_errDTM...")

# Build a tidy frame with all four cells (slope x frame).
make_frame <- function(d, frame_label, slope_col, slope_label) {
  d2 <- copy(d)
  d2[, slope_value := d2[[slope_col]]]
  d2[, frame_label := frame_label]
  d2[, slope_label := slope_label]
  d2[is.finite(slope_value) & is.finite(corr_errP3D_errDTM)]
}

all_pts <- rbindlist(list(
  make_frame(per_site,                          "18 sites (DTM-OK)",
             "b_rh98_uni",     "Univariate b_RH98"),
  make_frame(per_site,                          "18 sites (DTM-OK)",
             "b_rh98_partial", "Partial b_RH98 (control: cover, slope)"),
  make_frame(per_site[flag_status != "FLAGGED"], "15 non-flagged (DTM-OK)",
             "b_rh98_uni",     "Univariate b_RH98"),
  make_frame(per_site[flag_status != "FLAGGED"], "15 non-flagged (DTM-OK)",
             "b_rh98_partial", "Partial b_RH98 (control: cover, slope)")
), fill = TRUE)

# Per-cell correlation annotation.
ann_cells <- merge(corr_tbl, boot_tbl,
                   by = c("frame", "slope"), all.x = TRUE)
ann_cells[, frame_label := ifelse(frame == "18_sites_dtm_ok",
                                  "18 sites (DTM-OK)",
                                  "15 non-flagged (DTM-OK)")]
ann_cells[, slope_label := ifelse(slope == "b_rh98_uni",
                                  "Univariate b_RH98",
                                  "Partial b_RH98 (control: cover, slope)")]
ann_cells[, txt := sprintf("r = %+.2f (p = %.3f)\n95%% boot CI: [%+.2f, %+.2f]",
                           pearson_r, pearson_p, ci_lo_95, ci_hi_95)]

p_fig2 <- ggplot(all_pts, aes(slope_value, corr_errP3D_errDTM)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  geom_smooth(method = "lm", se = TRUE, colour = "grey20",
              linewidth = 0.5, alpha = 0.18, formula = y ~ x) +
  geom_point(aes(colour = flag_status, shape = dominant_lc), size = 2.6) +
  geom_text_repel(aes(label = site_short_name, colour = flag_status),
                  size = 2.5, min.segment.length = 0.2,
                  max.overlaps = 18, seed = 1L, show.legend = FALSE) +
  geom_text(data = ann_cells,
            aes(x = -Inf, y = Inf, label = txt),
            hjust = -0.05, vjust = 1.2, size = 2.8,
            colour = "grey25", inherit.aes = FALSE) +
  scale_colour_manual(values = c("OK" = "#1b6cb0",
                                 "FLAGGED" = "#c0392b",
                                 "DTM_excluded" = "#8e44ad")) +
  facet_grid(slope_label ~ frame_label, scales = "free_x") +
  labs(
    x = "Per-site b_RH98 (m of CNN-DTM error per m of GEDI canopy height)",
    y = "Per-site footprint-level r(err_p3d, err_dtm) from step 02",
    title = "Step 05 Fig 2: per-site CNN-DTM RH98 slope vs compensation strength",
    subtitle = "Predicted negative correlation: steeper b_RH98 (more pulling-into-canopy) -> more negative compensation r.",
    colour = "Flag status", shape = "Dominant LC"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(0.7, "lines"),
    plot.title = element_text(face = "bold")
  )

out_fig2 <- file.path(PLOT_DIR, "groundwork_task4b_cross_site_scatter.pdf")
ggsave(out_fig2, p_fig2, width = 12, height = 10)
log_progress(sprintf("Wrote %s", out_fig2))

# ---- 10. Checkpoint and summary -------------------------------------

save_checkpoint("groundwork_task4b", list(
  per_site   = per_site,
  corr_tbl   = corr_tbl,
  boot_tbl   = boot_tbl,
  loo_tbl    = loo_tbl,
  n_boot     = n_boot,
  fp_n       = nrow(fp_data),
  timestamp  = Sys.time()
))

cat("\n\n====================  SUMMARY  ====================\n\n")

cat("Per-site b_RH98 distribution (univariate):\n")
cat(sprintf("  18 sites: median %+0.3f, mean %+0.3f, range [%+0.3f, %+0.3f]\n",
            median(per_site$b_rh98_uni, na.rm = TRUE),
            mean(per_site$b_rh98_uni, na.rm = TRUE),
            min(per_site$b_rh98_uni, na.rm = TRUE),
            max(per_site$b_rh98_uni, na.rm = TRUE)))
nf <- per_site[flag_status != "FLAGGED"]
cat(sprintf("  15 non-flagged: median %+0.3f, mean %+0.3f, range [%+0.3f, %+0.3f]\n",
            median(nf$b_rh98_uni, na.rm = TRUE),
            mean(nf$b_rh98_uni, na.rm = TRUE),
            min(nf$b_rh98_uni, na.rm = TRUE),
            max(nf$b_rh98_uni, na.rm = TRUE)))

cat("\nCross-site test (b_RH98_uni vs corr_errP3D_errDTM):\n")
hl_corr <- corr_tbl[frame == "15_non_flagged_dtm_ok" & slope == "b_rh98_uni"]
hl_boot <- boot_tbl[frame == "15_non_flagged_dtm_ok" & slope == "b_rh98_uni"]
cat(sprintf("  15 non-flagged: r = %+0.3f (p = %.3f), 95%% boot CI [%+0.3f, %+0.3f]\n",
            hl_corr$pearson_r, hl_corr$pearson_p,
            hl_boot$ci_lo_95, hl_boot$ci_hi_95))
hl_corr_18 <- corr_tbl[frame == "18_sites_dtm_ok" & slope == "b_rh98_uni"]
hl_boot_18 <- boot_tbl[frame == "18_sites_dtm_ok" & slope == "b_rh98_uni"]
cat(sprintf("  18 sites:       r = %+0.3f (p = %.3f), 95%% boot CI [%+0.3f, %+0.3f]\n",
            hl_corr_18$pearson_r, hl_corr_18$pearson_p,
            hl_boot_18$ci_lo_95, hl_boot_18$ci_hi_95))

cat("\nPartial slope (controlling for cover and slope_mean):\n")
hl_corr_p <- corr_tbl[frame == "15_non_flagged_dtm_ok" & slope == "b_rh98_partial"]
hl_boot_p <- boot_tbl[frame == "15_non_flagged_dtm_ok" & slope == "b_rh98_partial"]
cat(sprintf("  15 non-flagged: r = %+0.3f (p = %.3f), 95%% boot CI [%+0.3f, %+0.3f]\n",
            hl_corr_p$pearson_r, hl_corr_p$pearson_p,
            hl_boot_p$ci_lo_95, hl_boot_p$ci_hi_95))

cat("\nLOO leverage (headline 15-site):\n")
top_lev <- loo_tbl[order(-abs(delta_r))][1:3]
for (i in seq_len(nrow(top_lev))) {
  cat(sprintf("  Site %02d %s: drop r from %+0.3f to %+0.3f (delta = %+0.3f)\n",
              top_lev$manuscript_site[i],
              top_lev$site_short_name[i],
              top_lev$r_full[i],
              top_lev$r_loo[i],
              top_lev$delta_r[i]))
}

cat("\nArtifacts:\n")
cat(sprintf("  %s\n", out_per_site))
cat(sprintf("  %s\n", out_corr))
cat(sprintf("  %s\n", out_boot))
cat(sprintf("  %s\n", out_loo))
cat(sprintf("  %s\n", out_fig1))
cat(sprintf("  %s\n", out_fig2))
cat(sprintf("  checkpoint: groundwork_task4b\n"))

cat("\n===========================================================\n\n")
log_progress("Step 05 complete.")
