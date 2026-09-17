# =====================================================================
# coupling_03b_site_correlates_metadata.R
#
# Patch for coupling_03_site_correlates.R.
#
# The first pass used guessed column names for three metadata covariates
# that turned out to be wrong in the enriched CSVs, yielding NA columns:
#   meta_offnad_avg        -> should be  meta_off_nadir_avg
#   meta_sunelev_avg       -> should be  meta_sun_elev_avg
#   meta_fwd_rev_ratio_avg -> doesn't exist; two separate columns instead:
#                             meta_fwd_ratio and meta_rev_ratio
#
# This patch:
#   - recomputes forest-footprint means for the three corrected covariates
#     plus a constructed meta_fwd_minus_rev = fwd_ratio - rev_ratio
#     plus meta_tot_ct (image count per footprint, useful proxy for
#     photogrammetric geometry quality)
#   - joins them into the existing phase25 per-site covariate CSV
#   - recomputes the full correlation table with the new covariates
#   - regenerates both PDFs with the filled-in panels
#   - updates the checkpoint
#
# Runtime: ~1 minute (reads only the 5 needed columns from 18 CSVs).
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
    stop("Could not locate analysis_config.R and analysis_utils.R.")
  }
  source(cfg)
  source(util)
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

MT_DIR     <- file.path(PROJECT_ROOT, "manuscript_tables")
ENRICH_DIR <- file.path(PROJECT_ROOT, "data", "enriched_by_site")
PLOT_DIR   <- file.path(PROJECT_ROOT, "plots", "groundwork")

log_progress(strrep("=", 70))
log_progress("Step 03 PATCH — add missing covariates, recompute")
log_progress(strrep("=", 70))

# ---- 2. Load existing step 03 per-site table ----------------------

prev_csv <- file.path(MT_DIR, "coupling_site_covariates.csv")
stopifnot(file.exists(prev_csv))
joined_prev <- fread(prev_csv)
log_progress(sprintf("Loaded existing step 03 per-site table: %d rows x %d cols",
                     nrow(joined_prev), ncol(joined_prev)))

# Drop the all-NA columns from the previous pass (we'll replace them).
drop_cols <- c("meta_offnad_avg_avg", "meta_sunelev_avg_avg",
               "meta_fwd_rev_ratio_avg_avg")
drop_cols <- intersect(drop_cols, names(joined_prev))
if (length(drop_cols)) {
  joined_prev[, (drop_cols) := NULL]
  log_progress(sprintf("Dropped %d all-NA columns from prior pass: %s",
                       length(drop_cols), paste(drop_cols, collapse = ", ")))
}

# ---- 3. Re-read enriched CSVs for the missing covariates ------------

FOREST_CLASSES <- c("EBF", "BDF", "ENF", "DNF")

# Columns to read this pass. LC column handled with a fallback below.
NEW_COLS <- c("meta_off_nadir_avg", "meta_sun_elev_avg",
              "meta_fwd_ratio", "meta_rev_ratio", "meta_tot_ct")

enriched_files <- list.files(ENRICH_DIR,
  pattern = "^site_[0-9]+_enriched\\.csv(\\.gz)?$",
  full.names = TRUE
)

get_tracker_id <- function(p) {
  m <- regmatches(basename(p), regexec("site_([0-9]+)_enriched", basename(p)))[[1]]
  if (length(m) == 2L) as.integer(m[2]) else NA_integer_
}

compute_new_means <- function(fp) {
  tid <- get_tracker_id(fp)
  hdr <- names(fread(fp, nrows = 0L))

  # Determine the LC column present in this file.
  lc_col <- intersect(c("lc_l1_code", "lc2022_mode_l1_code", "lc2022_l1_code"), hdr)[1]
  if (is.na(lc_col)) stop(sprintf("No LC column in %s", basename(fp)))

  keep <- unique(c(lc_col, intersect(NEW_COLS, hdr)))
  DT <- fread(fp, select = keep, showProgress = FALSE)
  if (lc_col != "lc_l1_code") setnames(DT, lc_col, "lc_l1_code")

  DT_forest <- DT[lc_l1_code %in% FOREST_CLASSES]

  # Means for the raw columns that are present.
  present <- intersect(NEW_COLS, names(DT_forest))
  out <- DT_forest[, lapply(.SD, function(x) mean(x, na.rm = TRUE)),
                   .SDcols = present]
  setnames(out, present, paste0(present, "_avg"))

  # Constructed covariate: forward-minus-reverse image fraction.
  # This is probably the closest raw analogue of the manuscript's
  # composite meta_fwdrev_z. Only makes sense if both are present.
  if (all(c("meta_fwd_ratio", "meta_rev_ratio") %in% present)) {
    fwd_mean <- out$meta_fwd_ratio_avg
    rev_mean <- out$meta_rev_ratio_avg
    out[, meta_fwd_minus_rev_avg := fwd_mean - rev_mean]
  } else {
    out[, meta_fwd_minus_rev_avg := NA_real_]
  }

  out[, tracker_site := tid]
  out[]
}

log_progress("Reading new covariates from 18 enriched CSVs...")
t0 <- Sys.time()
new_cov <- rbindlist(lapply(enriched_files, function(fp) {
  r <- tryCatch(compute_new_means(fp), error = function(e) {
    log_progress(sprintf("  WARN: %s — %s", basename(fp), e$message))
    NULL
  })
  if (!is.null(r))
    log_progress(sprintf("  site_%02d: OK", r$tracker_site))
  r
}), fill = TRUE)
t1 <- Sys.time()
log_progress(sprintf("  Done in %.1f sec (%d sites)",
                     as.numeric(difftime(t1, t0, units = "secs")),
                     nrow(new_cov)))

# ---- 4. Join new covariates into existing table ---------------------

joined <- merge(joined_prev, new_cov,
                by = "tracker_site", all.x = TRUE, sort = FALSE)
stopifnot(nrow(joined) == 18L)
log_progress(sprintf("Joined: %d rows x %d cols", nrow(joined), ncol(joined)))

out_cov_csv <- file.path(MT_DIR, "coupling_site_covariates.csv")
fwrite(joined, out_cov_csv)
log_progress(sprintf("Wrote updated per-site covariate table: %s", out_cov_csv))

# ---- 5. Rebuild covariate list and run correlations -----------------

# Full covariate set now (drop the _avg_avg columns that were replaced).
cov_avg_cols <- grep("_avg$", names(joined), value = TRUE)
# Exclude anything that clearly isn't a covariate mean (defensive).
cov_avg_cols <- setdiff(cov_avg_cols, c("n_forest_fp"))

# Display labels for the full set (old + new).
cov_display_map <- c(
  "cover_avg"                   = "Canopy cover (%)",
  "rh_98_avg"                   = "GEDI RH98 (m)",
  "wsci_avg"                    = "WSCI",
  "slope_mean_avg"              = "Slope mean (deg)",
  "meta_az_concentration_avg"   = "Azimuth concentration",
  "meta_stereo_ratio_avg"       = "Stereo ratio",
  "meta_leaf_on_ratio_avg"      = "Leaf-on fraction",
  "meta_abs_geoacc_avg_avg"     = "Abs. geolocation acc.",
  "meta_rel_geoacc_avg_avg"     = "Rel. geolocation acc.",
  # New this patch:
  "meta_off_nadir_avg_avg"      = "Off-nadir angle (deg)",
  "meta_sun_elev_avg_avg"       = "Sun elevation (deg)",
  "meta_fwd_ratio_avg"          = "Fwd image fraction",
  "meta_rev_ratio_avg"          = "Rev image fraction",
  "meta_fwd_minus_rev_avg"      = "Fwd minus rev fraction",
  "meta_tot_ct_avg"             = "Total image count"
)

# Keep only columns we have labels for (safety filter).
cov_avg_cols <- intersect(cov_avg_cols, names(cov_display_map))

log_progress(sprintf("Correlating r(err_p3d, err_dtm) against %d covariates", length(cov_avg_cols)))

cor_tbl <- lapply(cov_avg_cols, function(cc) {
  x <- joined[[cc]]
  y <- joined$corr_errP3D_errDTM
  ok_all <- is.finite(x) & is.finite(y)
  nf <- joined$flag_status != "FLAGGED"

  p_all  <- tryCatch(cor.test(x[ok_all],          y[ok_all],          method = "pearson"),
                     error = function(e) NULL)
  s_all  <- tryCatch(cor.test(x[ok_all],          y[ok_all],          method = "spearman", exact = FALSE),
                     error = function(e) NULL)
  p_nf   <- tryCatch(cor.test(x[ok_all & nf],     y[ok_all & nf],     method = "pearson"),
                     error = function(e) NULL)
  s_nf   <- tryCatch(cor.test(x[ok_all & nf],     y[ok_all & nf],     method = "spearman", exact = FALSE),
                     error = function(e) NULL)

  data.table(
    covariate       = cc,
    covariate_label = cov_display_map[cc],
    n_18            = sum(ok_all),
    pearson_r_18    = if (!is.null(p_all)) unname(p_all$estimate) else NA_real_,
    pearson_p_18    = if (!is.null(p_all)) p_all$p.value           else NA_real_,
    spearman_rho_18 = if (!is.null(s_all)) unname(s_all$estimate)  else NA_real_,
    spearman_p_18   = if (!is.null(s_all)) s_all$p.value           else NA_real_,
    n_15            = sum(ok_all & nf),
    pearson_r_15    = if (!is.null(p_nf)) unname(p_nf$estimate) else NA_real_,
    pearson_p_15    = if (!is.null(p_nf)) p_nf$p.value           else NA_real_,
    spearman_rho_15 = if (!is.null(s_nf)) unname(s_nf$estimate) else NA_real_,
    spearman_p_15   = if (!is.null(s_nf)) s_nf$p.value           else NA_real_
  )
}) |> rbindlist()

cor_tbl[, abs_r_15 := abs(pearson_r_15)]
setorder(cor_tbl, -abs_r_15)
cor_tbl[, abs_r_15 := NULL]

out_cor_csv <- file.path(MT_DIR, "coupling_site_correlations.csv")
fwrite(cor_tbl, out_cor_csv)
log_progress(sprintf("Wrote correlations table: %s", out_cor_csv))

log_progress("Top 8 correlations with r(err_p3d, err_dtm), 15-site:")
for (i in seq_len(min(8L, nrow(cor_tbl)))) {
  r15 <- cor_tbl$pearson_r_15[i]
  p15 <- cor_tbl$pearson_p_15[i]
  r18 <- cor_tbl$pearson_r_18[i]
  lab <- cor_tbl$covariate_label[i]
  log_progress(sprintf("  %-30s  r_15=%+.3f (p=%.3f) | r_18=%+.3f",
                       lab, r15, p15, r18))
}

# ---- 6. Rebuild plots ----------------------------------------------

plot_df <- melt(joined,
                id.vars = c("tracker_site", "manuscript_site",
                            "site_short_name", "flag_status", "dominant_lc",
                            "corr_errP3D_errDTM"),
                measure.vars = cov_avg_cols,
                variable.name = "covariate", value.name = "cov_value")
plot_df[, covariate_label := cov_display_map[as.character(covariate)]]

facet_order <- cor_tbl$covariate_label
plot_df[, covariate_label := factor(covariate_label, levels = facet_order)]
plot_df[, flag := ifelse(flag_status == "FLAGGED", "flagged",
                  ifelse(flag_status == "DTM_excluded", "DTM-excl", "ok"))]

annot <- cor_tbl[, .(covariate_label,
                     lbl = sprintf("r15 = %+0.2f (p=%.2f)\nr18 = %+0.2f",
                                   pearson_r_15, pearson_p_15, pearson_r_18))]
annot[, covariate_label := factor(covariate_label, levels = facet_order)]

p_matrix <- ggplot(plot_df, aes(cov_value, corr_errP3D_errDTM)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "grey20",
              linewidth = 0.4, alpha = 0.15, na.rm = TRUE, formula = y ~ x) +
  geom_point(aes(color = flag, shape = dominant_lc), size = 2.5) +
  geom_text_repel(aes(label = site_short_name, color = flag),
                  size = 2.3, min.segment.length = 0.2,
                  max.overlaps = 18, seed = 1L) +
  geom_text(data = annot, aes(x = -Inf, y = Inf, label = lbl),
            hjust = -0.05, vjust = 1.2, size = 2.7, color = "grey25",
            inherit.aes = FALSE) +
  scale_color_manual(values = c("ok" = "#1b6cb0",
                                "flagged" = "#c0392b",
                                "DTM-excl" = "#8e44ad")) +
  facet_wrap(~ covariate_label, scales = "free_x", ncol = 3) +
  labs(
    x = "Site-level mean of covariate (forest footprints)",
    y = "Footprint-level r(err_p3d, err_dtm) from step 02",
    title = "Cross-site predictors of CNN-DTM compensation strength",
    subtitle = sprintf("N=18 sites. Facets ordered by |15-site Pearson r|. (%d covariates)",
                       length(cov_avg_cols)),
    color = "Site status", shape = "Dominant LC"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    panel.spacing = unit(0.7, "lines"),
    plot.title = element_text(face = "bold")
  )

out_matrix_pdf <- file.path(PLOT_DIR, "groundwork_task4_phase25_scatter_matrix.pdf")
# Larger canvas because we now have more panels.
ggsave(out_matrix_pdf, p_matrix,
       width = 12, height = max(10, 3 * ceiling(length(cov_avg_cols) / 3)))
log_progress(sprintf("Wrote scatter matrix: %s", out_matrix_pdf))

# Top-panel figure (4 strongest).
top_n <- min(4L, nrow(cor_tbl))
top_covs <- cor_tbl$covariate_label[seq_len(top_n)]
plot_top <- plot_df[covariate_label %in% top_covs]
plot_top[, covariate_label := factor(covariate_label, levels = top_covs)]
annot_top <- annot[covariate_label %in% top_covs]
annot_top[, covariate_label := factor(covariate_label, levels = top_covs)]

p_top <- ggplot(plot_top, aes(cov_value, corr_errP3D_errDTM)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "grey20",
              linewidth = 0.5, alpha = 0.18, na.rm = TRUE, formula = y ~ x) +
  geom_point(aes(color = flag, shape = dominant_lc), size = 3.2) +
  geom_text_repel(aes(label = site_short_name, color = flag),
                  size = 2.9, min.segment.length = 0.2,
                  max.overlaps = 18, seed = 1L) +
  geom_text(data = annot_top, aes(x = -Inf, y = Inf, label = lbl),
            hjust = -0.05, vjust = 1.2, size = 3.2, color = "grey25",
            inherit.aes = FALSE) +
  scale_color_manual(values = c("ok" = "#1b6cb0",
                                "flagged" = "#c0392b",
                                "DTM-excl" = "#8e44ad")) +
  facet_wrap(~ covariate_label, scales = "free_x", ncol = 2) +
  labs(
    x = "Site-level mean of covariate (forest footprints)",
    y = "Footprint-level r(err_p3d, err_dtm) from step 02",
    title = "Step 03 — strongest cross-site correlates",
    subtitle = sprintf("Top %d covariates by |Pearson r| on 15 non-flagged sites", top_n),
    color = "Site status", shape = "Dominant LC"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(0.8, "lines"),
    plot.title = element_text(face = "bold")
  )

out_top_pdf <- file.path(PLOT_DIR, "groundwork_task4_phase25_top_panel.pdf")
ggsave(out_top_pdf, p_top, width = 10, height = 9)
log_progress(sprintf("Wrote top-panel figure: %s", out_top_pdf))

# ---- 7. Update checkpoint ------------------------------------------

save_checkpoint("groundwork_task4_phase25", list(
  joined        = joined,
  cor_tbl       = cor_tbl,
  cov_avg_cols  = cov_avg_cols,
  cov_display   = cov_display_map,
  patched       = TRUE,
  patch_time    = Sys.time()
))

log_progress(strrep("=", 70))
log_progress("PATCH COMPLETE")
log_progress(strrep("=", 70))
