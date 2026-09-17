# =====================================================================
# coupling_03_site_correlates.R
#
# Step 03: site-level correlates of CNN-compensation strength.
#
# For each of 18 sites, step 02 produced a scalar `corr_errP3D_errDTM`
# describing how inversely correlated DSM error and CNN-DTM error are
# at footprint level. This script asks: is that per-site scalar
# predictable from site-level means of the manuscript covariates
# (canopy cover, RH98, WSCI, slope, acquisition geometry, phenology)?
#
# Precedent: Supp Fig S21 regressed site-level covariate means against
# CHM and DTM site random intercepts and found no surviving Bonferroni-
# corrected effect across 17 tests. This script uses the same site-level
# covariate means but a *different* response — the compensation
# strength — and is framed exploratorily (no multiple-testing correction),
# one scalar per site, 18 points.
#
# Inputs (all on disk, no cluster-side extraction):
#   manuscript_tables/coupling_site_summary.csv       (step 02)
#   manuscript_tables/site_id_lookup.csv
#   data/enriched_by_site/site_NN_enriched.csv[.gz]           (18 files)
#
# Outputs:
#   manuscript_tables/coupling_site_covariates.csv
#   manuscript_tables/coupling_site_correlations.csv
#   plots/groundwork/groundwork_task4_phase25_scatter_matrix.pdf
#   plots/groundwork/groundwork_task4_phase25_top_panel.pdf
#   checkpoints/groundwork_task4_phase25.rds
#
# Runtime: <1 minute on 9% of one cluster node.
# =====================================================================

# ---- 1. Environment and dependencies --------------------------------

# Require analysis_config.R and analysis_utils.R sourced prior
# (provides PROJECT_ROOT, log_progress, save_checkpoint, etc.). If
# missing, source them from the conventional location.
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
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(stringr)
  library(purrr)
})

# ---- 2. Paths and banner --------------------------------------------

MT_DIR     <- file.path(PROJECT_ROOT, "manuscript_tables")
ENRICH_DIR <- file.path(PROJECT_ROOT, "data", "enriched_by_site")
PLOT_DIR   <- file.path(PROJECT_ROOT, "plots", "groundwork")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

log_progress(strrep("=", 70))
log_progress("Step 03 — correlates of CNN-compensation strength")
log_progress(strrep("=", 70))
log_progress(sprintf("PROJECT_ROOT:    %s", PROJECT_ROOT))
log_progress(sprintf("manuscript_tbls: %s", MT_DIR))
log_progress(sprintf("enriched dir:    %s", ENRICH_DIR))
log_progress(sprintf("plots dir:       %s", PLOT_DIR))

# ---- 3. Load step 02 per-site summary -------------------------------

site_sum_path <- file.path(MT_DIR, "coupling_site_summary.csv")
stopifnot(file.exists(site_sum_path))
site_sum <- fread(site_sum_path)

log_progress(sprintf("Loaded step 02 site summary: %d rows, %d cols",
                     nrow(site_sum), ncol(site_sum)))

# Basic sanity: expect 18 rows (all sites with 3DEP coverage) and the
# correlation column that names the response variable for this script.
stopifnot(nrow(site_sum) == 18L)
stopifnot("corr_errP3D_errDTM" %in% names(site_sum))
stopifnot("tracker_site"       %in% names(site_sum))
stopifnot("manuscript_site"    %in% names(site_sum))
stopifnot("flag_status"        %in% names(site_sum))

log_progress(sprintf("  r(err_p3d, err_dtm): median %.3f  mean %.3f  range [%.3f, %.3f]",
                     median(site_sum$corr_errP3D_errDTM),
                     mean(site_sum$corr_errP3D_errDTM),
                     min(site_sum$corr_errP3D_errDTM),
                     max(site_sum$corr_errP3D_errDTM)))
log_progress(sprintf("  Flag status counts: %s",
                     paste(sprintf("%s=%d", names(table(site_sum$flag_status)),
                                   as.integer(table(site_sum$flag_status))),
                           collapse = ", ")))

# ---- 4. Per-site forest-footprint covariate means -------------------
#
# Covariates, matched to Supp Fig S21's list (minus DTM-only vars that
# aren't relevant to the CHM-side compensation question):
#
#   Canopy / structure:        cover, rh_98, wsci
#   Terrain:                   slope_mean
#   Acquisition geometry:      meta_offnad_avg, meta_sunelev_avg,
#                              meta_az_concentration, meta_stereo_ratio,
#                              meta_fwd_rev_ratio, meta_leaf_on_ratio,
#                              meta_abs_geoacc_avg, meta_rel_geoacc_avg
#
# Forest footprints = lc_l1_code in c("EBF","BDF","ENF","DNF"),
# matching manuscript Set 1/3 convention and Supp Fig S21.

FOREST_CLASSES <- c("EBF", "BDF", "ENF", "DNF")

COV_COLS <- c(
  "cover", "rh_98", "wsci", "slope_mean",
  "meta_offnad_avg", "meta_sunelev_avg", "meta_az_concentration",
  "meta_stereo_ratio", "meta_fwd_rev_ratio", "meta_leaf_on_ratio",
  "meta_abs_geoacc_avg", "meta_rel_geoacc_avg"
)

# Discover per-site enriched CSVs (zero-padded filenames).
enriched_files <- list.files(ENRICH_DIR,
  pattern = "^site_[0-9]+_enriched\\.csv(\\.gz)?$",
  full.names = TRUE
)
log_progress(sprintf("Discovered %d enriched CSV files in %s",
                     length(enriched_files), ENRICH_DIR))

get_tracker_id <- function(p) {
  m <- regmatches(basename(p), regexec("site_([0-9]+)_enriched", basename(p)))[[1]]
  if (length(m) == 2L) as.integer(m[2]) else NA_integer_
}

# Compute site-level forest-footprint means.
compute_site_means <- function(fp) {
  tid <- get_tracker_id(fp)
  # Only keep the columns we need, reducing IO.
  hdr <- names(fread(fp, nrows = 0L))
  want <- unique(c("lc_l1_code", COV_COLS))
  have <- intersect(want, hdr)

  # Some column-name compatibility: older exports may have
  # lc2022_mode_l1_code instead of lc_l1_code.
  if (!("lc_l1_code" %in% have)) {
    if ("lc2022_mode_l1_code" %in% hdr) {
      DT <- fread(fp, select = c("lc2022_mode_l1_code", intersect(COV_COLS, hdr)),
                  showProgress = FALSE)
      setnames(DT, "lc2022_mode_l1_code", "lc_l1_code")
    } else if ("lc2022_l1_code" %in% hdr) {
      DT <- fread(fp, select = c("lc2022_l1_code", intersect(COV_COLS, hdr)),
                  showProgress = FALSE)
      setnames(DT, "lc2022_l1_code", "lc_l1_code")
    } else {
      stop(sprintf("No usable LC column in %s", basename(fp)))
    }
  } else {
    DT <- fread(fp, select = have, showProgress = FALSE)
  }

  # Forest-only filter matching manuscript Set 1/3 convention.
  DT_forest <- DT[lc_l1_code %in% FOREST_CLASSES]

  # Per-site means and SDs of each covariate (SDs useful to inspect
  # variation; main output uses means to match Supp Fig S21 precedent).
  present_covs <- intersect(COV_COLS, names(DT_forest))
  means <- DT_forest[, lapply(.SD, function(x) mean(x, na.rm = TRUE)),
                     .SDcols = present_covs]
  setnames(means, present_covs, paste0(present_covs, "_avg"))

  # Also report forest-footprint count per site as a sanity check.
  means[, n_forest_fp := nrow(DT_forest)]
  means[, tracker_site := tid]

  # Fill missing covariates with NA so rbindlist works across sites
  missing <- setdiff(paste0(COV_COLS, "_avg"), names(means))
  for (m in missing) means[, (m) := NA_real_]

  setcolorder(means, c("tracker_site", "n_forest_fp",
                       paste0(COV_COLS, "_avg")))
  means[]
}

log_progress("Computing forest-footprint covariate means per site...")
t0 <- Sys.time()
per_site_cov <- rbindlist(lapply(enriched_files, function(fp) {
  r <- tryCatch(compute_site_means(fp), error = function(e) {
    log_progress(sprintf("  WARN: %s — %s", basename(fp), e$message))
    NULL
  })
  if (!is.null(r))
    log_progress(sprintf("  site_%02d: %d forest footprints",
                         r$tracker_site, r$n_forest_fp))
  r
}), fill = TRUE)
t1 <- Sys.time()
log_progress(sprintf("  Done in %.1f sec (%d sites)",
                     as.numeric(difftime(t1, t0, units = "secs")),
                     nrow(per_site_cov)))

# ---- 5. Join to step 02 summary -------------------------------------

# The step 02 site_summary keys on tracker_site (1..20 except 16, 21).
joined <- merge(site_sum, per_site_cov,
                by = "tracker_site", all.x = TRUE, sort = FALSE)

stopifnot(nrow(joined) == 18L)
log_progress(sprintf("Joined table: %d rows x %d cols", nrow(joined), ncol(joined)))

# Write the joined per-site covariate table for reference.
out_cov_csv <- file.path(MT_DIR, "coupling_site_covariates.csv")
fwrite(joined, out_cov_csv)
log_progress(sprintf("Wrote per-site covariate table: %s", out_cov_csv))

# ---- 6. Cross-site correlations -------------------------------------

cov_avg_cols <- paste0(COV_COLS, "_avg")
cov_avg_cols <- intersect(cov_avg_cols, names(joined))

# Human-readable names for plotting and the output table.
cov_display_map <- c(
  "cover_avg"                  = "Canopy cover (%)",
  "rh_98_avg"                  = "GEDI RH98 (m)",
  "wsci_avg"                   = "WSCI",
  "slope_mean_avg"             = "Slope mean (deg)",
  "meta_offnad_avg_avg"        = "Off-nadir angle (mean deg)",
  "meta_sunelev_avg_avg"       = "Sun elevation (mean deg)",
  "meta_az_concentration_avg"  = "Azimuth concentration",
  "meta_stereo_ratio_avg"      = "Stereo ratio",
  "meta_fwd_rev_ratio_avg"     = "Fwd/rev image ratio",
  "meta_leaf_on_ratio_avg"     = "Leaf-on fraction",
  "meta_abs_geoacc_avg_avg"    = "Abs. geolocation acc.",
  "meta_rel_geoacc_avg_avg"    = "Rel. geolocation acc."
)

cor_tbl <- lapply(cov_avg_cols, function(cc) {
  # All 18 sites
  x <- joined[[cc]]
  y <- joined$corr_errP3D_errDTM
  ok_all <- is.finite(x) & is.finite(y)
  p_all  <- tryCatch(cor.test(x[ok_all], y[ok_all], method = "pearson"),
                     error = function(e) NULL)
  s_all  <- tryCatch(cor.test(x[ok_all], y[ok_all], method = "spearman",
                              exact = FALSE),
                     error = function(e) NULL)

  # Non-flagged 15 sites (sensitivity)
  nf <- joined$flag_status != "FLAGGED"
  p_nf <- tryCatch(cor.test(x[ok_all & nf], y[ok_all & nf], method = "pearson"),
                   error = function(e) NULL)
  s_nf <- tryCatch(cor.test(x[ok_all & nf], y[ok_all & nf], method = "spearman",
                            exact = FALSE),
                   error = function(e) NULL)

  data.table(
    covariate          = cc,
    covariate_label    = ifelse(cc %in% names(cov_display_map),
                                cov_display_map[cc], cc),
    n_18               = sum(ok_all),
    pearson_r_18       = if (!is.null(p_all)) unname(p_all$estimate) else NA_real_,
    pearson_p_18       = if (!is.null(p_all)) p_all$p.value           else NA_real_,
    spearman_rho_18    = if (!is.null(s_all)) unname(s_all$estimate)  else NA_real_,
    spearman_p_18      = if (!is.null(s_all)) s_all$p.value           else NA_real_,
    n_15               = sum(ok_all & nf),
    pearson_r_15       = if (!is.null(p_nf)) unname(p_nf$estimate) else NA_real_,
    pearson_p_15       = if (!is.null(p_nf)) p_nf$p.value           else NA_real_,
    spearman_rho_15    = if (!is.null(s_nf)) unname(s_nf$estimate) else NA_real_,
    spearman_p_15      = if (!is.null(s_nf)) s_nf$p.value           else NA_real_
  )
}) |> rbindlist()

# Rank by magnitude of 15-site Pearson r (cleaner read — 18-site is
# pulled by the flagged trio at positive r values).
cor_tbl[, abs_r_15 := abs(pearson_r_15)]
setorder(cor_tbl, -abs_r_15)
cor_tbl[, abs_r_15 := NULL]

out_cor_csv <- file.path(MT_DIR, "coupling_site_correlations.csv")
fwrite(cor_tbl, out_cor_csv)
log_progress(sprintf("Wrote correlations table: %s", out_cor_csv))

# Echo the top 5 strongest relationships to the console.
log_progress("Top cross-site correlations with r(err_p3d, err_dtm), 15-site:")
for (i in seq_len(min(5L, nrow(cor_tbl)))) {
  r15 <- cor_tbl$pearson_r_15[i]
  p15 <- cor_tbl$pearson_p_15[i]
  r18 <- cor_tbl$pearson_r_18[i]
  lab <- cor_tbl$covariate_label[i]
  log_progress(sprintf("  %-30s  r_15=%+.3f (p=%.3f) | r_18=%+.3f",
                       lab, r15, p15, r18))
}

# ---- 7. Scatter matrix (all covariates) -----------------------------

# Tidy shape for ggplot facets.
plot_df <- melt(joined,
                id.vars = c("tracker_site", "manuscript_site",
                            "site_short_name", "flag_status", "dominant_lc",
                            "corr_errP3D_errDTM"),
                measure.vars = cov_avg_cols,
                variable.name = "covariate", value.name = "cov_value")
plot_df[, covariate_label := ifelse(covariate %in% names(cov_display_map),
                                     cov_display_map[as.character(covariate)],
                                     as.character(covariate))]
# Preserve facet order by |r_15| ranking
facet_order <- cor_tbl$covariate_label
plot_df[, covariate_label := factor(covariate_label, levels = facet_order)]
plot_df[, flag := ifelse(flag_status == "FLAGGED", "flagged",
                  ifelse(flag_status == "DTM_excluded", "DTM-excl", "ok"))]

# Per-panel correlation annotation (15-site r).
annot <- cor_tbl[, .(covariate_label,
                     lbl = sprintf("r15 = %+0.2f (p=%.2f)\nr18 = %+0.2f",
                                   pearson_r_15, pearson_p_15, pearson_r_18))]
annot[, covariate_label := factor(covariate_label, levels = facet_order)]

# Free x across facets (covariates are on different scales).
p_matrix <- ggplot(plot_df, aes(cov_value, corr_errP3D_errDTM)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "grey20",
              linewidth = 0.4, alpha = 0.15, na.rm = TRUE,
              formula = y ~ x) +
  geom_point(aes(color = flag, shape = dominant_lc), size = 2.5) +
  geom_text_repel(aes(label = site_short_name, color = flag),
                  size = 2.4, min.segment.length = 0.2,
                  max.overlaps = 18, seed = 1L) +
  geom_text(data = annot, aes(x = -Inf, y = Inf, label = lbl),
            hjust = -0.05, vjust = 1.2, size = 2.8, color = "grey25",
            inherit.aes = FALSE) +
  scale_color_manual(values = c("ok" = "#1b6cb0",
                                "flagged" = "#c0392b",
                                "DTM-excl" = "#8e44ad")) +
  facet_wrap(~ covariate_label, scales = "free_x", ncol = 3) +
  labs(
    x = "Site-level mean of covariate (forest footprints)",
    y = "Footprint-level r(err_p3d, err_dtm) from step 02",
    title = "Cross-site predictors of CNN-DTM compensation strength",
    subtitle = "Each point is one site (N=18). Facets ordered by |15-site Pearson r|.",
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
ggsave(out_matrix_pdf, p_matrix, width = 12, height = 12)
log_progress(sprintf("Wrote scatter matrix: %s", out_matrix_pdf))

# ---- 8. Top-panel figure (4 strongest) ------------------------------

top_n <- min(4L, nrow(cor_tbl))
top_covs <- cor_tbl$covariate_label[seq_len(top_n)]
plot_top <- plot_df[covariate_label %in% top_covs]
plot_top[, covariate_label := factor(covariate_label, levels = top_covs)]
annot_top <- annot[covariate_label %in% top_covs]
annot_top[, covariate_label := factor(covariate_label, levels = top_covs)]

p_top <- ggplot(plot_top, aes(cov_value, corr_errP3D_errDTM)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "grey20",
              linewidth = 0.5, alpha = 0.18, na.rm = TRUE,
              formula = y ~ x) +
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

# ---- 9. Checkpoint --------------------------------------------------

save_checkpoint("groundwork_task4_phase25", list(
  joined        = joined,
  cor_tbl       = cor_tbl,
  per_site_cov  = per_site_cov,
  cov_avg_cols  = cov_avg_cols,
  cov_display   = cov_display_map
))

log_progress(strrep("=", 70))
log_progress("Step 03 COMPLETE")
log_progress(strrep("=", 70))
log_progress("Outputs:")
log_progress(sprintf("  %s", out_cov_csv))
log_progress(sprintf("  %s", out_cor_csv))
log_progress(sprintf("  %s", out_matrix_pdf))
log_progress(sprintf("  %s", out_top_pdf))
log_progress("  checkpoint: groundwork_task4_phase25")
log_progress("")
log_progress("Interpretation: see cor_tbl, particularly pearson_r_15 column.")
log_progress("Strong negative r (e.g. -0.5 or stronger) against a canopy")
log_progress("covariate would corroborate the CNN-contamination-drives-")
log_progress("compensation mechanism. Strong positive r against terrain or")
log_progress("geometry covariates would instead point to DSM-side drivers.")
