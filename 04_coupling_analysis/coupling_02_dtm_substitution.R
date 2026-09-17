# =====================================================================
# coupling_02_dtm_substitution.R
#
# Step 02: Substitute 3DEP DTM for the CNN-inferred P3D DTM in
# the CHM equation, and quantify how much CHM error "collapses" when
# the DTM is replaced with a ground-truthed terrain surface.
#
# Algebraic framing:
#   P3D_DSM = p3d_chm_mean + p3d_dtm_mean          (footprint-mean identity;
#                                                   verified by step 01)
#   chm_p3d    = p3d_chm_mean                        (original P3D CHM)
#   chm_alt    = P3D_DSM - dep_dtm_mean              (CNN DTM replaced by 3DEP)
#              = (p3d_chm_mean + p3d_dtm_mean) - dep_dtm_mean
#              = p3d_chm_mean + dtm_error_raw
#   chm_ref    = als_chm_mean                        (the CHM reference)
#
# Three footprint-level error series:
#   err_p3d = chm_p3d - chm_ref                      (original CHM error,
#                                                     same as Ch1 response)
#   err_alt = chm_alt - chm_ref                      (CHM error with ALS
#                                                     terrain instead of CNN)
#   err_dtm = p3d_dtm_mean - dep_dtm_mean            (DTM error; matches
#                                                     Ch1 dtm_error_mean)
#
# Identity:  err_alt = err_p3d - err_dtm
# So the variance decomposition is:
#   Var(err_alt) = Var(err_p3d) + Var(err_dtm) - 2*Cov(err_p3d, err_dtm)
# Interpretations:
#   - If err_dtm is a large component of err_p3d (positive correlation),
#     Var(err_alt) shrinks --> "the CNN DTM is the dominant CHM error
#     source; fixing the DTM largely fixes the CHM."
#   - If err_dtm is independent of err_p3d, Var(err_alt) INFLATES by
#     Var(err_dtm) --> "the DSM and the CNN DTM fail independently; the
#     CHM error is mostly DSM-side."
#
# Scope:
#   - 18 sites (19 manuscript sites minus Site 10 neon_clbj2022,
#     DTM-excluded; excluded because dep_dtm_mean is unavailable there).
#   - Forest footprints only (als_chm_p90 >= chm_forest_thresh_m = 2 m).
#   - Footprints must have all four of p3d_chm_mean, p3d_dtm_mean,
#     als_chm_mean, dep_dtm_mean finite with valid_frac >= 0.5 on each.
#
# Outputs (CSVs -> $PROJECT_ROOT/manuscript_tables/, figures ->
# $PROJECT_ROOT/plots/groundwork/):
#   coupling_footprint_counts.csv
#   coupling_site_summary.csv
#   coupling_overall_summary.csv
#   groundwork_task4_footprint_data.rds
#   groundwork_task4_err_dist_by_site.pdf
#   groundwork_task4_variance_reduction_bar.pdf
#   groundwork_task4_err_scatter_by_site.pdf
#   groundwork_task4_err_vs_dtm_error.pdf
#
# Checkpoint: groundwork_task4_phase2 (via save_checkpoint())
#
# Runtime estimate: ~15-30 minutes for ingest + analysis + figures
# (the data volume is similar to what section_01 handles).
# =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

if (!exists("log_progress"))   source("analysis_utils.R")
if (!exists("PROJECT_ROOT"))   source("analysis_config.R")

log_progress("=== Step 02: CHM alt partition analysis ===")

# ---- Environment sanity: verify we can find everything we need ----

cat(sprintf("  PROJECT_ROOT = %s\n", PROJECT_ROOT))
cat(sprintf("  ENRICHED_DIR = %s\n", ENRICHED_DIR))
cat(sprintf("  getwd()      = %s\n", getwd()))

MS_TABLES_DIR   <- file.path(PROJECT_ROOT, "manuscript_tables")
GROUNDWORK_PLOT <- file.path(PROJECT_ROOT, "plots", "groundwork")

.required_dirs <- c(
  "PROJECT_ROOT" = PROJECT_ROOT,
  "ENRICHED_DIR" = ENRICHED_DIR,
  "manuscript_tables" = MS_TABLES_DIR
)
.missing <- .required_dirs[!vapply(.required_dirs, dir.exists, logical(1))]
if (length(.missing)) {
  stop("Required directories are missing:\n  ",
       paste(names(.missing), "->", .missing, collapse = "\n  "))
}

dir.create(GROUNDWORK_PLOT, recursive = TRUE, showWarnings = FALSE)

# ---- Load the site ID lookup ----

lookup <- data.table::fread(file.path(MS_TABLES_DIR, "site_id_lookup.csv"))
log_progress(sprintf("Loaded site_id_lookup.csv (%d rows)", nrow(lookup)))

# Keep only manuscript sites (drop neon_nogp and nasa_sonoma_whole,
# which have NA manuscript_site in the lookup).
lookup_ms <- lookup[!is.na(manuscript_site)]

# Drop Site 10 (neon_clbj2022) -- no 3DEP coverage.
task4_sites_ms <- lookup_ms[manuscript_site != 10]
task4_trackers <- task4_sites_ms$tracker_site
log_progress(sprintf(
  "Scope: %d sites (manuscript %s), Site 10 excluded (no 3DEP)",
  nrow(task4_sites_ms),
  paste(range(task4_sites_ms$manuscript_site), collapse = "-")
))

# ---- Per-site ingest ----
#
# Read each enriched CSV, keeping only the columns step 02 needs, and
# track footprint counts through the filter stages.

.read_site <- function(site_id_tracker) {
  ms_row <- lookup_ms[tracker_site == site_id_tracker]
  ms_id  <- ms_row$manuscript_site[1]
  short  <- ms_row$site_short_name[1]
  flag   <- ms_row$flag_status[1]
  dlc    <- ms_row$dominant_lc[1]

  patterns <- c(
    file.path(ENRICHED_DIR, sprintf("site_%02d_enriched.csv.gz", site_id_tracker)),
    file.path(ENRICHED_DIR, sprintf("site_%d_enriched.csv.gz",   site_id_tracker)),
    file.path(ENRICHED_DIR, sprintf("site_%02d_enriched.csv",    site_id_tracker)),
    file.path(ENRICHED_DIR, sprintf("site_%d_enriched.csv",      site_id_tracker))
  )
  path <- patterns[file.exists(patterns)][1]
  if (is.na(path)) {
    log_progress(sprintf("  WARNING: no CSV for tracker site %d; skipping",
                         site_id_tracker))
    return(NULL)
  }

  hdr <- names(data.table::fread(path, nrows = 0))
  needed <- c(
    "shot_number", "site",
    "lc2022_l1_code", "lc2022_l1_name",
    "lc2022_mode_l1_code", "lc2022_mode_l1_name",
    "p3d_chm_mean", "p3d_dtm_mean",
    "als_chm_mean", "dep_dtm_mean",
    "als_chm_p90",
    "p3d_chm_valid_frac", "p3d_dtm_valid_frac",
    "als_chm_valid_frac", "dep_dtm_valid_frac",
    "slope_valid_frac"
  )
  have <- intersect(needed, hdr)
  DT <- data.table::fread(path, select = have,
                          colClasses = list(character = "shot_number"),
                          nThread = min(4L, CPU_BUDGET),
                          showProgress = FALSE)

  # LC column normalization, matching section_01_ingest.R
  if (!("lc2022_l1_code" %in% names(DT)) &&
      "lc2022_mode_l1_code" %in% names(DT)) {
    setnames(DT, "lc2022_mode_l1_code", "lc2022_l1_code")
  }
  if (!("lc2022_l1_name" %in% names(DT)) &&
      "lc2022_mode_l1_name" %in% names(DT)) {
    setnames(DT, "lc2022_mode_l1_name", "lc2022_l1_name")
  }

  n_read <- nrow(DT)

  # Count survivors at each filter stage.
  n_finite_p3d <- DT[is.finite(p3d_chm_mean) & is.finite(p3d_dtm_mean) &
                     is.finite(p3d_chm_valid_frac) & p3d_chm_valid_frac >= 0.5 &
                     is.finite(p3d_dtm_valid_frac) & p3d_dtm_valid_frac >= 0.5, .N]
  n_finite_als <- DT[is.finite(als_chm_mean) &
                     is.finite(als_chm_valid_frac) & als_chm_valid_frac >= 0.5, .N]
  n_finite_dep <- DT[is.finite(dep_dtm_mean) &
                     is.finite(dep_dtm_valid_frac) & dep_dtm_valid_frac >= 0.5, .N]

  # Full four-way filter + forest filter + slope valid.
  DT_clean <- DT[
    is.finite(p3d_chm_mean)       & is.finite(p3d_dtm_mean)       &
    is.finite(als_chm_mean)       & is.finite(dep_dtm_mean)       &
    is.finite(als_chm_p90)                                        &
    is.finite(p3d_chm_valid_frac) & p3d_chm_valid_frac >= 0.5     &
    is.finite(p3d_dtm_valid_frac) & p3d_dtm_valid_frac >= 0.5     &
    is.finite(als_chm_valid_frac) & als_chm_valid_frac >= 0.5     &
    is.finite(dep_dtm_valid_frac) & dep_dtm_valid_frac >= 0.5     &
    is.finite(slope_valid_frac)   & slope_valid_frac   >= 0.5
  ]
  n_all_four <- nrow(DT_clean)

  DT_clean <- DT_clean[als_chm_p90 >= chm_forest_thresh_m]
  n_forest <- nrow(DT_clean)

  if (n_forest == 0) {
    log_progress(sprintf(
      "  [tracker %d / manuscript %s] %s: 0 forest footprints after filter; skipping",
      site_id_tracker, as.character(ms_id), short
    ))
    return(NULL)
  }

  # Compute the three error series.
  # Match manuscript convention: positive = P3D overestimation.
  DT_clean[, `:=`(
    err_p3d = p3d_chm_mean - als_chm_mean,
    err_alt = (p3d_chm_mean + p3d_dtm_mean) - dep_dtm_mean - als_chm_mean,
    err_dtm = p3d_dtm_mean - dep_dtm_mean
  )]

  # |error| <= 100 m guard, matching section_01_ingest.R.
  DT_clean <- DT_clean[abs(err_p3d) <= 100 & abs(err_alt) <= 100 &
                       abs(err_dtm) <= 100]
  n_guarded <- nrow(DT_clean)

  # Attach metadata.
  DT_clean[, `:=`(
    tracker_site    = site_id_tracker,
    manuscript_site = ms_id,
    site_short_name = short,
    flag_status     = flag,
    dominant_lc     = dlc
  )]

  # Footprint-count record for this site.
  counts <- data.table(
    tracker_site           = site_id_tracker,
    manuscript_site        = ms_id,
    site_short_name        = short,
    flag_status            = flag,
    dominant_lc            = dlc,
    n_rows_read            = n_read,
    n_finite_p3d           = n_finite_p3d,
    n_finite_als_chm       = n_finite_als,
    n_finite_dep_dtm       = n_finite_dep,
    n_all_four_valid       = n_all_four,
    n_forest_p90_ge2       = n_forest,
    n_final_abs_err_le_100 = n_guarded
  )

  list(data = DT_clean, counts = counts)
}

log_progress("Reading per-site enriched CSVs...")

results <- lapply(task4_trackers, function(t) {
  log_progress(sprintf("  tracker site %d", t))
  .read_site(t)
})
results <- results[!sapply(results, is.null)]
log_progress(sprintf("Read %d sites successfully", length(results)))

fp_data <- rbindlist(lapply(results, `[[`, "data"), fill = TRUE, use.names = TRUE)
fp_counts <- rbindlist(lapply(results, `[[`, "counts"),
                       fill = TRUE, use.names = TRUE)
fp_counts <- fp_counts[order(manuscript_site)]

log_progress(sprintf("Total footprints across all sites: %s",
                     format(nrow(fp_data), big.mark = ",")))

# ---- Write footprint counts CSV ----

fwrite(fp_counts, file.path(MS_TABLES_DIR, "coupling_footprint_counts.csv"))
log_progress("Wrote coupling_footprint_counts.csv")

# ---- Per-site summary ----

.site_summary_stats <- function(d) {
  e1 <- d$err_p3d
  e2 <- d$err_alt
  e3 <- d$err_dtm

  v1 <- var(e1)
  v2 <- var(e2)

  # Correlation between err_p3d and err_dtm: a positive value means the
  # CNN DTM error systematically shares sign with the CHM error, i.e.
  # substitution should reduce variance.
  corr_p3d_dtm <- if (length(e1) >= 2 && sd(e1) > 0 && sd(e3) > 0) {
    cor(e1, e3)
  } else NA_real_

  data.table(
    n                 = length(e1),
    # err_p3d stats
    p3d_mean_bias     = mean(e1),
    p3d_median_bias   = median(e1),
    p3d_sd            = sd(e1),
    p3d_mae           = mean(abs(e1)),
    p3d_rmse          = sqrt(mean(e1 ^ 2)),
    # err_alt stats
    alt_mean_bias     = mean(e2),
    alt_median_bias   = median(e2),
    alt_sd            = sd(e2),
    alt_mae           = mean(abs(e2)),
    alt_rmse          = sqrt(mean(e2 ^ 2)),
    # err_dtm stats (for reference)
    dtm_mean_bias     = mean(e3),
    dtm_sd            = sd(e3),
    dtm_rmse          = sqrt(mean(e3 ^ 2)),
    # Partition metrics
    var_err_p3d       = v1,
    var_err_alt       = v2,
    var_reduction_frac = 1 - v2 / v1,   # Var(err_alt)/Var(err_p3d)
    rmse_reduction_frac = 1 - sqrt(mean(e2 ^ 2)) / sqrt(mean(e1 ^ 2)),
    corr_errP3D_errDTM = corr_p3d_dtm
  )
}

site_summary <- fp_data[, .site_summary_stats(.SD),
                        by = .(manuscript_site, tracker_site,
                               site_short_name, flag_status, dominant_lc)]
site_summary <- site_summary[order(manuscript_site)]

fwrite(site_summary, file.path(MS_TABLES_DIR, "coupling_site_summary.csv"))
log_progress("Wrote coupling_site_summary.csv")

# ---- Overall summaries (pooled) ----

.pooled_stats <- function(d, label) {
  s <- .site_summary_stats(d)
  s[, stratum := label]
  setcolorder(s, "stratum")
  s
}

overall_rows <- list(
  .pooled_stats(fp_data, "all_18_sites"),
  .pooled_stats(fp_data[flag_status != "FLAGGED"], "non_flagged_15_sites")
)

# Forest-type pooled rows (BDF, ENF, IWL, GRS, other).
# Use the site's dominant_lc for stratification since that's what the
# lookup encodes; footprint-level lc_l1_code exists too but per-site LC
# tagging from the lookup is what the §6 manuscript analysis used.
lc_classes <- unique(fp_data$dominant_lc)
for (lc in lc_classes) {
  sub <- fp_data[dominant_lc == lc]
  if (nrow(sub) >= 50) {
    overall_rows[[length(overall_rows) + 1L]] <-
      .pooled_stats(sub, sprintf("lc_%s", lc))
  }
}

overall_summary <- rbindlist(overall_rows, fill = TRUE, use.names = TRUE)
fwrite(overall_summary, file.path(MS_TABLES_DIR,
                                  "coupling_overall_summary.csv"))
log_progress("Wrote coupling_overall_summary.csv")

# ---- Save the footprint-level data for downstream refit (step 04) ----

save_checkpoint("groundwork_task4_phase2", list(
  fp_data         = fp_data,
  fp_counts       = fp_counts,
  site_summary    = site_summary,
  overall_summary = overall_summary,
  lookup_ms       = lookup_ms,
  timestamp       = Sys.time()
))
log_progress("Saved checkpoint groundwork_task4_phase2")

saveRDS(fp_data, file.path(MS_TABLES_DIR, "groundwork_task4_footprint_data.rds"),
        compress = "xz")
log_progress("Wrote groundwork_task4_footprint_data.rds")

# =====================================================================
# FIGURES
# =====================================================================

# Consistent site ordering and labeling for all figures.
site_order <- site_summary[order(manuscript_site)]
fp_data[, site_label := factor(
  sprintf("Site %02d\n%s", manuscript_site, site_short_name),
  levels = sprintf("Site %02d\n%s",
                   site_order$manuscript_site,
                   site_order$site_short_name)
)]

# Flag annotation for panels.
fp_data[, flag_annot := ifelse(flag_status == "FLAGGED", "★ FLAGGED", "")]

# Consistent color palette: P3D (orig) vs ALT (3DEP DTM swap).
err_pal <- c("err_p3d" = "#D45E5E", "err_alt" = "#2C6E9B")
err_lab <- c("err_p3d" = "Original CHM error (P3D DSM − P3D DTM)",
             "err_alt" = "Alt CHM error (P3D DSM − 3DEP DTM)")

# ---- Figure 1: per-site error distribution, err_p3d vs err_alt ----

log_progress("Figure 1: error distributions by site...")

long_err <- melt(
  fp_data[, .(manuscript_site, site_label, flag_status, dominant_lc,
              err_p3d, err_alt)],
  id.vars = c("manuscript_site", "site_label", "flag_status", "dominant_lc"),
  measure.vars = c("err_p3d", "err_alt"),
  variable.name = "series", value.name = "error"
)

# Per-site medians for annotation.
med_annot <- site_summary[, .(
  manuscript_site, site_label = sprintf("Site %02d\n%s", manuscript_site,
                                        site_short_name),
  med_p3d = p3d_median_bias,
  med_alt = alt_median_bias,
  n       = n,
  flag_status
)]
med_annot[, site_label := factor(site_label,
                                 levels = levels(fp_data$site_label))]

# Clamp for display (heavy tails compress the densities otherwise).
x_lim <- c(-25, 25)

p_fig1 <- ggplot(long_err, aes(x = error, fill = series, colour = series)) +
  geom_density(alpha = 0.35, linewidth = 0.4, adjust = 1.2) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "black") +
  facet_wrap(~ site_label, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = err_pal, labels = err_lab,
                    name = NULL, aesthetics = c("fill", "colour")) +
  coord_cartesian(xlim = x_lim) +
  labs(
    x = "CHM error (m), positive = P3D overestimates",
    y = "Density",
    title = "Figure 1: Per-site CHM error distributions — original vs 3DEP-DTM substitution",
    subtitle = "x-axis clamped to ±25 m for display; full tails included in summary statistics"
  ) +
  theme_cowplot(font_size = 10) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(size = 8, lineheight = 0.9),
    panel.spacing = unit(0.3, "lines")
  )

ggsave(file.path(GROUNDWORK_PLOT, "groundwork_task4_err_dist_by_site.pdf"),
       p_fig1, width = 14, height = 10)
log_progress("  wrote groundwork_task4_err_dist_by_site.pdf")

# ---- Figure 2: per-site variance reduction bar chart ----

log_progress("Figure 2: per-site variance reduction bar chart...")

bar_df <- site_summary[, .(
  manuscript_site, site_short_name, flag_status, dominant_lc,
  var_err_p3d, var_err_alt, var_reduction_frac,
  n
)]
bar_df[, label := sprintf("Site %02d %s%s",
                          manuscript_site,
                          site_short_name,
                          ifelse(flag_status == "FLAGGED", " ★", ""))]
bar_df[, label := factor(label,
                         levels = label[order(var_reduction_frac)])]
bar_df[, bar_fill := ifelse(flag_status == "FLAGGED", "FLAGGED",
                     ifelse(var_reduction_frac >= 0, "reduction", "inflation"))]

p_fig2 <- ggplot(bar_df, aes(x = var_reduction_frac, y = label,
                             fill = bar_fill)) +
  geom_col() +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.2f (n=%s)",
                                var_reduction_frac,
                                format(n, big.mark = ","))),
            hjust = ifelse(bar_df$var_reduction_frac >= 0, -0.05, 1.05),
            size = 3) +
  scale_fill_manual(values = c(FLAGGED   = "#B06DB3",
                               reduction = "#2C6E9B",
                               inflation = "#D45E5E")) +
  scale_x_continuous(
    labels = label_percent(),
    expand = expansion(mult = c(0.10, 0.20))
  ) +
  labs(
    x = "Variance reduction when CNN DTM is replaced with 3DEP DTM\n(1 − Var(err_alt) / Var(err_p3d))",
    y = NULL,
    title = "Figure 2: Per-site CHM variance reduction from DTM substitution",
    subtitle = "Positive = substituting the 3DEP DTM into the CHM equation reduces variance.\nNegative = substitution INFLATES variance (CNN DTM and DSM errors were cancelling each other).",
    fill = NULL
  ) +
  theme_cowplot(font_size = 10) +
  theme(legend.position = "top")

ggsave(file.path(GROUNDWORK_PLOT, "groundwork_task4_variance_reduction_bar.pdf"),
       p_fig2, width = 11, height = 9)
log_progress("  wrote groundwork_task4_variance_reduction_bar.pdf")

# ---- Figure 3: per-site err_p3d vs err_alt scatter ----

log_progress("Figure 3: per-site err_p3d vs err_alt scatter...")

# Subsample to keep file size manageable (20k points max per site).
set.seed(2025)
scatter_df <- fp_data[, .(manuscript_site, site_label, flag_status,
                          err_p3d, err_alt)]
scatter_df <- scatter_df[,
  .SD[sample(.N, min(.N, 20000L))],
  by = manuscript_site
]

p_fig3 <- ggplot(scatter_df, aes(x = err_p3d, y = err_alt)) +
  geom_hex(bins = 60) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey40") +
  facet_wrap(~ site_label, ncol = 5) +
  scale_fill_viridis_c(trans = "log10", name = "count") +
  coord_fixed(xlim = c(-30, 30), ylim = c(-30, 30)) +
  labs(
    x = "Original CHM error (err_p3d, m)",
    y = "Alt CHM error (err_alt, m)",
    title = "Figure 3: Footprint-level original vs alt CHM error",
    subtitle = "Points on the 1:1 line = DTM substitution had no effect at that footprint.\nPoints shifted toward y=0 = substitution reduced error magnitude. Clamped to ±30 m for display."
  ) +
  theme_cowplot(font_size = 10) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(size = 8, lineheight = 0.9),
    panel.spacing = unit(0.3, "lines")
  )

ggsave(file.path(GROUNDWORK_PLOT, "groundwork_task4_err_scatter_by_site.pdf"),
       p_fig3, width = 14, height = 10)
log_progress("  wrote groundwork_task4_err_scatter_by_site.pdf")

# ---- Figure 4: err_p3d vs err_dtm (does CNN DTM error drive CHM error?) ----

log_progress("Figure 4: err_p3d vs err_dtm scatter...")

scatter2_df <- fp_data[, .(manuscript_site, site_label, flag_status,
                           err_p3d, err_dtm)]
scatter2_df <- scatter2_df[,
  .SD[sample(.N, min(.N, 20000L))],
  by = manuscript_site
]

# Per-site Pearson r to annotate each panel.
corr_df <- site_summary[, .(
  manuscript_site,
  site_label = sprintf("Site %02d\n%s", manuscript_site, site_short_name),
  r = corr_errP3D_errDTM
)]
corr_df[, site_label := factor(site_label,
                               levels = levels(fp_data$site_label))]
corr_df[, label_txt := sprintf("r = %+.2f", r)]

p_fig4 <- ggplot(scatter2_df, aes(x = err_dtm, y = err_p3d)) +
  geom_hex(bins = 60) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey40") +
  geom_text(
    data = corr_df, inherit.aes = FALSE,
    aes(x = -Inf, y = Inf, label = label_txt),
    hjust = -0.1, vjust = 1.4, size = 3.2,
    fontface = "bold", colour = "black"
  ) +
  facet_wrap(~ site_label, ncol = 5) +
  scale_fill_viridis_c(trans = "log10", name = "count") +
  coord_fixed(xlim = c(-20, 20), ylim = c(-30, 30)) +
  labs(
    x = "DTM error (err_dtm = P3D DTM − 3DEP DTM, m)",
    y = "CHM error (err_p3d = P3D CHM − ALS CHM, m)",
    title = "Figure 4: Footprint-level CHM error vs DTM error",
    subtitle = "If CNN DTM error propagates through the DSM−DTM subtraction, we expect r > 0 (points along the dashed 1:1 line).\nIf the CHM error is DSM-driven, err_p3d is independent of err_dtm (r ≈ 0)."
  ) +
  theme_cowplot(font_size = 10) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(size = 8, lineheight = 0.9),
    panel.spacing = unit(0.3, "lines")
  )

ggsave(file.path(GROUNDWORK_PLOT, "groundwork_task4_err_vs_dtm_error.pdf"),
       p_fig4, width = 14, height = 10)
log_progress("  wrote groundwork_task4_err_vs_dtm_error.pdf")

# =====================================================================
# CONSOLE SUMMARY
# =====================================================================

cat("\n\n====================  SUMMARY  ====================\n\n")

cat("Per-site results (manuscript-numbered):\n\n")
print_cols <- c("manuscript_site", "site_short_name", "flag_status",
                "dominant_lc", "n",
                "p3d_rmse", "alt_rmse",
                "rmse_reduction_frac", "var_reduction_frac",
                "corr_errP3D_errDTM")
print(site_summary[, ..print_cols], row.names = FALSE, digits = 3)

cat("\n\nOverall / stratified results:\n\n")
print_cols2 <- c("stratum", "n", "p3d_rmse", "alt_rmse",
                 "rmse_reduction_frac", "var_reduction_frac",
                 "corr_errP3D_errDTM")
print(overall_summary[, ..print_cols2], row.names = FALSE, digits = 3)

# Quick qualitative verdict.
cat("\n\nInterpretation cheat-sheet:\n")
cat("  - var_reduction_frac > 0  ==> swapping in 3DEP DTM REDUCES CHM error variance.\n")
cat("    Magnitude = share of CHM variance attributable to CNN DTM error.\n")
cat("  - var_reduction_frac < 0  ==> swap INFLATES variance (err_dtm and err_p3d_dsm-only\n")
cat("    were negatively correlated; cancellation was happening).\n")
cat("  - corr_errP3D_errDTM      ==> footprint-level correlation of CHM and DTM errors.\n")
cat("    High positive r + large var_reduction_frac ==> CNN DTM drives CHM error.\n")
cat("    Low r                    ==> CHM error is DSM-driven.\n")

cat("\nArtifacts written:\n")
cat(sprintf("  %s\n", file.path(MS_TABLES_DIR, "coupling_footprint_counts.csv")))
cat(sprintf("  %s\n", file.path(MS_TABLES_DIR, "coupling_site_summary.csv")))
cat(sprintf("  %s\n", file.path(MS_TABLES_DIR, "coupling_overall_summary.csv")))
cat(sprintf("  %s\n", file.path(MS_TABLES_DIR, "groundwork_task4_footprint_data.rds")))
cat(sprintf("  %s (checkpoint)\n",
            file.path(Sys.getenv("PROJECT_ROOT"), "checkpoints",
                      "groundwork_task4_phase2.rds")))
cat(sprintf("  %s\n", file.path(GROUNDWORK_PLOT, "groundwork_task4_err_dist_by_site.pdf")))
cat(sprintf("  %s\n", file.path(GROUNDWORK_PLOT, "groundwork_task4_variance_reduction_bar.pdf")))
cat(sprintf("  %s\n", file.path(GROUNDWORK_PLOT, "groundwork_task4_err_scatter_by_site.pdf")))
cat(sprintf("  %s\n", file.path(GROUNDWORK_PLOT, "groundwork_task4_err_vs_dtm_error.pdf")))

cat("\n===========================================================\n\n")

log_progress("step 02 complete.")
