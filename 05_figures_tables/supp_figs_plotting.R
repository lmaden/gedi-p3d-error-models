#!/usr/bin/env Rscript
# =====================================================================
# supp_figs_plotting.R
#
# Renders the supplementary figures.
#
# CHANGES FROM v2:
#   S04 - drop 3 broken interaction/smooth panels (keep 17 mains);
#         add y-axis label; human-readable predictor labels
#   S06 - 2x2 facet layout (was 1x4); fixes annotation overlap;
#         caches subsample-with-fitted to avoid re-compute on re-runs
#   S08 - manuscript renumbering: tracker > 16 -> tracker - 1
#         (sites 17,18,19,20 -> manuscript 16,17,18,19)
#   S10 - manuscript renumbering (same as S08)
#   S17 - x-axis ordered by manuscript site number (not posterior mean);
#         bottom variance bar uses legend with embedded percentages
#         instead of in-bar text (fixes "lc_l1_cod..." truncation)
#   S18 - coord_cartesian xlim = c(-20, 20) for focus on density mass
#   S19 - MOST AGGRESSIVE: 2-panel layout. Top = main effects only,
#         bottom = land-cover interactions. Alphabetical sort within
#         each panel. Drops smooth basis function rows (sslope_*).
#         Uses geom_errorbar(orientation="y") to clear deprecation warn
#
# UNCHANGED: S12, S20, S21, S22
#
# Outputs: same 11 PNGs, overwriting v2 outputs.
# =====================================================================

suppressPackageStartupMessages({
  library(brms)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(scales)
})

Sys.setenv(DISPLAY = "")
options(device = pdf)

if (!exists("log_progress"))  source("analysis_utils.R")
if (!exists("PROJECT_ROOT"))  source("analysis_config.R")

manuscript_tables_dir <- file.path(PROJECT_ROOT, "manuscript_tables")
plots_dir <- file.path(PROJECT_ROOT, "plots", "section_L_supp_v05")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

COL_CHM  <- "#2D5F3F"
COL_DTM  <- "#C7522A"
COL_POOL <- "#555555"
COL_REF  <- "#A02020"
COL_EREG <- "#D8B044"
COL_LC   <- "#88B0D8"
DPI      <- 300
W_FULL   <- 6.5
H_HALF   <- 4.0
H_FULL   <- 6.5

theme_supp <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
          strip.background = element_rect(fill = "grey95", color = NA),
          strip.text       = element_text(size = base_size - 1, face = "bold"),
          plot.title       = element_text(size = base_size + 1, face = "bold"),
          legend.title     = element_text(size = base_size - 1),
          legend.text      = element_text(size = base_size - 1))
}

# Manuscript renumbering helper
# Tracker sites 17, 18, 19, 20 become manuscript sites 16, 17, 18, 19
# (site 16 was permanently omitted; canonical renumbering shifts everything
# after the gap)
to_ms_site <- function(tracker_id) {
  t <- suppressWarnings(as.integer(as.character(tracker_id)))
  ifelse(is.na(t), as.character(tracker_id),
         as.character(ifelse(t > 16, t - 1, t)))
}

# Human-readable predictor label map for S4
PREDICTOR_LABELS <- c(
  "slope_mean_z"   = "Slope mean (z)",
  "slope_sd_z"     = "Slope SD (z)",
  "wsci_z"         = "WSCI (z)",
  "rh_98_z"        = "RH98 (z)",
  "cover_z"        = "Canopy cover (z)",
  "aspect_sin_z"   = "Aspect sin (z)",
  "aspect_cos_z"   = "Aspect cos (z)",
  "meta_offnad_z"  = "Off-nadir angle (z)",
  "meta_sunel_z"   = "Sun elevation (z)",
  "meta_az_conc_z" = "Azimuth concentration (z)",
  "meta_stereo_z"  = "Stereo ratio (z)",
  "meta_leafon_z"  = "Leaf-on fraction (z)",
  "meta_fwdrev_z"  = "Fwd/rev ratio (z)",
  "meta_absgeo_z"  = "Abs geo accuracy (z)",
  "meta_relgeo_z"  = "Rel geo accuracy (z)",
  "view_az_sin_z"  = "View azimuth sin (z)",
  "view_az_cos_z"  = "View azimuth cos (z)"
)
prettify_pred <- function(x) {
  ifelse(x %in% names(PREDICTOR_LABELS), PREDICTOR_LABELS[x], x)
}

log_progress("================================================================")
log_progress("Re-rendering 11 figures")
log_progress(sprintf("Output directory: %s", plots_dir))
log_progress("================================================================")

# ---------------------------------------------------------------------
# Pre-flight (same as v2)
# ---------------------------------------------------------------------

log_subsection("Pre-flight: testing PDF save + gs PNG conversion")
GS_BIN <- Sys.which("gs")
if (!nzchar(GS_BIN)) stop("Ghostscript (gs) not found")
log_progress(sprintf("  gs: %s", GS_BIN))

test_path_pdf <- file.path(plots_dir, ".preflight_test.pdf")
test_path_png <- file.path(plots_dir, ".preflight_test.png")
ok_pdf <- tryCatch({
  p_test <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) +
    geom_point(color = COL_CHM, size = 3) +
    labs(title = "preflight") + theme_supp()
  ggsave(test_path_pdf, p_test, width = 3, height = 2, device = "pdf")
  file.exists(test_path_pdf) && file.info(test_path_pdf)$size > 200
}, error = function(e) FALSE)
if (!ok_pdf) stop("PDF save failed")
gs_cmd_test <- sprintf(
  "%s -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r%d -sOutputFile=%s %s",
  GS_BIN, DPI, shQuote(test_path_png), shQuote(test_path_pdf))
ok_gs <- system(gs_cmd_test, ignore.stdout = TRUE, ignore.stderr = TRUE) == 0 &&
         file.exists(test_path_png) && file.info(test_path_png)$size > 500
if (!ok_gs) stop("gs PNG conversion failed")
unlink(c(test_path_pdf, test_path_png))
log_progress("  Pre-flight passed.")

# ---------------------------------------------------------------------
# Save helpers (same as v2)
# ---------------------------------------------------------------------

save_fig <- function(p, tag, width = W_FULL, height = H_HALF) {
  pdf_path <- file.path(plots_dir, sprintf("%s.pdf", tag))
  png_path <- file.path(plots_dir, sprintf("%s.png", tag))
  ggsave(pdf_path, p, width = width, height = height,
         device = "pdf", bg = "white", units = "in")
  if (!file.exists(pdf_path) || file.info(pdf_path)$size < 200) return(FALSE)
  gs_cmd <- sprintf(
    "%s -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r%d -sOutputFile=%s %s",
    GS_BIN, DPI, shQuote(png_path), shQuote(pdf_path))
  ret <- system(gs_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  ret == 0 && file.exists(png_path) && file.info(png_path)$size > 500
}

render_fig <- function(tag, expr) {
  log_subsection(sprintf("Rendering %s", tag))
  tryCatch({
    eval(expr, envir = parent.frame())
    png_path <- file.path(plots_dir, sprintf("%s.png", tag))
    if (file.exists(png_path)) {
      log_progress(sprintf("  %s OK (%d bytes)", tag,
                           file.info(png_path)$size))
    } else {
      log_progress(sprintf("  %s: PNG missing post-render", tag))
    }
  }, error = function(e) {
    log_progress(sprintf("  ERROR in %s: %s", tag, e$message))
  })
}

# ---------------------------------------------------------------------
# Load data (same as v2)
# ---------------------------------------------------------------------

log_subsection("Loading data CSVs")
read_if <- function(fn) {
  path <- file.path(manuscript_tables_dir, fn)
  if (file.exists(path)) { log_progress(sprintf("  + %s", fn)); fread(path) }
  else { log_progress(sprintf("  - MISSING: %s", fn)); NULL }
}

chm_re        <- read_if("section_L_s17_chm_site_re_summary.csv")
dtm_re        <- read_if("section_L_s21_dtm_site_re_summary.csv")
chm_var_part  <- read_if("section_L_s22_chm_variance_partition.csv")
chm_resid     <- read_if("section_L_s20_chm_residuals.csv")
dtm_resid     <- read_if("section_L_s20_dtm_residuals.csv")
chm_var_emp   <- read_if("section_L_s20_chm_empirical_variogram.csv")
dtm_var_emp   <- read_if("section_L_s20_dtm_empirical_variogram.csv")
chm_var_par   <- read_if("section_L_s20_chm_variogram_params.csv")
dtm_var_par   <- read_if("section_L_s20_dtm_variogram_params.csv")
chm_ppc       <- read_if("section_L_s18_chm_ppc_draws.csv")
dtm_ppc       <- read_if("section_L_s18_dtm_ppc_draws.csv")
chm_mcmc      <- read_if("section_L_s12_chm_mcmc_diagnostics.csv")
dtm_mcmc      <- read_if("section_L_s12_dtm_mcmc_diagnostics.csv")
chm_biv       <- read_if("section_L_s21_chm_full17_predictors.csv")
dtm_biv       <- read_if("section_L_s21_dtm_full17_predictors.csv")
chm_fixef_cmp <- read_if("chm_sensitivity_fixef_comparison.csv")

# =====================================================================
# S4 - PATCHED: drop interactions/smooths, add y-axis label, pretty preds
# =====================================================================

render_fig("figS04", quote({
  ce_path <- file.path(manuscript_tables_dir,
                       "section_L_s4_chm_conditional_effects.rds")
  if (!file.exists(ce_path)) stop("conditional_effects rds not found")
  ce <- readRDS(ce_path)
  log_progress(sprintf("  Total effects in rds: %d", length(ce)))
  log_progress(sprintf("  Effect names: %s",
                       paste(names(ce), collapse = ", ")))

  # Filter: drop any whose name contains ":" (interactions and 2D smooths)
  keep_idx <- !grepl(":", names(ce))
  ce <- ce[keep_idx]
  log_progress(sprintf("  Keeping %d main-effect panels", length(ce)))

  plots <- list()
  for (i in seq_along(ce)) {
    df <- ce[[i]]
    if (!is.data.frame(df)) next
    xvar <- setdiff(names(df), c("estimate__", "se__", "lower__", "upper__",
                                 "cond__", "effect1__", "effect2__"))[1]
    if (is.na(xvar)) next
    p <- ggplot(df, aes(x = .data[[xvar]], y = estimate__)) +
      geom_ribbon(aes(ymin = lower__, ymax = upper__),
                  fill = COL_CHM, alpha = 0.25) +
      geom_line(color = COL_CHM, linewidth = 0.5) +
      labs(x = prettify_pred(xvar),
           y = "CHM error (m)") +
      theme_supp(7)
    plots[[length(plots) + 1]] <- p
  }
  if (length(plots) == 0) stop("no plottable panels")
  # 17 panels -> 5x4 grid (one cell empty bottom-right)
  save_fig(wrap_plots(plots, ncol = 4), "figS04",
           width = W_FULL, height = 9.0)
}))

# =====================================================================
# S6 - PATCHED: 2x2 facet layout, caches subsample for re-runs
# =====================================================================

render_fig("figS06", quote({
  s6_subs_cache <- file.path(manuscript_tables_dir,
                            "section_L_s6_subsample_with_fitted.csv")
  s6_r2_cache <- file.path(manuscript_tables_dir,
                          "section_L_s6_per_class_r2.csv")

  if (file.exists(s6_subs_cache) && file.exists(s6_r2_cache)) {
    log_progress("  Using cached S6 subsample-with-fitted CSV")
    subs <- fread(s6_subs_cache)
    r2_by_class <- fread(s6_r2_cache)
  } else {
    log_progress("  S6 cache miss - computing fitted() (slow)")
    if (!exists("cp_chm_l4")) {
      cp_chm_l4 <<- load_checkpoint("10b_chm_sensitivity")
    }
    fit_chm <- cp_chm_l4$fit_chm_16
    train <- as.data.table(fit_chm$data)
    forest_types <- c("BDF", "DNF", "EBF", "ENF")
    train_forest <- train[lc_l1_code %in% forest_types]

    set.seed(20260512)
    per_class_n <- min(12500L, floor(nrow(train_forest) /
                                     length(forest_types)))
    subs <- train_forest[, .SD[sample(.N, min(.N, per_class_n))],
                         by = lc_l1_code]

    log_progress("  Computing fitted() on full forest training data")
    fitted_full <- fitted(fit_chm, newdata = train_forest, summary = TRUE)
    train_forest[, pred_mean := fitted_full[, "Estimate"]]
    r2_by_class <- train_forest[, .(
      n  = .N,
      r2 = 1 - sum((chm_error_mean - pred_mean)^2) /
               sum((chm_error_mean - mean(chm_error_mean))^2)
    ), by = lc_l1_code]
    fwrite(r2_by_class, s6_r2_cache)

    log_progress("  Computing fitted() on subsample")
    subs_pred <- fitted(fit_chm, newdata = subs, summary = TRUE)
    subs[, pred_mean := subs_pred[, "Estimate"]]
    fwrite(subs, s6_subs_cache)
  }

  log_progress("  Per-class R^2 anchors:")
  for (rr in seq_len(nrow(r2_by_class))) {
    log_progress(sprintf("    %s: R^2 = %.3f, n = %d",
                         r2_by_class$lc_l1_code[rr],
                         r2_by_class$r2[rr], r2_by_class$n[rr]))
  }

  axis_lim <- range(c(subs$chm_error_mean, subs$pred_mean), na.rm = TRUE)
  axis_lim <- c(floor(axis_lim[1]), ceiling(axis_lim[2]))
  label_df <- r2_by_class[, .(
    lc_l1_code,
    lab = sprintf("R^2 = %.3f\nn = %s", r2, format(n, big.mark = ","))
  )]
  label_df[, x := axis_lim[1] + 0.05 * diff(axis_lim)]
  label_df[, y := axis_lim[2] - 0.10 * diff(axis_lim)]

  p <- ggplot(subs, aes(x = pred_mean, y = chm_error_mean)) +
    geom_hex(bins = 40) +
    scale_fill_viridis_c(trans = "log10", name = "Count") +
    geom_abline(slope = 1, intercept = 0, color = COL_REF,
                linetype = "dashed", linewidth = 0.5) +
    geom_text(data = label_df, aes(x = x, y = y, label = lab),
              hjust = 0, vjust = 1, size = 3, inherit.aes = FALSE) +
    # 2x2 layout
    facet_wrap(~ lc_l1_code, nrow = 2, ncol = 2) +
    coord_fixed(xlim = axis_lim, ylim = axis_lim) +
    labs(x = "Predicted CHM error (m)", y = "Observed CHM error (m)") +
    theme_supp()
  save_fig(p, "figS06", width = W_FULL, height = 6.5)
}))

# =====================================================================
# S8 - PATCHED: manuscript renumbering
# =====================================================================

render_fig("figS08", quote({
  if (is.null(chm_resid)) stop("S20 CHM residuals not loaded")
  chm_resid_ms <- copy(chm_resid)
  chm_resid_ms[, ms_site := to_ms_site(site)]
  chm_resid_ms[, ms_site := factor(ms_site,
                                   levels = sort(unique(suppressWarnings(
                                     as.integer(ms_site)))))]
  lim <- max(abs(quantile(chm_resid_ms$resid_mean, c(0.02, 0.98),
                          na.rm = TRUE)))
  p <- ggplot(chm_resid_ms, aes(x = x_alb, y = y_alb,
                                color = resid_mean)) +
    geom_point(size = 0.2, alpha = 0.6) +
    scale_color_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                          midpoint = 0, limits = c(-lim, lim),
                          oob = squish, name = "Residual (m)") +
    facet_wrap(~ ms_site, scales = "free", ncol = 4,
              labeller = labeller(ms_site = function(x) paste("Site", x))) +
    labs(x = "Easting (m, EPSG:5070)", y = "Northing (m, EPSG:5070)") +
    theme_supp(8) +
    theme(axis.text = element_text(size = 5), legend.position = "bottom")
  save_fig(p, "figS08", width = W_FULL, height = 8.5)
}))

# =====================================================================
# S10 - PATCHED: manuscript renumbering
# =====================================================================

render_fig("figS10", quote({
  if (is.null(chm_resid)) stop("S20 CHM residuals not loaded")
  chm_resid_ms <- copy(chm_resid)
  chm_resid_ms[, ms_site := to_ms_site(site)]
  chm_resid_ms[, ms_site := factor(ms_site,
                                   levels = sort(unique(suppressWarnings(
                                     as.integer(ms_site)))))]
  p <- ggplot(chm_resid_ms, aes(x = resid_mean)) +
    geom_histogram(bins = 50, fill = COL_CHM, color = NA, alpha = 0.8) +
    geom_vline(xintercept = 0, color = COL_REF,
               linetype = "dashed", linewidth = 0.4) +
    facet_wrap(~ ms_site, scales = "free_y", ncol = 4,
              labeller = labeller(ms_site = function(x) paste("Site", x))) +
    labs(x = "CHM model residual (m)", y = "Count") + theme_supp(8)
  save_fig(p, "figS10", width = W_FULL, height = 7.5)
}))

# =====================================================================
# S12 - UNCHANGED
# =====================================================================

render_fig("figS12", quote({
  if (is.null(chm_mcmc) || is.null(dtm_mcmc)) stop("MCMC not loaded")
  chm_mcmc[, product := "CHM (16-site, n=187)"]
  dtm_mcmc[, product := "DTM (18-site, n=191)"]
  mcmc_all <- rbind(chm_mcmc, dtm_mcmc, fill = TRUE)
  p_rhat <- ggplot(mcmc_all, aes(x = rhat, fill = product)) +
    geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = c(COL_CHM, COL_DTM)) +
    geom_vline(xintercept = 1.01, color = COL_REF,
               linetype = "dashed", linewidth = 0.4) +
    facet_wrap(~ product, ncol = 1, scales = "free_y") +
    labs(x = expression(hat(R)), y = "Count") +
    theme_supp() + theme(legend.position = "none")
  p_ess <- ggplot(mcmc_all, aes(x = neff_ratio, fill = product)) +
    geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = c(COL_CHM, COL_DTM)) +
    geom_vline(xintercept = 0.1, color = COL_REF,
               linetype = "dashed", linewidth = 0.4) +
    facet_wrap(~ product, ncol = 1, scales = "free_y") +
    labs(x = "Effective sample size ratio", y = "Count") +
    theme_supp() + theme(legend.position = "none")
  save_fig(p_rhat | p_ess, "figS12", width = W_FULL, height = 5.0)
}))

# =====================================================================
# S17 - PATCHED: manuscript site x-axis order, legend instead of in-bar
# =====================================================================

render_fig("figS17", quote({
  if (is.null(chm_re)) stop("CHM site REs not loaded")
  if (is.null(chm_var_part)) stop("variance partition not loaded")

  # Top panel: order by manuscript site number ascending
  re_df <- copy(chm_re)
  re_df[, ms_site := to_ms_site(site_label)]
  re_df[, ms_site_num := suppressWarnings(as.integer(ms_site))]
  setorder(re_df, ms_site_num)
  re_df[, site_lab := factor(ms_site, levels = ms_site)]

  p_top <- ggplot(re_df, aes(x = site_lab, y = re_mean)) +
    geom_errorbar(aes(ymin = re_q025, ymax = re_q975),
                  width = 0, color = COL_POOL, linewidth = 0.4) +
    geom_point(color = COL_CHM, size = 2) +
    geom_hline(yintercept = 0, color = COL_REF,
               linetype = "dashed", linewidth = 0.4) +
    labs(x = "Site (manuscript numbering)",
         y = "Site random intercept (m)") +
    theme_supp() +
    theme(axis.text.x = element_text(size = 8))

  # Bottom panel: variance partition with legend (no in-bar text)
  vp <- chm_var_part[, .(grouping, mean = share_mean)]
  vp <- vp[order(-mean)]
  nice <- c("site" = "Site", "ecoregion" = "Ecoregion",
            "lc_l1_code" = "Land cover")
  vp[, nice_label := nice[grouping]]
  vp[, legend_label := sprintf("%s (%.1f%%)", nice_label, mean * 100)]
  vp[, legend_label := factor(legend_label, levels = legend_label)]
  vp[, x := 1]

  color_map <- setNames(c(COL_CHM, COL_EREG, COL_LC),
                       vp$legend_label)

  p_bot <- ggplot(vp, aes(x = x, y = mean, fill = legend_label)) +
    geom_col(width = 0.4, color = "white", linewidth = 0.5) +
    scale_fill_manual(values = color_map, name = NULL) +
    coord_flip() +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = c(0, 0)) +
    labs(x = NULL, y = "Share of total variance") +
    theme_supp() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 9))

  save_fig(p_top / p_bot + plot_layout(heights = c(3, 1.2)),
           "figS17", width = W_FULL, height = 6.5)
}))

# =====================================================================
# S18 - PATCHED: zoom to [-20, 20] for both panels
# =====================================================================

render_fig("figS18", quote({
  if (is.null(chm_ppc) || is.null(dtm_ppc)) stop("PPC not loaded")
  chm_ppc[, product := "CHM (16-site)"]
  dtm_ppc[, product := "DTM (18-site)"]
  plot_ppc <- function(d, color_main) {
    draws_to_plot <- sample(unique(d$draw),
                           min(30, length(unique(d$draw))))
    d_sub <- d[draw %in% draws_to_plot]
    obs_density <- d[draw == draw[1], .(observed)]
    ggplot() +
      geom_density(data = d_sub, aes(x = yrep, group = draw),
                   color = alpha(color_main, 0.15), linewidth = 0.3) +
      geom_density(data = obs_density, aes(x = observed),
                   color = "black", linewidth = 0.8) +
      coord_cartesian(xlim = c(-20, 20)) +
      labs(x = "Error (m)", y = "Density", title = unique(d$product)) +
      theme_supp() +
      theme(plot.title = element_text(hjust = 0.5))
  }
  save_fig(plot_ppc(chm_ppc, COL_CHM) / plot_ppc(dtm_ppc, COL_DTM),
           "figS18", width = W_FULL, height = 5.5)
}))

# =====================================================================
# S19 - PATCHED: 2-panel split (main effects + interactions),
#                drop smooth basis rows, alphabetical sort,
#                geom_errorbar(orientation="y") to clear warning
# =====================================================================

render_fig("figS19", quote({
  if (is.null(chm_fixef_cmp)) stop("fixef cmp not loaded")
  df <- copy(chm_fixef_cmp)

  # Categorize predictors
  df[, is_smooth      := grepl("^sslope_", predictor)]
  df[, is_interaction := grepl(":", predictor) & !is_smooth]
  df[, is_main        := !is_smooth & !is_interaction]

  n_smooth <- sum(df$is_smooth)
  n_interaction <- sum(df$is_interaction)
  n_main <- sum(df$is_main)
  log_progress(sprintf("  Categorized: %d main, %d interactions, %d smooths (dropped)",
                       n_main, n_interaction, n_smooth))

  # Drop smooth basis function rows
  df <- df[!is_smooth]

  build_long <- function(d) {
    rbind(
      d[, .(predictor, estimate = est_full,
            lower = lo_full, upper = hi_full,
            frame = "19-site (full)")],
      d[, .(predictor, estimate = est_16,
            lower = lo_16, upper = hi_16,
            frame = "16-site (primary)")]
    )
  }

  # Shared x-axis limits across both panels
  all_long <- build_long(df)
  x_lims <- range(c(all_long$lower, all_long$upper), na.rm = TRUE)

  build_panel <- function(d_subset, title) {
    plot_df <- build_long(d_subset)
    plot_df[, predictor := factor(predictor,
                                  levels = sort(unique(predictor),
                                                decreasing = TRUE))]
    ggplot(plot_df, aes(x = estimate, y = predictor, color = frame)) +
      geom_vline(xintercept = 0, color = COL_REF,
                 linetype = "dashed", linewidth = 0.4) +
      geom_errorbar(aes(xmin = lower, xmax = upper),
                    width = 0, orientation = "y",
                    position = position_dodge(width = 0.5),
                    linewidth = 0.5) +
      geom_point(position = position_dodge(width = 0.5), size = 2) +
      scale_color_manual(values = c("19-site (full)" = COL_POOL,
                                    "16-site (primary)" = COL_CHM)) +
      coord_cartesian(xlim = x_lims) +
      labs(x = NULL, y = NULL, color = "Model frame", title = title) +
      theme_supp() +
      theme(plot.title = element_text(hjust = 0, size = 10))
  }

  p_main <- build_panel(df[is_main],
                        sprintf("Main effects (n = %d)", n_main))
  p_int  <- build_panel(df[is_interaction],
                        sprintf("Land-cover interactions (n = %d)",
                                n_interaction))

  # Share legend by collecting from one
  p_combined <- (p_main / p_int) +
    plot_layout(guides = "collect", heights = c(1, 1.6)) +
    plot_annotation(caption = "Standardized coefficient (95% CI)") &
    theme(legend.position = "right")

  save_fig(p_combined, "figS19", width = W_FULL, height = 9.0)
}))

# =====================================================================
# S20, S21, S22 - UNCHANGED from v2
# =====================================================================

render_fig("figS20", quote({
  if (is.null(chm_var_emp) || is.null(dtm_var_emp))
    stop("variogram CSVs not loaded")
  build_var_panel <- function(emp, par, color_main, title) {
    dist_seq <- seq(0, max(emp$dist), length.out = 200)
    nugget <- par$nugget; psill <- par$partial_sill; rng <- par$range_m
    sph_model <- ifelse(dist_seq <= rng,
                       nugget + psill * (1.5 * (dist_seq / rng) -
                                        0.5 * (dist_seq / rng)^3),
                       nugget + psill)
    model_df <- data.table(dist = dist_seq, gamma = sph_model)
    label <- sprintf("Nugget/sill = %.2f\nRange = %.2f km",
                     par$nugget_to_sill, par$range_km)
    ggplot(emp, aes(x = dist / 1000, y = gamma)) +
      geom_point(color = color_main, size = 1.5) +
      geom_line(data = model_df, aes(x = dist / 1000, y = gamma),
                color = color_main, linewidth = 0.6) +
      annotate("text", x = Inf, y = -Inf, label = label,
               hjust = 1.1, vjust = -0.5, size = 3.0) +
      labs(x = "Distance (km)", y = expression(gamma~"(m"^2*")"),
           title = title) +
      theme_supp() +
      theme(plot.title = element_text(hjust = 0.5))
  }
  save_fig(build_var_panel(chm_var_emp, chm_var_par, COL_CHM, "CHM (16-site)") |
           build_var_panel(dtm_var_emp, dtm_var_par, COL_DTM, "DTM (18-site)"),
           "figS20", width = W_FULL, height = 3.5)
}))

render_fig("figS21", quote({
  if (is.null(chm_biv) || is.null(dtm_biv)) stop("S21 CSVs not loaded")
  build_biv_panel <- function(biv, color_main, title, top_n = 17) {
    d <- biv[order(-r_squared)][seq_len(min(top_n, .N))]
    d[, predictor := factor(predictor, levels = rev(predictor))]
    d[, sig := ifelse(bonferroni_significant, "*", "")]
    ggplot(d, aes(x = r_squared, y = predictor)) +
      geom_col(fill = color_main, alpha = 0.8, width = 0.7) +
      geom_text(aes(label = sprintf("%.3f %s", r_squared, sig),
                    x = r_squared),
                hjust = -0.15, size = 2.8) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.25)),
                         limits = c(0, max(d$r_squared) * 1.25)) +
      labs(x = expression("Bivariate r"^2), y = NULL, title = title) +
      theme_supp() +
      theme(plot.title = element_text(hjust = 0.5))
  }
  save_fig(build_biv_panel(chm_biv, COL_CHM, "CHM (16-site)") |
           build_biv_panel(dtm_biv, COL_DTM, "DTM (18-site)"),
           "figS21", width = W_FULL, height = 5.5)
}))

render_fig("figS22", quote({
  if (is.null(chm_var_part)) stop("variance partition not loaded")
  vp <- chm_var_part[, .(grouping, mean = share_mean,
                         q025 = share_q025, q975 = share_q975)]
  vp[, grouping := factor(grouping,
                          levels = c("lc_l1_code", "ecoregion", "site"),
                          labels = c("Land cover (l1)", "Ecoregion", "Site"))]
  vp[, label := sprintf("%.1f%% [%.1f%%, %.1f%%]",
                       mean * 100, q025 * 100, q975 * 100)]
  p <- ggplot(vp, aes(x = grouping, y = mean, fill = grouping)) +
    geom_col(width = 0.55, alpha = 0.85) +
    geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.15,
                  linewidth = 0.5) +
    geom_text(aes(label = label, y = q975), vjust = -0.6, size = 3.0) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = c("Land cover (l1)" = COL_LC,
                                 "Ecoregion" = COL_EREG,
                                 "Site" = COL_CHM)) +
    labs(x = NULL, y = "Share of explained variance (95% CI)") +
    theme_supp() + theme(legend.position = "none")
  save_fig(p, "figS22", width = W_FULL, height = 4.0)
}))

# ---------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------

log_progress("================================================================")
log_progress("SUMMARY")
log_progress("================================================================")
expected <- sprintf("figS%02d", c(4, 6, 8, 10, 12, 17, 18, 19, 20, 21, 22))
n_ok <- 0
for (tag in expected) {
  png_p <- file.path(plots_dir, sprintf("%s.png", tag))
  if (file.exists(png_p)) {
    log_progress(sprintf("  %s.png  OK (%s)", tag,
                         format(file.info(png_p)$size, big.mark = ",")))
    n_ok <- n_ok + 1
  } else {
    log_progress(sprintf("  %s.png  MISSING", tag))
  }
}
log_progress("")
log_progress(sprintf("%d / %d figures saved", n_ok, length(expected)))
log_progress("================================================================")
