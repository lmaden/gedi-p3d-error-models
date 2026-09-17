#!/usr/bin/env Rscript
# =====================================================================
# supp_figs_rerender.R
#
# PURPOSE
#   Regenerates four supplementary figures:
#
#     S6  — CHM Predicted vs Observed Error by Forest Type (4-panel hex)
#           wrong R^2 formula (1 - SS_res/SS_tot instead
#           of cor(pred,obs)^2), flipped axes, missing RMSE and title.
#
#     S12 — MCMC Convergence Diagnostics (6-panel: R-hat, ESS, by type,
#           trace plots). only produced 2-panel histogram
#           comparison; missing 4 of v04's 6 panels entirely.
#
#     S18 — Posterior Predictive Checks (2-panel CHM/DTM density overlay
#           with 95% PI coverage annotations). missing
#           title, missing 95% PI coverage subtitles.
#
#     S20 — Semi-variograms of CHM/DTM model residuals (2-panel side-by-
#           side with fitted Sph/Exp/Gau/Mat model selection, GEDI swath
#           reference line, sill/range markers). used
#           20 km cutoff instead of 50 km; binning differed; no model
#           selection. Numerical results didn't match v04's DTM panel
#           even though DTM model is unchanged.
#
# SCRIPT SOURCES (which the four blocks below mirror)
#   S6  : publication_figures.R lines 149-340
#   S12 : generate_convergence_diagnostics.R (whole)
#   S18 : generate_ppc_overall.R (whole)
#   S20 : section_12_spatial.R lines 463-700, with DTM-side block added
#
# MODELS
#   CHM: 16-site refit from checkpoints/10b_chm_sensitivity.rds
#        (stored as `fit_chm_16` — sites 1, 2, 3 excluded for reference-
#        CHM underestimation per §10b script comment)
#   DTM: unchanged 18-site Stage 2 from checkpoints/10_models_stage2.rds
#        (stored as `fit_dtm_s2`)
#
# OUTPUTS
#   plots/groundwork/figS06_pred_vs_obs_by_forest_16site.pdf
#   plots/groundwork/figS12_convergence_16site.pdf
#   plots/groundwork/figS18_ppc_overall_16site.pdf
#   plots/groundwork/figS20_variogram_2panel_16site.pdf
#   manuscript_tables/section_L5_figS06_r2_anchors.csv
#   manuscript_tables/section_L5_figS18_coverage.csv
#   manuscript_tables/section_L5_figS20_variogram_params.csv
#
# RUN: from /gpfs/data1/vclgp/lmaden/chpt1, with PROJECT_ROOT exported
#   Rscript supp_figs_rerender.R
#
# =====================================================================

Sys.setenv(DISPLAY = "")
options(device = pdf)

suppressPackageStartupMessages({
  library(brms)
  library(bayesplot)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(scales)
  library(viridis)
  library(data.table)
  library(sf)
  library(gstat)
})

bayesplot::bayesplot_theme_set(bayesplot::theme_default())

# ---------------------------------------------------------------------
# Source project infrastructure
# ---------------------------------------------------------------------

if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("PROJECT_ROOT")) source("analysis_config.R")

stopifnot(exists("PROJECT_ROOT"), exists("CHECKPOINT_DIR"))

out_dir <- file.path(PROJECT_ROOT, "plots", "groundwork")
manuscript_tables_dir <- file.path(PROJECT_ROOT, "manuscript_tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(manuscript_tables_dir, showWarnings = FALSE, recursive = TRUE)

log_progress(strrep("=", 70))
log_progress("Regenerating figS06, figS12, figS18, figS20")
log_progress(strrep("=", 70))

# ---------------------------------------------------------------------
# Load checkpoints
# ---------------------------------------------------------------------

log_progress("Loading 16-site CHM refit (10b_chm_sensitivity)...")
cp_chm <- load_checkpoint("10b_chm_sensitivity")
stopifnot("fit_chm_16" %in% names(cp_chm))
fit_chm <- cp_chm$fit_chm_16

log_progress(sprintf("  CHM model: %d obs across %d sites",
                     nrow(fit_chm$data),
                     length(unique(as.character(fit_chm$data$site)))))

log_progress("Loading 18-site DTM (10_models_stage2)...")
cp_s2 <- load_checkpoint("10_models_stage2")
stopifnot("fit_dtm_s2" %in% names(cp_s2))
fit_dtm <- cp_s2$fit_dtm_s2

log_progress(sprintf("  DTM model: %d obs across %d sites",
                     nrow(fit_dtm$data),
                     length(unique(as.character(fit_dtm$data$site)))))

# ---------------------------------------------------------------------
# Constants used across figures
# ---------------------------------------------------------------------

FOREST_CLASSES <- c("BDF", "DNF", "EBF", "ENF")
FOREST_LABELS <- c(
  "BDF" = "Broadleaf Deciduous",
  "DNF" = "Deciduous Needleleaf",
  "EBF" = "Evergreen Broadleaf",
  "ENF" = "Evergreen Needleleaf"
)

# Publication theme (matches section_15 / generate_convergence)
theme_pub <- function(base_size = 10) {
  theme_cowplot(font_size = base_size) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.line  = element_line(color = "black", linewidth = 0.5),
      axis.text  = element_text(color = "black", size = base_size - 1),
      axis.title = element_text(color = "black", size = base_size),
      plot.title    = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0, color = "gray40")
    )
}

# =====================================================================
# figS06 — CHM Predicted vs Observed Error by Forest Type
# Adapted from publication_figures.R lines 149-340
# Two differences from publication_figures.R:
#   (1) R^2 = cor(pred, obs)^2, not 1 - SS_res/SS_tot. The sum-of-squares R^2 is also written to the
#       anchors CSV, not displayed on the figure.
#   (2) Axes: observed on x, predicted on y.
# =====================================================================

log_progress("")
log_progress(strrep("-", 70))
log_progress("figS06 — Predicted vs Observed by Forest Type")
log_progress(strrep("-", 70))

# Get fitted values on the model's training data restricted to forest
chm_train <- as.data.table(fit_chm$data)
chm_forest <- chm_train[lc_l1_code %in% FOREST_CLASSES]
log_progress(sprintf("  Forest-class training rows: %s",
                     format(nrow(chm_forest), big.mark = ",")))

log_progress("  Computing fitted() on full forest training data...")
# fitted() defaults to re_formula = NULL = include all random effects
fitted_full <- fitted(fit_chm, newdata = chm_forest, summary = TRUE, ndraws = 100)
chm_forest[, predicted := fitted_full[, "Estimate"]]
chm_forest[, observed  := chm_error_mean]

# Per-forest-type anchors: BOTH R^2 definitions for transparency
r2_anchors <- chm_forest[, .(
  n       = .N,
  r2_cor  = cor(observed, predicted, use = "complete.obs")^2,
  r2_ss   = 1 - sum((observed - predicted)^2, na.rm = TRUE) /
                sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE),
  rmse    = sqrt(mean((observed - predicted)^2, na.rm = TRUE)),
  mae     = mean(abs(observed - predicted), na.rm = TRUE),
  bias    = mean(predicted - observed, na.rm = TRUE)
), by = lc_l1_code]
setorder(r2_anchors, lc_l1_code)

log_progress("  Per-forest-type anchors:")
for (rr in seq_len(nrow(r2_anchors))) {
  log_progress(sprintf(
    "    %s: cor^2 = %.3f | SS R^2 = %.3f | RMSE = %.2f m | n = %s",
    r2_anchors$lc_l1_code[rr],
    r2_anchors$r2_cor[rr],
    r2_anchors$r2_ss[rr],
    r2_anchors$rmse[rr],
    format(r2_anchors$n[rr], big.mark = ",")))
}

fwrite(r2_anchors,
       file.path(manuscript_tables_dir, "section_L5_figS06_r2_anchors.csv"))

# Per-panel label using v04's format: R² = X.XX / RMSE = X.X m / n = X,XXX
label_df <- r2_anchors[, .(
  lc_l1_code,
  lab = sprintf("R\u00b2 = %.2f\nRMSE = %.1f m\nn = %s",
                r2_cor, rmse, format(n, big.mark = ","))
)]

# Symmetric axis limits at 99.5%ile of |obs ∪ pred|
lim_val <- quantile(abs(c(chm_forest$observed, chm_forest$predicted)),
                    0.995, na.rm = TRUE)
lim_val <- ceiling(lim_val)
log_progress(sprintf("  Axis limits: [-%d, %d] m", lim_val, lim_val))

label_df[, x := -lim_val * 0.9]
label_df[, y :=  lim_val * 0.85]

p_s06 <- ggplot(chm_forest, aes(x = observed, y = predicted)) +
  geom_hex(bins = 50, aes(fill = after_stat(count))) +
  scale_fill_viridis_c(trans = "log10", name = "Count",
                       labels = scales::label_comma()) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "red", linewidth = 0.6) +
  geom_text(data = label_df,
            aes(x = x, y = y, label = lab),
            hjust = 0, vjust = 1, size = 3, inherit.aes = FALSE) +
  facet_wrap(~ lc_l1_code, ncol = 2,
             labeller = labeller(lc_l1_code = FOREST_LABELS)) +
  coord_fixed(xlim = c(-lim_val, lim_val), ylim = c(-lim_val, lim_val)) +
  labs(
    title = "CHM Predicted vs Observed Error by Forest Type",
    x     = "Observed error (m)",
    y     = "Predicted error (m)"
  ) +
  theme_pub() +
  theme(legend.position = "right")

ggsave(file.path(out_dir, "figS06_pred_vs_obs_by_forest_16site.pdf"),
       p_s06, width = 10, height = 10, bg = "white", device = pdf)
log_progress(sprintf("  Saved: %s/figS06_pred_vs_obs_by_forest_16site.pdf", out_dir))

# =====================================================================
# figS12 — MCMC Convergence Diagnostics
# Adapted from generate_convergence_diagnostics.R (whole script)
# 6-panel layout: (a) R-hat distribution overlay, (b) ESS ratio overlay,
# (c) CHM R-hat by parameter type, (d) DTM same, (e) CHM traces, (f) DTM
# =====================================================================

log_progress("")
log_progress(strrep("-", 70))
log_progress("figS12 — MCMC Convergence Diagnostics")
log_progress(strrep("-", 70))

extract_diagnostics <- function(fit, product_name) {
  rhat_vals <- rhat(fit)
  neff_vals <- neff_ratio(fit)
  data.frame(
    parameter  = names(rhat_vals),
    Rhat       = as.numeric(rhat_vals),
    neff_ratio = as.numeric(neff_vals),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      type = case_when(
        grepl("^b_sigma",         parameter) ~ "Sigma Fixed",
        grepl("^b_",              parameter) ~ "Location Fixed",
        grepl("^r_site__sigma",   parameter) ~ "Sigma Site RE",
        grepl("^r_site",          parameter) ~ "Site RE",
        grepl("^r_ecoregion",     parameter) ~ "Ecoregion RE",
        grepl("^r_lc_l1_code",    parameter) ~ "Land Cover RE",
        grepl("^sd_",             parameter) ~ "Group-level SD",
        grepl("^cor_",            parameter) ~ "Correlations",
        grepl("nu",               parameter) ~ "Degrees of Freedom",
        grepl("sslope|^s_",       parameter) ~ "Spline Terms",
        TRUE                                 ~ "Other"
      ),
      product = product_name
    )
}

diag_chm <- extract_diagnostics(fit_chm, "CHM")
diag_dtm <- extract_diagnostics(fit_dtm, "DTM")
diag_all <- bind_rows(diag_chm, diag_dtm)

log_progress(sprintf("  CHM: %d parameters, max R-hat = %.4f, min ESS ratio = %.3f",
                     nrow(diag_chm),
                     max(diag_chm$Rhat, na.rm = TRUE),
                     min(diag_chm$neff_ratio, na.rm = TRUE)))
log_progress(sprintf("  DTM: %d parameters, max R-hat = %.4f, min ESS ratio = %.3f",
                     nrow(diag_dtm),
                     max(diag_dtm$Rhat, na.rm = TRUE),
                     min(diag_dtm$neff_ratio, na.rm = TRUE)))

# Panel (a): R-hat distribution overlay (CHM blue, DTM orange)
p_rhat <- ggplot(diag_all, aes(x = Rhat, fill = product)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity",
                 color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 1.01, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  annotate("text", x = 1.01, y = Inf, label = "  1.01 threshold",
           hjust = 0, vjust = 1.5, size = 3, color = "red") +
  scale_fill_manual(values = c("CHM" = "#3182bd", "DTM" = "#e6550d"), name = "Model") +
  scale_x_continuous(limits = c(0.999, NA),
                     breaks = seq(1.000, 1.010, by = 0.002)) +
  labs(title = "(a) R-hat distribution",
       x = expression(hat(R)), y = "Number of parameters") +
  theme_pub() +
  theme(legend.position = c(0.82, 0.82))

# Panel (b): ESS ratio distribution overlay
p_neff <- ggplot(diag_all, aes(x = neff_ratio, fill = product)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity",
                 color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.1, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  annotate("text", x = 0.1, y = Inf, label = "  0.1 threshold",
           hjust = 0, vjust = 1.5, size = 3, color = "red") +
  scale_fill_manual(values = c("CHM" = "#3182bd", "DTM" = "#e6550d"), name = "Model") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  labs(title = "(b) Effective sample size ratio",
       x = "ESS / Total draws", y = "Number of parameters") +
  theme_pub() +
  theme(legend.position = "none")

# Panel (c): CHM R-hat by parameter type
type_order_chm <- diag_chm %>%
  group_by(type) %>%
  summarise(med = median(Rhat, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(type)
diag_chm$type_f <- factor(diag_chm$type, levels = rev(type_order_chm))

p_rhat_chm <- ggplot(diag_chm %>% filter(!is.na(type_f)),
                     aes(x = Rhat, y = type_f)) +
  geom_jitter(alpha = 0.4, size = 1, color = "#3182bd", height = 0.2) +
  geom_vline(xintercept = 1.01, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  scale_x_continuous(limits = c(0.999, 1.012),
                     breaks = seq(1.000, 1.010, by = 0.002)) +
  labs(title = "(c) CHM: R-hat by parameter type",
       x = expression(hat(R)), y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 8))

# Panel (d): DTM R-hat by parameter type
type_order_dtm <- diag_dtm %>%
  group_by(type) %>%
  summarise(med = median(Rhat, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(type)
diag_dtm$type_f <- factor(diag_dtm$type, levels = rev(type_order_dtm))

p_rhat_dtm <- ggplot(diag_dtm %>% filter(!is.na(type_f)),
                     aes(x = Rhat, y = type_f)) +
  geom_jitter(alpha = 0.4, size = 1, color = "#e6550d", height = 0.2) +
  geom_vline(xintercept = 1.01, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  scale_x_continuous(limits = c(0.999, 1.012),
                     breaks = seq(1.000, 1.010, by = 0.002)) +
  labs(title = "(d) DTM: R-hat by parameter type",
       x = expression(hat(R)), y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 8))

# Panels (e, f): trace plots
key_chm <- c("b_Intercept", "b_wsci_z", "b_slope_mean_z",
             "b_sigma_Intercept", "b_sigma_slope_mean_z", "nu")
avail_chm <- intersect(key_chm, variables(fit_chm))
log_progress(sprintf("  CHM trace params (%d): %s",
                     length(avail_chm), paste(avail_chm, collapse = ", ")))

chm_labels <- c(
  "b_Intercept"          = "Intercept",
  "b_wsci_z"             = "WSCI",
  "b_slope_mean_z"       = "Slope",
  "b_sigma_Intercept"    = "Sigma Intercept",
  "b_sigma_slope_mean_z" = "Sigma Slope",
  "nu"                   = "Degrees of Freedom (v)"
)

p_trace_chm <- mcmc_trace(
  fit_chm, pars = avail_chm,
  facet_args = list(ncol = 1, strip.position = "right",
                    labeller = as_labeller(chm_labels))) +
  labs(title = "(e) CHM: Trace plots") +
  theme_pub(base_size = 9) +
  theme(strip.text = element_text(size = 6),
        strip.text.y.right = element_text(angle = 0, hjust = 0,
                                          margin = margin(l = 2, r = 2)),
        plot.margin   = margin(5, 40, 5, 5),
        legend.position = "none")

key_dtm <- c("b_Intercept", "b_rh_98_z", "b_slope_mean_z",
             "b_sigma_Intercept", "b_sigma_slope_mean_z", "nu")
avail_dtm <- intersect(key_dtm, variables(fit_dtm))
log_progress(sprintf("  DTM trace params (%d): %s",
                     length(avail_dtm), paste(avail_dtm, collapse = ", ")))

dtm_labels <- c(
  "b_Intercept"          = "Intercept",
  "b_rh_98_z"            = "RH98",
  "b_slope_mean_z"       = "Slope",
  "b_sigma_Intercept"    = "Sigma Intercept",
  "b_sigma_slope_mean_z" = "Sigma Slope",
  "nu"                   = "Degrees of Freedom (v)"
)

p_trace_dtm <- mcmc_trace(
  fit_dtm, pars = avail_dtm,
  facet_args = list(ncol = 1, strip.position = "right",
                    labeller = as_labeller(dtm_labels))) +
  labs(title = "(f) DTM: Trace plots") +
  theme_pub(base_size = 9) +
  theme(strip.text = element_text(size = 6),
        strip.text.y.right = element_text(angle = 0, hjust = 0,
                                          margin = margin(l = 2, r = 2)),
        plot.margin   = margin(5, 40, 5, 5),
        legend.position = "none")

top_row <- p_rhat     | p_neff
mid_row <- p_rhat_chm | p_rhat_dtm
bot_row <- p_trace_chm | p_trace_dtm

p_s12 <- top_row / mid_row / bot_row +
  plot_layout(heights = c(1, 1.2, 1.5)) +
  plot_annotation(
    title = "MCMC Convergence Diagnostics",
    subtitle = sprintf(
      "CHM: %d parameters, max R-hat = %.4f | DTM: %d parameters, max R-hat = %.4f",
      nrow(diag_chm), max(diag_chm$Rhat, na.rm = TRUE),
      nrow(diag_dtm), max(diag_dtm$Rhat, na.rm = TRUE)
    ),
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )
  )

ggsave(file.path(out_dir, "figS12_convergence_16site.pdf"),
       p_s12, width = 16, height = 16, bg = "white", device = pdf)
log_progress(sprintf("  Saved: %s/figS12_convergence_16site.pdf", out_dir))

# =====================================================================
# figS18 — Posterior Predictive Checks (CHM/DTM density overlay)
# Adapted from generate_ppc_overall.R (whole script)
# Subtitle updated to "95% PI coverage: X.X% (full holdout)" to match v04
# Coverage computed on the n_samples = 10000 stratified subsample, which
# is the earlier convention. The step 01 preflight reported 94.41%/94.44% which is in
# the same regime as v04's 94.0%/93.4% — DTM differs from v04 because the
# v04 subsample seed differs from this script's seed.
# =====================================================================

log_progress("")
log_progress(strrep("-", 70))
log_progress("figS18 — Posterior Predictive Checks")
log_progress(strrep("-", 70))

n_samples <- 10000
n_draws_ppc <- 500
n_draws_plot <- 50
stage2_frac <- as.numeric(Sys.getenv("STAGE2_FRAC", "0.33"))
set.seed(2025)

# Reconstruct holdout sets (same recipe as generate_ppc_overall.R)
log_progress("  Loading full dataset for holdout reconstruction...")
data_chk <- load_checkpoint("01_data_ingest")
chm_full_raw <- data_chk$chm_df
dtm_full_raw <- data_chk$dtm_df

chm_full <- chm_full_raw %>%
  filter(is.finite(chm_error_mean), !is.na(lc_l1_code),
         !is.na(site), !is.na(ecoregion)) %>%
  mutate(lc_l1_code = factor(lc_l1_code),
         ecoregion  = factor(ecoregion),
         site       = factor(site))
dtm_full <- dtm_full_raw %>%
  filter(is.finite(dtm_error_mean), !is.na(lc_l1_code),
         !is.na(site), !is.na(ecoregion)) %>%
  mutate(lc_l1_code = factor(lc_l1_code),
         ecoregion  = factor(ecoregion),
         site       = factor(site))

# For 16-site CHM model: filter to the 16 retained sites BEFORE training
# split, mirroring how section_10b_chm_sensitivity_refit.R handled it
sites_excluded_chm <- cp_chm$sites_excluded
log_progress(sprintf("  CHM 16-site refit excludes: %s",
                     paste(sites_excluded_chm, collapse = ", ")))
chm_full <- chm_full %>% filter(!as.character(site) %in% sites_excluded_chm)
chm_full <- chm_full %>% mutate(site = droplevels(factor(site)))

chm_full$original_row <- seq_len(nrow(chm_full))
dtm_full$original_row <- seq_len(nrow(dtm_full))

set.seed(2025)
chm_train_smp <- chm_full %>% group_by(site) %>%
  sample_frac(stage2_frac) %>% ungroup()
dtm_train_smp <- dtm_full %>% group_by(site) %>%
  sample_frac(stage2_frac) %>% ungroup()

train_lc_chm   <- levels(droplevels(chm_train_smp$lc_l1_code))
train_site_chm <- levels(droplevels(chm_train_smp$site))
train_eco_chm  <- levels(droplevels(chm_train_smp$ecoregion))
train_lc_dtm   <- levels(droplevels(dtm_train_smp$lc_l1_code))
train_site_dtm <- levels(droplevels(dtm_train_smp$site))
train_eco_dtm  <- levels(droplevels(dtm_train_smp$ecoregion))

chm_holdout <- chm_full %>%
  filter(!original_row %in% chm_train_smp$original_row) %>%
  select(-original_row) %>%
  filter(lc_l1_code %in% train_lc_chm,
         site       %in% train_site_chm,
         ecoregion  %in% train_eco_chm) %>%
  mutate(lc_l1_code = factor(lc_l1_code, levels = train_lc_chm),
         site       = factor(site,       levels = train_site_chm),
         ecoregion  = factor(ecoregion,  levels = train_eco_chm))

dtm_holdout <- dtm_full %>%
  filter(!original_row %in% dtm_train_smp$original_row) %>%
  select(-original_row) %>%
  filter(lc_l1_code %in% train_lc_dtm,
         site       %in% train_site_dtm,
         ecoregion  %in% train_eco_dtm) %>%
  mutate(lc_l1_code = factor(lc_l1_code, levels = train_lc_dtm),
         site       = factor(site,       levels = train_site_dtm),
         ecoregion  = factor(ecoregion,  levels = train_eco_dtm))

# Site 10 still excluded from DTM (per project convention)
dtm_holdout <- dtm_holdout %>% filter(!site %in% "10")

log_progress(sprintf("  CHM holdout: %s obs | DTM holdout: %s obs",
                     format(nrow(chm_holdout), big.mark = ","),
                     format(nrow(dtm_holdout), big.mark = ",")))

# Stratified subsample for PPC
chm_ppc_sample <- chm_holdout %>%
  group_by(site) %>%
  slice_sample(n = ceiling(n_samples / n_distinct(chm_holdout$site)),
               replace = FALSE) %>% ungroup() %>%
  slice_sample(n = min(n_samples, nrow(.)))
dtm_ppc_sample <- dtm_holdout %>%
  group_by(site) %>%
  slice_sample(n = ceiling(n_samples / n_distinct(dtm_holdout$site)),
               replace = FALSE) %>% ungroup() %>%
  slice_sample(n = min(n_samples, nrow(.)))

log_progress(sprintf("  CHM PPC sample: %d obs | DTM PPC sample: %d obs",
                     nrow(chm_ppc_sample), nrow(dtm_ppc_sample)))

log_progress("  Generating CHM posterior predictions...")
chm_yrep <- posterior_predict(fit_chm, newdata = chm_ppc_sample,
                              ndraws = n_draws_ppc, allow_new_levels = TRUE)
log_progress("  Generating DTM posterior predictions...")
dtm_yrep <- posterior_predict(fit_dtm, newdata = dtm_ppc_sample,
                              ndraws = n_draws_ppc, allow_new_levels = TRUE)

compute_coverage <- function(y_obs, y_rep, alpha = 0.95) {
  lower <- (1 - alpha) / 2; upper <- 1 - lower
  pi_lo <- apply(y_rep, 2, quantile, probs = lower, na.rm = TRUE)
  pi_hi <- apply(y_rep, 2, quantile, probs = upper, na.rm = TRUE)
  mean((y_obs >= pi_lo) & (y_obs <= pi_hi), na.rm = TRUE)
}

chm_y_obs <- chm_ppc_sample$chm_error_mean
dtm_y_obs <- dtm_ppc_sample$dtm_error_mean

cov_chm <- compute_coverage(chm_y_obs, chm_yrep)
cov_dtm <- compute_coverage(dtm_y_obs, dtm_yrep)

log_progress(sprintf("  CHM 95%% PI coverage: %.2f%%", cov_chm * 100))
log_progress(sprintf("  DTM 95%% PI coverage: %.2f%%", cov_dtm * 100))

# Persist for reference
fwrite(data.table(
  product       = c("CHM", "DTM"),
  n_sample      = c(length(chm_y_obs), length(dtm_y_obs)),
  coverage_95pi = c(cov_chm, cov_dtm)),
  file.path(manuscript_tables_dir, "section_L5_figS18_coverage.csv"))

# Drop non-finite cols before plotting
chm_finite <- is.finite(chm_y_obs) &
  apply(chm_yrep, 2, function(col) all(is.finite(col)))
dtm_finite <- is.finite(dtm_y_obs) &
  apply(dtm_yrep, 2, function(col) all(is.finite(col)))
chm_y_clean    <- chm_y_obs[chm_finite]
chm_yrep_clean <- chm_yrep[, chm_finite, drop = FALSE]
dtm_y_clean    <- dtm_y_obs[dtm_finite]
dtm_yrep_clean <- dtm_yrep[, dtm_finite, drop = FALSE]

n_plot_chm <- min(n_draws_plot, nrow(chm_yrep_clean))
n_plot_dtm <- min(n_draws_plot, nrow(dtm_yrep_clean))

p_ppc_chm <- ppc_dens_overlay(chm_y_clean,
                              chm_yrep_clean[seq_len(n_plot_chm), ]) +
  coord_cartesian(xlim = c(-30, 30)) +
  labs(
    title    = sprintf("CHM (n = %s)", format(length(chm_y_clean), big.mark = ",")),
    subtitle = sprintf("95%% PI coverage: %.1f%% (full holdout)", cov_chm * 100),
    x = "CHM Error (m)", y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle      = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.title         = element_text(size = 11),
    axis.text          = element_text(size = 10),
    panel.grid.minor   = element_blank()
  )

p_ppc_dtm <- ppc_dens_overlay(dtm_y_clean,
                              dtm_yrep_clean[seq_len(n_plot_dtm), ]) +
  coord_cartesian(xlim = c(-30, 30)) +
  labs(
    title    = sprintf("DTM (n = %s)", format(length(dtm_y_clean), big.mark = ",")),
    subtitle = sprintf("95%% PI coverage: %.1f%% (full holdout)", cov_dtm * 100),
    x = "DTM Error (m)", y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle      = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.title         = element_text(size = 11),
    axis.text          = element_text(size = 10),
    panel.grid.minor   = element_blank()
  )

ppc_combined <- plot_grid(
  p_ppc_chm, p_ppc_dtm,
  nrow = 1,
  labels = c("(a)", "(b)"),
  label_size = 12, label_fontface = "bold"
)

p_s18 <- plot_grid(
  ggdraw() +
    draw_label("Posterior Predictive Checks",
               fontface = "bold", size = 14) +
    draw_label("Observed error distribution (dark line) vs. posterior predictive draws (light lines)",
               size = 10, y = 0.2, color = "gray40"),
  ppc_combined,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

ggsave(file.path(out_dir, "figS18_ppc_overall_16site.pdf"),
       p_s18, width = 10, height = 5, bg = "white", device = pdf)
log_progress(sprintf("  Saved: %s/figS18_ppc_overall_16site.pdf", out_dir))

# =====================================================================
# figS20 — Semi-variograms of CHM/DTM residuals (2-panel)
# Adapted from section_12_spatial.R lines 463-700, with the DTM panel
# block written to mirror the CHM block exactly.
#
# Methodology (from section_12_spatial.R):
#   - TARGET_CRS = 5070 (NAD83 Albers Equal Area)
#   - MAX_LAG_DISTANCE = 50000 m (50 km)
#   - MIN_LAG_DISTANCE = 100 m
#   - n_lags = 20  =>  lag_width = 2495 m
#   - Try Exp / Sph / Gau / Mat models, pick lowest weighted SSE
#   - Stratified sampling: 5000 obs per site, set.seed(42)
# =====================================================================

log_progress("")
log_progress(strrep("-", 70))
log_progress("figS20 — Semi-variograms (CHM 16-site + DTM 18-site)")
log_progress(strrep("-", 70))

TARGET_CRS         <- 5070
SOURCE_CRS         <- 4326
GEDI_SWATH_WIDTH   <- 6000
MIN_LAG_DISTANCE   <- 100
MAX_LAG_DISTANCE   <- 50000
N_LAGS             <- 20

# Both checkpoints should have x, y attached on the data
chm_df_all <- data_chk$chm_df
dtm_df_all <- data_chk$dtm_df
stopifnot(all(c("x", "y") %in% names(chm_df_all)))
stopifnot(all(c("x", "y") %in% names(dtm_df_all)))

# Filter CHM to 16 sites to match the 16-site refit
chm_df_all <- chm_df_all %>%
  filter(!as.character(site) %in% sites_excluded_chm)

# Sample for spatial analysis (matches section_12)
set.seed(42)
chm_s <- chm_df_all %>%
  tidyr::drop_na(x, y, chm_error_mean) %>%
  group_by(site) %>%
  group_modify(~ slice_sample(.x, n = min(5000, nrow(.x)), replace = FALSE)) %>%
  ungroup()
dtm_s <- dtm_df_all %>%
  tidyr::drop_na(x, y, dtm_error_mean) %>%
  group_by(site) %>%
  group_modify(~ slice_sample(.x, n = min(5000, nrow(.x)), replace = FALSE)) %>%
  ungroup()
log_progress(sprintf("  CHM spatial sample: %s | DTM spatial sample: %s",
                     format(nrow(chm_s), big.mark = ","),
                     format(nrow(dtm_s), big.mark = ",")))

# Project to Albers (meters)
chm_sf <- st_transform(
  st_as_sf(chm_s, coords = c("x", "y"), crs = SOURCE_CRS, remove = FALSE),
  TARGET_CRS)
dtm_sf <- st_transform(
  st_as_sf(dtm_s, coords = c("x", "y"), crs = SOURCE_CRS, remove = FALSE),
  TARGET_CRS)

chm_s <- chm_s %>%
  mutate(x_proj = st_coordinates(chm_sf)[, 1],
         y_proj = st_coordinates(chm_sf)[, 2])
dtm_s <- dtm_s %>%
  mutate(x_proj = st_coordinates(dtm_sf)[, 1],
         y_proj = st_coordinates(dtm_sf)[, 2])

# Extract model residuals using compatible-factor-levels approach
extract_brms_residuals <- function(model, newdata, response_col, name) {
  model_data <- model$data
  factor_vars <- names(model_data)[sapply(model_data, is.factor)]
  valid_idx <- rep(TRUE, nrow(newdata))
  for (fvar in factor_vars) {
    if (fvar %in% names(newdata)) {
      model_levels <- levels(model_data[[fvar]])
      valid_idx <- valid_idx & (as.character(newdata[[fvar]]) %in% model_levels)
    }
  }
  if (sum(valid_idx) < 100) {
    log_progress(sprintf("    %s: too few valid obs; using raw error", name))
    return(list(residuals = newdata[[response_col]],
                valid_idx = seq_len(nrow(newdata))))
  }
  newdata_valid <- newdata[valid_idx, ]
  for (fvar in factor_vars) {
    if (fvar %in% names(newdata_valid)) {
      newdata_valid[[fvar]] <- factor(newdata_valid[[fvar]],
                                       levels = levels(model_data[[fvar]]))
    }
  }
  log_progress(sprintf("    %s: predicting on %d obs (factor-compatible)...",
                       name, nrow(newdata_valid)))
  pred <- fitted(model, newdata = newdata_valid,
                 summary = TRUE, ndraws = 100)
  resid_vec <- newdata_valid[[response_col]] - pred[, "Estimate"]
  list(residuals = resid_vec, valid_idx = which(valid_idx))
}

log_progress("  CHM residual extraction:")
chm_res_obj <- extract_brms_residuals(fit_chm, chm_s, "chm_error_mean", "CHM")
chm_s_v <- chm_s[chm_res_obj$valid_idx, ]
chm_s_v$residual <- chm_res_obj$residuals

log_progress("  DTM residual extraction:")
dtm_res_obj <- extract_brms_residuals(fit_dtm, dtm_s, "dtm_error_mean", "DTM")
dtm_s_v <- dtm_s[dtm_res_obj$valid_idx, ]
dtm_s_v$residual <- dtm_res_obj$residuals

# Empirical variograms
vg_calc_projected <- function(df_with_coords, res_col,
                              max_dist = MAX_LAG_DISTANCE,
                              n_lags = N_LAGS,
                              min_lag = MIN_LAG_DISTANCE) {
  df <- df_with_coords %>%
    mutate(res = .data[[res_col]]) %>%
    filter(is.finite(res))
  sf_obj <- st_as_sf(df, coords = c("x_proj", "y_proj"), crs = TARGET_CRS)
  g <- gstat::gstat(id = "res", formula = res ~ 1, data = sf_obj)
  lag_width <- (max_dist - min_lag) / n_lags
  gstat::variogram(g, cutoff = max_dist, width = lag_width)
}

log_progress("  Computing empirical variograms (50 km cutoff, 20 bins)...")
vg_chm <- vg_calc_projected(chm_s_v, "residual")
vg_dtm <- vg_calc_projected(dtm_s_v, "residual")

# Model fitting (Exp / Sph / Gau / Mat with selection by weighted SSE)
fit_variogram_models <- function(vg_empirical, name) {
  if (is.null(vg_empirical)) return(NULL)
  results <- list()
  max_gamma <- max(vg_empirical$gamma)
  max_dist  <- max(vg_empirical$dist)
  init_nugget <- vg_empirical$gamma[1] * 0.8
  init_sill   <- max_gamma - init_nugget
  gamma_target <- init_nugget + 0.63 * init_sill
  init_range  <- vg_empirical$dist[which.min(abs(vg_empirical$gamma - gamma_target))]
  if (init_range < MIN_LAG_DISTANCE) init_range <- max_dist / 4

  for (m in c("Exp", "Sph", "Gau", "Mat")) {
    tryCatch({
      vg_model <- if (m == "Mat") {
        gstat::vgm(psill = init_sill, model = m, range = init_range,
                   nugget = init_nugget, kappa = 1.5)
      } else {
        gstat::vgm(psill = init_sill, model = m, range = init_range,
                   nugget = init_nugget)
      }
      fit <- gstat::fit.variogram(vg_empirical, vg_model)
      predicted <- gstat::variogramLine(fit, maxdist = max_dist,
                                        n = nrow(vg_empirical))
      pred_at_emp <- approx(predicted$dist, predicted$gamma,
                            vg_empirical$dist)$y
      weights <- vg_empirical$np / sum(vg_empirical$np)
      sse <- sum(weights * (vg_empirical$gamma - pred_at_emp)^2, na.rm = TRUE)
      nugget       <- fit$psill[1]
      partial_sill <- fit$psill[2]
      total_sill   <- sum(fit$psill)
      range_param  <- fit$range[2]
      results[[m]] <- list(
        model = fit, sse = sse, nugget = nugget,
        partial_sill = partial_sill, sill = total_sill,
        range = range_param, nugget_sill_ratio = nugget / total_sill
      )
      log_progress(sprintf("    %s %s: nugget=%.2f, sill=%.2f, range=%.0f m, SSE=%.4f",
                           name, m, nugget, total_sill, range_param, sse))
    }, error = function(e) {
      log_progress(sprintf("    %s %s: fit failed (%s)", name, m, e$message))
    })
  }
  if (length(results) == 0) return(NULL)
  sse_vals <- sapply(results, function(x) x$sse)
  best <- names(which.min(sse_vals))
  results$best     <- best
  results$best_fit <- results[[best]]
  log_progress(sprintf("    %s best model: %s", name, best))
  results
}

vg_fits_chm <- fit_variogram_models(vg_chm, "CHM")
vg_fits_dtm <- fit_variogram_models(vg_dtm, "DTM")

# Persist parameters
params_df <- data.table(
  product           = c("CHM", "DTM"),
  best_model        = c(vg_fits_chm$best,                 vg_fits_dtm$best),
  nugget            = c(vg_fits_chm$best_fit$nugget,      vg_fits_dtm$best_fit$nugget),
  sill              = c(vg_fits_chm$best_fit$sill,        vg_fits_dtm$best_fit$sill),
  range_m           = c(vg_fits_chm$best_fit$range,       vg_fits_dtm$best_fit$range),
  nugget_sill_ratio = c(vg_fits_chm$best_fit$nugget_sill_ratio,
                        vg_fits_dtm$best_fit$nugget_sill_ratio)
)
fwrite(params_df,
       file.path(manuscript_tables_dir, "section_L5_figS20_variogram_params.csv"))

# Build CHM panel (mirrors section_12_spatial.R lines 679-700)
build_variogram_panel <- function(vg_data, vg_fit, title_prefix) {
  vg_plot_data <- vg_data
  fitted_line <- gstat::variogramLine(vg_fit$best_fit$model,
                                       maxdist = max(vg_data$dist), n = 100)
  ggplot() +
    geom_point(data = vg_plot_data,
               aes(x = dist / 1000, y = gamma, size = np), alpha = 0.7) +
    geom_line(data = fitted_line,
              aes(x = dist / 1000, y = gamma),
              color = "red", linewidth = 1) +
    geom_hline(yintercept = vg_fit$best_fit$sill,
               linetype = "dashed", color = "blue", alpha = 0.7) +
    geom_vline(xintercept = vg_fit$best_fit$range / 1000,
               linetype = "dashed", color = "forestgreen", alpha = 0.7) +
    geom_vline(xintercept = GEDI_SWATH_WIDTH / 1000,
               linetype = "dotted", color = "orange", alpha = 0.5) +
    annotate("text", x = GEDI_SWATH_WIDTH / 1000 + 1,
             y = max(vg_plot_data$gamma) * 0.95,
             label = "GEDI swath", color = "orange",
             hjust = 0, size = 3) +
    labs(
      title    = sprintf("%s: Semi-variogram with Fitted Model", title_prefix),
      subtitle = sprintf("%s model: nugget=%.2f m\u00b2, sill=%.2f m\u00b2, range=%.1f km\nNugget/Sill ratio=%.3f",
                         vg_fit$best, vg_fit$best_fit$nugget,
                         vg_fit$best_fit$sill, vg_fit$best_fit$range / 1000,
                         vg_fit$best_fit$nugget_sill_ratio),
      x = "Distance (km)", y = "Semivariance (m\u00b2)", size = "N pairs"
    ) +
    theme_cowplot() +
    theme(plot.subtitle = element_text(size = 9))
}

p_vg_chm <- build_variogram_panel(vg_chm, vg_fits_chm, "CHM")
p_vg_dtm <- build_variogram_panel(vg_dtm, vg_fits_dtm, "DTM")

p_s20 <- plot_grid(p_vg_chm, p_vg_dtm, ncol = 2,
                   labels = c("(a)", "(b)"),
                   label_size = 12, label_fontface = "bold")

ggsave(file.path(out_dir, "figS20_variogram_2panel_16site.pdf"),
       p_s20, width = 16, height = 7, bg = "white", device = pdf)
log_progress(sprintf("  Saved: %s/figS20_variogram_2panel_16site.pdf", out_dir))

# ---------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------

log_progress("")
log_progress(strrep("=", 70))
log_progress("Supplementary figure re-render complete. Outputs:")
log_progress(sprintf("  %s/figS06_pred_vs_obs_by_forest_16site.pdf", out_dir))
log_progress(sprintf("  %s/figS12_convergence_16site.pdf",           out_dir))
log_progress(sprintf("  %s/figS18_ppc_overall_16site.pdf",           out_dir))
log_progress(sprintf("  %s/figS20_variogram_2panel_16site.pdf",      out_dir))
log_progress("")
log_progress("Next: convert each PDF to PNG at 300 DPI via gs, e.g.:")
log_progress("  cd plots/groundwork && for f in figS06_*.pdf figS12_*.pdf \\")
log_progress("    figS18_*.pdf figS20_*.pdf; do")
log_progress("    gs -dQUIET -dBATCH -dNOPAUSE -sDEVICE=png16m -r300 \\")
log_progress("       -sOutputFile=\"${f%.pdf}.png\" \"$f\"; done")
log_progress(strrep("=", 70))
