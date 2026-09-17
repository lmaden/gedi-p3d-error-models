# =====================================================================
# fig5_pred_vs_obs.R   (CANONICAL HOLDOUT RECIPE)
# Figure 5: predicted-versus-observed scatter on holdout data.
#
# v04 -> v10 changes:
#   - CHM panel (A) regenerated from the 16-site CHM fit:
#       Story-lock §2A targets: R2 = 0.112, RMSE = 4.479 m, MAE 2.974 m,
#       bias +0.485 m, PI coverage 0.940, PI width median 18.20 m
#       (n = 253,689 holdout footprints).
#   - DTM panel (B) unchanged (DTM not refit under 16-site).
#   - APA34 (id=140): aspect ratio rebalanced to square panels;
#     headline stats moved INSIDE the panel in a left-justified upper
#     block at base_size + 1.
#   - APA35 (id=139): hexbin bin count drops the lower-count exclusion
#     threshold so all populated hexagons render. Cap count via
#     log-scale fill, not via filtering.
#
# HOLDOUT RECIPE (canonical; matches verify_predictive_accuracy.R and
# ppc_by_landcover.R, the scripts that produced the +0.485 m
# bias and 253,689 anchor in story-lock §2A):
#   1. Load FULL chm_df / dtm_df from 01_data_ingest.rds.
#   2. Filter is.finite(<response>), !is.na(lc_l1_code, site, ecoregion).
#   3. Factor() conversions matching section_08.
#   4. set.seed(2025), add original_row, group_by(site) %>%
#      sample_frac(0.33) -> training set (matches section_08 stage 2).
#   5. Anti-join via original_row -> 19-site holdout.
#   6. Filter to training factor levels (drops novel LCs/sites/ecoregions
#      that didn't make it into training -- prevents brm prediction failures).
#   7. For CHM 16-site frame: drop sites 1, 2, 3 from holdout
#      (matches section_10b SITES_TO_EXCLUDE).
#   8. For DTM 18-site frame: drop site 10 from holdout
#      (matches SITES_EXCLUDE_DTM in analysis_config.R).
#   9. posterior_predict(fit, newdata = holdout, re_formula = NULL).
#
# This recipe should reproduce the 253,689 / +0.485 m anchor exactly,
# bit-for-bit, modulo posterior_predict Monte Carlo noise on n_draws.
# Use N_DRAWS = 1000 here for tight bias agreement.
#
# Response columns are chm_error_mean / dtm_error_mean (per section_08).
#
# Sources:
#   - load_checkpoint("10b_chm_sensitivity")$fit_chm_16
#   - load_checkpoint("10_models_stage2")$fit_dtm_s2
#   - load_checkpoint("01_data_ingest")  (raw chm_df, dtm_df pre-split)
# =====================================================================

source("fig_common.R")
fig_banner("Figure 5",
                 "Predicted vs observed on canonical holdout (CHM 16-site, DTM 18-site)")

suppressPackageStartupMessages({
  library(brms); library(posterior); library(ggplot2); library(scales)
  library(data.table); library(dplyr); library(patchwork)
})

# ---------------------------------------------------------------------
# 1. Load fits
# ---------------------------------------------------------------------
log_subsection("Loading fits")

chm_ck <- load_checkpoint("10b_chm_sensitivity")
fit_chm_16 <- chm_ck$fit_chm_16

s2 <- load_checkpoint("10_models_stage2")
fit_dtm <- s2$fit_dtm_s2

# ---------------------------------------------------------------------
# 2. Reconstruct canonical holdout
# ---------------------------------------------------------------------
log_subsection("Reconstructing canonical holdout (seed 2025, frac 0.33 per site)")

stage2_frac <- as.numeric(Sys.getenv("STAGE2_FRAC", "0.33"))
log_progress(sprintf("  stage2_frac = %.2f", stage2_frac))

# Sites excluded at fit time. CHM 16-site: 1, 2, 3. DTM 18-site: 10.
SITES_EXCLUDE_CHM_16 <- c("1", "2", "3")
SITES_EXCLUDE_DTM_18 <- c("10")

# Load full unsampled data from data-ingest checkpoint.
full_ck <- load_checkpoint("01_data_ingest")
if (!all(c("chm_df", "dtm_df") %in% names(full_ck))) {
  stop("01_data_ingest checkpoint missing chm_df/dtm_df. Names: ",
       paste(names(full_ck), collapse = ", "))
}

chm_df_full <- as.data.frame(full_ck$chm_df) %>%
  filter(is.finite(chm_error_mean),
         !is.na(lc_l1_code), !is.na(site), !is.na(ecoregion)) %>%
  mutate(lc_l1_code = factor(lc_l1_code),
         ecoregion  = factor(ecoregion),
         site       = factor(site))

dtm_df_full <- as.data.frame(full_ck$dtm_df) %>%
  filter(is.finite(dtm_error_mean),
         !is.na(lc_l1_code), !is.na(site), !is.na(ecoregion)) %>%
  mutate(lc_l1_code = factor(lc_l1_code),
         ecoregion  = factor(ecoregion),
         site       = factor(site))

log_progress(sprintf("  Full filtered: %s CHM, %s DTM",
                     formatC(nrow(chm_df_full), big.mark = ",", format = "d"),
                     formatC(nrow(dtm_df_full), big.mark = ",", format = "d")))

# Add original_row, then training split using the canonical seed.
chm_df_full$original_row <- seq_len(nrow(chm_df_full))
dtm_df_full$original_row <- seq_len(nrow(dtm_df_full))

set.seed(2025)
chm_train <- chm_df_full %>% group_by(site) %>%
  sample_frac(stage2_frac) %>% ungroup()
dtm_train <- dtm_df_full %>% group_by(site) %>%
  sample_frac(stage2_frac) %>% ungroup()

log_progress(sprintf("  Training (recreated): %s CHM, %s DTM",
                     formatC(nrow(chm_train), big.mark = ",", format = "d"),
                     formatC(nrow(dtm_train), big.mark = ",", format = "d")))

# Holdout = anti-join via original_row, then restrict to training factor levels.
train_lc_chm   <- levels(droplevels(chm_train$lc_l1_code))
train_site_chm <- levels(droplevels(chm_train$site))
train_eco_chm  <- levels(droplevels(chm_train$ecoregion))

train_lc_dtm   <- levels(droplevels(dtm_train$lc_l1_code))
train_site_dtm <- levels(droplevels(dtm_train$site))
train_eco_dtm  <- levels(droplevels(dtm_train$ecoregion))

chm_holdout <- chm_df_full %>%
  filter(!original_row %in% chm_train$original_row) %>%
  select(-original_row) %>%
  filter(lc_l1_code %in% train_lc_chm,
         site       %in% train_site_chm,
         ecoregion  %in% train_eco_chm) %>%
  mutate(lc_l1_code = factor(lc_l1_code, levels = train_lc_chm),
         site       = factor(site,       levels = train_site_chm),
         ecoregion  = factor(ecoregion,  levels = train_eco_chm))

dtm_holdout <- dtm_df_full %>%
  filter(!original_row %in% dtm_train$original_row) %>%
  select(-original_row) %>%
  filter(lc_l1_code %in% train_lc_dtm,
         site       %in% train_site_dtm,
         ecoregion  %in% train_eco_dtm) %>%
  mutate(lc_l1_code = factor(lc_l1_code, levels = train_lc_dtm),
         site       = factor(site,       levels = train_site_dtm),
         ecoregion  = factor(ecoregion,  levels = train_eco_dtm))

log_progress(sprintf("  19-site CHM holdout (pre-site-exclude): %s",
                     formatC(nrow(chm_holdout), big.mark = ",", format = "d")))
log_progress(sprintf("  19-site DTM holdout (pre-site-exclude): %s",
                     formatC(nrow(dtm_holdout), big.mark = ",", format = "d")))

# Apply CHM 16-site exclusion (drop sites 1, 2, 3 from holdout).
chm_holdout_16 <- chm_holdout %>%
  filter(!site %in% SITES_EXCLUDE_CHM_16) %>%
  mutate(site       = droplevels(site),
         lc_l1_code = droplevels(lc_l1_code),
         ecoregion  = droplevels(ecoregion))

# Apply DTM exclusion (drop site 10).
dtm_holdout_18 <- dtm_holdout %>%
  filter(!site %in% SITES_EXCLUDE_DTM_18) %>%
  mutate(site       = droplevels(site),
         lc_l1_code = droplevels(lc_l1_code),
         ecoregion  = droplevels(ecoregion))

# Add eco_on if the model formula expects it.
fit_chm_vars <- all.vars(formula(fit_chm_16))
fit_dtm_vars <- all.vars(formula(fit_dtm))
if ("eco_on" %in% fit_chm_vars) chm_holdout_16$eco_on <- 1
if ("eco_on" %in% fit_dtm_vars) dtm_holdout_18$eco_on <- 1

log_progress(sprintf("  CHM 16-site holdout (target 253,689): %s",
                     formatC(nrow(chm_holdout_16), big.mark = ",", format = "d")))
log_progress(sprintf("  DTM 18-site holdout: %s",
                     formatC(nrow(dtm_holdout_18), big.mark = ",", format = "d")))

# ---------------------------------------------------------------------
# 3. Posterior predictive means on holdout
# Use a relatively dense N_DRAWS to keep Monte Carlo noise on the
# bias well below the third decimal place (so we cleanly see +0.485).
# ---------------------------------------------------------------------
log_subsection("Generating posterior predictive means")

N_DRAWS <- as.integer(Sys.getenv("FIG5_N_DRAWS", "1000"))
log_progress(sprintf("  N_DRAWS = %d", N_DRAWS))

predict_y <- function(fit, df, ndraws = N_DRAWS) {
  pm <- posterior_predict(fit, newdata = df, ndraws = ndraws,
                          allow_new_levels = FALSE,
                          re_formula = NULL)
  colMeans(pm)
}

chm_holdout_16$y_obs  <- chm_holdout_16$chm_error_mean
chm_holdout_16$y_pred <- predict_y(fit_chm_16, chm_holdout_16)

dtm_holdout_18$y_obs  <- dtm_holdout_18$dtm_error_mean
dtm_holdout_18$y_pred <- predict_y(fit_dtm,    dtm_holdout_18)

# ---------------------------------------------------------------------
# 4. Panel stats
# Bias convention: pred - obs (matches canonical verify_predictive_accuracy.R
# line 139: Bias = mean(pred - obs); story-lock §2A anchor +0.485 m).
# Positive bias = model over-predicts observed error.
# R^2, RMSE, MAE are sign-insensitive (squared/absolute).
# ---------------------------------------------------------------------
panel_stats <- function(df) {
  r_signed <- df$y_pred - df$y_obs     # pred - obs (lab convention)
  r_resid  <- df$y_obs - df$y_pred     # statistical residual (for R^2 / RMSE / MAE)
  data.frame(
    n     = nrow(df),
    R2    = 1 - var(r_resid, na.rm = TRUE) / var(df$y_obs, na.rm = TRUE),
    RMSE  = sqrt(mean(r_resid^2, na.rm = TRUE)),
    MAE   = mean(abs(r_resid), na.rm = TRUE),
    bias  = mean(r_signed, na.rm = TRUE)   # pred - obs
  )
}
chm_stats <- panel_stats(chm_holdout_16)
dtm_stats <- panel_stats(dtm_holdout_18)

log_progress("Panel stats:")
print(rbind(cbind(product = "CHM", chm_stats),
            cbind(product = "DTM", dtm_stats)))

# ---------------------------------------------------------------------
# 5. Hexbin scatter per panel
# ---------------------------------------------------------------------
log_subsection("Composing hexbin panels")

build_hex_panel <- function(df, stats, title_lab, axis_lab_product) {
  xlim <- c(-20, 20)
  ylim <- c(-15, 15)
  ggplot(df, aes(x = y_obs, y = y_pred)) +
    geom_hex(bins = 80) +
    geom_abline(slope = 1, intercept = 0,
                color = COLOR_REF_LINE, linetype = "dashed",
                linewidth = 0.5) +
    geom_hline(yintercept = 0, color = COLOR_ZERO_LINE,
               linewidth = 0.2, linetype = "dotted") +
    geom_vline(xintercept = 0, color = COLOR_ZERO_LINE,
               linewidth = 0.2, linetype = "dotted") +
    scale_fill_viridis_c(trans = "log10", option = "viridis",
                         name = "Count") +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    annotate("text",
             x = xlim[1] + 0.5, y = ylim[2] - 0.5,
             label = sprintf("n = %s\nR^2 = %.3f\nRMSE = %.2f m\nMAE = %.2f m\nBias = %s%.2f m",
                             formatC(stats$n, big.mark = ",", format = "d"),
                             stats$R2,
                             stats$RMSE,
                             stats$MAE,
                             ifelse(stats$bias < 0, MINUS, "+"),
                             abs(stats$bias)),
             hjust = 0, vjust = 1, size = 3.4,
             color = "black", lineheight = 1.05,
             family = "mono") +
    labs(
      title = title_lab,
      x = sprintf("Observed %s error (m)", axis_lab_product),
      y = sprintf("Predicted %s error (m)", axis_lab_product)
    ) +
    theme_section_G(base_size = 11) +
    theme(
      legend.position    = "right",
      panel.grid.major   = element_line(color = "gray95", linewidth = 0.2),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(size = 12, face = "bold")
    )
}

p_chm <- build_hex_panel(chm_holdout_16, chm_stats,
                         "(a) CHM (16-site primary)", "CHM")
p_dtm <- build_hex_panel(dtm_holdout_18, dtm_stats,
                         "(b) DTM (18-site)", "DTM")

final_plot <- p_chm + p_dtm + patchwork::plot_layout(ncol = 2)

# ---------------------------------------------------------------------
# 6. Save and verify
# ---------------------------------------------------------------------
save_figure(final_plot, "fig05_pred_vs_obs",
            width_in = 10.5, height_in = 5.0)

log_subsection("VERIFICATION (story-lock anchors)")
cat("\nMust match story-lock §2A (CHM 16-site primary) exactly:\n")
cat("  n_holdout = 253,689\n")
cat("  R2        = 0.112\n")
cat("  RMSE      = 4.479 m\n")
cat("  MAE       = 2.974 m\n")
cat("  Bias      = +0.485 m  (pred - obs, per canonical verify_predictive_accuracy.R)\n")
cat("\nMust match story-lock §2B (DTM 18-site):\n")
cat("  R2        = 0.24\n")
cat("  RMSE      = 3.01 m\n")
cat("  MAE       = 2.07 m\n\n")

cat("Observed:\n")
cat(sprintf("  CHM 16-site:  n = %s, R2 = %.3f, RMSE = %.3f m, MAE = %.3f m, bias = %+.3f m\n",
            formatC(chm_stats$n, big.mark = ",", format = "d"),
            chm_stats$R2, chm_stats$RMSE, chm_stats$MAE, chm_stats$bias))
cat(sprintf("  DTM 18-site:  n = %s, R2 = %.3f, RMSE = %.3f m, MAE = %.3f m, bias = %+.3f m\n",
            formatC(dtm_stats$n, big.mark = ",", format = "d"),
            dtm_stats$R2, dtm_stats$RMSE, dtm_stats$MAE, dtm_stats$bias))

cat("\nDeltas vs story-lock:\n")
cat(sprintf("  CHM 16-site:  delta_n = %+d, delta_R2 = %+.4f, delta_RMSE = %+.4f m, delta_bias = %+.4f m\n",
            chm_stats$n - 253689L,
            chm_stats$R2 - 0.112,
            chm_stats$RMSE - 4.479,
            chm_stats$bias - 0.485))

cat("\nTolerance guidance:\n")
cat("  |delta_n|   < ~500 (0.2%%): expected, due to factor-level filter / sample_frac variance.\n")
cat("  |delta_R2|  < ~0.01:        expected MC + sample-difference noise.\n")
cat("  |delta_RMSE|< ~0.02 m:      excellent match.\n")
cat("  |delta_bias|< ~0.02 m:      excellent match (raise FIG5_N_DRAWS if larger).\n")

log_progress("fig5_pred_vs_obs.R complete.")
