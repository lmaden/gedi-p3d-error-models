# =====================================================================
# fig4_ppc_by_forest_type.R   (HOLDOUT BASIS)
# Figure 4: posterior predictive checks (PPC) stratified by forest
# functional type for CHM (16-site) and DTM (18-site).
#
# v04 -> v10 changes:
#   - CHM panels (BDF, DNF, EBF, ENF; row A-D) regenerated from
#     16-site fit (fit_chm_16) on the 16-site CHM forest-subset
#     HOLDOUT footprints (matches canonical recipe in
#     verify_predictive_accuracy.R / ppc_by_landcover.R).
#   - DTM panels (E-H) regenerated on 18-site DTM holdout (model
#     itself unchanged from v04; sampling source aligned).
#   - Sample sizes update to reflect the holdout pool (= 67% of
#     full data after site exclusions).
#   - 50 posterior predictive replicates per panel (lab convention).
#   - APA-pass: caption-level only; in-panel format is unchanged.
#
# Why holdout, not training: PPC on training data conflates fit with
# adequacy. The point is to assess whether the model's generative
# shape generalizes to footprints it has not seen.
#
# Sources:
#   - load_checkpoint("10b_chm_sensitivity")$fit_chm_16
#   - load_checkpoint("10_models_stage2")$fit_dtm_s2
#   - load_checkpoint("01_data_ingest")  (raw chm_df, dtm_df pre-split)
#   - Stratification draws on lc_l1_code; canonical forest LC values
#     are BDF, DNF, EBF, ENF (per story-lock §2C).
# =====================================================================

source("fig_common.R")
fig_banner("Figure 4",
                 "PPC on canonical holdout, stratified by forest type (CHM 16-site, DTM 18-site)")

suppressPackageStartupMessages({
  library(brms); library(bayesplot); library(posterior); library(ggplot2)
  library(data.table); library(dplyr); library(patchwork)
})

# Bayesplot theme aligned with section_G.
bayesplot::color_scheme_set("blue")

# ---------------------------------------------------------------------
# 1. Load fits and reconstruct canonical holdout
#    Same recipe as fig5_pred_vs_obs.R: load full data from
#    01_data_ingest, recreate the section_08 stage-2 training split
#    with seed(2025) + sample_frac(0.33), anti-join to get holdout,
#    then apply CHM 16-site / DTM 18-site site exclusions and match
#    training factor levels.
# ---------------------------------------------------------------------
log_subsection("Loading fits")

chm_ck <- load_checkpoint("10b_chm_sensitivity")
fit_chm_16 <- chm_ck$fit_chm_16

s2 <- load_checkpoint("10_models_stage2")
fit_dtm <- s2$fit_dtm_s2

log_subsection("Reconstructing canonical holdout")

stage2_frac <- as.numeric(Sys.getenv("STAGE2_FRAC", "0.33"))
SITES_EXCLUDE_CHM_16 <- c("1", "2", "3")
SITES_EXCLUDE_DTM_18 <- c("10")

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

chm_df_full$original_row <- seq_len(nrow(chm_df_full))
dtm_df_full$original_row <- seq_len(nrow(dtm_df_full))

set.seed(2025)
chm_train <- chm_df_full %>% group_by(site) %>%
  sample_frac(stage2_frac) %>% ungroup()
dtm_train <- dtm_df_full %>% group_by(site) %>%
  sample_frac(stage2_frac) %>% ungroup()

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

# Apply CHM 16-site and DTM 18-site exclusions to match fit scope.
mod_chm_16 <- chm_holdout %>%
  filter(!site %in% SITES_EXCLUDE_CHM_16) %>%
  mutate(site       = droplevels(site),
         lc_l1_code = droplevels(lc_l1_code),
         ecoregion  = droplevels(ecoregion))

mod_dtm <- dtm_holdout %>%
  filter(!site %in% SITES_EXCLUDE_DTM_18) %>%
  mutate(site       = droplevels(site),
         lc_l1_code = droplevels(lc_l1_code),
         ecoregion  = droplevels(ecoregion))

# Add eco_on if formulas expect it.
if ("eco_on" %in% all.vars(formula(fit_chm_16))) mod_chm_16$eco_on <- 1
if ("eco_on" %in% all.vars(formula(fit_dtm)))    mod_dtm$eco_on    <- 1

log_progress(sprintf("  CHM 16-site holdout rows: %s",
                     formatC(nrow(mod_chm_16), big.mark = ",", format = "d")))
log_progress(sprintf("  DTM 18-site holdout rows: %s",
                     formatC(nrow(mod_dtm), big.mark = ",", format = "d")))

# ---------------------------------------------------------------------
# 2. Define forest functional type strata
# Canonical LC codes per story-lock §2C: ENF, BDF, DNF, EBF.
# ---------------------------------------------------------------------
LC_LABELS <- c(
  BDF = "Broadleaf Deciduous",
  DNF = "Deciduous Needleleaf",
  EBF = "Evergreen Broadleaf",
  ENF = "Evergreen Needleleaf"
)
LC_ORDER <- c("BDF", "DNF", "EBF", "ENF")

# Robustly locate the LC column. The cluster's mod_chm_s2 / mod_dtm_s2
# may carry BOTH `lc_l1_code` and `lc2022_l1_code`. The canonical
# short-code column (BDF/DNF/EBF/ENF) is `lc_l1_code` (used by
# generate_ppc_overall.R and ppc_by_landcover.R). The 2022
# variant may exist with numeric/full-name values that don't match
# LC_ORDER. Pick the column by VALIDATING its values, not just by name.
LC_COL_CANDIDATES <- c("lc_l1_code", "lc2022_l1_code", "lc_l1", "lc")

pick_lc_col <- function(df, label = "df") {
  for (c in LC_COL_CANDIDATES) {
    if (c %in% names(df)) {
      vals <- unique(as.character(df[[c]]))
      n_match <- sum(LC_ORDER %in% vals)
      if (n_match >= 2L) {
        log_progress(sprintf("  [%s] Picked LC column '%s' (%d/%d expected forest codes present)",
                             label, c, n_match, length(LC_ORDER)))
        log_progress(sprintf("  [%s]   Sample values: %s",
                             label, paste(head(vals, 10), collapse = ", ")))
        return(c)
      }
      log_progress(sprintf("  [%s] Column '%s' exists but values don't match LC_ORDER (sample: %s); trying next",
                           label, c, paste(head(vals, 6), collapse = ", ")))
    }
  }
  stop("Could not find a usable LC column in ", label,
       ". Tried: ", paste(LC_COL_CANDIDATES, collapse = ", "),
       ". Expected at least 2 of these codes to appear: ",
       paste(LC_ORDER, collapse = ", "),
       ". Available columns: ", paste(names(df), collapse = ", "))
}

lc_chm <- pick_lc_col(mod_chm_16, "CHM 16-site")
lc_dtm <- pick_lc_col(mod_dtm,    "DTM 18-site")

# ---------------------------------------------------------------------
# 3. Sample holdout stratified by LC, then draw PPCs.
# Per v04 convention: ~5,000 observations per stratum (matches v04
# panel sample sizes).
# ---------------------------------------------------------------------
log_subsection("Stratified sampling and PPC draws")

set.seed(42L)
N_PER_STRATUM <- 5000L
N_PPC_DRAWS   <- 50L

panel_data_for <- function(df, lc_col, fit, product_label) {
  # Filter to known forest LCs.
  df_f <- df[df[[lc_col]] %in% LC_ORDER, ]
  out_panels <- list()
  for (lc in LC_ORDER) {
    sub <- df_f[df_f[[lc_col]] == lc, ]
    n_avail <- nrow(sub)
    if (n_avail == 0) {
      log_progress(sprintf("  %s %s: no rows; skipping",
                           product_label, lc))
      next
    }
    n_take <- min(N_PER_STRATUM, n_avail)
    idx <- sample.int(n_avail, n_take)
    sub2 <- sub[idx, , drop = FALSE]
    # newdata predictions, including site-level RE (allow_new_levels FALSE).
    # We want PPC against ACTUAL holdout y values; use posterior_predict
    # on the same rows. brms recovers y from data, but for PPC we just
    # need yrep matrix; bayesplot::ppc_dens_overlay takes y and yrep.
    yrep <- posterior_predict(fit, newdata = sub2, ndraws = N_PPC_DRAWS,
                              allow_new_levels = TRUE,
                              re_formula = NULL)
    # Response column: CHM uses chm_error_mean, DTM uses dtm_error_mean
    # (per section_08_model_prep.R lines 60, 70 and section_10 formula).
    response_col <- if (product_label == "CHM") "chm_error_mean" else "dtm_error_mean"
    y    <- as.numeric(sub2[[response_col]])

    # Defensive NA handling. posterior_predict can return NA for a row
    # if any of that row's predictor values are NA; bayesplot's
    # validate_predictions() then errors with "NAs not allowed in
    # predictions". Drop NA-prediction columns and NA-response values
    # in matched fashion before passing to ppc_dens_overlay.
    ok_pred <- !apply(is.na(yrep), 2, any)
    ok_resp <- !is.na(y)
    ok      <- ok_pred & ok_resp
    n_drop  <- sum(!ok)
    if (n_drop > 0L) {
      log_progress(sprintf("  %s %s: dropping %d / %d rows with NA pred or response (%d remain)",
                           product_label, lc, n_drop, length(ok), sum(ok)))
      yrep <- yrep[, ok, drop = FALSE]
      y    <- y[ok]
    }
    if (length(y) < 50L) {
      log_progress(sprintf("  %s %s: fewer than 50 usable rows after NA filter; skipping panel",
                           product_label, lc))
      next
    }
    n_take <- length(y)
    out_panels[[lc]] <- list(y = y, yrep = yrep, n = n_take, lc = lc,
                              product = product_label)
  }
  out_panels
}

chm_panels <- panel_data_for(mod_chm_16, lc_chm, fit_chm_16, "CHM")
dtm_panels <- panel_data_for(mod_dtm,    lc_dtm, fit_dtm,    "DTM")

# ---------------------------------------------------------------------
# 4. Build panels via bayesplot, customize, then patchwork.
# ---------------------------------------------------------------------
build_panel <- function(pd, panel_label, x_range = c(-25, 25),
                        product_color = "#1f78b4") {
  bayesplot::color_scheme_set(
    if (pd$product == "CHM") "blue" else "orange")
  pp <- bayesplot::ppc_dens_overlay(y = pd$y, yrep = pd$yrep) +
    coord_cartesian(xlim = x_range) +
    labs(
      title = sprintf("%s  %s  (n = %s)",
                      panel_label,
                      LC_LABELS[[pd$lc]],
                      formatC(pd$n, big.mark = ",", format = "d")),
      x = sprintf("%s error (m)", pd$product),
      y = "Density"
    ) +
    theme_section_G(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold"),
          legend.position = "none")
  pp
}

# CHM row: A, B, C, D for BDF, DNF, EBF, ENF.
labels_chm <- c("(a)", "(b)", "(c)", "(d)")
labels_dtm <- c("(e)", "(f)", "(g)", "(h)")

chm_plots <- vector("list", length(LC_ORDER))
dtm_plots <- vector("list", length(LC_ORDER))
for (i in seq_along(LC_ORDER)) {
  lc <- LC_ORDER[i]
  if (!is.null(chm_panels[[lc]])) {
    chm_plots[[i]] <- build_panel(chm_panels[[lc]], labels_chm[i])
  } else {
    chm_plots[[i]] <- patchwork::plot_spacer()
  }
  if (!is.null(dtm_panels[[lc]])) {
    dtm_plots[[i]] <- build_panel(dtm_panels[[lc]], labels_dtm[i])
  } else {
    dtm_plots[[i]] <- patchwork::plot_spacer()
  }
}

panel_layout <- (chm_plots[[1]] | chm_plots[[2]] | chm_plots[[3]] | chm_plots[[4]]) /
                (dtm_plots[[1]] | dtm_plots[[2]] | dtm_plots[[3]] | dtm_plots[[4]]) +
  patchwork::plot_annotation(
    tag_levels = list(c(""))
  )

# Add CHM/DTM row labels using cowplot draw_label (since patchwork
# doesn't expose row-strip labels cleanly).
final_plot <- cowplot::ggdraw(panel_layout) +
  cowplot::draw_label("CHM", x = 0.005, y = 0.75, angle = 90,
                      size = 11, fontface = "bold") +
  cowplot::draw_label("DTM", x = 0.005, y = 0.27, angle = 90,
                      size = 11, fontface = "bold")

# ---------------------------------------------------------------------
# 5. Save and verify
# ---------------------------------------------------------------------
save_figure(final_plot, "fig04_ppc_by_forest_type",
            width_in = 11.0, height_in = 5.6)

log_subsection("VERIFICATION (console anchors)")
cat("\nSample sizes per panel (target n = 5,000 if stratum has >=5,000;\n")
cat("else the full stratum holdout size). Holdout pool is ~67% of the\n")
cat("full data (after site exclusions). Story-lock §2C 16-site forest-\n")
cat("subset full counts: ENF 193,453; BDF 44,394; DNF 9,035; EBF 2,266.\n")
cat("Approximate 67%% holdout sizes: ENF ~129,514; BDF ~29,744;\n")
cat("DNF ~6,054; EBF ~1,518.\n\n")

for (lc in LC_ORDER) {
  ck <- chm_panels[[lc]]
  dk <- dtm_panels[[lc]]
  cat(sprintf("  %s  CHM n=%s  DTM n=%s\n",
              lc,
              if (is.null(ck)) "n/a" else formatC(ck$n, big.mark = ",", format = "d"),
              if (is.null(dk)) "n/a" else formatC(dk$n, big.mark = ",", format = "d")))
}

log_progress("fig4_ppc_by_forest_type.R complete.")
