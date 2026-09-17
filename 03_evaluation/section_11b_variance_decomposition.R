# =====================================================================
# section_11b_variance_decomposition.R
# Three-level random effect variance decomposition
#
# PURPOSE:
#   Decompose random-effect variance across site, ecoregion, and land
#   cover grouping levels for both CHM and DTM models. Produces:
#     1. Hyperparameter-based decomposition (intercept SDs only)
#     2. Observation-level decomposition accounting for random slopes
#     3. Posterior summaries with credible intervals
#     4. Export tables for manuscript (Table 5 revision)
#
# REQUIRES:
#   - Stage 2 model objects (fit_chm_s2, fit_dtm_s2)
#   - Loaded via checkpoint 10_models_stage2 or already in environment
#
# OUTPUTS:
#   tables/variance_decomposition_hyperparameter.csv
#   tables/variance_decomposition_observation_level.csv
#   tables/variance_decomposition_full_posterior.csv
#   tables/variance_decomposition_manuscript.csv
#   plots/11b_variance_decomposition.pdf
#   plots/11b_variance_decomposition_bar.pdf
# =====================================================================

# --- Dependencies ----------------------------------------------------
if (!exists("log_progress")) source("analysis_utils.R")
if (!exists("PROJECT_ROOT")) source("analysis_config.R")

suppressPackageStartupMessages({
  library(brms)
  library(tidybayes)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(cowplot)
})

# --- Load models -----------------------------------------------------
if (!exists("fit_chm_s2") || !exists("fit_dtm_s2")) {
  log_progress("Loading Stage 2 models from checkpoint...")
  if (checkpoint_exists("10_models_stage2")) {
    cp <- load_checkpoint("10_models_stage2")
    fit_chm_s2 <- cp$fit_chm_s2
    fit_dtm_s2 <- cp$fit_dtm_s2
  } else {
    stop("Checkpoint 10_models_stage2 not found. Run section_10 first.")
  }
}

log_progress("=" %>% strrep(60))
log_progress("THREE-LEVEL VARIANCE DECOMPOSITION")
log_progress("=" %>% strrep(60))

# =====================================================================
# APPROACH 1: Hyperparameter-based (intercept SDs from posterior)
# =====================================================================
#
# For each posterior draw, extract:
#   sigma_site      = sd(Intercept) for site grouping
#   sigma_ecoregion = sd(Intercept) for ecoregion grouping
#   sigma_lc        = sd(Intercept) for lc_l1_code grouping
#
# Then:  prop_k = sigma_k^2 / sum(sigma_k^2)
#
# NOTE: This uses intercept SDs only. Ecoregion has a random slope
# on slope_mean_z and lc_l1_code has a random slope on wsci_z;
# those slope contributions are captured in Approach 2.
# =====================================================================

log_progress("--- Approach 1: Hyperparameter-based (intercept SDs) ---")

extract_intercept_sd_draws <- function(fit, label) {
  # Get all SD parameter draws from the mu submodel
  # brms names these as: sd_<group>__Intercept
  draws <- as_draws_df(fit)
  
  # Identify SD columns for intercepts in the mu submodel
  sd_cols <- grep("^sd_.*__Intercept$", names(draws), value = TRUE)
  # Exclude sigma-submodel SDs (prefixed with sigma_)
  sd_cols <- sd_cols[!grepl("^sd_sigma_", sd_cols)]
  
  log_progress(sprintf("  [%s] SD intercept parameters found: %s",
                        label, paste(sd_cols, collapse = ", ")))
  
  # Extract and rename
  sd_draws <- draws[, sd_cols, drop = FALSE]
  
  # Parse group names from column names
  # Pattern: sd_<groupname>__Intercept
  group_names <- sub("^sd_", "", sub("__Intercept$", "", sd_cols))
  names(sd_draws) <- group_names
  
  as_tibble(sd_draws)
}

compute_hyperpar_decomp <- function(sd_draws) {
  # For each posterior draw, compute variance proportions
  var_draws <- sd_draws^2
  total_var <- rowSums(var_draws)
  
  prop_draws <- var_draws / total_var
  names(prop_draws) <- paste0("prop_", names(prop_draws))
  
  # Combine variance and proportion draws
  out <- bind_cols(
    var_draws %>% rename_with(~ paste0("var_", .x)),
    prop_draws,
    tibble(total_re_var = total_var)
  )
  out
}

summarize_decomp <- function(decomp_draws, label) {
  # Posterior summary: median [95% CI]
  cols <- names(decomp_draws)
  
  summary_list <- lapply(cols, function(col) {
    x <- decomp_draws[[col]]
    tibble(
      parameter = col,
      median = median(x),
      mean = mean(x),
      q025 = quantile(x, 0.025),
      q975 = quantile(x, 0.975),
      sd = sd(x)
    )
  })
  
  summary_df <- bind_rows(summary_list)
  summary_df$product <- label
  summary_df
}

# CHM
sd_draws_chm <- extract_intercept_sd_draws(fit_chm_s2, "CHM")
decomp_chm <- compute_hyperpar_decomp(sd_draws_chm)
summary_hyp_chm <- summarize_decomp(decomp_chm, "CHM")

# DTM
sd_draws_dtm <- extract_intercept_sd_draws(fit_dtm_s2, "DTM")
decomp_dtm <- compute_hyperpar_decomp(sd_draws_dtm)
summary_hyp_dtm <- summarize_decomp(decomp_dtm, "DTM")

# Combine and export
summary_hyp <- bind_rows(summary_hyp_chm, summary_hyp_dtm)

log_progress("  Hyperparameter-based variance proportions (intercepts only):")
summary_hyp %>%
  filter(grepl("^prop_", parameter)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  select(product, parameter, median, q025, q975) %>%
  print(n = 20)

write_csv(summary_hyp, file.path(out_tables, "variance_decomposition_hyperparameter.csv"))
log_progress("  -> Exported: variance_decomposition_hyperparameter.csv")


# =====================================================================
# APPROACH 2: Observation-level decomposition (accounts for slopes)
# =====================================================================
#
# For each posterior draw s and each observation i, compute the
# group-level contribution to the linear predictor:
#
#   re_site[s,i]      = u_site_intercept[s, site[i]]
#   re_ecoregion[s,i] = u_eco_intercept[s, eco[i]]
#                        + u_eco_slope[s, eco[i]] * slope_mean_z[i]
#   re_lc[s,i]        = u_lc_intercept[s, lc[i]]
#                        + u_lc_slope[s, lc[i]] * wsci_z[i]
#
# Then for draw s:
#   V_site[s]      = Var_i( re_site[s,i] )
#   V_ecoregion[s] = Var_i( re_ecoregion[s,i] )
#   V_lc[s]        = Var_i( re_lc[s,i] )
#
# Proportion:  prop_k[s] = V_k[s] / (V_site + V_eco + V_lc)[s]
#
# This is the most complete decomposition because it accounts for
# how random slopes interact with the actual predictor distribution.
# =====================================================================

log_progress("--- Approach 2: Observation-level (with random slopes) ---")

obs_level_decomp <- function(fit, label, n_draws = 500) {
  
  # --- Extract model data ---
  dat <- fit$data
  n_obs <- nrow(dat)
  
  # Group membership vectors
  sites <- as.character(dat$site)
  ecoregions <- as.character(dat$ecoregion)
  lc_codes <- as.character(dat$lc_l1_code)
  
  # Predictor values for random slopes
  slope_z <- dat$slope_mean_z
  wsci_z <- dat$wsci_z
  
  unique_sites <- sort(unique(sites))
  unique_ecos <- sort(unique(ecoregions))
  unique_lcs <- sort(unique(lc_codes))
  
  log_progress(sprintf("  [%s] n_obs=%d, n_sites=%d, n_eco=%d, n_lc=%d",
                        label, n_obs, length(unique_sites),
                        length(unique_ecos), length(unique_lcs)))
  
  # --- Extract posterior draws of group-level effects ---
  # Convert to plain data.frame immediately — draws_df class has special

  # subsetting semantics that break standard [,] and [[]] indexing below
  draws <- as.data.frame(as_draws_df(fit))
  n_total_draws <- nrow(draws)
  
  # Subsample draws for computational tractability
  if (n_total_draws > n_draws) {
    draw_idx <- sort(sample(n_total_draws, n_draws))
  } else {
    draw_idx <- seq_len(n_total_draws)
    n_draws <- n_total_draws
  }
  
  log_progress(sprintf("  Using %d posterior draws (of %d total)",
                        n_draws, n_total_draws))
  
  # --- Helper: extract group-level draws ---
  # brms stores group effects as r_<group>[<level>,<coef>]
  get_re_matrix <- function(draws, group, coef, levels) {
    # Returns n_draws x n_levels matrix
    cols <- grep(sprintf("^r_%s\\[", group), names(draws), value = TRUE)
    cols <- cols[grepl(paste0(",", coef, "\\]$"), cols)]
    
    if (length(cols) == 0) {
      log_progress(sprintf("    WARNING: No columns found for r_%s[*,%s]", group, coef))
      return(matrix(0, nrow = n_draws, ncol = length(levels)))
    }
    
    # Extract level names from column names
    col_levels <- sub(sprintf("^r_%s\\[", group), "",
                       sub(sprintf(",%s\\]$", coef), "", cols))
    
    # Build matrix aligned to 'levels' ordering (use numeric indexing
    # to avoid character-subscript edge cases with special level names)
    mat <- matrix(0, nrow = n_draws, ncol = length(levels))
    colnames(mat) <- levels
    
    for (j in seq_along(cols)) {
      lev <- col_levels[j]
      col_idx <- match(lev, levels)
      if (!is.na(col_idx)) {
        vals <- as.numeric(draws[[cols[j]]])
        mat[, col_idx] <- vals[draw_idx]
      }
    }
    mat
  }
  
  # Site: intercept only
  re_site_int <- get_re_matrix(draws, "site", "Intercept", unique_sites)
  
  # Ecoregion: intercept + slope on slope_mean_z
  re_eco_int <- get_re_matrix(draws, "ecoregion", "Intercept", unique_ecos)
  re_eco_slope <- get_re_matrix(draws, "ecoregion", "slope_mean_z", unique_ecos)
  
  # Land cover: intercept + slope on wsci_z
  re_lc_int <- get_re_matrix(draws, "lc_l1_code", "Intercept", unique_lcs)
  re_lc_slope <- get_re_matrix(draws, "lc_l1_code", "wsci_z", unique_lcs)
  
  # --- Compute observation-level RE contributions per draw ---
  # Pre-compute index vectors (1-based into unique_* vectors)
  site_idx <- match(sites, unique_sites)
  eco_idx <- match(ecoregions, unique_ecos)
  lc_idx <- match(lc_codes, unique_lcs)
  
  log_progress("  Computing observation-level variance decomposition...")
  
  # Preallocate results
  V_site <- numeric(n_draws)
  V_eco <- numeric(n_draws)
  V_lc <- numeric(n_draws)
  
  pb <- progress::progress_bar$new(
    format = "  Draw :current/:total [:bar] :percent eta: :eta",
    total = n_draws, clear = FALSE, width = 60
  )
  
  for (s in seq_len(n_draws)) {
    pb$tick()
    
    # Site contribution: intercept only
    contrib_site <- re_site_int[s, site_idx]
    
    # Ecoregion contribution: intercept + slope * slope_z
    contrib_eco <- re_eco_int[s, eco_idx] + re_eco_slope[s, eco_idx] * slope_z
    
    # Land cover contribution: intercept + slope * wsci_z
    contrib_lc <- re_lc_int[s, lc_idx] + re_lc_slope[s, lc_idx] * wsci_z
    
    # Compute variances across observations
    V_site[s] <- var(contrib_site)
    V_eco[s] <- var(contrib_eco)
    V_lc[s] <- var(contrib_lc)
  }
  
  # --- Assemble results ---
  total_V <- V_site + V_eco + V_lc
  
  tibble(
    product = label,
    draw = seq_len(n_draws),
    var_site = V_site,
    var_ecoregion = V_eco,
    var_lc = V_lc,
    var_total_re = total_V,
    prop_site = V_site / total_V,
    prop_ecoregion = V_eco / total_V,
    prop_lc = V_lc / total_V
  )
}

# Run decomposition (500 draws is usually plenty for stable summaries)
N_DRAWS <- as.integer(Sys.getenv("VARDECOMP_NDRAWS", "500"))

obs_chm <- obs_level_decomp(fit_chm_s2, "CHM", n_draws = N_DRAWS)
obs_dtm <- obs_level_decomp(fit_dtm_s2, "DTM", n_draws = N_DRAWS)

obs_all <- bind_rows(obs_chm, obs_dtm)

# Summarize
obs_summary <- obs_all %>%
  group_by(product) %>%
  summarize(
    # Variances
    var_site_med     = median(var_site),
    var_site_q025    = quantile(var_site, 0.025),
    var_site_q975    = quantile(var_site, 0.975),
    
    var_eco_med      = median(var_ecoregion),
    var_eco_q025     = quantile(var_ecoregion, 0.025),
    var_eco_q975     = quantile(var_ecoregion, 0.975),
    
    var_lc_med       = median(var_lc),
    var_lc_q025      = quantile(var_lc, 0.025),
    var_lc_q975      = quantile(var_lc, 0.975),
    
    var_total_med    = median(var_total_re),
    var_total_q025   = quantile(var_total_re, 0.025),
    var_total_q975   = quantile(var_total_re, 0.975),
    
    # Proportions
    prop_site_med    = median(prop_site),
    prop_site_q025   = quantile(prop_site, 0.025),
    prop_site_q975   = quantile(prop_site, 0.975),
    
    prop_eco_med     = median(prop_ecoregion),
    prop_eco_q025    = quantile(prop_ecoregion, 0.025),
    prop_eco_q975    = quantile(prop_ecoregion, 0.975),
    
    prop_lc_med      = median(prop_lc),
    prop_lc_q025     = quantile(prop_lc, 0.025),
    prop_lc_q975     = quantile(prop_lc, 0.975),
    
    .groups = "drop"
  )

log_progress("\n  Observation-level variance proportions (with random slopes):")
obs_summary %>%
  select(product, starts_with("prop_")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  print()

write_csv(obs_summary, file.path(out_tables, "variance_decomposition_observation_level.csv"))
log_progress("  -> Exported: variance_decomposition_observation_level.csv")

# Full posterior draws for downstream use
write_csv(obs_all, file.path(out_tables, "variance_decomposition_full_posterior.csv"))
log_progress("  -> Exported: variance_decomposition_full_posterior.csv")


# =====================================================================
# MANUSCRIPT-READY SUMMARY TABLE
# =====================================================================

log_progress("--- Manuscript summary ---")

manuscript_table <- obs_summary %>%
  transmute(
    Product = product,
    
    `Site var (m^2)` = sprintf("%.2f [%.2f, %.2f]",
                               var_site_med, var_site_q025, var_site_q975),
    `Site %` = sprintf("%.1f [%.1f, %.1f]",
                        prop_site_med * 100, prop_site_q025 * 100, prop_site_q975 * 100),
    
    `Ecoregion var (m^2)` = sprintf("%.2f [%.2f, %.2f]",
                                    var_eco_med, var_eco_q025, var_eco_q975),
    `Ecoregion %` = sprintf("%.1f [%.1f, %.1f]",
                             prop_eco_med * 100, prop_eco_q025 * 100, prop_eco_q975 * 100),
    
    `Land cover var (m^2)` = sprintf("%.2f [%.2f, %.2f]",
                                     var_lc_med, var_lc_q025, var_lc_q975),
    `Land cover %` = sprintf("%.1f [%.1f, %.1f]",
                              prop_lc_med * 100, prop_lc_q025 * 100, prop_lc_q975 * 100),
    
    `Total RE var (m^2)` = sprintf("%.2f [%.2f, %.2f]",
                                   var_total_med, var_total_q025, var_total_q975)
  )

cat("\n")
log_progress("Manuscript-ready variance decomposition:")
print(manuscript_table, width = 200)

write_csv(manuscript_table, file.path(out_tables, "variance_decomposition_manuscript.csv"))
log_progress("  -> Exported: variance_decomposition_manuscript.csv")


# =====================================================================
# FIGURE: Posterior distributions of variance proportions
# =====================================================================

log_progress("--- Generating variance decomposition figure ---")

plot_data <- obs_all %>%
  select(product, draw, prop_site, prop_ecoregion, prop_lc) %>%
  pivot_longer(
    cols = starts_with("prop_"),
    names_to = "level",
    values_to = "proportion"
  ) %>%
  mutate(
    level = recode(level,
      prop_site = "Site",
      prop_ecoregion = "Ecoregion",
      prop_lc = "Land Cover"
    ),
    level = factor(level, levels = c("Land Cover", "Ecoregion", "Site"))
  )

# Summarize: median + 50% CI + 95% CI
plot_summary <- plot_data %>%
  group_by(product, level) %>%
  summarize(
    median = median(proportion),
    q025   = quantile(proportion, 0.025),
    q975   = quantile(proportion, 0.975),
    q25    = quantile(proportion, 0.25),
    q75    = quantile(proportion, 0.75),
    .groups = "drop"
  )

p_decomp <- ggplot(plot_summary, aes(y = level, x = median, color = level)) +
  # 95% CI (thin line)
  geom_linerange(aes(xmin = q025, xmax = q975), linewidth = 0.6) +
  # 50% CI (thick line)
  geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 2) +
  # Median point
  geom_point(size = 3, shape = 21, fill = "white", stroke = 1.2) +
  facet_wrap(~ product, ncol = 1) +
  scale_x_log10(
    labels = function(x) ifelse(x >= 0.01,
                                 scales::percent(x, accuracy = 1),
                                 scales::percent(x, accuracy = 0.001)),
    breaks = c(0.0001, 0.001, 0.01, 0.1, 1),
    limits = c(1e-8, 1)
  ) +
  scale_color_manual(
    values = c("Site" = "#2166AC", "Ecoregion" = "#B2182B", "Land Cover" = "#4DAF4A"),
    guide = "none"
  ) +
  annotation_logticks(sides = "b", size = 0.3) +
  labs(
    x = "Proportion of random-effect variance (log scale)",
    y = NULL,
    title = "Variance decomposition across hierarchical levels",
    subtitle = "Posterior median with 50% (thick) and 95% (thin) credible intervals"
  ) +
  theme_cowplot() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
  )

ggsave(
  file.path(out_plots, "11b_variance_decomposition.pdf"),
  p_decomp, width = 8, height = 4.5, bg = "white"
)
log_progress("  -> Saved: 11b_variance_decomposition.pdf")

# Stacked bar version for quick reference
median_props <- obs_all %>%
  group_by(product) %>%
  summarize(
    Site = median(prop_site),
    Ecoregion = median(prop_ecoregion),
    `Land Cover` = median(prop_lc),
    .groups = "drop"
  ) %>%
  pivot_longer(-product, names_to = "level", values_to = "proportion") %>%
  mutate(level = factor(level, levels = c("Land Cover", "Ecoregion", "Site")))

p_bar <- ggplot(median_props, aes(x = product, y = proportion, fill = level)) +
  geom_col(width = 0.5) +
  geom_text(
    aes(label = sprintf("%.1f%%", proportion * 100)),
    position = position_stack(vjust = 0.5),
    size = 4, fontface = "bold"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = c("Site" = "#2166AC", "Ecoregion" = "#B2182B", "Land Cover" = "#4DAF4A"),
    name = "Grouping Level"
  ) +
  labs(
    x = NULL, y = "Proportion of random-effect variance",
    title = "Random-effect variance partitioning (posterior medians)"
  ) +
  theme_cowplot() +
  theme(legend.position = "right")

ggsave(
  file.path(out_plots, "11b_variance_decomposition_bar.pdf"),
  p_bar, width = 7, height = 5, bg = "white"
)
log_progress("  -> Saved: 11b_variance_decomposition_bar.pdf")


# =====================================================================
# CONSISTENCY CHECK: Compare approaches
# =====================================================================

log_progress("--- Approach comparison ---")

for (prod in c("CHM", "DTM")) {
  hyp <- summary_hyp %>% filter(product == prod, grepl("^prop_", parameter))
  obs <- obs_summary %>% filter(product == prod)
  
  cat(sprintf("\n  [%s] Hyperparameter (intercept-only) vs Observation-level:\n", prod))
  
  for (grp in c("site", "ecoregion", "lc_l1_code")) {
    grp_label <- switch(grp,
      site = "site", ecoregion = "ecoregion", lc_l1_code = "lc"
    )
    hyp_val <- hyp %>% filter(parameter == paste0("prop_", grp)) %>% pull(median)
    obs_val <- obs[[paste0("prop_", grp_label, "_med")]]
    
    if (length(hyp_val) > 0 && length(obs_val) > 0) {
      cat(sprintf("    %-12s  hyperpar: %.1f%%   obs-level: %.1f%%\n",
                  grp, hyp_val * 100, obs_val * 100))
    }
  }
}

cat("\n")
log_progress("NOTE: Differences between approaches reflect random slope contributions.")
log_progress("      The observation-level approach (Approach 2) is recommended for the")
log_progress("      manuscript because it accounts for slope_mean_z and wsci_z random")
log_progress("      slopes in the ecoregion and land cover groupings, respectively.")


# =====================================================================
# CHECKPOINT
# =====================================================================

if (!exists("BATCH_MODE") || !BATCH_MODE) {
  save_checkpoint("11b_variance_decomposition", list(
    obs_all = obs_all,
    obs_summary = obs_summary,
    summary_hyp = summary_hyp,
    manuscript_table = manuscript_table
  ))
}

log_progress("=" %>% strrep(60))
log_progress("VARIANCE DECOMPOSITION COMPLETE")
log_progress("=" %>% strrep(60))
