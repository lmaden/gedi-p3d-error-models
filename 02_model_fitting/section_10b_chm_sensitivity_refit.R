#!/usr/bin/env Rscript
# =====================================================================
# section_10b_chm_sensitivity_refit.R
#
# CHM-only sensitivity refit excluding the three sites excluded for
# systematic reference-CHM underestimation:
#   Site 1 (usda_me)       - Feb 2015 leaf-off flight; CHM RE +5.06 m
#   Site 2 (nasa_howland)  - Oct 2013 flight;          CHM RE +7.26 m
#   Site 3 (neon_sawb)     - 2013 early-NEON flight;   CHM RE +10.75 m
#
# All three show large negative als_chm_p90 - rh_98 offsets (-8, -8,
# -16 m) that co-rank with the CHM site random intercepts at r = -0.90
# across the 19 sites, indicating the reference CHM rasters
# systematically underestimate canopy at these sites and inflate the
# measured P3D - ALS error.
#
# This script refits only the CHM Stage 2 model on the 16 non-flagged
# sites, reusing the exact formula and priors stored in the stage-2
# checkpoint, and hardcoded sampler settings matched to section_10.
# It then produces side-by-side comparison tables (global parameters
# and per-predictor fixed effects) with the full 19-site fit.
#
# Inputs:
#   /gpfs/data1/vclgp/lmaden/chpt1/checkpoints/08_model_prep.rds
#   /gpfs/data1/vclgp/lmaden/chpt1/checkpoints/10_models_stage2.rds
#
# Outputs:
#   checkpoints/10b_chm_sensitivity.rds
#   manuscript_tables/chm_sensitivity_comparison.csv
#   manuscript_tables/chm_sensitivity_fixef_comparison.csv
#
# Config:
#   Sources analysis_config.R (matches rest of pipeline on cluster)
# =====================================================================

suppressPackageStartupMessages({
  library(brms)
  library(data.table)
  library(posterior)
  library(loo)
})

# Source project infrastructure (provides CHECKPOINT_DIR, MODELS_DIR,
# CPU_BUDGET, log_progress, save_checkpoint, load_checkpoint, etc.).
# Uses the _FIXED config that matches the rest of the pipeline on the
# cluster (section_15 also sources _FIXED; upstream checkpoints we load
# were produced by scripts running against this config).
if (!exists("log_progress"))  source("analysis_utils.R")
if (!exists("PROJECT_ROOT"))  source("analysis_config.R")

# Sanity-check required config variables up front. Failing here with a
# clear message beats failing deep inside brm() or fwrite() later.
.required_vars <- c("PROJECT_ROOT", "CPU_BUDGET", "RAM_FRAC", "MODELS_DIR")
.missing_vars  <- .required_vars[!vapply(.required_vars, exists, logical(1))]
if (length(.missing_vars) > 0) {
  stop("Missing required config variables: ",
       paste(.missing_vars, collapse = ", "),
       ". Check analysis_config.R.")
}

# 16-site outputs go to manuscript_tables/ (per 16-site checkpoint),
# NOT to the pipeline's out_tables (which resolves to tables/).
manuscript_tables_dir <- file.path(PROJECT_ROOT, "manuscript_tables")

SITES_TO_EXCLUDE <- c("1", "2", "3")

log_progress("================================================================")
log_progress("CHM sensitivity refit (16 sites, excluding the three flagged)")
log_progress(sprintf("  Excluded: %s", paste(SITES_TO_EXCLUDE, collapse = ", ")))
log_progress("================================================================")

# ---------------------------------------------------------------------
# Load upstream checkpoints
# ---------------------------------------------------------------------

log_progress("Loading stage-2 checkpoint (for formula, priors, full fit)...")
s2 <- load_checkpoint("10_models_stage2")
stopifnot(
  "fit_chm_s2"     %in% names(s2),
  "chm_formula_s2" %in% names(s2)
)
fit_chm_full   <- s2$fit_chm_s2
chm_formula_s2 <- s2$chm_formula_s2

log_progress("Loading model-prep checkpoint (for mod_chm_s2)...")
mp <- load_checkpoint("08_model_prep")
stopifnot("mod_chm_s2" %in% names(mp))
mod_chm_s2 <- mp$mod_chm_s2

# ---------------------------------------------------------------------
# Filter to 16 sites
# ---------------------------------------------------------------------
# The `site` column in mod_chm_s2 holds character site codes (confirmed
# in section_08 line 66: `site = factor(site)` where `site` is the
# character column from ingest). Filter by name, not by tracker ID.

keep_rows <- !as.character(mod_chm_s2$site) %in% SITES_TO_EXCLUDE
mod_chm_s2_16 <- mod_chm_s2[keep_rows, ]

# Drop unused factor levels so the model does not allocate RE slots
# for zero-row sites.
mod_chm_s2_16$site <- droplevels(factor(mod_chm_s2_16$site))

# Hard sanity check: the filter must have dropped exactly the three
# target sites. This is the guard against the v1 silent-no-op bug.
remaining_sites <- sort(unique(as.character(mod_chm_s2_16$site)))
stopifnot(
  length(remaining_sites) == 16L,
  !any(SITES_TO_EXCLUDE %in% remaining_sites)
)

log_progress(sprintf("  Full dataset  : %s rows across %d sites",
                     format(nrow(mod_chm_s2),   big.mark = ","),
                     length(unique(mod_chm_s2$site))))
log_progress(sprintf("  16-site subset: %s rows across %d sites",
                     format(nrow(mod_chm_s2_16), big.mark = ","),
                     length(unique(mod_chm_s2_16$site))))
log_progress(sprintf("  Rows removed  : %s",
                     format(nrow(mod_chm_s2) - nrow(mod_chm_s2_16),
                            big.mark = ",")))

# ---------------------------------------------------------------------
# RAM estimation (mirrors section_10 lines 73-91)
# ---------------------------------------------------------------------
# Warns if the refit would exceed the share of system RAM allocated by
# RAM_FRAC in analysis_config.R. The 16-site subset is smaller than the
# full dataset, so this should pass unconditionally if section_10 did.

estimate_ram_needed <- function(n_rows, n_chains = 4) {
  base_mb <- 500
  per_row_mb <- 0.15
  per_chain_mb <- base_mb + (n_rows * per_row_mb)
  (per_chain_mb * n_chains) / 1024
}

ram_needed_gb <- estimate_ram_needed(nrow(mod_chm_s2_16), n_chains = 4)
ram_avail_gb  <- tryCatch({
  as.numeric(system("free -g | awk '/^Mem:/ {print $2}'", intern = TRUE)) * RAM_FRAC
}, error = function(e) 16)

log_progress(sprintf("  RAM estimate: need ~%.1f GB, budget ~%.1f GB",
                     ram_needed_gb, ram_avail_gb))
if (ram_needed_gb > ram_avail_gb * 0.8) {
  warning(sprintf("Estimated RAM %.1f GB exceeds 80%% of %.1f GB budget",
                  ram_needed_gb, ram_avail_gb))
}

# ---------------------------------------------------------------------
# Priors and sampler settings (matched to section_10)
# ---------------------------------------------------------------------
# Priors pull reliably from the fitted brmsfit object across backends.
# Sampler control settings DO NOT pull reliably from `fit@sim` on
# cmdstanr-backed fits, so they are hardcoded here to match section_10
# lines 218-229 exactly.

log_progress("Defining priors and control settings (matched to section_10)...")

priors_chm <- fit_chm_full$prior

ctrl_chm <- list(
  adapt_delta   = 0.98,
  max_treedepth = 15,
  step_size     = 0.0005
)

n_chains   <- 4
n_iter     <- 3000
n_warmup   <- 1500

parallel_chains   <- max(1L, min(n_chains, CPU_BUDGET))
threads_per_chain <- if (CPU_BUDGET >= 8) 2L else 1L

use_cmdstanr   <- requireNamespace("cmdstanr", quietly = TRUE)
backend_choice <- if (use_cmdstanr) "cmdstanr" else "rstan"

full_backend <- tryCatch(fit_chm_full$backend, error = function(e) "unknown")
log_progress(sprintf("  Backend: %s (full fit used: %s)",
                     backend_choice, full_backend))
if (!is.na(full_backend) && full_backend != "unknown" &&
    full_backend != backend_choice) {
  warning(sprintf(
    "Backend mismatch: full fit used '%s' but refit will use '%s'. Results are still valid but introduce a minor nuisance source of variation.",
    full_backend, backend_choice))
}
log_progress(sprintf("  Chains: %d parallel, %d thread(s) per chain",
                     parallel_chains, threads_per_chain))
log_progress(sprintf("  Iter: %d (warmup %d)", n_iter, n_warmup))
log_progress(sprintf("  Control: adapt_delta=%.2f, max_treedepth=%d, step_size=%.4f",
                     ctrl_chm$adapt_delta, ctrl_chm$max_treedepth,
                     ctrl_chm$step_size))

# Ensure cache directory exists (analysis_config.R defines the path but
# doesn't necessarily create it).
dir.create(MODELS_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------
# Refit CHM on the 16-site subset
# ---------------------------------------------------------------------

log_progress("Refitting CHM Stage 2 on 16-site subset...")

fit_chm_16 <- brm(
  formula = chm_formula_s2,
  data    = mod_chm_s2_16,
  prior   = priors_chm,
  family  = student(),
  warmup  = n_warmup,
  iter    = n_iter,
  chains  = n_chains,
  cores   = parallel_chains,
  threads = threading(threads_per_chain),
  control = ctrl_chm,
  backend = backend_choice,
  file    = file.path(MODELS_DIR, "fit_chm_s2_16site"),
  seed    = 2025
)

log_progress("  Refit complete")

# Post-fit convergence check (mirrors section_10 lines 262-267)
max_rhat_16 <- max(brms::rhat(fit_chm_16),       na.rm = TRUE)
min_ess_16  <- min(brms::neff_ratio(fit_chm_16),  na.rm = TRUE)
log_progress(sprintf("  Max Rhat: %.3f, Min ESS ratio: %.3f",
                     max_rhat_16, min_ess_16))
if (max_rhat_16 > 1.01) {
  warning("Some 16-site CHM parameters have Rhat > 1.01")
}

# ---------------------------------------------------------------------
# LOO on the refit's own data (not comparable to the full-fit LOO)
# ---------------------------------------------------------------------

log_progress("Computing LOO for the 16-site refit...")
loo_chm_16 <- tryCatch({
  log_progress("  Attempting LOO with moment matching...")
  loo(fit_chm_16, moment_match = TRUE, cores = parallel_chains)
}, error = function(e) {
  log_progress(sprintf("  Moment matching failed: %s", e$message))
  tryCatch({
    log_progress("  Falling back to basic LOO...")
    loo(fit_chm_16, cores = parallel_chains)
  }, error = function(e2) {
    log_progress(sprintf("  LOO failed entirely: %s", e2$message))
    NULL
  })
})

if (!is.null(loo_chm_16)) {
  log_progress(sprintf("  16-site LOO ELPD: %.1f (SE: %.1f)",
                       loo_chm_16$estimates["elpd_loo", "Estimate"],
                       loo_chm_16$estimates["elpd_loo", "SE"]))
}

# ---------------------------------------------------------------------
# Side-by-side global-parameter comparison
# ---------------------------------------------------------------------

extract_summary <- function(fit, label) {
  fx    <- fixef(fit, probs = c(0.025, 0.5, 0.975))
  b0    <- fx["Intercept", ]

  vc    <- VarCorr(fit, probs = c(0.025, 0.5, 0.975))
  sd_site <- vc$site$sd["Intercept", ]

  re       <- ranef(fit)$site[, "Estimate", "Intercept"]
  re_range <- range(re)
  re_sd    <- sd(re)

  # Pull only the two scalar-summary variables we need (memory-efficient
  # vs. pulling all draws into a data frame).
  pars      <- as_draws_df(fit, variable = c("nu", "b_sigma_Intercept"))
  nu_vec    <- if ("nu" %in% names(pars)) pars$nu else NA_real_
  # sigma submodel uses log link: b_sigma_Intercept is on log scale
  sigma_vec <- if ("b_sigma_Intercept" %in% names(pars))
                 exp(pars$b_sigma_Intercept) else NA_real_

  data.table(
    fit                     = label,
    n_obs                   = nrow(fit$data),
    n_sites                 = length(unique(as.character(fit$data$site))),
    intercept_mean          = unname(b0["Estimate"]),
    intercept_lo            = unname(b0["Q2.5"]),
    intercept_hi            = unname(b0["Q97.5"]),
    sd_site_mean            = unname(sd_site["Estimate"]),
    sd_site_lo              = unname(sd_site["Q2.5"]),
    sd_site_hi              = unname(sd_site["Q97.5"]),
    re_site_min             = re_range[1],
    re_site_max             = re_range[2],
    re_site_sd_posteriormean = re_sd,
    nu_median               = if (is.numeric(nu_vec))    median(nu_vec)    else NA_real_,
    sigma_median            = if (is.numeric(sigma_vec)) median(sigma_vec) else NA_real_
  )
}

log_progress("Building global-parameter comparison table...")
cmp <- rbind(
  extract_summary(fit_chm_full, "CHM 19-site (full)"),
  extract_summary(fit_chm_16,   "CHM 16-site (sensitivity)")
)
print(cmp)

# ---------------------------------------------------------------------
# Fixed-effect coefficient comparison
# ---------------------------------------------------------------------

log_progress("Comparing fixed-effect coefficients...")

fx_full <- as.data.table(fixef(fit_chm_full, probs = c(0.025, 0.975)),
                         keep.rownames = "predictor")
fx_16   <- as.data.table(fixef(fit_chm_16,   probs = c(0.025, 0.975)),
                         keep.rownames = "predictor")

fx_join <- merge(
  fx_full[, .(predictor, est_full = Estimate, lo_full = Q2.5, hi_full = Q97.5)],
  fx_16[,   .(predictor, est_16   = Estimate, lo_16   = Q2.5, hi_16   = Q97.5)],
  by = "predictor", all = TRUE
)

fx_join[, shift := est_16 - est_full]

# Per-predictor stability metrics.
#   includes_zero_*     : CI straddles zero (i.e., predictor "not significant")
#   status_preserved    : zero-exclusion status identical across both fits
#                         (this is the metric 16-site's decision rule hinges
#                         on; "stable" = status_preserved for all predictors)
#   cis_overlap_bw_fits : auxiliary, whether the two fits' 95% CIs overlap
#                         each other (secondary diagnostic)
fx_join[, includes_zero_full   := (lo_full <= 0) & (hi_full >= 0)]
fx_join[, includes_zero_16     := (lo_16   <= 0) & (hi_16   >= 0)]
fx_join[, status_preserved     := includes_zero_full == includes_zero_16]
fx_join[, cis_overlap_bw_fits  := !(hi_full < lo_16 | hi_16 < lo_full)]

# Handle predictors present in one fit but not the other (structural
# change from dropping the three flagged sites). The merge with
# all = TRUE leaves NAs on one side; without this block those rows
# would be silently excluded from the decision-rule flip count, which
# is the opposite of what 16-site wants. Treat NA as status-not-preserved.
n_missing_rows <- sum(is.na(fx_join$status_preserved))
if (n_missing_rows > 0) {
  missing_pred <- fx_join[is.na(status_preserved), predictor]
  warning(sprintf(
    "%d predictor(s) present in only one fit: %s. Treating as status-not-preserved.",
    n_missing_rows, paste(missing_pred, collapse = ", ")))
  fx_join[is.na(status_preserved), status_preserved := FALSE]
}

fx_join <- fx_join[order(-abs(shift))]

print(fx_join[seq_len(min(.N, 12))])

# ---------------------------------------------------------------------
# Write CSV outputs
# ---------------------------------------------------------------------
# 16-site outputs go to manuscript_tables/ (sibling of checkpoints/,
# alongside other Track A/B outputs like als_point_density_full.csv and
# site_effects_vs_offset.csv). NOT to the pipeline's out_tables, which
# resolves to tables/ and is used for different, pipeline-internal files.

dir.create(manuscript_tables_dir, showWarnings = FALSE, recursive = TRUE)

fwrite(cmp,     file.path(manuscript_tables_dir, "chm_sensitivity_comparison.csv"))
fwrite(fx_join, file.path(manuscript_tables_dir, "chm_sensitivity_fixef_comparison.csv"))
log_progress(sprintf("  Wrote: %s/chm_sensitivity_comparison.csv",       manuscript_tables_dir))
log_progress(sprintf("  Wrote: %s/chm_sensitivity_fixef_comparison.csv", manuscript_tables_dir))

# ---------------------------------------------------------------------
# Apply 16-site decision rule
# ---------------------------------------------------------------------
# Decision rule (from 16-site checkpoint):
#   - bias shift < 0.3 m AND all predictors preserve zero-exclusion status
#       -> minimal update: Supp Text S2 + Supp Tables S8/S9a only
#   - bias shift >= 0.3 m OR any predictor flips zero-exclusion status
#       -> full update: also Abstract, Results CHM section, Fig 2, Table 6
#
# Computed before save_checkpoint so the verdict is persisted for
# downstream chats and programmatic access.

bias_full    <- cmp[fit == "CHM 19-site (full)",         intercept_mean]
bias_16      <- cmp[fit == "CHM 16-site (sensitivity)",  intercept_mean]
bias_shift   <- bias_16 - bias_full
n_flipped    <- sum(!fx_join$status_preserved, na.rm = TRUE)
flipped_pred <- fx_join[status_preserved == FALSE, predictor]
minimal_path <- (abs(bias_shift) < 0.3) && (n_flipped == 0)
verdict      <- if (minimal_path) "MINIMAL UPDATE PATH" else "FULL UPDATE PATH"

# ---------------------------------------------------------------------
# Save checkpoint (including verdict)
# ---------------------------------------------------------------------

checkpoint_data <- list(
  fit_chm_16       = fit_chm_16,
  chm_formula_s2   = chm_formula_s2,
  sites_excluded   = SITES_TO_EXCLUDE,
  comparison_table = cmp,
  fixef_comparison = fx_join,
  loo_chm_16       = loo_chm_16,
  max_rhat_16      = max_rhat_16,
  min_ess_16       = min_ess_16,
  bias_full        = bias_full,
  bias_16          = bias_16,
  bias_shift       = bias_shift,
  n_flipped        = n_flipped,
  flipped_pred     = flipped_pred,
  minimal_path     = minimal_path,
  verdict          = verdict
)
save_checkpoint("10b_chm_sensitivity", checkpoint_data)

# ---------------------------------------------------------------------
# Print verdict and next steps
# ---------------------------------------------------------------------

log_progress("")
log_progress("================================================================")
log_progress("16-site decision-rule verdict")
log_progress("================================================================")
log_progress(sprintf("  Global CHM bias:      full = %+.3f m,  16-site = %+.3f m",
                     bias_full, bias_16))
log_progress(sprintf("  Bias shift:           %+.3f m  (threshold: 0.3 m)",
                     bias_shift))
log_progress(sprintf("  Predictors flipped:   %d of %d",
                     n_flipped, nrow(fx_join)))
if (n_flipped > 0) {
  log_progress(sprintf("    Flipped: %s",
                       paste(flipped_pred, collapse = ", ")))
}

log_progress("")
log_progress(sprintf("  VERDICT: %s", verdict))
if (minimal_path) {
  log_progress("    - Fill Supp Text S2 placeholders with refit values")
  log_progress("    - Update Supp Tables S8 and S9a with CHM sensitivity values")
  log_progress("    - Abstract, Results, Figure 2, Table 6: unchanged")
} else {
  log_progress("    - All minimal-path updates, PLUS:")
  log_progress("    - Update Abstract quantitative CHM bias")
  log_progress("    - Update Results CHM summary paragraphs")
  log_progress("    - Update Figure 2 and caption if values are quoted")
  log_progress("    - Recompute Table 6 (CHM holdout accuracy) on 16-site refit")
  log_progress("    - Add Discussion paragraph on reference-data variance contribution")
}

log_progress("")
log_progress("Next steps:")
log_progress("  1. Inspect chm_sensitivity_comparison.csv. Note changes in")
log_progress("     intercept_mean, sd_site_mean, and re_site_min/max.")
log_progress("  2. Inspect chm_sensitivity_fixef_comparison.csv. The")
log_progress("     status_preserved column (FALSE = predictor flipped its")
log_progress("     zero-exclusion status between fits) drives 16-site's")
log_progress("     stability call. cis_overlap_bw_fits is the auxiliary")
log_progress("     two-fit CI-overlap diagnostic.")
log_progress("  3. Follow the verdict above: either minimal or full update")
log_progress("     path into Supp Text S2 and downstream manuscript sections.")
log_progress("================================================================")
