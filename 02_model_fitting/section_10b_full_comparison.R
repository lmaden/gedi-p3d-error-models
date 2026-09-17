# ============================================================================
# section_10b_full_comparison.R
#
# Comprehensive 19-site vs 16-site CHM model comparison.
# Full-versus-restricted predictor comparison for the CHM sensitivity refit.
#
# Inputs (under CHKPT_DIR):
#   10_models_stage2.rds      19-site CHM/DTM fits
#   10b_chm_sensitivity.rds   16-site CHM sensitivity fit
#   08_model_prep.rds         training/holdout split
#   01_data_ingest.rds        full data (fallback for holdout derivation)
#
# Outputs (under OUT_DIR, manuscript_tables/):
#   sensitivity_variance_partition.csv
#   sensitivity_nu_sigma.csv
#   sensitivity_holdout_overall.csv
#   sensitivity_holdout_by_forest_type.csv
#   sensitivity_pi_widths.csv
#   sensitivity_re_decomposition.csv
#   sensitivity_site_re_comparison.csv
#   sensitivity_conditional_effects.csv
#
# Intermediate cache (re-used across runs, deletable):
#   checkpoints/_cache_holdout_predictions.rds
#
# Runtime: ~1-3 hours dominated by posterior_predict on the 19-site holdout.
# Each section is wrapped so a failure midway does not lose earlier outputs.
#
# Usage from a session that already has fits in memory (skips reload):
#   source("section_10b_full_comparison.R")
#
# Usage from a fresh R session:
#   Rscript section_10b_full_comparison.R
# ============================================================================


# ---- CONFIG -----------------------------------------------------------------
# Adjust here, not below. CHM data column names are inferred from the project
# conventions in v5 and Methods §2.3; verify the assertions in Section 1 if
# this script is run for the first time on a checkpoint.

CHKPT_DIR <- "checkpoints"
OUT_DIR   <- "manuscript_tables"
CACHE_DIR <- "checkpoints"

# Column names in the holdout data
COL_RESPONSE    <- "err_chm"          # CHM error response
COL_SITE        <- "site"             # site identifier (1-19, character or factor)
COL_FOREST_TYPE <- "lc_l1_code"       # ENF / BDF / DNF / EBF (and other classes)
FOREST_CODES    <- c("ENF", "BDF", "DNF", "EBF")  # functional types reported in Table 6

# Sites excluded from the 16-site fit
SITES_FLAGGED <- c("1", "2", "3")

# Posterior draws used for posterior_predict. Lower this for testing.
N_DRAWS_PP <- 500

# Predictors for which to write conditional-effects tables
COND_PREDICTORS <- c("slope_mean_z", "wsci_z", "rh_98_z", "cover_z",
                     "meta_offnad_z", "meta_sunel_z", "meta_leafon_z",
                     "meta_stereo_z")

# DEBUG_MODE: subset holdout to this many rows for fast end-to-end test.
# Set to NULL or 0 for the production run.
DEBUG_MODE <- NULL

SEED <- 20260428


# ---- 0. Setup ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(brms)
  library(data.table)
  library(posterior)
})

if (!dir.exists(OUT_DIR))   dir.create(OUT_DIR,   recursive = TRUE)
if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, recursive = TRUE)

set.seed(SEED)

.t0 <- Sys.time()
log_msg <- function(...) {
  elapsed <- format(round(as.numeric(difftime(Sys.time(), .t0, units = "mins")), 1))
  message(sprintf("[%s | +%s min] %s",
                  format(Sys.time(), "%H:%M:%S"), elapsed,
                  paste0(..., collapse = "")))
}

run_section <- function(name, expr) {
  log_msg(sprintf("=== START: %s", name))
  out <- tryCatch(force(expr),
                  error = function(e) {
                    log_msg(sprintf("!!! %s FAILED: %s", name, conditionMessage(e)))
                    NULL
                  })
  log_msg(sprintf("=== END:   %s", name))
  invisible(out)
}

write_csv_log <- function(x, fname) {
  path <- file.path(OUT_DIR, fname)
  fwrite(x, path)
  log_msg(sprintf("  wrote %s (%d rows)", path, nrow(x)))
}


# ---- 1. Load fits and data --------------------------------------------------

# Checkpoint payloads are wrapped: list(data = ..., timestamp = ..., session_info = ...)
load_payload <- function(rds_path, candidates) {
  if (!file.exists(rds_path)) stop(sprintf("checkpoint not found: %s", rds_path))
  payload <- readRDS(rds_path)
  data <- if (is.list(payload) && !is.null(payload$data)) payload$data else payload
  for (n in candidates) {
    if (!is.null(data[[n]])) return(list(name = n, obj = data[[n]]))
  }
  stop(sprintf("none of (%s) found in %s. available: %s",
               paste(candidates, collapse = ", "), rds_path,
               paste(names(data), collapse = ", ")))
}

if (!exists("fit_chm_full")) {
  log_msg("loading fit_chm_full from 10_models_stage2.rds")
  hit <- load_payload(file.path(CHKPT_DIR, "10_models_stage2.rds"),
                      c("fit_chm_full", "fit_chm_s2", "fit_chm"))
  fit_chm_full <- hit$obj
  log_msg(sprintf("  loaded as %s", hit$name))
}

if (!exists("fit_chm_16")) {
  log_msg("loading fit_chm_16 from 10b_chm_sensitivity.rds")
  hit <- load_payload(file.path(CHKPT_DIR, "10b_chm_sensitivity.rds"),
                      c("fit_chm_16", "fit_chm_sensitivity", "fit_chm_s2_16"))
  fit_chm_16 <- hit$obj
  log_msg(sprintf("  loaded as %s", hit$name))
}

# Holdout: try the standard prep checkpoint, then fall back to deriving from full data
load_holdout <- function() {
  prep_path <- file.path(CHKPT_DIR, "08_model_prep.rds")
  full_path <- file.path(CHKPT_DIR, "01_data_ingest.rds")

  if (file.exists(prep_path)) {
    payload <- readRDS(prep_path)
    data <- if (is.list(payload) && !is.null(payload$data)) payload$data else payload
    candidates <- c("holdout_chm_s2", "holdout_chm", "chm_holdout",
                    "test_chm_s2", "test_chm", "mod_chm_s2_holdout")
    for (n in candidates) {
      obj <- data[[n]]
      if (!is.null(obj) && is.data.frame(obj)) {
        log_msg(sprintf("  found holdout as %s in 08_model_prep.rds (n = %d)", n, nrow(obj)))
        return(as.data.table(obj))
      }
    }
    log_msg("  no named holdout in 08_model_prep.rds; falling back to anti-join")
  }

  if (!file.exists(full_path) || !file.exists(prep_path)) {
    stop("cannot derive holdout: missing 01_data_ingest.rds or 08_model_prep.rds")
  }

  full_d <- readRDS(full_path)
  full_d <- if (is.list(full_d) && !is.null(full_d$data)) full_d$data else full_d
  chm_full <- full_d[["chm_df"]]
  prep_d <- readRDS(prep_path)
  prep_d <- if (is.list(prep_d) && !is.null(prep_d$data)) prep_d$data else prep_d
  train  <- prep_d[["mod_chm_s2"]]

  if (is.null(chm_full) || is.null(train)) {
    stop("anti-join fallback: chm_df or mod_chm_s2 not found")
  }
  if (!"shot_number" %in% names(chm_full) || !"shot_number" %in% names(train)) {
    stop("anti-join fallback requires shot_number in both full and training data")
  }
  hd <- as.data.table(chm_full)[!shot_number %in% train$shot_number]
  log_msg(sprintf("  derived holdout via shot_number anti-join (n = %d)", nrow(hd)))
  hd
}

if (!exists("holdout_chm")) {
  log_msg("loading holdout_chm")
  holdout_chm <- load_holdout()
}

# Verify column assumptions before any heavy lifting
required_cols <- c(COL_RESPONSE, COL_SITE, COL_FOREST_TYPE)
missing_cols  <- setdiff(required_cols, names(holdout_chm))
if (length(missing_cols)) {
  stop(sprintf("holdout missing required columns: %s. Available: %s",
               paste(missing_cols, collapse = ", "),
               paste(names(holdout_chm), collapse = ", ")))
}

# Force site to character to match SITES_FLAGGED
holdout_chm[, (COL_SITE) := as.character(get(COL_SITE))]

if (!is.null(DEBUG_MODE) && DEBUG_MODE > 0) {
  log_msg(sprintf("DEBUG_MODE: subsetting holdout to %d rows", DEBUG_MODE))
  holdout_chm <- holdout_chm[sample(.N, min(.N, DEBUG_MODE))]
}

holdout_chm_16 <- holdout_chm[!get(COL_SITE) %in% SITES_FLAGGED]
log_msg(sprintf("holdout sizes: 19-site n = %d, 16-site n = %d",
                nrow(holdout_chm), nrow(holdout_chm_16)))


# ---- 2. Variance partitioning -----------------------------------------------
# Conditional R² (FE+RE) and marginal R² (FE only). Shares are FE_only / total
# and (total - FE_only) / total per posterior draw.

run_section("variance_partition", {
  vp_one <- function(fit, label) {
    R2_full   <- bayes_R2(fit, summary = FALSE)[, 1]
    R2_femarg <- bayes_R2(fit, summary = FALSE, re_formula = NA)[, 1]
    re_share  <- pmax(0, R2_full - R2_femarg) / pmax(R2_full, .Machine$double.eps)
    fe_share  <- 1 - re_share

    qsum <- function(x) {
      c(mean = mean(x), median = median(x),
        q025 = unname(quantile(x, 0.025)),
        q975 = unname(quantile(x, 0.975)))
    }
    rbind(
      data.table(metric = "R2_total",   t(qsum(R2_full))),
      data.table(metric = "R2_FE_only", t(qsum(R2_femarg))),
      data.table(metric = "RE_share",   t(qsum(re_share))),
      data.table(metric = "FE_share",   t(qsum(fe_share)))
    )[, model := label]
  }

  out <- rbindlist(list(vp_one(fit_chm_full, "19_site"),
                        vp_one(fit_chm_16,   "16_site")), use.names = TRUE)
  setcolorder(out, c("model", "metric", "mean", "median", "q025", "q975"))
  write_csv_log(out, "sensitivity_variance_partition.csv")
})


# ---- 3. Student-t nu and residual sigma -------------------------------------

run_section("nu_sigma", {
  ns_one <- function(fit, label) {
    dr <- as_draws_df(fit)
    nu_col <- intersect(c("nu", "b_nu"), names(dr))
    if (!length(nu_col)) nu_col <- grep("^nu(_Intercept)?$", names(dr), value = TRUE)
    if (!length(nu_col)) {
      log_msg(sprintf("  no nu in %s; available: %s",
                      label, paste(head(names(dr), 30), collapse = ", ")))
      return(data.table())
    }
    nu <- dr[[nu_col[1]]]

    sig_col <- grep("^b_sigma_Intercept$|^sigma$", names(dr), value = TRUE)
    sigma_row <- if (length(sig_col)) {
      s <- dr[[sig_col[1]]]
      if (sig_col[1] == "b_sigma_Intercept") s <- exp(s)  # log-link
      data.table(parameter = "sigma_intercept",
                 mean = mean(s), median = median(s),
                 q025 = unname(quantile(s, 0.025)),
                 q975 = unname(quantile(s, 0.975)))
    } else data.table()

    rbindlist(list(
      data.table(parameter = "nu",
                 mean = mean(nu), median = median(nu),
                 q025 = unname(quantile(nu, 0.025)),
                 q975 = unname(quantile(nu, 0.975))),
      sigma_row
    ), fill = TRUE)[, model := label]
  }

  out <- rbindlist(list(ns_one(fit_chm_full, "19_site"),
                        ns_one(fit_chm_16,   "16_site")),
                   use.names = TRUE, fill = TRUE)
  setcolorder(out, c("model", "parameter", "mean", "median", "q025", "q975"))
  write_csv_log(out, "sensitivity_nu_sigma.csv")
})


# ---- 4. Posterior predictions on holdout (cached) ---------------------------

CACHE_PATH <- file.path(CACHE_DIR, "_cache_holdout_predictions.rds")
pp_cache <- if (file.exists(CACHE_PATH)) readRDS(CACHE_PATH) else list()

ensure_pp <- function(fit, newdata, key) {
  if (!is.null(pp_cache[[key]]) && ncol(pp_cache[[key]]) == nrow(newdata)) {
    log_msg(sprintf("  cache hit: %s", key))
    return(pp_cache[[key]])
  }
  log_msg(sprintf("  posterior_predict: %s (n = %d, ndraws = %d)",
                  key, nrow(newdata), N_DRAWS_PP))
  pp <- posterior_predict(fit, newdata = newdata, ndraws = N_DRAWS_PP,
                          allow_new_levels = FALSE)
  pp_cache[[key]] <<- pp
  saveRDS(pp_cache, CACHE_PATH)
  pp
}

run_section("posterior_predict", {
  ensure_pp(fit_chm_full, holdout_chm,    "pp_19_full")    # 19-site model on full holdout
  ensure_pp(fit_chm_full, holdout_chm_16, "pp_19_on_16")   # 19-site model on 16-site holdout
  ensure_pp(fit_chm_16,   holdout_chm_16, "pp_16_full")    # 16-site model on 16-site holdout
  invisible(NULL)
})


# ---- 5. Holdout overall accuracy --------------------------------------------

acc_metrics <- function(pp, y) {
  pred_mean <- colMeans(pp)
  pi_lo <- apply(pp, 2, quantile, 0.025)
  pi_hi <- apply(pp, 2, quantile, 0.975)
  resid <- y - pred_mean
  list(
    n             = length(y),
    bias          = mean(pred_mean - y),
    rmse          = sqrt(mean(resid^2)),
    mae           = mean(abs(resid)),
    r2            = 1 - sum(resid^2) / sum((y - mean(y))^2),
    coverage95    = mean(y >= pi_lo & y <= pi_hi),
    pi_width_med  = median(pi_hi - pi_lo),
    pi_width_mean = mean(pi_hi - pi_lo)
  )
}

run_section("holdout_overall", {
  rows <- list(
    cbind(model = "19_site", evaluated_on = "19_site_holdout",
          as.data.table(acc_metrics(pp_cache$pp_19_full,  holdout_chm[[COL_RESPONSE]]))),
    cbind(model = "19_site", evaluated_on = "16_site_holdout",
          as.data.table(acc_metrics(pp_cache$pp_19_on_16, holdout_chm_16[[COL_RESPONSE]]))),
    cbind(model = "16_site", evaluated_on = "16_site_holdout",
          as.data.table(acc_metrics(pp_cache$pp_16_full,  holdout_chm_16[[COL_RESPONSE]])))
  )
  out <- rbindlist(rows, use.names = TRUE)
  write_csv_log(out, "sensitivity_holdout_overall.csv")
})


# ---- 6. Holdout by forest type ----------------------------------------------

acc_by_group <- function(pp, dat, grp_col, codes) {
  y <- dat[[COL_RESPONSE]]
  g <- as.character(dat[[grp_col]])
  out <- rbindlist(lapply(codes, function(cd) {
    idx <- which(g == cd)
    if (!length(idx)) return(NULL)
    cbind(forest_type = cd,
          as.data.table(acc_metrics(pp[, idx, drop = FALSE], y[idx])))
  }), use.names = TRUE, fill = TRUE)
  out
}

run_section("holdout_by_forest_type", {
  rows <- list(
    cbind(model = "19_site", evaluated_on = "19_site_holdout",
          acc_by_group(pp_cache$pp_19_full, holdout_chm,
                       COL_FOREST_TYPE, FOREST_CODES)),
    cbind(model = "19_site", evaluated_on = "16_site_holdout",
          acc_by_group(pp_cache$pp_19_on_16, holdout_chm_16,
                       COL_FOREST_TYPE, FOREST_CODES)),
    cbind(model = "16_site", evaluated_on = "16_site_holdout",
          acc_by_group(pp_cache$pp_16_full, holdout_chm_16,
                       COL_FOREST_TYPE, FOREST_CODES))
  )
  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  write_csv_log(out, "sensitivity_holdout_by_forest_type.csv")
})


# ---- 7. Prediction interval widths (overall + by site) ----------------------

pi_widths_by_site <- function(pp, dat, label) {
  pi_lo <- apply(pp, 2, quantile, 0.025)
  pi_hi <- apply(pp, 2, quantile, 0.975)
  width <- pi_hi - pi_lo
  d <- data.table(site = as.character(dat[[COL_SITE]]), width = width)
  out <- d[, .(n = .N,
               width_mean = mean(width),
               width_median = median(width),
               width_q025 = quantile(width, 0.025),
               width_q975 = quantile(width, 0.975)), by = site]
  out[, model := label][]
}

run_section("pi_widths", {
  # Overall: site-aggregated medians/min/max as reported in §3.4 (3.5-29.7 m for 19-site CHM)
  by_site <- rbind(
    pi_widths_by_site(pp_cache$pp_19_full, holdout_chm,    "19_site"),
    pi_widths_by_site(pp_cache$pp_16_full, holdout_chm_16, "16_site")
  )

  overall <- by_site[, .(
    site_count       = .N,
    width_overall_med = median(width_median),
    width_min_site    = min(width_median),
    width_max_site    = max(width_median)
  ), by = model]

  setcolorder(by_site, c("model", "site", "n", "width_mean", "width_median",
                         "width_q025", "width_q975"))
  write_csv_log(by_site, "sensitivity_pi_widths.csv")
  # also write the rolled-up overall for convenience
  write_csv_log(overall, "sensitivity_pi_widths_overall.csv")
})


# ---- 8. Three-level RE decomposition (site / ecoregion / land cover) --------
# Uses sd_*__Intercept draws to compute variance shares per draw, summarized
# over the posterior. Mirrors Table S12 in the manuscript.

run_section("re_decomposition", {
  re_decomp_one <- function(fit, label) {
    dr <- as_draws_df(fit)
    # brms names random-effect SDs as sd_<group>__Intercept (or with slope)
    sd_cols <- grep("^sd_.*__Intercept$", names(dr), value = TRUE)
    if (!length(sd_cols)) {
      log_msg(sprintf("  no sd_*__Intercept columns in %s", label))
      return(data.table())
    }
    # Map column names back to grouping factor names: sd_<group>__Intercept
    group_names <- sub("^sd_(.*)__Intercept$", "\\1", sd_cols)

    var_mat <- as.matrix(dr[, sd_cols, drop = FALSE])^2
    colnames(var_mat) <- group_names
    total <- rowSums(var_mat)
    share <- sweep(var_mat, 1, pmax(total, .Machine$double.eps), "/")

    qsum <- function(x) c(mean = mean(x), median = median(x),
                          q025 = unname(quantile(x, 0.025)),
                          q975 = unname(quantile(x, 0.975)))

    rbindlist(lapply(group_names, function(g) {
      rbindlist(list(
        data.table(group = g, quantity = "variance",
                   t(qsum(var_mat[, g]))),
        data.table(group = g, quantity = "share",
                   t(qsum(share[, g])))
      ))
    }))[, model := label][]
  }

  out <- rbindlist(list(re_decomp_one(fit_chm_full, "19_site"),
                        re_decomp_one(fit_chm_16,   "16_site")), use.names = TRUE)
  setcolorder(out, c("model", "group", "quantity", "mean", "median", "q025", "q975"))
  write_csv_log(out, "sensitivity_re_decomposition.csv")
})


# ---- 9. Site-level RE comparison (16 shared sites) --------------------------
# Per-site posterior median and 95% CI from each fit, plus the diff and the
# correlation across sites. The sensitivity model only has site REs for sites
# 4-19; for the 19-site fit we also report sites 1-3 for completeness.

run_section("site_re_comparison", {
  ranef_summary <- function(fit, label) {
    re <- ranef(fit, summary = TRUE)
    # ranef() returns a list per grouping factor; take "site"
    if (is.null(re$site)) {
      log_msg(sprintf("  no $site in ranef(%s); available: %s",
                      label, paste(names(re), collapse = ", ")))
      return(data.table())
    }
    arr <- re$site  # [n_sites, 4 (Estimate, Est.Error, Q2.5, Q97.5), n_terms]
    # We want the Intercept term
    if ("Intercept" %in% dimnames(arr)[[3]]) {
      m <- arr[, , "Intercept"]
    } else {
      m <- arr[, , 1]  # take the first term if Intercept missing
    }
    data.table(
      site = rownames(m),
      re_mean   = m[, "Estimate"],
      re_se     = m[, "Est.Error"],
      re_q025   = m[, "Q2.5"],
      re_q975   = m[, "Q97.5"]
    )[, model := label][]
  }

  re_full <- ranef_summary(fit_chm_full, "19_site")
  re_16   <- ranef_summary(fit_chm_16,   "16_site")

  # Side-by-side for the 16 shared sites, and full table including flagged sites
  shared <- merge(
    re_full[, .(site, re_mean_19 = re_mean, re_q025_19 = re_q025, re_q975_19 = re_q975)],
    re_16[,   .(site, re_mean_16 = re_mean, re_q025_16 = re_q025, re_q975_16 = re_q975)],
    by = "site", all = TRUE
  )
  shared[, diff := re_mean_19 - re_mean_16]

  # Pearson correlation across the shared 16 sites (where both fits have estimates)
  shared_ok <- shared[!is.na(re_mean_19) & !is.na(re_mean_16)]
  rho <- if (nrow(shared_ok) >= 3) cor(shared_ok$re_mean_19, shared_ok$re_mean_16) else NA_real_
  log_msg(sprintf("  site RE correlation (shared sites, n=%d): r = %.3f",
                  nrow(shared_ok), rho))

  setcolorder(shared, c("site", "re_mean_19", "re_q025_19", "re_q975_19",
                        "re_mean_16", "re_q025_16", "re_q975_16", "diff"))
  write_csv_log(shared, "sensitivity_site_re_comparison.csv")
})


# ---- 10. Conditional effects -------------------------------------------------
# conditional_effects() returns a list of data.frames; we stack them with a
# predictor column and a model label.

run_section("conditional_effects", {
  ce_one <- function(fit, label, predictors) {
    pieces <- list()
    for (p in predictors) {
      ce <- tryCatch(conditional_effects(fit, effects = p)[[1]],
                     error = function(e) {
                       log_msg(sprintf("  conditional_effects failed for %s in %s: %s",
                                       p, label, conditionMessage(e)))
                       NULL
                     })
      if (is.null(ce)) next
      ce <- as.data.table(ce)
      keep <- intersect(c(p, "estimate__", "se__", "lower__", "upper__"), names(ce))
      ce <- ce[, ..keep]
      setnames(ce, p, "x")
      ce[, predictor := p]
      pieces[[p]] <- ce
    }
    if (!length(pieces)) return(data.table())
    rbindlist(pieces, use.names = TRUE, fill = TRUE)[, model := label][]
  }

  out <- rbind(
    ce_one(fit_chm_full, "19_site", COND_PREDICTORS),
    ce_one(fit_chm_16,   "16_site", COND_PREDICTORS),
    fill = TRUE
  )
  if (nrow(out)) {
    setcolorder(out, c("model", "predictor", "x"))
    write_csv_log(out, "sensitivity_conditional_effects.csv")
  } else {
    log_msg("  no conditional effects produced")
  }
})


# ---- Summary -----------------------------------------------------------------

log_msg("=== run summary ===")
expected <- c("sensitivity_variance_partition.csv",
              "sensitivity_nu_sigma.csv",
              "sensitivity_holdout_overall.csv",
              "sensitivity_holdout_by_forest_type.csv",
              "sensitivity_pi_widths.csv",
              "sensitivity_pi_widths_overall.csv",
              "sensitivity_re_decomposition.csv",
              "sensitivity_site_re_comparison.csv",
              "sensitivity_conditional_effects.csv")
for (f in expected) {
  path <- file.path(OUT_DIR, f)
  status <- if (file.exists(path)) "OK" else "MISSING"
  log_msg(sprintf("  [%s] %s", status, path))
}
log_msg(sprintf("total elapsed: %.1f min",
                as.numeric(difftime(Sys.time(), .t0, units = "mins"))))
