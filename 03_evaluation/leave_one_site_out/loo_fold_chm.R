#!/usr/bin/env Rscript
# =====================================================================
# loo_fold_chm.R  --  one leave-one-site-out fold for the CHM model
#
# Canopy-height counterpart of loo_fold_dtm.R. One site per run, so folds
# can be spread across nodes.
#
# FIVE THINGS THAT DIFFER FROM THE DTM VERSION, ALL VERIFIED
#
#   1. SITES 1, 2 AND 3 ARE EXCLUDED, as in revision_confirm_fullfit.R,
#      because their ALS reference is unreliable. That leaves 16 sites:
#      4 5 6 7 8 9 10 11 12 13 14 15 17 18 19 20, and 125,086 rows,
#      which becomes the published 124,966 after brms drops 120
#      incomplete cases. So there are 16 folds here, not 18.
#
#   2. Response is chm_error_mean, data is mod_chm_s2 (204,792 rows
#      across 19 sites before the exclusion).
#
#   3. NO lkj_corr_cholesky prior. The CHM confirmation fit sets seven
#      priors and leaves the random-effect correlation at the brms
#      default; the DTM fit set lkj(2). Mirroring each product's own fit
#      matters more than making the two scripts look alike.
#
#   4. adapt_delta 0.99 and max_treedepth 14, NOT the confirmation fit's
#      0.95 and 12. Reasoning, since this is a departure:
#        - The CHM confirmation fit already touched its cap of 12 on
#          0.18% of draws. It passed, but with no headroom.
#        - Holding out a site removes a random-effect level, which
#          sharpened the funnel on DTM: the DTM full fit was clean at
#          adapt_delta 0.98 yet 6 of 18 folds diverged at the same value.
#        - Raising adapt_delta lowers the step size, which raises tree
#          depth, so a cap that is already binding would bind harder.
#        - On DTM, moving 0.98 to 0.99 changed fold RMSE by a median of
#          0.05%, so this affects sampling reliability, not the estimates.
#      Both are overridable: --adapt-delta and --max-treedepth.
#
#   5. Oracle predictions come from fit_chm_16_tuned.rds, which was fitted
#      on the same 16 sites the folds train on.
#
# THE THREE RUNGS  (identical in meaning to the DTM version)
#   1. ORACLE, site offset known. Predictions from the full 16-site fit,
#      in sample. Deliberately optimistic and NOT a validation: the site
#      offset was estimated from the very footprints being predicted. It
#      exists to express the site effect in metres rather than variance
#      percentages, and must be labelled an upper bound wherever it
#      appears.
#   2. NEW SITE, ecoregion and land cover known. The operational case.
#   3. NEW SITE AND NEW ECOREGION, land cover known. Tests transfer
#      beyond the sampled ecoregions. Land cover is never treated as
#      unknown: it is a fixed-effect covariate readable from GLC-FCS30
#      for any pixel, and a new level could not enter the design matrix.
#
#   If a held-out site is the only one in its ecoregion, rung 2 cannot use
#   a known ecoregion effect and collapses into rung 3. The script detects
#   this and flags eco_known = FALSE rather than mislabelling it.
#
# USAGE
#   export PROJECT_ROOT=/gpfs/data1/vclgp/lmaden/chpt1
#   export TMPDIR=$PROJECT_ROOT/tmp
#   Rscript loo_fold_chm.R --site 18
#
#   --site           required
#   --iter           total iterations, half warmup. Default 3000.
#   --threads        per chain. Default 16, so 64 cores for 4 chains.
#   --adapt-delta    default 0.99
#   --max-treedepth  default 14
#   --init           "random" (Stan default, uniform(-2,2) unconstrained)
#                    or a number such as 0. Use 0 when a single chain gets
#                    stuck at the treedepth ceiling while the others are
#                    fine, which shows up as exactly 25.00% saturation.
#   --smoke          tiny subsample, checks the script runs
#
#   A non-default adapt_delta or max_treedepth adds a suffix to every
#   output name, so nothing overwrites anything.
#
# EXPECTED TIME
#   The CHM confirmation fit took 5.81 h at 0.95 and max_treedepth 12 on
#   all 16 sites. Fold cost scales with rows retained, so holding out a
#   small site costs nearly a full fit while site 5, which is 56.3% of the
#   data, costs less than half. At 0.99 and depth 14 expect meaningfully
#   longer, but the multiplier is unmeasured, which is why the first fold
#   should be run alone before committing to all sixteen.
# =====================================================================

suppressPackageStartupMessages({
  library(brms); library(data.table); library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
SITE          <- getarg("--site")
ITER          <- as.integer(getarg("--iter", "3000"))
THREADS       <- as.integer(getarg("--threads", "16"))
ADAPT_DELTA   <- as.numeric(getarg("--adapt-delta", "0.99"))
MAX_TREEDEPTH <- as.integer(getarg("--max-treedepth", "14"))
INIT          <- getarg("--init", "random")
SMOKE         <- "--smoke" %in% args

if (is.null(SITE)) {
  cat("\nYou must name a site to hold out, for example:\n")
  cat("  Rscript loo_fold_chm.R --site 18\n\n")
  quit(status = 2)
}
WARMUP <- ITER %/% 2

SITES_FLAGGED <- c("1", "2", "3")   # unreliable ALS reference, per the CHM fit
NDRAWS_EVAL   <- 1000
CHUNK_ROWS    <- 5000

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
CHKPT   <- file.path(PROJECT_ROOT, "checkpoints")
FULLFIT <- file.path(PROJECT_ROOT, "models", "confirm", "fit_chm_16_tuned.rds")
OUTMOD  <- file.path(PROJECT_ROOT, "models", "loo")
OUTTBL  <- file.path(PROJECT_ROOT, "manuscript_tables", "loo")
dir.create(OUTMOD, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTTBL, recursive = TRUE, showWarnings = FALSE)

SUF <- ""
if (abs(ADAPT_DELTA - 0.99) > 1e-9)
  SUF <- paste0(SUF, "_ad", sub("\\.", "", format(ADAPT_DELTA, nsmall = 2)))
if (MAX_TREEDEPTH != 14) SUF <- paste0(SUF, "_td", MAX_TREEDEPTH)
if (INIT != "random")    SUF <- paste0(SUF, "_init", INIT)
TAG <- sprintf("chm_site%s%s%s", SITE, SUF, if (SMOKE) "_smoke" else "")

LOG <- file.path(OUTTBL, sprintf("loo_%s_%s.log", TAG, format(Sys.Date(), "%Y%m%d")))
con <- file(LOG, open = "wt"); sink(con, split = TRUE)

hr  <- function() cat(strrep("=", 92), "\n")
sec <- function(x) { cat("\n"); hr(); cat("## ", x, "\n"); hr() }
unwrap <- function(p) if (is.list(p) && !is.null(p$data)) p$data else p

cat("loo_fold_chm.R  --  holding out site ", SITE, if (SMOKE) "  [SMOKE]" else "", "\n", sep = "")
cat("started: ", format(Sys.time()), "\n")
cat(sprintf("iterations %d (%d warmup) | 4 chains x %d threads = %d cores\n",
            ITER, WARMUP, THREADS, 4 * THREADS))
cat(sprintf("adapt_delta %.3f | max_treedepth %d | init %s | output tag: %s\n",
            ADAPT_DELTA, MAX_TREEDEPTH, INIT, TAG))


# =====================================================================
sec("STEP 0  --  WRITE-LOCATION GUARD")
# =====================================================================
ALLOWED <- "/gpfs/data1/vclgp/lmaden"
cat("tempdir(): ", tempdir(), "\n")
if (!startsWith(normalizePath(tempdir(), mustWork = FALSE), ALLOWED)) {
  cat("\nFATAL: tempdir() is outside ", ALLOWED, "\n", sep = "")
  cat("Run  export TMPDIR=$PROJECT_ROOT/tmp  then relaunch.\n")
  sink(); close(con); quit(status = 1)
}
cat("OK\n")


# =====================================================================
sec("STEP 1  --  DATA, SPLIT INTO TRAINING AND HELD-OUT")
# =====================================================================

mp  <- unwrap(readRDS(file.path(CHKPT, "08_model_prep.rds")))
dat <- as.data.table(mp$mod_chm_s2)
rm(mp); gc()
cat("mod_chm_s2 rows:", format(nrow(dat), big.mark = ","),
    " sites:", uniqueN(as.character(dat$site)), "\n")

dat <- dat[!as.character(site) %in% SITES_FLAGGED]
cat("after excluding sites", paste(SITES_FLAGGED, collapse = "/"), ":",
    format(nrow(dat), big.mark = ","), "rows,",
    uniqueN(as.character(dat$site)), "sites\n")
cat("(the published CHM fit used 124,966 rows after dropping incomplete cases)\n")

KEEP <- c("ENF", "IWL", "BDF", "RCP", "GRS", "UNK", "DNF", "WTR", "IMP", "EBF")
dat[, lc_grp := factor(fifelse(as.character(lc_l1_code) %in% KEEP,
                               as.character(lc_l1_code), "OTHER"))]
dat[, lc_grp := relevel(lc_grp, ref = "ENF")]
dat[, `:=`(site = factor(site), ecoregion = factor(ecoregion))]
cat("land-cover levels:", nlevels(dat$lc_grp), "->",
    paste(levels(dat$lc_grp), collapse = " "), "\n")

# Drop rows with NAs in any model variable. brms excludes them from the
# fit silently, but posterior_predict returns NaN for them.
MODEL_VARS <- c("chm_error_mean", "slope_mean_z", "slope_sd_z", "wsci_z",
                "rh_98_z", "cover_z", "aspect_sin_z", "aspect_cos_z",
                "meta_offnad_z", "meta_sunel_z", "meta_az_conc_z",
                "meta_stereo_z", "meta_leafon_z", "meta_fwdrev_z",
                "meta_absgeo_z", "meta_relgeo_z", "view_az_sin_z",
                "view_az_cos_z", "site", "ecoregion", "lc_grp")
missing_var <- setdiff(MODEL_VARS, names(dat))
if (length(missing_var)) {
  cat("\nFATAL: columns absent from mod_chm_s2:",
      paste(missing_var, collapse = " "), "\n")
  sink(); close(con); quit(status = 1)
}
n_raw <- nrow(dat)
dat <- dat[complete.cases(dat[, ..MODEL_VARS])]
cat("rows dropped for incomplete model variables:", n_raw - nrow(dat), "\n")

if (!(SITE %in% levels(dat$site))) {
  cat("\nFATAL: site ", SITE, " is not available. Sites are: ",
      paste(levels(droplevels(dat$site)), collapse = " "), "\n", sep = "")
  sink(); close(con); quit(status = 1)
}

train <- dat[site != SITE]
hold  <- dat[site == SITE]
cat(sprintf("\ntraining  %s rows, %d sites\n",
            format(nrow(train), big.mark = ","),
            uniqueN(as.character(train$site))))
cat(sprintf("held out  %s rows (%.2f%% of the CHM analysis data)\n",
            format(nrow(hold), big.mark = ","), 100 * nrow(hold) / nrow(dat)))

eco_held <- unique(as.character(hold$ecoregion))
ECO_KNOWN <- all(eco_held %in% unique(as.character(train$ecoregion)))
cat("held-out ecoregion(s):", paste(eco_held, collapse = " "),
    " also in training data:", ECO_KNOWN, "\n")
if (!ECO_KNOWN) {
  cat("\n>> This site is the only one in its ecoregion, so rung 2 cannot use a\n")
  cat(">> known ecoregion effect and collapses into rung 3. Both are still\n")
  cat(">> computed; rung 2 is flagged eco_known = FALSE.\n")
}

train[, `:=`(site = droplevels(site), ecoregion = droplevels(ecoregion),
             lc_grp = droplevels(lc_grp))]

if (SMOKE) {
  set.seed(2026)
  idx <- integer(0)
  for (s in levels(train$site)) {
    rows <- which(train$site == s)
    idx <- c(idx, rows[sample.int(length(rows), min(length(rows), 250))])
  }
  for (l in levels(train$lc_grp)) {
    rows <- which(train$lc_grp == l)
    if (length(rows) && sum(idx %in% rows) < 30)
      idx <- c(idx, rows[sample.int(length(rows), min(length(rows), 30))])
  }
  train <- train[sort(unique(idx))]
  train[, `:=`(site = droplevels(site), ecoregion = droplevels(ecoregion),
               lc_grp = droplevels(lc_grp))]
  if (nrow(hold) > 500) hold <- hold[sample.int(nrow(hold), 500)]
  ITER <- 200; WARMUP <- 100
  cat(sprintf("\nSMOKE: training %s rows, held out %d, 200 iterations.\n",
              format(nrow(train), big.mark = ","), nrow(hold)))
  cat("Diagnostics will look bad. The question is only whether it finishes.\n")
}


# =====================================================================
sec("STEP 2  --  FIT WITHOUT THE HELD-OUT SITE")
# =====================================================================

main_effects <- paste(
  "slope_mean_z + slope_sd_z + wsci_z + rh_98_z + cover_z +",
  "aspect_sin_z + aspect_cos_z + meta_offnad_z + meta_sunel_z +",
  "meta_az_conc_z + meta_stereo_z + meta_leafon_z + meta_fwdrev_z +",
  "meta_absgeo_z + meta_relgeo_z + view_az_sin_z + view_az_cos_z"
)
f_fold <- bf(
  as.formula(paste(
    "chm_error_mean ~", main_effects,
    "+ s(slope_mean_z, wsci_z, k = 20)",
    "+ slope_mean_z:lc_grp + wsci_z:lc_grp",
    "+ (1 + slope_mean_z | ecoregion) + (1 | site) + (1 + wsci_z | lc_grp)"
  )),
  sigma ~ 1 + slope_mean_z + (1 | site) + (1 | lc_grp)
)

# Seven priors, matching revision_confirm_fullfit.R. No lkj_corr_cholesky:
# the CHM fit leaves the random-effect correlation at the brms default.
priors <- c(
  prior(normal(0, 2),       class = "b"),
  prior(normal(0, 10),      class = "Intercept"),
  prior(gamma(2, 0.1),      class = "nu"),
  prior(student_t(3, 0, 2), class = "sd"),
  prior(normal(0, 1.5),     class = "sds"),
  prior(normal(0, 1),       class = "b",         dpar = "sigma"),
  prior(student_t(3, 0, 5), class = "Intercept", dpar = "sigma")
)
cat("\nFormula:\n"); print(f_fold)
cat("\nPriors:", nrow(priors), "rows (no lkj_corr_cholesky, matching the CHM fit)\n")

fitfile <- file.path(OUTMOD, sprintf("fit_%s", TAG))
if (file.exists(paste0(fitfile, ".rds"))) {
  cat("\n*** A fit for this fold already exists and will be REUSED.\n")
  cat("*** Delete ", paste0(fitfile, ".rds"), " to refit.\n\n", sep = "")
}

t0 <- Sys.time()
fit <- brm(
  formula = f_fold, data = train, family = student(), prior = priors,
  chains = 4, iter = ITER, warmup = WARMUP,
  cores = 4, threads = threading(THREADS),
  control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
  init = if (INIT == "random") "random" else as.numeric(INIT),
  backend = "cmdstanr", file = fitfile, seed = 2026, refresh = 100
)
elapsed_h <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
cat(sprintf("\nfold fit took %.2f hours\n", elapsed_h))
cat("(the 16-site CHM confirmation fit took 5.81 h at adapt_delta 0.95, depth 12)\n")


# =====================================================================
sec("STEP 3  --  FOLD DIAGNOSTICS")
# =====================================================================

np <- as.data.table(brms::nuts_params(fit))
td <- np[Parameter == "treedepth__", Value]
d  <- sum(np[Parameter == "divergent__", Value])
mr <- max(brms::rhat(fit), na.rm = TRUE)
me <- min(brms::neff_ratio(fit), na.rm = TRUE)
ceil <- mean(td >= MAX_TREEDEPTH)

cat(sprintf("  divergences     %6d   %s  (need 0)\n", d, ifelse(d == 0, "PASS", "FAIL")))
cat(sprintf("  max Rhat        %6.4f   %s  (need < 1.01)\n", mr, ifelse(mr < 1.01, "PASS", "FAIL")))
cat(sprintf("  min ESS ratio   %6.3f   %s  (need > 0.10)\n", me, ifelse(me > 0.10, "PASS", "FAIL")))
cat(sprintf("  mean treedepth  %6.2f          (confirmation fit was 10.24)\n", mean(td)))
cat(sprintf("  %% at ceiling     %6.2f%%   %s  (need < 1%%)\n",
            100 * ceil, ifelse(ceil < 0.01, "PASS", "CHECK")))
cat(sprintf("  %% at 12 or more  %6.2f%%          (confirmation fit was 0.18%%)\n",
            100 * mean(td >= 12)))
# Per-chain, because a single stuck chain looks like a global failure in
# the pooled numbers: 25.00% at the ceiling is exactly one chain of four.
tdc <- np[Parameter == "treedepth__", .(mean_td = round(mean(Value), 2),
                                        pct_ceiling = round(100 * mean(Value >= MAX_TREEDEPTH), 1)),
          by = Chain][order(Chain)]
cat("\n  per-chain treedepth:\n")
for (i in seq_len(nrow(tdc)))
  cat(sprintf("    chain %d: mean %.2f, %.1f%% at ceiling%s\n",
              tdc$Chain[i], tdc$mean_td[i], tdc$pct_ceiling[i],
              if (tdc$pct_ceiling[i] > 50) "   <-- STUCK" else ""))
if (any(tdc$pct_ceiling > 50) && !all(tdc$pct_ceiling > 50)) {
  cat("\n  ** One or more chains is stuck at the treedepth ceiling while others\n")
  cat("     are healthy. This is an initialisation problem, not model geometry.\n")
  cat("     Rerun with --init 0. Do NOT raise adapt_delta: that lowers the step\n")
  cat("     size, deepens trees, and makes a pinned chain worse. **\n")
}
if (d > 0) {
  cat("\n  ** Divergences at adapt_delta ", ADAPT_DELTA, ". Rerun this fold with\n", sep = "")
  cat("     --adapt-delta 0.995 rather than accepting them. **\n")
}
if (ceil >= 0.01) {
  cat("\n  ** Treedepth saturating at ", MAX_TREEDEPTH, ". Rerun with\n", sep = "")
  cat("     --max-treedepth 15. A truncated trajectory can bias the posterior. **\n")
}
fold_ok <- (d == 0) && (mr < 1.01) && (me > 0.10) && (ceil < 0.01)


# =====================================================================
sec("STEP 4  --  PREDICT THE HELD-OUT SITE UNDER THREE RUNGS")
# =====================================================================

y <- hold$chm_error_mean

crps_cols <- function(draws, obs) {
  m <- nrow(draws)
  vapply(seq_along(obs), function(j) {
    x <- sort(draws[, j]); yy <- obs[j]; i <- seq_len(m)
    (2 / m^2) * sum((x - yy) * (m * (yy < x) - i + 0.5))
  }, numeric(1))
}

score <- function(draws, obs) {
  bad <- which(!is.finite(colSums(draws)))
  if (length(bad)) {
    cat("    WARNING: dropping", length(bad), "of", ncol(draws),
        "rows with non-finite predictive draws\n")
    draws <- draws[, -bad, drop = FALSE]; obs <- obs[-bad]
  }
  pm <- colMeans(draws)
  q  <- apply(draws, 2, quantile, probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
  res <- pm - obs
  list(summary = data.table(
         n = length(obs), rmse = sqrt(mean(res^2)), bias = mean(res),
         mae = mean(abs(res)), crps = mean(crps_cols(draws, obs)),
         cov50 = mean(obs >= q[3, ] & obs <= q[4, ]),
         cov90 = mean(obs >= q[2, ] & obs <= q[5, ]),
         cov95 = mean(obs >= q[1, ] & obs <= q[6, ]),
         width90 = mean(q[5, ] - q[2, ])),
       rows = data.table(pred_mean = pm, lo90 = q[2, ], hi90 = q[5, ]))
}

predict_blocks <- function(model, newdata, ...) {
  n <- nrow(newdata)
  nd <- min(NDRAWS_EVAL, brms::ndraws(model))
  cat("    scoring with", nd, "posterior draws\n")
  starts <- seq(1, n, by = CHUNK_ROWS)
  out <- vector("list", length(starts))
  for (k in seq_along(starts)) {
    rng <- starts[k]:min(starts[k] + CHUNK_ROWS - 1, n)
    out[[k]] <- posterior_predict(model, newdata = newdata[rng], ndraws = nd, ...)
    if (k %% 5 == 0) cat("    block", k, "of", length(starts), "\n")
  }
  do.call(cbind, out)
}

results <- list(); rowpreds <- list()

cat("\n[rung 1] oracle: full 16-site fit, in sample. UPPER BOUND, NOT VALIDATION.\n")
if (file.exists(FULLFIT)) {
  full <- readRDS(FULLFIT)
  set.seed(1)
  pp <- predict_blocks(full, hold)
  s <- score(pp, y); rm(pp, full); gc()
  results[["1_oracle_site_known"]] <- s$summary
  rowpreds[["1_oracle_site_known"]] <- s$rows
  cat(sprintf("  RMSE %.3f  bias %+.3f  CRPS %.3f  cov90 %.3f\n",
              s$summary$rmse, s$summary$bias, s$summary$crps, s$summary$cov90))
} else {
  cat("  SKIPPED: ", FULLFIT, " not found\n", sep = "")
}

cat("\n[rung 2] new site; ecoregion and land cover known\n")
set.seed(2)
pp <- predict_blocks(fit, hold, allow_new_levels = TRUE, sample_new_levels = "gaussian")
s <- score(pp, y); rm(pp); gc()
s$summary[, eco_known := ECO_KNOWN]
results[["2_new_site_eco_lc_known"]] <- s$summary
rowpreds[["2_new_site_eco_lc_known"]] <- s$rows
cat(sprintf("  RMSE %.3f  bias %+.3f  CRPS %.3f  cov90 %.3f%s\n",
            s$summary$rmse, s$summary$bias, s$summary$crps, s$summary$cov90,
            if (ECO_KNOWN) "" else "   [ecoregion NOT actually known]"))

cat("\n[rung 3] new site AND new ecoregion; land cover known\n")
hold3 <- copy(hold)
# lc_grp deliberately left alone: it enters the fixed effects through
# slope_mean_z:lc_grp and wsci_z:lc_grp, so a new level cannot enter the
# design matrix, and land cover is always knowable from GLC-FCS30.
hold3[, `:=`(site = factor("__new_site__"), ecoregion = factor("__new_eco__"))]
set.seed(3)
pp <- predict_blocks(fit, hold3, allow_new_levels = TRUE, sample_new_levels = "gaussian")
s <- score(pp, y); rm(pp); gc()
results[["3_new_site_new_ecoregion"]] <- s$summary
rowpreds[["3_new_site_new_ecoregion"]] <- s$rows
cat(sprintf("  RMSE %.3f  bias %+.3f  CRPS %.3f  cov90 %.3f\n",
            s$summary$rmse, s$summary$bias, s$summary$crps, s$summary$cov90))


# =====================================================================
sec("STEP 5  --  WRITE RESULTS")
# =====================================================================

tab <- rbindlist(lapply(names(results), function(k)
  cbind(data.table(product = "CHM", site = SITE, rung = k), results[[k]])), fill = TRUE)
tab[, `:=`(n_train = nrow(train), n_held = nrow(hold), iter = ITER,
           threads = THREADS, adapt_delta = ADAPT_DELTA,
           max_treedepth = MAX_TREEDEPTH, fit_hours = elapsed_h,
           divergences = d, max_rhat = mr, min_ess_ratio = me,
           pct_at_ceiling = 100 * ceil, fold_diagnostics_ok = fold_ok)]
print(tab)
fwrite(tab, file.path(OUTTBL, sprintf("loo_metrics_%s.csv", TAG)))

rp <- rbindlist(lapply(names(rowpreds), function(k)
  cbind(data.table(rung = k, site = SITE,
                   shot_number = if ("shot_number" %in% names(hold)) hold$shot_number else NA,
                   observed = y), rowpreds[[k]])), fill = TRUE)
fwrite(rp, file.path(OUTTBL, sprintf("loo_rowpreds_%s.csv", TAG)))

sec("SUMMARY")
cat(sprintf("
  site held out .............. %s  (%s rows, %.2f%% of data)
  fold fit hours ............. %.2f   at %d iter, %d threads, adapt_delta %.3f
  fold diagnostics ........... %s
  divergences ................ %d
  %%%% at treedepth ceiling ..... %.2f%%
  RMSE, oracle / rung2 / rung3 %s
  cov90, oracle / rung2 / rung3 %s

  The gap between the oracle and rung 2 is what knowing the site is worth,
  in metres. For DTM that gap was 0.29 m, about 12%%. If it is much larger
  here, that contrast between products is the paper's central claim stated
  operationally rather than as a variance percentage.
",
SITE, format(nrow(hold), big.mark = ","), 100 * nrow(hold) / nrow(dat),
elapsed_h, ITER, THREADS, ADAPT_DELTA, ifelse(fold_ok, "PASS", "FAIL"), d,
100 * ceil,
paste(sprintf("%.3f", tab$rmse), collapse = " / "),
paste(sprintf("%.3f", tab$cov90), collapse = " / ")))

hr(); cat("finished: ", format(Sys.time()), "\n"); cat("log: ", LOG, "\n"); hr()
sink(); close(con)
