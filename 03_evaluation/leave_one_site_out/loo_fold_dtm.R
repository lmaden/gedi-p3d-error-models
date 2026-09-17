#!/usr/bin/env Rscript
# =====================================================================
# loo_fold_dtm.R  --  one leave-one-site-out fold for the DTM model
#
# WHAT IT DOES
#   Refits the tuned DTM model with one site removed, then predicts that
#   site's footprints under three levels of prior knowledge, and scores
#   the predictions. One site per run, so folds can be spread over nodes.
#
# THE THREE RUNGS
#   Each answers "how well can we predict a place we have not measured?",
#   with progressively less knowledge available.
#
#   1. ORACLE. Predictions come from the full 18-site fit, in sample.
#      This is what you would get if the site's offset were known in
#      advance. It is deliberately optimistic and NOT a validation: the
#      offset was estimated from the very footprints being predicted.
#      It exists to put the site effect in metres instead of variance
#      percentages, and must be labelled as an upper bound wherever it
#      appears. Taken from the full fit rather than grafted onto the fold
#      fit, so that everything comes from one posterior.
#
#   2. NEW SITE, ECOREGION AND LAND COVER KNOWN. Fold fit. The site is a
#      new level, drawn from the fitted between-site distribution; the
#      ecoregion and land-cover effects are the estimated ones.
#
#   3. NEW SITE AND NEW ECOREGION, LAND COVER KNOWN. Fold fit, with new
#      levels drawn for both site and ecoregion. This asks whether the
#      model transfers outside the nine sampled ecoregions.
#      Land cover is deliberately never treated as unknown: it is a
#      fixed-effect covariate read from GLC-FCS30, available for any pixel,
#      and a new level could not enter the design matrix in any case.
#
#   If the held-out site is the only site in its ecoregion, rung 2 is not
#   available, since ecoregion is then also a new level. The script
#   detects this and marks the rung accordingly rather than mislabelling.
#
# MODEL SPECIFICATION
#   Identical to revision_confirm_fullfit_dtm.R, which passed on 26 July:
#   19.61 h, zero divergences, max Rhat 1.0050, zero sign flips.
#   adapt_delta 0.98, max_treedepth 14, student(), same priors.
#
# USAGE
#   export PROJECT_ROOT=/gpfs/data1/vclgp/lmaden/chpt1
#   export TMPDIR=$PROJECT_ROOT/tmp
#
#   Rscript loo_fold_dtm.R --site 18 --iter 1500 --threads 8
#
#   --site     required, the site to hold out
#   --iter     total iterations, half are warmup. Default 3000.
#   --threads  threads per chain. Default 16 (= 64 cores for 4 chains).
#   --adapt-delta  sampler target acceptance. Default 0.98, matching the
#              confirmation fit. A non-default value adds a suffix to every
#              output name (0.99 -> _ad099), so reruns never overwrite the
#              originals and both remain on disk for comparison.
#   --smoke    tiny subsample, for checking the script runs
#
# EXPECTED TIME (from the 19.61 h full fit, scaled by rows retained)
#   site 18 (23 rows held out):     ~9.8 h at 1500 iter, 16 threads
#   site 1  (76,631 rows held out): ~6.1 h at 1500 iter, 16 threads
#   Add roughly a third at 8 threads.
# =====================================================================

suppressPackageStartupMessages({
  library(brms); library(data.table); library(posterior)
})

# ---- arguments ------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
SITE    <- getarg("--site")
ITER    <- as.integer(getarg("--iter", "3000"))
THREADS <- as.integer(getarg("--threads", "16"))
ADAPT_DELTA <- as.numeric(getarg("--adapt-delta", "0.98"))
SMOKE   <- "--smoke" %in% args

if (is.null(SITE)) {
  cat("\nYou must name a site to hold out, for example:\n")
  cat("  Rscript loo_fold_dtm.R --site 18 --iter 1500 --threads 8\n\n")
  quit(status = 2)
}
WARMUP <- ITER %/% 2

MAX_TREEDEPTH <- 14
NDRAWS_EVAL   <- 1000    # draws used for scoring; 1000 is ample and keeps memory sane
CHUNK_ROWS    <- 5000    # predict in blocks so a 76k-row site does not blow memory

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
CHKPT   <- file.path(PROJECT_ROOT, "checkpoints")
FULLFIT <- file.path(PROJECT_ROOT, "models", "confirm", "fit_dtm_18_tuned.rds")
OUTMOD  <- file.path(PROJECT_ROOT, "models", "loo")
OUTTBL  <- file.path(PROJECT_ROOT, "manuscript_tables", "loo")
dir.create(OUTMOD, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTTBL, recursive = TRUE, showWarnings = FALSE)

AD_TAG <- if (abs(ADAPT_DELTA - 0.98) < 1e-9) "" else
            paste0("_ad", sub("\\.", "", format(ADAPT_DELTA, nsmall = 2)))
TAG <- sprintf("dtm_site%s%s%s", SITE, AD_TAG, if (SMOKE) "_smoke" else "")
LOG <- file.path(OUTTBL, sprintf("loo_%s_%s.log", TAG, format(Sys.Date(), "%Y%m%d")))
con <- file(LOG, open = "wt"); sink(con, split = TRUE)

hr  <- function() cat(strrep("=", 92), "\n")
sec <- function(x) { cat("\n"); hr(); cat("## ", x, "\n"); hr() }
unwrap <- function(p) if (is.list(p) && !is.null(p$data)) p$data else p

cat("loo_fold_dtm.R  --  holding out site ", SITE, if (SMOKE) "  [SMOKE]" else "", "\n", sep = "")
cat("started: ", format(Sys.time()), "\n")
cat(sprintf("iterations %d (%d warmup) | 4 chains x %d threads = %d cores\n",
            ITER, WARMUP, THREADS, 4 * THREADS))
cat(sprintf("adapt_delta %.3f | max_treedepth %d | output tag: %s\n",
            ADAPT_DELTA, MAX_TREEDEPTH, TAG))


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
dat <- as.data.table(mp$mod_dtm_s2)
rm(mp); gc()

KEEP <- c("ENF", "IWL", "BDF", "RCP", "GRS", "UNK", "DNF", "WTR", "IMP", "EBF")
dat[, lc_grp := factor(fifelse(as.character(lc_l1_code) %in% KEEP,
                               as.character(lc_l1_code), "OTHER"))]
dat[, lc_grp := relevel(lc_grp, ref = "ENF")]
dat[, `:=`(site = factor(site), ecoregion = factor(ecoregion))]

if (!(SITE %in% levels(dat$site))) {
  cat("\nFATAL: site ", SITE, " is not in mod_dtm_s2. Available: ",
      paste(levels(dat$site), collapse = " "), "\n", sep = "")
  sink(); close(con); quit(status = 1)
}

# Drop rows with NAs in any model variable. brms silently excludes them
# from the fit, but posterior_predict returns NaN for them, which is what
# killed the site 1 pilot during scoring.
MODEL_VARS <- c("dtm_error_mean", "slope_mean_z", "slope_sd_z", "wsci_z",
                "rh_98_z", "cover_z", "aspect_sin_z", "aspect_cos_z",
                "meta_offnad_z", "meta_sunel_z", "meta_az_conc_z",
                "meta_stereo_z", "meta_leafon_z", "meta_fwdrev_z",
                "meta_absgeo_z", "meta_relgeo_z", "view_az_sin_z",
                "view_az_cos_z", "site", "ecoregion", "lc_grp")
n_raw <- nrow(dat)
dat <- dat[complete.cases(dat[, ..MODEL_VARS])]
cat("rows dropped for incomplete model variables: ", n_raw - nrow(dat), "\n")

train <- dat[site != SITE]
hold  <- dat[site == SITE]
cat(sprintf("training  %s rows, %d sites\n",
            format(nrow(train), big.mark = ","),
            length(unique(as.character(train$site)))))
cat(sprintf("held out  %s rows (%.2f%% of the data)\n",
            format(nrow(hold), big.mark = ","), 100 * nrow(hold) / nrow(dat)))

# Is the held-out site alone in its ecoregion? If so rung 2 is unavailable.
eco_held <- unique(as.character(hold$ecoregion))
eco_in_train <- eco_held %in% unique(as.character(train$ecoregion))
ECO_KNOWN <- all(eco_in_train)
cat("held-out site's ecoregion(s): ", paste(eco_held, collapse = " "), "\n")
cat("also present in training data: ", ECO_KNOWN, "\n")
if (!ECO_KNOWN) {
  cat("\n>> This site is the only one in its ecoregion. Rung 2 cannot use a\n")
  cat(">> known ecoregion effect, so it collapses into rung 3. Both are still\n")
  cat(">> computed; rung 2 is flagged eco_known = FALSE in the output.\n")
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
  cat(sprintf("\nSMOKE: training %s rows, held out %d rows, 200 iterations.\n",
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
    "dtm_error_mean ~", main_effects,
    "+ s(slope_mean_z, wsci_z, k = 20)",
    "+ slope_mean_z:lc_grp + wsci_z:lc_grp",
    "+ (1 + slope_mean_z | ecoregion) + (1 | site) + (1 + wsci_z | lc_grp)"
  )),
  sigma ~ 1 + slope_mean_z + (1 | site) + (1 | lc_grp)
)
priors <- c(
  prior(normal(0, 2),         class = "b"),
  prior(normal(0, 10),        class = "Intercept"),
  prior(gamma(2, 0.1),        class = "nu"),
  prior(student_t(3, 0, 2),   class = "sd"),
  prior(normal(0, 1.5),       class = "sds"),
  prior(lkj_corr_cholesky(2), class = "L"),
  prior(normal(0, 1),         class = "b",         dpar = "sigma"),
  prior(student_t(3, 0, 5),   class = "Intercept", dpar = "sigma")
)

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
  backend = "cmdstanr", file = fitfile, seed = 2026, refresh = 100
)
elapsed_h <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
cat(sprintf("\nfold fit took %.2f hours\n", elapsed_h))


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
cat(sprintf("  mean treedepth  %6.2f          (full fit was 11.3 at 3000 iter)\n", mean(td)))
cat(sprintf("  %% at ceiling     %6.1f%%   %s  (need < 1%%)\n",
            100 * ceil, ifelse(ceil < 0.01, "PASS", "CHECK")))
if (d > 0 && ITER < 3000) {
  cat("\n  ** Divergences at reduced iterations. If this is the pilot, that is\n")
  cat("     the answer: warmup is too short and the full run needs 3000. **\n")
}
fold_ok <- (d == 0) && (mr < 1.01) && (me > 0.10) && (ceil < 0.01)


# =====================================================================
sec("STEP 4  --  PREDICT THE HELD-OUT SITE UNDER THREE RUNGS")
# =====================================================================

y <- hold$dtm_error_mean

# CRPS from posterior draws, via the sorted-sample identity. O(m log m)
# per row, so it stays cheap. draws: matrix of m draws x n rows.
crps_cols <- function(draws, obs) {
  m <- nrow(draws)
  vapply(seq_along(obs), function(j) {
    x <- sort(draws[, j]); yy <- obs[j]
    i <- seq_len(m)
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
  pm  <- colMeans(draws)
  q   <- apply(draws, 2, quantile, probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
  res <- pm - obs
  list(
    summary = data.table(
      n          = length(obs),
      rmse       = sqrt(mean(res^2)),
      bias       = mean(res),
      mae        = mean(abs(res)),
      crps       = mean(crps_cols(draws, obs)),
      cov50      = mean(obs >= q[3, ] & obs <= q[4, ]),
      cov90      = mean(obs >= q[2, ] & obs <= q[5, ]),
      cov95      = mean(obs >= q[1, ] & obs <= q[6, ]),
      width90    = mean(q[5, ] - q[2, ])
    ),
    rows = data.table(pred_mean = pm, lo90 = q[2, ], hi90 = q[5, ])
  )
}

# Predict in row-blocks and stitch, so a 76k-row site does not exhaust RAM.
predict_blocks <- function(model, newdata, ...) {
  n <- nrow(newdata)
  nd <- min(NDRAWS_EVAL, brms::ndraws(model))   # a short fit has fewer draws
  cat("    scoring with", nd, "posterior draws\n")
  starts <- seq(1, n, by = CHUNK_ROWS)
  out <- vector("list", length(starts))
  for (k in seq_along(starts)) {
    rng <- starts[k]:min(starts[k] + CHUNK_ROWS - 1, n)
    out[[k]] <- posterior_predict(model, newdata = newdata[rng],
                                  ndraws = nd, ...)
    if (k %% 5 == 0) cat("    block", k, "of", length(starts), "\n")
  }
  do.call(cbind, out)
}

results <- list(); rowpreds <- list()

# ---- rung 1: oracle, from the full 18-site fit ----------------------
cat("\n[rung 1] oracle: full 18-site fit, in sample. UPPER BOUND, NOT VALIDATION.\n")
if (file.exists(FULLFIT)) {
  full <- readRDS(FULLFIT)
  set.seed(1)
  pp <- predict_blocks(full, hold)
  s  <- score(pp, y); rm(pp, full); gc()
  results[["1_oracle_site_known"]] <- s$summary
  rowpreds[["1_oracle_site_known"]] <- s$rows
  cat(sprintf("  RMSE %.3f  bias %+.3f  CRPS %.3f  cov90 %.3f\n",
              s$summary$rmse, s$summary$bias, s$summary$crps, s$summary$cov90))
} else {
  cat("  SKIPPED: ", FULLFIT, " not found\n", sep = "")
}

# ---- rung 2: new site, ecoregion and land cover known ---------------
cat("\n[rung 2] new site; ecoregion and land cover known\n")
set.seed(2)
pp <- predict_blocks(fit, hold, allow_new_levels = TRUE,
                     sample_new_levels = "gaussian")
s <- score(pp, y); rm(pp); gc()
s$summary[, eco_known := ECO_KNOWN]
results[["2_new_site_eco_lc_known"]] <- s$summary
rowpreds[["2_new_site_eco_lc_known"]] <- s$rows
cat(sprintf("  RMSE %.3f  bias %+.3f  CRPS %.3f  cov90 %.3f%s\n",
            s$summary$rmse, s$summary$bias, s$summary$crps, s$summary$cov90,
            if (ECO_KNOWN) "" else "   [ecoregion NOT actually known]"))

# ---- rung 3: nothing known ------------------------------------------
cat("\n[rung 3] new site AND new ecoregion; land cover known\n")
hold3 <- copy(hold)
# lc_grp is deliberately left alone: it enters the FIXED effects through
# slope_mean_z:lc_grp and wsci_z:lc_grp, so a new level cannot be put in
# the design matrix. Land cover is also always knowable from GLC-FCS30, so
# treating it as unknown was never the right question.
hold3[, `:=`(site      = factor("__new_site__"),
             ecoregion = factor("__new_eco__"))]
set.seed(3)
pp <- predict_blocks(fit, hold3, allow_new_levels = TRUE,
                     sample_new_levels = "gaussian")
s <- score(pp, y); rm(pp); gc()
results[["3_new_site_new_ecoregion"]] <- s$summary
rowpreds[["3_new_site_new_ecoregion"]] <- s$rows
cat(sprintf("  RMSE %.3f  bias %+.3f  CRPS %.3f  cov90 %.3f\n",
            s$summary$rmse, s$summary$bias, s$summary$crps, s$summary$cov90))
cat("  This asks whether the model transfers beyond the nine sampled\n")
cat("  ecoregions. It differs from rung 2 only in the ecoregion random\n")
cat("  effect, since land cover stays known in both.\n")


# =====================================================================
sec("STEP 5  --  WRITE RESULTS")
# =====================================================================

tab <- rbindlist(lapply(names(results), function(k)
  cbind(data.table(product = "DTM", site = SITE, rung = k), results[[k]])),
  fill = TRUE)
tab[, `:=`(n_train = nrow(train), n_held = nrow(hold), iter = ITER,
           threads = THREADS, fit_hours = elapsed_h, divergences = d,
           max_rhat = mr, min_ess_ratio = me, pct_at_ceiling = 100 * ceil,
           fold_diagnostics_ok = fold_ok)]
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
  fold fit hours ............. %.2f   at %d iterations, %d threads
  fold diagnostics ........... %s
  divergences ................ %d
  RMSE, oracle / rung2 / rung3 %s
  cov90, oracle / rung2 / rung3 %s

  The gap between the oracle and rung 2 is what knowing the site is worth,
  in metres. That gap is the paper's argument, stated operationally.
  Coverage well below 0.90 at rungs 2 and 3 is the direct evidence on R3.
",
SITE, format(nrow(hold), big.mark = ","), 100 * nrow(hold) / nrow(dat),
elapsed_h, ITER, THREADS, ifelse(fold_ok, "PASS", "FAIL"), d,
paste(sprintf("%.3f", tab$rmse), collapse = " / "),
paste(sprintf("%.3f", tab$cov90), collapse = " / ")))

hr(); cat("finished: ", format(Sys.time()), "\n"); cat("log: ", LOG, "\n"); hr()
sink(); close(con)
