#!/usr/bin/env Rscript
# =====================================================================
# rework_fold.R  --  one transfer fold: fit without a fold's sites, predict them
# Paper: Table 3b and Table S2 (held-out sites), fold diagnostics of Text S2.
#
# Six folds of three sites (folds.csv from rework_model_prep.R), the same
# folds for the CHM, DTM and DSM:
#   --response chm | dtm | dsm   --fold 1..K (folds.csv)
#
# Fits the model of rework_fullfit.R (fold sampler settings adapt_delta
# 0.99, max_treedepth 14) on the other fifteen sites and scores the three
# held-out sites:
#   new_site         the model with each held-out site as a new site (its
#                    offset drawn from the fitted between-site spread;
#                    ecoregion and land cover used where training has them)
#   constant_global  one number, the training sites' mean error
#   constant_site    the held-out site's own mean error (in-sample bound)
# The site-known row comes from rework_fullfit.R (oracle_by_fold.csv).
# Metrics are written pooled over the fold and per held-out site.
#
# OUTPUT
#   models/rework/folds/fit_<resp>_fold<k>.rds
#   manuscript_tables/rework/folds/metrics_<resp>_fold<k>.csv
#   manuscript_tables/rework/folds/rowpreds_<resp>_fold<k>.csv.gz
#   manuscript_tables/rework/folds/fold_<resp>_fold<k>_<date>.log
# =====================================================================

suppressPackageStartupMessages({
  library(brms); library(data.table); library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1]
}
RESP          <- getarg("--response")
FOLD          <- as.integer(getarg("--fold", NA))
ITER          <- as.integer(getarg("--iter", "3000"))
THREADS       <- as.integer(getarg("--threads", "8"))
ADAPT_DELTA   <- as.numeric(getarg("--adapt-delta", "0.99"))
MAX_TREEDEPTH <- as.integer(getarg("--max-treedepth", "14"))
INIT          <- getarg("--init", "random")
SMOKE         <- "--smoke" %in% args
if (is.null(RESP) || !RESP %in% c("chm", "dtm", "dsm") || is.na(FOLD)) {
  cat("usage: Rscript rework_fold.R --response chm|dtm|dsm --fold k\n"); quit(status = 2)
}
YVAR <- paste0(RESP, "_error_mean")
DEFAULT_AD <- 0.99
SUF <- ""
if (abs(ADAPT_DELTA - DEFAULT_AD) > 1e-9) SUF <- paste0(SUF, "_ad", sub("\\.", "", format(ADAPT_DELTA, nsmall = 2)))
if (MAX_TREEDEPTH != 14) SUF <- paste0(SUF, "_td", MAX_TREEDEPTH)
if (INIT != "random")    SUF <- paste0(SUF, "_init", INIT)
if (ITER != 3000L)       SUF <- paste0(SUF, "_it", ITER)
TAG  <- sprintf("%s_fold%d%s%s", RESP, FOLD, SUF, if (SMOKE) "_smoke" else "")
NDRAWS_EVAL <- 1000; CHUNK_ROWS <- 5000

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "<CLUSTER_ROOT>/chpt1")
CHKPT  <- file.path(PROJECT_ROOT, "checkpoints")
OUTMOD <- file.path(PROJECT_ROOT, "models", "rework", "folds")
OUTTBL <- file.path(PROJECT_ROOT, "manuscript_tables", "rework", "folds")
dir.create(OUTMOD, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTTBL, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(OUTTBL, sprintf("fold_%s_%s.log", TAG, format(Sys.Date(), "%Y%m%d")))
con <- file(LOG, open = "wt"); sink(con, split = TRUE)
hr  <- function() cat(strrep("=", 92), "\n")
sec <- function(x) { cat("\n"); hr(); cat("## ", x, "\n"); hr() }
fatal <- function(...) { cat("\nFATAL: ", ..., "\n", sep = ""); sink(); close(con); quit(status = 1) }

cat("rework_fold.R  ", TAG, "\nstarted: ", format(Sys.time()), "\n", sep = "")
cat(sprintf("iterations %d | 4 chains x %d threads | adapt_delta %.3f | max_treedepth %d | init %s\n",
            ITER, THREADS, ADAPT_DELTA, MAX_TREEDEPTH, INIT))

sec("STEP 0  --  WRITE-LOCATION GUARD")
ALLOWED <- "<CLUSTER_ROOT>"   # the user storage root that tempdir() must lie in
cat("tempdir(): ", tempdir(), "\n")
if (!startsWith(normalizePath(tempdir(), mustWork = FALSE), ALLOWED))
  fatal("tempdir() is outside ", ALLOWED, ". export TMPDIR=$PROJECT_ROOT/tmp/$(hostname -s) and relaunch.")

sec("STEP 1  --  DATA AND SPLIT")
f09 <- file.path(CHKPT, "09_rework_prep.rds")
if (!file.exists(f09)) fatal("missing ", f09, ". Run rework_model_prep.R first.")
prep  <- readRDS(f09)
dat   <- as.data.table(prep$dat); folds <- as.data.table(prep$folds); rm(prep)
KEEP <- c("ENF", "IWL", "BDF", "RCP", "GRS", "UNK", "DNF", "WTR", "IMP", "EBF")
dat[, lc_grp := factor(fifelse(as.character(lc_l1_code) %in% KEEP, as.character(lc_l1_code), "OTHER"))]
dat[, lc_grp := relevel(lc_grp, ref = "ENF")]
dat[, `:=`(site = factor(site), ecoregion = factor(ecoregion))]
dat[, y := get(YVAR)]
MODEL_VARS <- c("y", "slope_mean_z", "slope_sd_z", "wsci_z", "rh_98_z", "cover_z",
                "aspect_sin_z", "aspect_cos_z", "meta_offnad_z", "meta_sunel_z", "meta_az_conc_z",
                "meta_stereo_z", "meta_leafon_z", "meta_fwdrev_z", "meta_absgeo_z", "meta_relgeo_z",
                "view_az_sin_z", "view_az_cos_z", "site", "ecoregion", "lc_grp")
dat <- dat[complete.cases(dat[, ..MODEL_VARS])]
HS <- as.character(folds[fold == FOLD]$site)
if (!length(HS) %in% 2:4) fatal("fold ", FOLD, " holds ", length(HS), " sites (expected 2 to 4)")
train <- dat[!as.character(site) %in% HS]; hold <- dat[as.character(site) %in% HS]
train[, `:=`(site = droplevels(site), ecoregion = droplevels(ecoregion), lc_grp = droplevels(lc_grp))]
cat(sprintf("held out sites %s: %s rows | training %s rows on %d sites\n",
            paste(HS, collapse = " "), format(nrow(hold), big.mark = ","),
            format(nrow(train), big.mark = ","), uniqueN(train$site)))
eco_new <- setdiff(unique(as.character(hold$ecoregion)), unique(as.character(train$ecoregion)))
cat("held-out ecoregions absent from training (drawn as new): ",
    if (length(eco_new)) paste(eco_new, collapse = " ") else "none", "\n", sep = "")
if (SMOKE) {
  set.seed(2026)
  idx <- train[, .I[sample.int(.N, min(.N, 250L))], by = site]$V1
  for (l in levels(train$lc_grp)) {
    rows <- which(train$lc_grp == l)
    if (length(rows) && sum(idx %in% rows) < 30) idx <- c(idx, rows[sample.int(length(rows), min(length(rows), 30))])
  }
  train <- train[sort(unique(idx))]
  train[, `:=`(site = droplevels(site), ecoregion = droplevels(ecoregion), lc_grp = droplevels(lc_grp))]
  hold <- hold[sample.int(.N, min(.N, 500L))]
  ITER <- 200L
  cat("SMOKE: ", nrow(train), " training rows, ", nrow(hold), " held out, 200 iterations\n", sep = "")
}

nl <- setdiff(as.character(unique(hold$lc_grp)), levels(train$lc_grp))
if (length(nl)) fatal("held-out land-cover levels absent from training (fixed-effect interactions cannot predict them): ",
                      paste(nl, collapse = " "))

sec("STEP 2  --  FIT WITHOUT THE THREE SITES")
main_effects <- paste(
  "slope_mean_z + slope_sd_z + wsci_z + rh_98_z + cover_z +",
  "aspect_sin_z + aspect_cos_z + meta_offnad_z + meta_sunel_z +",
  "meta_az_conc_z + meta_stereo_z + meta_leafon_z + meta_fwdrev_z +",
  "meta_absgeo_z + meta_relgeo_z + view_az_sin_z + view_az_cos_z")
f_mod <- bf(
  as.formula(paste("y ~", main_effects,
                   "+ s(slope_mean_z, wsci_z, k = 20)",
                   "+ slope_mean_z:lc_grp + wsci_z:lc_grp",
                   "+ (1 + slope_mean_z | ecoregion) + (1 | site) + (1 + wsci_z | lc_grp)")),
  sigma ~ 1 + slope_mean_z + (1 | site) + (1 | lc_grp))
priors <- c(
  prior(normal(0, 2),       class = "b"),
  prior(normal(0, 10),      class = "Intercept"),
  prior(gamma(2, 0.1),      class = "nu"),
  prior(student_t(3, 0, 2), class = "sd"),
  prior(normal(0, 1.5),     class = "sds"),
  prior(normal(0, 1),       class = "b",         dpar = "sigma"),
  prior(student_t(3, 0, 5), class = "Intercept", dpar = "sigma"))
fitfile <- file.path(OUTMOD, paste0("fit_", TAG))
if (file.exists(paste0(fitfile, ".rds"))) cat("*** ", fitfile, ".rds exists and will be REUSED\n", sep = "")
t0 <- Sys.time()
fit <- brm(f_mod, data = train, family = student(), prior = priors,
           chains = 4, iter = ITER, warmup = ITER %/% 2, cores = 4, threads = threading(THREADS),
           control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
           init = if (INIT == "random") "random" else as.numeric(INIT),
           backend = "cmdstanr", file = fitfile, file_refit = "on_change", seed = 2026, refresh = 100)
if (nrow(fit$data) != nrow(train)) fatal("stale fit file ", fitfile, ".rds (rows differ); delete it and rerun")
elapsed_h <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
cat(sprintf("\nfold fit took %.2f hours\n", elapsed_h))

sec("STEP 3  --  DIAGNOSTICS")
np <- as.data.table(brms::nuts_params(fit))
td <- np[Parameter == "treedepth__", Value]
d  <- sum(np[Parameter == "divergent__", Value])
mr <- max(brms::rhat(fit), na.rm = TRUE)
me <- min(brms::neff_ratio(fit), na.rm = TRUE)
ceil <- mean(td >= MAX_TREEDEPTH)
cat(sprintf("  divergences     %6d   %s\n", d, ifelse(d == 0, "PASS", "FAIL")))
cat(sprintf("  max Rhat        %6.4f   %s\n", mr, ifelse(mr < 1.01, "PASS", "FAIL")))
cat(sprintf("  min ESS ratio   %6.3f   %s\n", me, ifelse(me > 0.10, "PASS", "FAIL")))
cat(sprintf("  %% at ceiling    %6.2f%%   %s\n", 100 * ceil, ifelse(ceil < 0.01, "PASS", "CHECK")))
tdc <- np[Parameter == "treedepth__", .(pct = round(100 * mean(Value >= MAX_TREEDEPTH), 1)), by = Chain]
if (any(tdc$pct > 50) && !all(tdc$pct > 50)) cat("  ** a chain is stuck at the ceiling: rerun with --init 0 **\n")
if (d > 0) cat("  ** divergences: rerun with --adapt-delta 0.995 **\n")
fold_ok <- (d == 0) && (mr < 1.01) && (me > 0.10) && (ceil < 0.01)
cat("\n  FOLD DIAGNOSTICS: ", ifelse(fold_ok, "PASS", "FAIL"), "\n")

sec("STEP 4  --  PREDICT THE HELD-OUT SITES")
y <- hold$y
set.seed(FOLD)
starts <- seq(1, nrow(hold), by = CHUNK_ROWS)
pp <- do.call(cbind, lapply(starts, function(a) posterior_predict(
  fit, newdata = hold[a:min(a + CHUNK_ROWS - 1, nrow(hold))], ndraws = min(NDRAWS_EVAL, ndraws(fit)),
  allow_new_levels = TRUE, sample_new_levels = "gaussian")))
bad <- which(!is.finite(colSums(pp)))
if (length(bad)) {
  cat("  WARNING: dropping", length(bad), "held-out rows with non-finite predictive draws\n")
  pp <- pp[, -bad, drop = FALSE]; hold <- hold[-bad]; y <- y[-bad]
}
pm <- colMeans(pp); q <- apply(pp, 2, quantile, probs = c(0.05, 0.95)); rm(pp); gc()
g  <- mean(train$y)
rows <- data.table(site = as.character(hold$site), shot_number = as.character(hold$shot_number),
                   eco_known = !(as.character(hold$ecoregion) %in% eco_new),
                   observed = y, pred_new_site = pm, lo90 = q[1, ], hi90 = q[2, ],
                   pred_const_global = g)
rows[, pred_const_site := mean(observed), by = site]
m1 <- function(x) {
  f <- function(p) sqrt(mean((p - x$observed)^2))
  rbind(
    data.table(rung = "new_site", n = nrow(x), rmse = f(x$pred_new_site),
               bias = mean(x$pred_new_site - x$observed), mae = mean(abs(x$pred_new_site - x$observed)),
               cov90 = mean(x$observed >= x$lo90 & x$observed <= x$hi90), width90 = mean(x$hi90 - x$lo90)),
    data.table(rung = "constant_global", n = nrow(x), rmse = f(x$pred_const_global),
               bias = mean(x$pred_const_global - x$observed),
               mae = mean(abs(x$pred_const_global - x$observed)), cov90 = NA_real_, width90 = NA_real_),
    data.table(rung = "constant_site", n = nrow(x), rmse = f(x$pred_const_site),
               bias = mean(x$pred_const_site - x$observed),
               mae = mean(abs(x$pred_const_site - x$observed)), cov90 = NA_real_, width90 = NA_real_))
}
tab <- rbind(cbind(data.table(scope = "fold"), m1(rows)),
             rbindlist(lapply(split(rows, rows$site), function(x) cbind(data.table(scope = x$site[1]), m1(x)))))
tab[, `:=`(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH, init = INIT, iter = ITER)]
eco_by_site <- rows[, .(eco_known = all(eco_known)), by = site]
tab[eco_by_site, eco_known := i.eco_known, on = c(scope = "site")]
tab[scope == "fold", eco_known := all(rows$eco_known)]
tab[, `:=`(product = toupper(RESP), fold = FOLD, held_sites = paste(HS, collapse = " "),
           n_train = nrow(train), fit_hours = elapsed_h, divergences = d, max_rhat = mr,
           min_ess_ratio = me, pct_at_ceiling = 100 * ceil, fold_diagnostics_ok = fold_ok)]
print(tab[, .(scope, rung, n, rmse = round(rmse, 3), bias = round(bias, 3), cov90 = round(cov90, 3))])
fwrite(tab, file.path(OUTTBL, sprintf("metrics_%s.csv", TAG)))
fwrite(rows, file.path(OUTTBL, sprintf("rowpreds_%s.csv.gz", TAG)))

f <- tab[scope == "fold"]
cat(sprintf("\n  %s fold %d (sites %s): model %.3f vs constant %.3f RMSE | cov90 %.3f | diagnostics %s\n",
            toupper(RESP), FOLD, paste(HS, collapse = " "), f[rung == "new_site"]$rmse,
            f[rung == "constant_global"]$rmse, f[rung == "new_site"]$cov90, ifelse(fold_ok, "PASS", "FAIL")))
hr(); cat("finished: ", format(Sys.time()), "\nlog: ", LOG, "\n", sep = ""); hr()
sink(); close(con)
if (!fold_ok && !SMOKE) quit(status = 3)
