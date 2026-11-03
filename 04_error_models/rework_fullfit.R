#!/usr/bin/env Rscript
# =====================================================================
# rework_fullfit.R  --  full error model for one product, all sites
# Paper: Table 3a, Table S5 and Text S2 (the full fits; also read by rework_divcheck.R and rework_ppc.R).
#
# The error model (formula, seven priors, land-cover regrouping) fitted
# to one response on the shared training set of rework_model_prep.R:
#   --response chm | dtm | dsm
#
# At the end it also scores the "site known" row of the transfer table:
# the full fit predicts the footprints of each fold's three sites with
# their own fitted site offsets (in sample, an upper bound), so the fold
# runs do not wait on this fit.
#
# OUTPUT
#   models/rework/fit_<resp>_15.rds   (the name rework_divcheck.R and rework_ppc.R read)
#   manuscript_tables/rework/<resp>/  params, fixef, ranef_site, r2, priors,
#                                     oracle_by_fold.csv, fullfit_<date>.log
#
# RUN (one product per run, 4 chains x 8 threads = 32 cores)
#   export PROJECT_ROOT=<CLUSTER_ROOT>/chpt1 R_LIBS=<R_LIBRARY>
#   export CMDSTAN=<CLUSTER_ROOT>/cmdstan/cmdstan-2.37.0
#   export TMPDIR=$PROJECT_ROOT/tmp/$(hostname -s); mkdir -p $TMPDIR
#   Rscript rework_fullfit.R --response dtm     (normally launched by run_rework_queue.sh)
#   --smoke: tiny subsample and 200 iterations, checks the script runs
# =====================================================================

suppressPackageStartupMessages({
  library(brms); library(data.table); library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1]
}
RESP          <- getarg("--response")
ITER          <- as.integer(getarg("--iter", "3000"))
THREADS       <- as.integer(getarg("--threads", "8"))
ADAPT_DELTA   <- as.numeric(getarg("--adapt-delta", "0.98"))
MAX_TREEDEPTH <- as.integer(getarg("--max-treedepth", "14"))
INIT          <- getarg("--init", "random")
SMOKE         <- "--smoke" %in% args
if (is.null(RESP) || !RESP %in% c("chm", "dtm", "dsm")) {
  cat("usage: Rscript rework_fullfit.R --response chm|dtm|dsm\n"); quit(status = 2)
}
YVAR <- paste0(RESP, "_error_mean")
DEFAULT_AD <- 0.98
SUF <- ""
if (abs(ADAPT_DELTA - DEFAULT_AD) > 1e-9) SUF <- paste0(SUF, "_ad", sub("\\.", "", format(ADAPT_DELTA, nsmall = 2)))
if (MAX_TREEDEPTH != 14) SUF <- paste0(SUF, "_td", MAX_TREEDEPTH)
if (INIT != "random")    SUF <- paste0(SUF, "_init", INIT)
if (ITER != 3000L)       SUF <- paste0(SUF, "_it", ITER)
TAG  <- paste0(RESP, SUF, if (SMOKE) "_smoke" else "")
NDRAWS_EVAL <- 1000; CHUNK_ROWS <- 5000

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "<CLUSTER_ROOT>/chpt1")
CHKPT  <- file.path(PROJECT_ROOT, "checkpoints")
OUTMOD <- file.path(PROJECT_ROOT, "models", "rework")
OUTTBL <- file.path(PROJECT_ROOT, "manuscript_tables", "rework", TAG)
dir.create(OUTMOD, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTTBL, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(OUTTBL, sprintf("fullfit_%s.log", format(Sys.Date(), "%Y%m%d")))
con <- file(LOG, open = "wt"); sink(con, split = TRUE)
hr  <- function() cat(strrep("=", 92), "\n")
sec <- function(x) { cat("\n"); hr(); cat("## ", x, "\n"); hr() }
fatal <- function(...) { cat("\nFATAL: ", ..., "\n", sep = ""); sink(); close(con); quit(status = 1) }

cat("rework_fullfit.R  --response ", RESP, if (SMOKE) "  [SMOKE]" else "", "\nstarted: ",
    format(Sys.time()), "\n", sep = "")
cat(sprintf("iterations %d | 4 chains x %d threads | adapt_delta %.3f | max_treedepth %d | init %s\n",
            ITER, THREADS, ADAPT_DELTA, MAX_TREEDEPTH, INIT))

sec("STEP 0  --  WRITE-LOCATION GUARD")
ALLOWED <- "<CLUSTER_ROOT>"   # the user storage root that tempdir() must lie in
cat("tempdir(): ", tempdir(), "\n")
if (!startsWith(normalizePath(tempdir(), mustWork = FALSE), ALLOWED))
  fatal("tempdir() is outside ", ALLOWED, ". export TMPDIR=$PROJECT_ROOT/tmp/$(hostname -s) and relaunch.")

sec("STEP 1  --  DATA (rework_model_prep.R)")
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
miss <- setdiff(MODEL_VARS, names(dat)); if (length(miss)) fatal("columns absent: ", paste(miss, collapse = " "))
n_raw <- nrow(dat); dat <- dat[complete.cases(dat[, ..MODEL_VARS])]
cat("rows ", format(nrow(dat), big.mark = ","), " (", n_raw - nrow(dat), " dropped for incomplete variables), ",
    uniqueN(dat$site), " sites\n", sep = "")
if (!setequal(as.integer(as.character(unique(dat$site))), folds$site)) fatal("training sites and fold sites differ")
if (SMOKE) {
  set.seed(2026)
  dat <- dat[dat[, .I[sample.int(.N, min(.N, 250L))], by = site]$V1]
  ITER <- 200L
  cat("SMOKE: ", nrow(dat), " rows, 200 iterations\n", sep = "")
}
cat(sprintf("response %s: mean %+.3f  sd %.3f\n", YVAR, mean(dat$y), sd(dat$y)))

sec("STEP 2  --  MODEL (the submitted specification)")
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
print(f_mod)

sec("STEP 3  --  FITTING")
fitfile <- file.path(OUTMOD, sprintf("fit_%s_15", TAG))
if (file.exists(paste0(fitfile, ".rds"))) cat("*** ", fitfile, ".rds exists and will be REUSED\n", sep = "")
t0 <- Sys.time()
fit <- brm(f_mod, data = dat, family = student(), prior = priors,
           chains = 4, iter = ITER, warmup = ITER %/% 2, cores = 4, threads = threading(THREADS),
           control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
           init = if (INIT == "random") "random" else as.numeric(INIT),
           backend = "cmdstanr", file = fitfile, file_refit = "on_change", seed = 2026, refresh = 100)
if (nrow(fit$data) != nrow(dat)) fatal("stale fit file ", fitfile, ".rds (rows differ); delete it and rerun")
elapsed_h <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
cat(sprintf("\nWALL CLOCK: %.2f hours\n", elapsed_h))

sec("STEP 4  --  DIAGNOSTICS")
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
fit_ok <- (d == 0) && (mr < 1.01) && (me > 0.10) && (ceil < 0.01)
cat("\n  FULL FIT DIAGNOSTICS: ", ifelse(fit_ok, "PASS", "FAIL"), "\n")

sec("STEP 5  --  SUMMARIES")
s <- as.data.table(posterior::summarise_draws(
  posterior::as_draws_array(fit), median, mean, sd,
  ~posterior::quantile2(.x, probs = c(0.025, 0.975)), posterior::rhat, posterior::ess_bulk))
setnames(s, "variable", "parameter")
fwrite(s, file.path(OUTTBL, "params.csv"))
fwrite(s[grepl("^b_|^bs_|^sds_|^sd_|^cor_|^nu$", parameter)], file.path(OUTTBL, "fixef.csv"))
re <- tryCatch({
  r <- ranef(fit, probs = c(0.025, 0.975))$site[, , "Intercept"]
  data.table(site = rownames(r), median = r[, "Estimate"], q2.5 = r[, "Q2.5"], q97.5 = r[, "Q97.5"])
}, error = function(e) NULL)
if (!is.null(re)) fwrite(re, file.path(OUTTBL, "ranef_site.csv"))
set.seed(2026)
r2t <- bayes_R2(fit, ndraws = min(1000, ndraws(fit)), summary = FALSE)
r2f <- bayes_R2(fit, ndraws = min(1000, ndraws(fit)), summary = FALSE, re_formula = NA)
qq <- function(x) quantile(x, c(0.5, 0.025, 0.975))
r2 <- data.table(quantity = c("R2_total", "R2_fixed_only"),
                 median = c(qq(r2t)[1], qq(r2f)[1]), q2.5 = c(qq(r2t)[2], qq(r2f)[2]),
                 q97.5 = c(qq(r2t)[3], qq(r2f)[3]))
fwrite(r2, file.path(OUTTBL, "r2.csv"))
fwrite(as.data.table(prior_summary(fit)), file.path(OUTTBL, "priors.csv"))
sds <- s[parameter == "sd_site__Intercept"]; rh <- s[parameter == "b_rh_98_z"]
cat(sprintf("  R2 total %.3f | R2 fixed-only %.3f | site SD %.3f | RH98 %+.3f m per SD | nu %.2f\n",
            r2$median[1], r2$median[2], sds$median, rh$median, s[parameter == "nu"]$median))

sec("STEP 6  --  SITE-KNOWN ROW OF THE TRANSFER TABLE, BY FOLD")
score <- function(draws, obs) {
  pm <- colMeans(draws); q <- apply(draws, 2, quantile, probs = c(0.05, 0.95)); res <- pm - obs
  data.table(n = length(obs), rmse = sqrt(mean(res^2)), bias = mean(res), mae = mean(abs(res)),
             cov90 = mean(obs >= q[1, ] & obs <= q[2, ]), width90 = mean(q[2, ] - q[1, ]))
}
orc <- list()
for (k in sort(unique(folds$fold))) {
  hs <- as.character(folds[fold == k]$site)
  h <- dat[as.character(site) %in% hs]
  if (!nrow(h)) next
  set.seed(k)
  starts <- seq(1, nrow(h), by = CHUNK_ROWS)
  pp <- do.call(cbind, lapply(starts, function(a) posterior_predict(
    fit, newdata = h[a:min(a + CHUNK_ROWS - 1, nrow(h))], ndraws = min(NDRAWS_EVAL, ndraws(fit)))))
  sc <- score(pp, h$y); rm(pp); gc()
  orc[[length(orc) + 1L]] <- cbind(data.table(product = toupper(RESP), fold = k,
                                              sites = paste(hs, collapse = " ")), sc,
                                   data.table(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH,
                                              init = INIT, iter = ITER, fit_diagnostics_ok = fit_ok))
  cat(sprintf("  fold %d (sites %s): RMSE %.3f  cov90 %.3f\n", k, paste(hs, collapse = " "), sc$rmse, sc$cov90))
}
fwrite(rbindlist(orc), file.path(OUTTBL, "oracle_by_fold.csv"))
cat(sprintf("\n  diagnostics %s | %.2f h\n", ifelse(fit_ok, "PASS", "FAIL"), elapsed_h))
hr(); cat("finished: ", format(Sys.time()), "\nlog: ", LOG, "\n", sep = ""); hr()
sink(); close(con)
if (!fit_ok && !SMOKE) quit(status = 3)
