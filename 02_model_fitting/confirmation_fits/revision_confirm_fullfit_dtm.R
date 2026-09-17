#!/usr/bin/env Rscript
# =====================================================================
# revision_confirm_fullfit_dtm.R
#
# Terrain counterpart of revision_confirm_fullfit.R. Refits the DTM model
# at full size (202,412 rows before brms drops incomplete cases, 18 sites)
# with the same three changes, none of which alter the science:
#
#   1. Drop the redundant slope_mean_z:wsci_z term.
#      Redundant with s(slope_mean_z, wsci_z, k = 20) over the same pair.
#   2. Group the near-empty land-cover classes into "OTHER".
#      Kept:    ENF IWL BDF RCP GRS UNK DNF WTR IMP EBF
#      Grouped: CWL SHR SVG MFT BAL ICP        (LMS does not occur in DTM)
#      16 classes -> 11.
#   3. Fix the sampler settings: no hardcoded step size, adapt_delta 0.95
#      (published was 0.98).
#
#   Everything else matches the published fit: same rows, same predictors,
#   same priors, same iterations, same family.
#
# THREE PLACES THIS DIFFERS FROM THE CHM SCRIPT, ALL DELIBERATE
#
#   a. NO SITE FILTERING. The CHM script drops sites 1, 2 and 3 because
#      their ALS reference is bad. The DTM analysis keeps them: sites in
#      mod_dtm_s2 are 1 2 3 4 5 6 7 8 9 11 12 13 14 15 17 18 19 20, and
#      202,412 minus 134 incomplete rows equals the published fit's
#      202,278 exactly. That arithmetic only works with zero filtering.
#
#   b. max_treedepth = 14, not 12. The CHM cap of 12 was chosen because
#      the published CHM never reached 13. The published DTM has mean
#      treedepth 11.273, max 13, with 29.5% of draws at 12 or more and
#      0.617% at 13. A cap of 12 would truncate about a third of the
#      trajectories and fail the "<1% at ceiling" check. 14 keeps the cap
#      non-binding, which is the reasoning the CHM choice rested on. The
#      speed comes from adapt_delta and the free step size, not the cap.
#
#   c. One extra prior: lkj_corr_cholesky(2) on class L. The published DTM
#      fit sets this (its prior table shows source = user). Omitting it
#      would be a fourth, undocumented change.
#
# NOTE ON THE BASELINE
#   The published DTM reference land-cover level is BAL, which this script
#   pools into OTHER, so the baseline moves to ENF as it does for CHM.
#   Intercept, slope_mean_z and wsci_z are therefore measured against a
#   different reference class and are reported separately in Step 6 rather
#   than counted as coefficient movement.
#
# RESOURCES
#   4 chains x 16 threads = 64 of 192 cores. Expect 30 to 60 GB RAM.
#
# HOW TO RUN
#   cd $PROJECT_ROOT/scripts/reviewed
#   export TMPDIR=$PROJECT_ROOT/tmp          # REQUIRED, see the guard below
#
#   Smoke test first, 10 to 25 minutes:
#     Rscript revision_confirm_fullfit_dtm.R --smoke
#   Then the real run:
#     nohup nice -n 10 Rscript revision_confirm_fullfit_dtm.R > confirm_dtm.out 2>&1 &
#
#   Check on it:   tail -30 confirm_dtm.out
#   Is it alive?   ps -u $USER -o pid,etime,time,%cpu,comm --sort=-%cpu | head
#   To stop it:    pkill -f revision_confirm_fullfit_dtm
#
# EXPECTED TIME
#   Published slowest chain was 161.04 h. CHM gained 10.0x from the same
#   changes, so expect roughly 16 hours, and 12 to 25 h is unsurprising.
#   Past 40 hours, stop it and look at the treedepth line.
#
# IF IT CRASHES AND YOU RERUN IT
#   It saves progress. Rerunning reuses a completed fit rather than
#   starting over, and says so loudly if it does, because a reused fit
#   makes the timing number meaningless. To force a fresh start, delete:
#   $PROJECT_ROOT/models/confirm/fit_dtm_18_tuned.rds
# =====================================================================

suppressPackageStartupMessages({
  library(brms); library(data.table); library(posterior)
})

SMOKE <- "--smoke" %in% commandArgs(trailingOnly = TRUE)

# ---- the one line to change if you want the CHM cap instead ----------
MAX_TREEDEPTH <- 14

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
CHKPT   <- file.path(PROJECT_ROOT, "checkpoints")
TBL     <- file.path(PROJECT_ROOT, "manuscript_tables")
OUTDIR  <- file.path(PROJECT_ROOT, "models", "confirm")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TBL,    recursive = TRUE, showWarnings = FALSE)

TAG <- if (SMOKE) "_smoke" else ""
LOG <- file.path(TBL, sprintf("confirm_fullfit_dtm%s_%s.log", TAG, format(Sys.Date(), "%Y%m%d")))
con <- file(LOG, open = "wt"); sink(con, split = TRUE)

hr  <- function() cat(strrep("=", 92), "\n")
sec <- function(x) { cat("\n"); hr(); cat("## ", x, "\n"); hr() }
unwrap <- function(p) if (is.list(p) && !is.null(p$data)) p$data else p

cat("revision_confirm_fullfit_dtm.R", if (SMOKE) "  [SMOKE TEST]" else "", "\n")
cat("started: ", format(Sys.time()), "\n")
cat("cores on node: ", parallel::detectCores(), " | using 4 chains x 16 threads = 64\n")


# =====================================================================
sec("STEP 0  --  WRITE-LOCATION GUARD")
# =====================================================================
# R fixes tempdir() at startup from TMPDIR. Setting TMPDIR from inside R
# is too late, so if it was not exported before launching, brms would
# compile the Stan model into /tmp. Nothing here may write outside
# /gpfs/data1/vclgp/lmaden, so this aborts rather than risk it.

ALLOWED <- "/gpfs/data1/vclgp/lmaden"
cat("TMPDIR:    ", Sys.getenv("TMPDIR"), "\n")
cat("tempdir(): ", tempdir(), "\n")
if (!startsWith(normalizePath(tempdir(), mustWork = FALSE), ALLOWED)) {
  cat("\nFATAL: tempdir() is outside", ALLOWED, "\n")
  cat("Run this first, then relaunch:\n")
  cat("  export TMPDIR=$PROJECT_ROOT/tmp\n")
  cat("The bottom-right tmux pane does not set TMPDIR; only the R pane does.\n")
  sink(); close(con); quit(status = 1)
}
cat("OK, all writes stay inside", ALLOWED, "\n")

free_gb <- tryCatch({
  k <- system2("df", c("-Pk", shQuote(PROJECT_ROOT)), stdout = TRUE)
  as.numeric(strsplit(trimws(gsub(" +", " ", k[2])), " ")[[1]][4]) / 1024^2
}, error = function(e) NA_real_)
cat("free space on PROJECT_ROOT: ",
    if (is.na(free_gb)) "unknown" else sprintf("%.1f GB", free_gb), "\n")
if (!SMOKE && !is.na(free_gb) && free_gb < 20) {
  cat("\nFATAL: under 20 GB free. Clear space before an overnight run.\n")
  sink(); close(con); quit(status = 1)
}


# =====================================================================
sec("STEP 1  --  LOAD THE SAME DATA THE PUBLISHED FIT USED")
# =====================================================================

mp <- unwrap(readRDS(file.path(CHKPT, "08_model_prep.rds")))
dat <- as.data.table(mp$mod_dtm_s2)
rm(mp); gc()

# No site exclusions. See note (a) in the header.
cat("rows: ", format(nrow(dat), big.mark = ","), "\n")
cat("sites: ", length(unique(as.character(dat$site))), " -> ",
    paste(sort(unique(as.character(dat$site))), collapse = " "), "\n")
cat("(the published DTM fit used 202,278 rows after brms dropped incomplete cases)\n")
cat("w_dtm present in the data but NOT used: the published fit has no weights() term\n")

stopifnot("dtm_error_mean" %in% names(dat))


# =====================================================================
sec("STEP 2  --  REGROUP LAND COVER")
# =====================================================================

KEEP <- c("ENF", "IWL", "BDF", "RCP", "GRS", "UNK", "DNF", "WTR", "IMP", "EBF")

dat[, lc_grp := factor(fifelse(as.character(lc_l1_code) %in% KEEP,
                               as.character(lc_l1_code), "OTHER"))]
dat[, lc_grp := relevel(lc_grp, ref = "ENF")]   # largest class as baseline

cat("\nBefore and after:\n")
before <- dat[, .N, by = lc_l1_code][order(-N)]
before[, pct := round(100 * N / sum(N), 3)]
before[, kept := fifelse(as.character(lc_l1_code) %in% KEEP, "kept", "-> OTHER")]
print(before)

cat("\nFinal categories used in the model:\n")
after <- dat[, .N, by = lc_grp][order(-N)]
after[, pct := round(100 * N / sum(N), 2)]
print(after)
cat(sprintf("\n%d categories -> %d categories\n", nrow(before), nrow(after)))
fwrite(before, file.path(TBL, "landcover_regrouping_dtm.csv"))

cat("\nPublished DTM baseline was BAL, which is pooled here, so the baseline\n")
cat("moves to ENF. Step 6 reports the affected coefficients separately.\n")

dat[, `:=`(site = factor(site), ecoregion = factor(ecoregion))]


# =====================================================================
sec("STEP 3  --  THE MODEL")
# =====================================================================

main_effects <- paste(
  "slope_mean_z + slope_sd_z + wsci_z + rh_98_z + cover_z +",
  "aspect_sin_z + aspect_cos_z + meta_offnad_z + meta_sunel_z +",
  "meta_az_conc_z + meta_stereo_z + meta_leafon_z + meta_fwdrev_z +",
  "meta_absgeo_z + meta_relgeo_z + view_az_sin_z + view_az_cos_z"
)
missing_pred <- setdiff(trimws(strsplit(gsub("\\+", " ", main_effects), " +")[[1]]), names(dat))
missing_pred <- missing_pred[nzchar(missing_pred)]
if (length(missing_pred)) stop("predictors absent from mod_dtm_s2: ",
                               paste(missing_pred, collapse = " "))

f_tuned <- bf(
  as.formula(paste(
    "dtm_error_mean ~", main_effects,
    "+ s(slope_mean_z, wsci_z, k = 20)",
    "+ slope_mean_z:lc_grp + wsci_z:lc_grp",           # slope:wsci REMOVED
    "+ (1 + slope_mean_z | ecoregion) + (1 | site) + (1 + wsci_z | lc_grp)"
  )),
  sigma ~ 1 + slope_mean_z + (1 | site) + (1 | lc_grp)
)

cat("\nFormula:\n"); print(f_tuned)

priors <- c(
  prior(normal(0, 2),            class = "b"),
  prior(normal(0, 10),           class = "Intercept"),
  prior(gamma(2, 0.1),           class = "nu"),
  prior(student_t(3, 0, 2),      class = "sd"),
  prior(normal(0, 1.5),          class = "sds"),
  prior(lkj_corr_cholesky(2),    class = "L"),          # published DTM sets this
  prior(normal(0, 1),            class = "b",         dpar = "sigma"),
  prior(student_t(3, 0, 5),      class = "Intercept", dpar = "sigma")
)
cat("\nPriors:\n"); print(priors)
cat("brms defaults are deliberately left in place for the variance-submodel\n")
cat("random-effect SDs, matching Supplementary Text S1.1.\n")

cat("\nSampler settings:\n")
cat("  chains 4, iterations 3000 (1500 warmup)  <- same as published\n")
cat("  adapt_delta 0.95        <- published was 0.98\n")
cat(sprintf("  max_treedepth %-2d        <- published was 15; DTM reached 13\n", MAX_TREEDEPTH))
cat("  step_size: not set      <- published forced 0.0005\n")
cat("  threads 16 per chain    <- published used 2\n")

if (SMOKE) {
  set.seed(2026)
  idx <- integer(0)
  for (s in unique(dat$site)) {
    rows <- which(dat$site == s)
    idx <- c(idx, rows[sample.int(length(rows), min(length(rows), 350))])
  }
  for (l in levels(dat$lc_grp)) {                 # keep every level alive
    rows <- which(dat$lc_grp == l)
    if (length(rows) && sum(idx %in% rows) < 30)
      idx <- c(idx, rows[sample.int(length(rows), min(length(rows), 30))])
  }
  dat <- dat[sort(unique(idx))]
  dat[, `:=`(site = droplevels(site), ecoregion = droplevels(ecoregion),
             lc_grp = droplevels(lc_grp))]
  cat(sprintf("\nSMOKE TEST: %s rows, %d sites, %d land-cover levels, 2 chains x 200 iter.\n",
              format(nrow(dat), big.mark = ","), length(levels(dat$site)),
              length(levels(dat$lc_grp))))
  cat("The quality checks below WILL fail at 200 iterations. That is expected.\n")
  cat("The only question is whether the script reaches the end and writes its files.\n")
}


# =====================================================================
sec("STEP 4  --  FITTING (this is the long part)")
# =====================================================================

fitfile <- file.path(OUTDIR, if (SMOKE) "fit_dtm_18_smoke" else "fit_dtm_18_tuned")
was_cached <- file.exists(paste0(fitfile, ".rds"))
if (was_cached) {
  cat("\n*** A saved fit already exists at", paste0(fitfile, ".rds"), "\n")
  cat("*** brms will LOAD it instead of refitting, so the wall-clock number\n")
  cat("*** below will be meaningless. Delete that file for a real timing run.\n\n")
}

t0 <- Sys.time()
cat("fitting started: ", format(t0), "\n\n")

fit <- brm(
  formula = f_tuned,
  data    = dat,
  family  = student(),
  prior   = priors,
  chains  = if (SMOKE) 2    else 4,
  iter    = if (SMOKE) 200  else 3000,
  warmup  = if (SMOKE) 100  else 1500,
  cores   = if (SMOKE) 2    else 4,
  threads = threading(if (SMOKE) 2 else 16),
  control = list(adapt_delta = 0.95, max_treedepth = MAX_TREEDEPTH),
  backend = "cmdstanr",
  file    = fitfile,
  seed    = 2026, refresh = 100
)

elapsed_h <- as.numeric(difftime(Sys.time(), t0, units = "hours"))


# =====================================================================
sec("STEP 5  --  RESULTS")
# =====================================================================

np <- as.data.table(brms::nuts_params(fit))
td <- np[Parameter == "treedepth__", Value]
dv <- np[Parameter == "divergent__", Value]

# The published 161.04 h is a slowest-chain figure, so report both.
et <- tryCatch(rstan::get_elapsed_time(fit$fit), error = function(e) NULL)
warm_h <- samp_h <- slow_h <- NA_real_
if (!is.null(et)) {
  slow_h <- max(rowSums(et)) / 3600
  warm_h <- max(et[, "warmup"]) / 3600
  samp_h <- max(et[, "sample"]) / 3600
  cat("elapsed hours per chain:\n")
  print(round(cbind(et / 3600, total = rowSums(et) / 3600), 3))
}

cat(sprintf("\nWALL CLOCK:            %.2f hours   <-- THE NUMBER WE NEEDED\n", elapsed_h))
if (!is.na(slow_h)) cat(sprintf("slowest chain:         %.2f hours\n", slow_h))
cat(sprintf("published fit took:    161.04 hours (slowest chain)\n"))
cat(sprintf("SPEEDUP:               %.1f x\n\n",
            161.04 / max(if (is.na(slow_h)) elapsed_h else slow_h, 1e-6)))

cat("QUALITY CHECKS (all must pass):\n")
d  <- sum(dv)
mr <- max(brms::rhat(fit), na.rm = TRUE)
me <- min(brms::neff_ratio(fit), na.rm = TRUE)
ceil <- mean(td >= MAX_TREEDEPTH)
cat(sprintf("  divergences        %6d      %s  (need 0)\n",
            d, ifelse(d == 0, "PASS", "FAIL")))
cat(sprintf("  max Rhat           %6.4f      %s  (need < 1.01)\n",
            mr, ifelse(mr < 1.01, "PASS", "FAIL")))
cat(sprintf("  min ESS ratio      %6.3f      %s  (need > 0.10)\n",
            me, ifelse(me > 0.10, "PASS", "FAIL")))
cat(sprintf("  mean treedepth     %6.2f            (published was 11.27)\n", mean(td)))
cat(sprintf("  %% at ceiling        %6.1f%%      %s  (need < 1%%)\n",
            100 * ceil, ifelse(ceil < 0.01, "PASS", "CHECK")))
cat(sprintf("  %% at 12 or more     %6.1f%%            (published was 29.5%%)\n",
            100 * mean(td >= 12)))

if (d > 0) {
  cat("\n  ** Divergences present. Rerun with adapt_delta = 0.98.\n")
  cat("     It will be slower but still far faster than 161 hours. **\n")
}
if (ceil >= 0.01) {
  cat(sprintf("\n  ** Treedepth saturating at %d. Raise MAX_TREEDEPTH near the top\n",
              MAX_TREEDEPTH))
  cat("     of this script and rerun. Do not accept it: a truncated\n")
  cat("     trajectory can bias the posterior. **\n")
}


# =====================================================================
sec("STEP 6  --  DID SIMPLIFYING CHANGE ANY CONCLUSIONS?")
# =====================================================================

pub_path <- file.path(CHKPT, "10_models_stage2.rds")
if (file.exists(pub_path)) {
  pub <- unwrap(readRDS(pub_path))$fit_dtm_s2

  gx <- function(f) {
    s <- as.data.table(posterior::summarise_draws(
      posterior::as_draws_df(f), mean, ~quantile(.x, c(0.025, 0.975))))
    setnames(s, c("variable", "mean", "lo", "hi"))
    s[grepl("^b_(?!sigma)", variable, perl = TRUE)]
  }
  a <- gx(pub); b <- gx(fit)
  cmp <- merge(a, b, by = "variable", suffixes = c("_pub", "_new"))
  cmp[, shift_SD := round((mean_new - mean_pub) / ((hi_pub - lo_pub) / 3.92), 2)]
  cmp[, sign_flip := sign(mean_pub) != sign(mean_new) &
        (lo_pub > 0 | hi_pub < 0) & (lo_new > 0 | hi_new < 0)]

  # Baseline moved from BAL to ENF, so these three are not like-for-like.
  BASELINE_DEP <- c("b_Intercept", "b_slope_mean_z", "b_wsci_z")
  cmp[, baseline_shifted := variable %in% BASELINE_DEP]
  cmp <- cmp[order(-baseline_shifted, -abs(shift_SD))]

  cat("\nCoefficients present in both models, largest movement first.\n")
  cat("'shift_SD' = how far the estimate moved, in units of its own\n")
  cat("uncertainty. Under 0.5 is negligible. Over 1.0 needs a look.\n")
  cat("Rows marked baseline_shifted are measured against a different\n")
  cat("reference land-cover class and are NOT evidence of drift.\n\n")
  print(cmp[, .(variable, mean_pub = round(mean_pub, 3),
                mean_new = round(mean_new, 3), shift_SD, sign_flip,
                baseline_shifted)][1:30])

  core <- cmp[baseline_shifted == FALSE]
  cat(sprintf("\n  comparable coefficients ............... %d\n", nrow(core)))
  cat(sprintf("  of those, moved more than 1 SD ....... %d\n",
              sum(abs(core$shift_SD) > 1, na.rm = TRUE)))
  cat(sprintf("  CREDIBLE SIGN FLIPS .................. %d   <-- must be 0\n",
              sum(core$sign_flip, na.rm = TRUE)))
  cat(sprintf("  baseline-shifted, reported separately  %d\n", nrow(cmp) - nrow(core)))

  rh <- cmp[variable == "b_rh_98_z"]
  if (nrow(rh)) cat(sprintf("\n  RH98 effect (the paper's key DTM number): published %.3f -> new %.3f\n",
                            rh$mean_pub, rh$mean_new))
  fwrite(cmp, file.path(TBL, sprintf("confirm_coefficient_comparison_dtm%s.csv", TAG)))
  rm(pub); gc()
} else {
  cat("published fit not found, skipping comparison\n")
}


# =====================================================================
sec("STEP 7  --  WHAT THIS MEANS FOR THE LOO BUDGET")
# =====================================================================

if (!SMOKE && !is.na(warm_h)) {
  n_tot <- nrow(dat)
  tab <- dat[, .N, by = site][order(-N)]
  tab[, share := round(N / n_tot, 4)]
  cat("Site shares of the DTM fitting data. A fold that holds out a large\n")
  cat("site is cheap; one that holds out a small site costs nearly a full fit,\n")
  cat("so the wall clock is set by the most expensive folds, not the average.\n\n")
  print(tab)

  fold_frac <- 1 - tab$N / n_tot
  hoursA <- fold_frac * (warm_h * 0.5 + samp_h * 0.5)   # 1500 total: 750 + 750
  hoursB <- fold_frac * (warm_h * 1.0 + samp_h * 1.0)   # 1500 post-warmup, 3000 total

  makespan <- function(h, slots) {                      # longest-job-first packing
    h <- sort(h, decreasing = TRUE); load <- numeric(slots)
    for (x in h) { i <- which.min(load); load[i] <- load[i] + x }
    max(load)
  }

  cat(sprintf("\nmeasured on the slowest chain: warmup %.2f h, sampling %.2f h\n",
              warm_h, samp_h))
  cat("\n'1500 iterations per fold' is ambiguous, so both readings are shown:\n")
  cat("  A = 1500 total (750 warmup + 750 post-warmup)\n")
  cat("  B = 1500 post-warmup on top of 1500 warmup, i.e. 3000 as in this fit\n")
  cat(sprintf("\nDTM, %d folds:\n", length(fold_frac)))
  cat(sprintf("  A: total %.1f h, longest fold %.2f h\n", sum(hoursA), max(hoursA)))
  cat(sprintf("  B: total %.1f h, longest fold %.2f h\n", sum(hoursB), max(hoursB)))
  cat(sprintf("\nEach fold needs 64 cores, so one fold per node.\n"))

  out <- rbindlist(lapply(c(1, 4, 5, 8), function(slots) {
    mA <- makespan(hoursA, slots); mB <- makespan(hoursB, slots)
    chmA <- 5.81 * 0.5 * 16 / slots; chmB <- 5.81 * 16 / slots
    cat(sprintf("  %d slot(s): DTM %.1f h (A) / %.1f h (B); CHM approx %.1f / %.1f h; both %.1f / %.1f days\n",
                slots, mA, mB, chmA, chmB, (mA + chmA) / 24, (mB + chmB) / 24))
    data.table(parallel_slots = slots, dtm_h_A = mA, dtm_h_B = mB,
               chm_h_A = chmA, chm_h_B = chmB,
               total_days_A = (mA + chmA) / 24, total_days_B = (mB + chmB) / 24)
  }))
  fwrite(out, file.path(TBL, "loo_budget_projection.csv"))

  cat("\nCaveats. Cost is taken as proportional to rows, which holds per\n")
  cat("leapfrog step. The CHM figures scale the known 5.81 h fit by fold\n")
  cat("count and divide by slots; because site 5 is 56.3% of CHM training\n")
  cat("data, real CHM fold costs are very unequal and the true CHM makespan\n")
  cat("will be worse than a flat divide. Recompute with CHM site counts\n")
  cat("before quoting a number to anyone.\n")
} else if (SMOKE) {
  cat("skipped in smoke mode\n")
}


sec("SUMMARY FOR THE NEXT CONVERSATION")
cat(sprintf("
  hours for one full DTM fit ............ %.2f  (slowest chain %.2f)
  divergences ........................... %d
  max Rhat .............................. %.4f
  %% of draws at the treedepth ceiling ... %.1f%%
  comparable coefficients moved > 1 SD .. %s
  credible sign flips ................... %s

  If sign flips is 0, divergences is 0, and the ceiling percentage is
  under 1, this specification is safe and the LOO plan can be sized from
  the hours figure in Step 7.
",
elapsed_h, ifelse(is.na(slow_h), elapsed_h, slow_h), d, mr, 100 * ceil,
ifelse(exists("core"), sum(abs(core$shift_SD) > 1, na.rm = TRUE), "n/a"),
ifelse(exists("core"), sum(core$sign_flip, na.rm = TRUE), "n/a")))

hr(); cat("finished: ", format(Sys.time()), "\n")
cat("log: ", LOG, "\n"); hr()
sink(); close(con)
