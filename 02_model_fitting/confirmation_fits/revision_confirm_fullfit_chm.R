#!/usr/bin/env Rscript
# =====================================================================
# revision_confirm_fullfit_chm.R
#
# WHAT THIS DOES
#   Refits the CHM model at FULL SIZE (all 124,966 rows, 16 sites) with
#   three changes, none of which alter the science:
#
#     1. Drop the redundant slope:wsci term.
#        (The speed experiment showed it contributes nothing: LOO score
#         identical to the published model.)
#     2. Group the 7 near-empty land-cover classes into "OTHER".
#        Kept separately: ENF IWL BDF RCP GRS UNK DNF WTR IMP EBF
#        Grouped:         CWL SHR SVG MFT BAL ICP LMS  (~2% of data)
#        (LMS currently has 2 observations and 5 parameters.)
#     3. Fix the sampler settings: remove the hardcoded tiny step size,
#        set adapt_delta to 0.95 (0.90 caused a few errors in testing;
#        the published value was 0.98).
#
#   Everything else is identical to the published fit: same data, same
#   predictors, same priors, same number of iterations.
#
# WHY
#   To get a real wall-clock number instead of an extrapolation, and to
#   confirm the simplified model gives the SAME ANSWERS as the published
#   one. The script prints a side-by-side coefficient comparison at the
#   end so we can see if anything moved.
#
# NOTE ON THE DATA
#   This uses the existing Stage 2 sample, the same one behind the
#   published fit, so the comparison is apples to apples. The known
#   19.9% overlap with Stage 1 gets fixed later, when we build the clean
#   three-way split. Do not worry about it here.
#
# RESOURCES
#   4 chains x 16 threads = 64 of 192 cores. Roughly 20-40 GB RAM.
#
# HOW TO RUN
#   cd $PROJECT_ROOT/scripts/reviewed
#   export TMPDIR=$PROJECT_ROOT/tmp
#   nohup nice -n 10 Rscript revision_confirm_fullfit.R > confirm.out 2>&1 &
#
#   Check on it:   tail -30 confirm.out
#   Is it alive?   ps -u $USER -o pid,etime,time,%cpu,comm --sort=-%cpu | head
#
# EXPECTED TIME
#   2 to 6 hours. If it passes 10 hours, stop it and tell me.
#   To stop it:  pkill -f revision_confirm_fullfit
#
# IF IT CRASHES AND YOU RERUN IT
#   It saves its progress. Rerunning will reuse a completed fit rather
#   than starting over. To force a fresh start, delete:
#   $PROJECT_ROOT/models/confirm/fit_chm_16_tuned.rds
# =====================================================================

suppressPackageStartupMessages({
  library(brms); library(data.table); library(posterior)
})

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
CHKPT   <- file.path(PROJECT_ROOT, "checkpoints")
TBL     <- file.path(PROJECT_ROOT, "manuscript_tables")
OUTDIR  <- file.path(PROJECT_ROOT, "models", "confirm")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TBL,    recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(TBL, sprintf("confirm_fullfit_%s.log", format(Sys.Date(), "%Y%m%d")))
con <- file(LOG, open = "wt"); sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

hr  <- function() cat(strrep("=", 92), "\n")
sec <- function(x) { cat("\n"); hr(); cat("## ", x, "\n"); hr() }
unwrap <- function(p) if (is.list(p) && !is.null(p$data)) p$data else p

cat("revision_confirm_fullfit.R\n")
cat("started: ", format(Sys.time()), "\n")
cat("cores on node: ", parallel::detectCores(), " | using 4 chains x 16 threads = 64\n")


# =====================================================================
sec("STEP 1  --  LOAD THE SAME DATA THE PUBLISHED FIT USED")
# =====================================================================

mp <- unwrap(readRDS(file.path(CHKPT, "08_model_prep.rds")))
dat <- as.data.table(mp$mod_chm_s2)
rm(mp); gc()

SITES_FLAGGED <- c("1", "2", "3")          # dropped from the CHM analysis
dat <- dat[!as.character(site) %in% SITES_FLAGGED]

cat("rows after excluding the 3 flagged sites: ", format(nrow(dat), big.mark = ","), "\n")
cat("sites: ", length(unique(as.character(dat$site))), "\n")
cat("(the published fit used 124,966 rows after dropping incomplete rows)\n")


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
fwrite(before, file.path(TBL, "landcover_regrouping.csv"))

# free diagnostic: what is UNK?
cat("\nWHAT IS 'UNK'? (unknown land cover, ~2.6% of the data)\n")
unk <- dat[lc_l1_code == "UNK"]
if (nrow(unk) > 0) {
  cat("  which sites it appears at:\n")
  print(unk[, .N, by = site][order(-N)])
  cat(sprintf("  median canopy cover: %.3f   (all data: %.3f)\n",
              median(unk$cover, na.rm = TRUE), median(dat$cover, na.rm = TRUE)))
  cat(sprintf("  median GEDI RH98:    %.2f m  (all data: %.2f m)\n",
              median(unk$rh_98, na.rm = TRUE), median(dat$rh_98, na.rm = TRUE)))
  cat(sprintf("  median CHM error:    %.2f m  (all data: %.2f m)\n",
              median(unk$chm_error_mean, na.rm = TRUE),
              median(dat$chm_error_mean, na.rm = TRUE)))
}

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

f_tuned <- bf(
  as.formula(paste(
    "chm_error_mean ~", main_effects,
    "+ s(slope_mean_z, wsci_z, k = 20)",
    "+ slope_mean_z:lc_grp + wsci_z:lc_grp",           # slope:wsci REMOVED
    "+ (1 + slope_mean_z | ecoregion) + (1 | site) + (1 + wsci_z | lc_grp)"
  )),
  sigma ~ 1 + slope_mean_z + (1 | site) + (1 | lc_grp)
)

cat("\nFormula:\n"); print(f_tuned)

priors <- c(
  prior(normal(0, 2),       class = "b"),
  prior(normal(0, 10),      class = "Intercept"),
  prior(gamma(2, 0.1),      class = "nu"),
  prior(student_t(3, 0, 2), class = "sd"),
  prior(normal(0, 1.5),     class = "sds"),
  prior(normal(0, 1),       class = "b",         dpar = "sigma"),
  prior(student_t(3, 0, 5), class = "Intercept", dpar = "sigma")
)

cat("\nSampler settings:\n")
cat("  chains 4, iterations 3000 (1500 warmup)  <- same as published\n")
cat("  adapt_delta 0.95        <- published was 0.98\n")
cat("  max_treedepth 12        <- published was 15 (never reached 13)\n")
cat("  step_size: not set      <- published forced 0.0005\n")
cat("  threads 16 per chain    <- published used 2\n")


# =====================================================================
sec("STEP 4  --  FITTING (this is the long part)")
# =====================================================================

t0 <- Sys.time()
cat("fitting started: ", format(t0), "\n\n")

fit <- brm(
  formula = f_tuned,
  data    = dat,
  family  = student(),
  prior   = priors,
  chains  = 4, iter = 3000, warmup = 1500,
  cores   = 4, threads = threading(16),
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  backend = "cmdstanr",
  file    = file.path(OUTDIR, "fit_chm_16_tuned"),
  seed    = 2026, refresh = 100
)

elapsed_h <- as.numeric(difftime(Sys.time(), t0, units = "hours"))


# =====================================================================
sec("STEP 5  --  RESULTS")
# =====================================================================

np <- as.data.table(brms::nuts_params(fit))
td <- np[Parameter == "treedepth__", Value]
dv <- np[Parameter == "divergent__", Value]

cat(sprintf("\nWALL CLOCK:            %.2f hours   <-- THE NUMBER WE NEEDED\n", elapsed_h))
cat(sprintf("published fit took:    58.10 hours\n"))
cat(sprintf("SPEEDUP:               %.1f x\n\n", 58.1 / max(elapsed_h, 1e-6)))

cat("QUALITY CHECKS (all must pass):\n")
d  <- sum(dv)
mr <- max(brms::rhat(fit), na.rm = TRUE)
me <- min(brms::neff_ratio(fit), na.rm = TRUE)
cat(sprintf("  divergences        %6d      %s  (need 0)\n",
            d, ifelse(d == 0, "PASS", "FAIL")))
cat(sprintf("  max Rhat           %6.4f      %s  (need < 1.01)\n",
            mr, ifelse(mr < 1.01, "PASS", "FAIL")))
cat(sprintf("  min ESS ratio      %6.3f      %s  (need > 0.10)\n",
            me, ifelse(me > 0.10, "PASS", "FAIL")))
cat(sprintf("  mean treedepth     %6.2f            (published was 10.73)\n", mean(td)))
cat(sprintf("  %% at ceiling        %6.1f%%      %s  (need < 1%%)\n",
            100 * mean(td >= 12), ifelse(mean(td >= 12) < 0.01, "PASS", "CHECK")))

if (d > 0) {
  cat("\n  ** Divergences present. Rerun with adapt_delta = 0.98 (line ~200).\n")
  cat("     It will be slower but still far faster than 58 hours. **\n")
}


# =====================================================================
sec("STEP 6  --  DID SIMPLIFYING CHANGE ANY CONCLUSIONS?")
# =====================================================================

pub_path <- file.path(CHKPT, "10b_chm_sensitivity.rds")
if (file.exists(pub_path)) {
  pub <- unwrap(readRDS(pub_path))$fit_chm_16

  gx <- function(f) {
    s <- as.data.table(posterior::summarise_draws(
      posterior::as_draws_df(f), mean, ~quantile(.x, c(0.025, 0.975))))
    setnames(s, c("variable", "mean", "lo", "hi"))
    s[grepl("^b_(?!sigma)", variable, perl = TRUE)]
  }
  a <- gx(pub); b <- gx(fit)
  cmp <- merge(a, b, by = "variable", suffixes = c("_pub", "_new"))
  cmp[, shift_SD := round((mean_new - mean_pub) /
                            ((hi_pub - lo_pub) / 3.92), 2)]
  cmp[, sign_flip := sign(mean_pub) != sign(mean_new) &
        (lo_pub > 0 | hi_pub < 0) & (lo_new > 0 | hi_new < 0)]
  cmp <- cmp[order(-abs(shift_SD))]

  cat("\nCoefficients present in both models, largest movement first.\n")
  cat("'shift_SD' = how far the estimate moved, in units of its own\n")
  cat("uncertainty. Under 0.5 is negligible. Over 1.0 needs a look.\n\n")
  print(cmp[, .(variable, mean_pub = round(mean_pub, 3),
                mean_new = round(mean_new, 3), shift_SD, sign_flip)][1:25])

  cat(sprintf("\n  coefficients that moved more than 1 SD: %d of %d\n",
              sum(abs(cmp$shift_SD) > 1, na.rm = TRUE), nrow(cmp)))
  cat(sprintf("  CREDIBLE SIGN FLIPS: %d   <-- must be 0\n",
              sum(cmp$sign_flip, na.rm = TRUE)))

  # the paper's headline coefficient
  rh <- cmp[variable == "b_rh_98_z"]
  if (nrow(rh)) cat(sprintf("\n  RH98 effect (the paper's key number): published %.3f -> new %.3f\n",
                            rh$mean_pub, rh$mean_new))
  fwrite(cmp, file.path(TBL, "confirm_coefficient_comparison.csv"))
  rm(pub); gc()
} else {
  cat("published fit not found, skipping comparison\n")
}


sec("SUMMARY FOR THE NEXT CONVERSATION")
cat(sprintf("
  hours for one full CHM fit ............ %.2f
  divergences ........................... %d
  max Rhat .............................. %.4f
  coefficients moved > 1 SD ............. %s
  credible sign flips ................... %s

  If sign flips is 0 and divergences is 0, this specification is safe
  and we can size the cross-validation plan from the hours figure.
",
elapsed_h, d, mr,
ifelse(exists("cmp"), sum(abs(cmp$shift_SD) > 1, na.rm = TRUE), "n/a"),
ifelse(exists("cmp"), sum(cmp$sign_flip, na.rm = TRUE), "n/a")))

hr(); cat("finished: ", format(Sys.time()), "\n")
cat("log: ", LOG, "\n"); hr()
