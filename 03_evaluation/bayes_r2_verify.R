# =====================================================================
# bayes_r2_verify.R
# PURPOSE: Verify the Table 6 in-sample Bayesian R2 rows for the CHM
#          16-site fit, against the published values:
#              Bayesian R2 (full model) = 0.131
#              R2 from fixed effects     = 0.017  (12.8% of explained)
#              R2 from random effects    = 0.114  (87.2% of explained)
#          Method per Table 4 caption: Bayesian R2 / variance partition
#          computed on TRAINING data (33% stratified subsample) following
#          Gelman et al. (2019); observation-level decomposition that
#          accounts for random slopes (Nakagawa & Schielzeth 2013).
#
# NO FABRICATION: this script only RECOMPUTES and COMPARES. It prints a
# PASS/FLAG gate. It does not write into the manuscript and does not alter
# any committed object. Run on gsapp22; paste the console block back.
#
# Mechanics:
#   conditional R2 (re_formula = NULL, all group effects) -> R2_full
#   marginal    R2 (re_formula = NA,   fixed only)        -> R2_fixed
#   R2_random  = R2_full - R2_fixed                        (per draw)
#   fixed_share  = R2_fixed  / R2_full                     (per draw)
#   random_share = R2_random / R2_full                     (per draw)
#   All summarized as POSTERIOR MEDIANS (matches "reported as posterior
#   medians" in the Table 4 caption).
# =====================================================================

suppressPackageStartupMessages({
  library(brms)
})

PROJECT_ROOT <- "/gpfs/data1/vclgp/lmaden/chpt1"
CKPT <- file.path(PROJECT_ROOT, "checkpoints", "10b_chm_sensitivity.rds")

cat("=== Table 6 Bayesian R2 verification (CHM 16-site PRIMARY) ===\n\n")

stopifnot(file.exists(CKPT))
obj <- readRDS(CKPT)
fit <- obj$data$fit_chm_16
stopifnot(inherits(fit, "brmsfit"))

# ---- G-train gate: confirm this is the primary 16-site fit ----
nd  <- fit$data
n   <- nrow(nd)
ns  <- length(unique(nd$site))
cat(sprintf("G-train: N=%d  sites=%d  family=%s  response=%s\n",
            n, ns, fit$family$family, all.vars(fit$formula$formula)[1]))
if (n != 124966) cat("  ** WARNING: N != 124,966 (expected primary 16-site training subsample)\n")
if (ns != 16)    cat("  ** WARNING: sites != 16\n")
cat("\n")

# ---- per-draw Bayesian R2: conditional (full) and marginal (fixed) ----
set.seed(1)  # bayes_R2 draws from the existing posterior; seed only affects any internal sampling
cat("Computing conditional Bayesian R2 (re_formula = NULL, all REs)...\n")
r2_full_draws  <- as.numeric(brms::bayes_R2(fit, re_formula = NULL, summary = FALSE))
cat("Computing marginal Bayesian R2 (re_formula = NA, fixed only)...\n")
r2_fixed_draws <- as.numeric(brms::bayes_R2(fit, re_formula = NA,   summary = FALSE))

# align draw counts defensively
L <- min(length(r2_full_draws), length(r2_fixed_draws))
r2_full_draws  <- r2_full_draws[seq_len(L)]
r2_fixed_draws <- r2_fixed_draws[seq_len(L)]

r2_random_draws  <- r2_full_draws - r2_fixed_draws
fixed_share_draws  <- r2_fixed_draws  / r2_full_draws
random_share_draws <- r2_random_draws / r2_full_draws

med <- function(x) median(x)
R2_full   <- med(r2_full_draws)
R2_fixed  <- med(r2_fixed_draws)
R2_random <- med(r2_random_draws)
share_fixed  <- med(fixed_share_draws)
share_random <- med(random_share_draws)

cat(sprintf("\n  draws used = %d\n\n", L))
cat(sprintf("  R2 full   (conditional) = %.4f   [published 0.131]\n", R2_full))
cat(sprintf("  R2 fixed  (marginal)    = %.4f   [published 0.017]\n", R2_fixed))
cat(sprintf("  R2 random (full-fixed)  = %.4f   [published 0.114]\n", R2_random))
cat(sprintf("  fixed share  = %.1f%%            [published 12.8%%]\n", 100*share_fixed))
cat(sprintf("  random share = %.1f%%            [published 87.2%%]\n", 100*share_random))

# also the ratio-of-medians shares (alternative summary), for transparency
cat(sprintf("\n  (ratio-of-medians shares: fixed %.1f%%, random %.1f%%)\n",
            100*R2_fixed/R2_full, 100*R2_random/R2_full))

# ---- GATE ----
ok_full   <- abs(R2_full   - 0.131) < 0.006
ok_fixed  <- abs(R2_fixed  - 0.017) < 0.006
ok_random <- abs(R2_random - 0.114) < 0.006
ok_share  <- abs(100*share_fixed - 12.8) < 1.5 && abs(100*share_random - 87.2) < 1.5

cat("\n=== GATE ===\n")
cat(sprintf("  R2 full   reproduces 0.131 : %s\n", if (ok_full)   "PASS" else "FLAG"))
cat(sprintf("  R2 fixed  reproduces 0.017 : %s\n", if (ok_fixed)  "PASS" else "FLAG"))
cat(sprintf("  R2 random reproduces 0.114 : %s\n", if (ok_random) "PASS" else "FLAG"))
cat(sprintf("  shares reproduce 12.8/87.2 : %s\n", if (ok_share)  "PASS" else "FLAG"))
if (ok_full && ok_fixed && ok_random && ok_share) {
  cat("\nRESULT: PASS - Table 6 Bayesian R2 rows are confirmed 16-site. No edit.\n")
} else {
  cat("\nRESULT: FLAG - one or more rows do NOT reproduce from fit_chm_16.\n")
  cat("        Review these before changing any published value.\n")
}
cat("\nDone.\n")
