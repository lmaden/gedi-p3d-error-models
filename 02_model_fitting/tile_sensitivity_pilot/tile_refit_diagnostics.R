#!/usr/bin/env Rscript
# =====================================================================
# tile_refit_diagnostics.R  --  READ-ONLY. Attribute m2_refit.R's STEP 3
# Rhat / ESS failures to specific parameters.
#
#   Rscript tile_refit_diagnostics.R --product dtm
#   Rscript tile_refit_diagnostics.R --product chm
#   Rscript tile_refit_diagnostics.R --product dtm --tag dtm_tile_init0   # a rerun
#
# WHY THIS EXISTS
# m2_refit.R STEP 3 computes
#     max(brms::rhat(fit))        and     min(brms::neff_ratio(fit))
# over EVERY parameter in the fit. Adding (1 | tile) grew that set by
# ~19,331 r_tile[...] intercepts, and STEP 1 reports a median of 3
# footprints per tile. One badly-mixed tile intercept therefore fails
# the whole gate. The 1.01 / 0.10 thresholds were set in checkpoint
# section 6.5 for the LOO fold fits, which carried no tile level.
#
# This script splits the same two numbers by parameter class, names the
# offenders, and tests whether the failures track tile size. It writes
# nothing except its own stdout and it DECIDES NOTHING -- the rerun call
# is ours to make from what it reports.
#
# USAGE (exports required before launch; HOME is node-local)
#   export PROJECT_ROOT=/gpfs/data1/vclgp/lmaden/chpt1
#   export TMPDIR=$PROJECT_ROOT/tmp
#   export R_LIBS=/gpfs/data1/vclgp/lmaden/Rlib
#   export CMDSTAN=/gpfs/data1/vclgp/lmaden/cmdstan/cmdstan-2.37.0
#   nohup nice -n 10 Rscript $PROJECT_ROOT/scripts/reviewed/tile_refit_diagnostics.R \
#     --product dtm > $PROJECT_ROOT/manuscript_tables/tile_refit_diagnostics_dtm.log 2>&1 &
#
# Expect 5-15 minutes: most of it is reading a multi-GB fit off gpfs and
# recomputing Rhat over ~19,400 parameters. One core, no fitting.
# =====================================================================

suppressPackageStartupMessages({
  library(data.table); library(brms); library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) args[i + 1] else default
}
PRODUCT <- tolower(getarg("--product", ""))
if (!PRODUCT %in% c("chm", "dtm"))
  stop("Usage: Rscript tile_refit_diagnostics.R --product chm|dtm [--tag <fit tag>]")
TAG <- getarg("--tag", paste0(PRODUCT, "_tile"))

THRESH_RHAT <- 1.01     # checkpoint 6.5 pass criterion
THRESH_ESS  <- 0.10     # checkpoint 6.5 pass criterion

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
FITFILE <- file.path(PROJECT_ROOT, "models", "m2_refit", paste0("fit_", TAG, ".rds"))

hr  <- function() cat(strrep("=", 92), "\n")
sec <- function(x) { cat("\n"); hr(); cat("## ", x, "\n"); hr() }

cat("tile_refit_diagnostics.R -- product ", toupper(PRODUCT), " | tag ", TAG, "\n", sep = "")
cat("started:  ", format(Sys.time()), "\n", sep = "")
cat("fit file: ", FITFILE, "\n", sep = "")
if (!file.exists(FITFILE)) stop("fit not found: ", FITFILE)
cat("on disk:  ", round(file.size(FITFILE) / 2^30, 2), " GB\n", sep = "")

t0 <- Sys.time()
fit <- readRDS(FITFILE)
cat(sprintf("loaded in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---------------------------------------------------------------------
sec("STEP 1 -- REPRODUCE THE TWO NUMBERS m2_refit.R PRINTED")
# Identical estimators to m2_refit.R STEP 3, deliberately: if these do
# not match the log, something other than the diagnostic is wrong.
t0 <- Sys.time()
rh <- brms::rhat(fit)
nr <- brms::neff_ratio(fit)
cat(sprintf("computed in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

if (is.null(names(rh)))
  stop("brms::rhat() returned an unnamed vector -- cannot attribute by parameter")

nm <- names(rh)
cat(sprintf("  parameters entering the gate: %s\n", format(length(nm), big.mark = ",")))
cat(sprintf("  max Rhat        %.4f     <- m2_refit.R STEP 3 printed this\n",
            max(rh, na.rm = TRUE)))
cat(sprintf("  min ESS ratio   %.3f     <- m2_refit.R STEP 3 printed this\n",
            min(nr, na.rm = TRUE)))
cat(sprintf("  worst Rhat is on: %s\n", nm[which.max(rh)]))
cat(sprintf("  worst ESS  is on: %s\n", names(nr)[which.min(nr)]))

# ---------------------------------------------------------------------
sec("STEP 2 -- WHICH PARAMETER CLASS FAILS")
# General to specific; later assignments override earlier ones.
klass <- rep("other", length(nm))
klass[grepl("^b_",          nm)] <- "b_ fixed effects"
klass[grepl("^b_sigma_",    nm)] <- "b_sigma_ sigma model"
klass[grepl("^bs_",         nm)] <- "bs_ spline fixed"
klass[grepl("^sds_",        nm)] <- "sds_ spline sd"
klass[grepl("^s_",          nm)] <- "s_ spline coefs"
klass[grepl("^sd_",         nm)] <- "sd_ group sds"
klass[grepl("^cor_|^L_",    nm)] <- "cor_ correlations"
klass[grepl("^r_ecoregion", nm)] <- "r_ecoregion"
klass[grepl("^r_lc_grp",    nm)] <- "r_lc_grp"
klass[grepl("^r_site",      nm)] <- "r_site"
klass[grepl("^r_tile",      nm)] <- "r_tile"
klass[nm == "nu"]                <- "nu"
klass[nm == "sigma"]             <- "sigma"
klass[nm %in% c("lp__", "lprior")] <- "lp__/lprior"

DT <- data.table(param = nm, klass = klass,
                 rhat = as.numeric(rh), ess = as.numeric(nr[nm]))

tab <- DT[, .(n            = .N,
              max_rhat     = round(max(rhat, na.rm = TRUE), 4),
              n_rhat_fail  = sum(rhat >= THRESH_RHAT, na.rm = TRUE),
              min_ess      = round(min(ess, na.rm = TRUE), 3),
              n_ess_fail   = sum(ess <= THRESH_ESS, na.rm = TRUE)),
          by = klass][order(-max_rhat)]
print(tab, nrows = 100)

# ---------------------------------------------------------------------
sec("STEP 3 -- THE QUANTITIES THE MANUSCRIPT ACTUALLY REPORTS")
# CORE = Table 5 (fixed effects), the sigma model, the variance
# partition (group sds), correlations and nu. Excludes every r_*
# (individual random-effect levels), the spline basis coefs, lp__.
CORE <- DT[!grepl("^r_|^s_", param) & !param %in% c("lp__", "lprior")]
core_rhat <- max(CORE$rhat, na.rm = TRUE)
core_ess  <- min(CORE$ess,  na.rm = TRUE)

cat(sprintf("  CORE set: %d parameters\n", nrow(CORE)))
cat(sprintf("  max Rhat        %.4f   %s  (need < %.2f)\n",
            core_rhat, ifelse(core_rhat < THRESH_RHAT, "PASS", "FAIL"), THRESH_RHAT))
cat(sprintf("  min ESS ratio   %.3f   %s  (need > %.2f)\n",
            core_ess, ifelse(core_ess > THRESH_ESS, "PASS", "FAIL"), THRESH_ESS))

cat("\n  worst 15 CORE parameters by Rhat:\n")
print(head(CORE[order(-rhat), .(param, rhat = round(rhat, 4), ess = round(ess, 3))], 15))
cat("\n  worst 15 CORE parameters by ESS ratio:\n")
print(head(CORE[order(ess), .(param, rhat = round(rhat, 4), ess = round(ess, 3))], 15))

SPL <- DT[grepl("^s_", param)]
if (nrow(SPL))
  cat(sprintf("\n  spline basis coefs (s_*, not reported directly): n %d, max Rhat %.4f, min ESS %.3f\n",
              nrow(SPL), max(SPL$rhat, na.rm = TRUE), min(SPL$ess, na.rm = TRUE)))

# ---------------------------------------------------------------------
sec("STEP 4 -- DO THE FAILURES TRACK TILE SIZE?")
d <- fit$data
if (!"tile" %in% names(d)) {
  cat("  fit$data carries no 'tile' column -- skipping this step\n")
} else {
  tn <- table(as.character(d$tile))
  RT <- DT[grepl("^r_tile\\[", param)]
  if (!nrow(RT)) {
    cat("  no r_tile parameters found -- skipping\n")
  } else {
    RT[, id    := sub("^r_tile\\[([^,]+),.*$", "\\1", param)]
    RT[, nfoot := as.integer(tn[id])]
    cat(sprintf("  tiles: %s | holding 1 footprint: %s | <= 3 footprints: %s | median %.0f\n",
                format(length(tn), big.mark = ","),
                format(sum(tn == 1), big.mark = ","),
                format(sum(tn <= 3), big.mark = ","),
                median(as.numeric(tn))))
    RT[, fail := (rhat >= THRESH_RHAT) | (ess <= THRESH_ESS)]
    cat(sprintf("  r_tile failing either gate: %s of %s (%.2f%%)\n",
                format(sum(RT$fail, na.rm = TRUE), big.mark = ","),
                format(nrow(RT), big.mark = ","),
                100 * mean(RT$fail, na.rm = TRUE)))
    cat("\n  footprints per tile, failing vs passing:\n")
    print(RT[, .(n_params = .N,
                 min_n    = min(nfoot, na.rm = TRUE),
                 median_n = as.numeric(median(nfoot, na.rm = TRUE)),
                 mean_n   = round(mean(nfoot, na.rm = TRUE), 1),
                 max_n    = max(nfoot, na.rm = TRUE)), by = fail])
    cat("\n  worst 10 r_tile by Rhat:\n")
    print(head(RT[order(-rhat),
                  .(param, nfoot, rhat = round(rhat, 4), ess = round(ess, 3))], 10))
  }
}

# ---------------------------------------------------------------------
sec("STEP 5 -- IS ONE CHAIN SITTING SOMEWHERE ELSE?")
# Chain 4 ran one treedepth level deeper than chains 1-3 (12.00 vs
# 11.00) at 0% ceiling, so it is not pinned -- but a chain exploring a
# different region would show a consistent offset across many CORE
# parameters. Checkpoint 6.6 is the precedent for this test AND for the
# rule that it is diagnostic only: selecting chains after inspecting
# results is indefensible. All four leave-one-out values are printed so
# nothing is cherry-picked.
dr <- posterior::as_draws_array(fit, variable = CORE$param)
nch <- dim(dr)[2]; ndraws <- prod(dim(dr)[1:2])
cat(sprintf("  CORE draws array: %d iterations x %d chains x %d parameters (%d draws)\n",
            dim(dr)[1], nch, dim(dr)[3], ndraws))

pc    <- apply(dr, c(2, 3), mean)                 # [chain, parameter]
sdall <- apply(dr, 3, sd)
z     <- sweep(sweep(pc, 2, colMeans(pc), "-"), 2, sdall, "/")
cat("\n  per-chain mean minus overall mean, in posterior SD units:\n")
for (i in seq_len(nch))
  cat(sprintf("    chain %d:  mean |z| %.3f   max |z| %.2f   (on %s)\n",
              i, mean(abs(z[i, ])), max(abs(z[i, ])),
              colnames(pc)[which.max(abs(z[i, ]))]))

cat("\n  max CORE Rhat with each chain dropped in turn (diagnostic only --\n")
cat("  NOT a licence to report a subset of chains, see checkpoint 6.6):\n")
for (cc in seq_len(nch)) {
  sub <- dr[, setdiff(seq_len(nch), cc), , drop = FALSE]
  r3  <- apply(sub, 3, function(m) posterior::rhat(as.matrix(m)))
  cat(sprintf("    drop chain %d  ->  max CORE Rhat %.4f\n", cc, max(r3, na.rm = TRUE)))
}

worst <- head(CORE[order(-rhat)]$param, 8)
cat("\n  per-chain posterior means, 8 worst CORE parameters:\n")
pcw <- t(round(pc[, worst, drop = FALSE], 4))
colnames(pcw) <- paste0("chain", seq_len(nch))
print(pcw)

# ---------------------------------------------------------------------
sec("STEP 6 -- BULK AND TAIL ESS FOR THE CORE SET")
sm <- as.data.table(posterior::summarise_draws(dr, "rhat", "ess_bulk", "ess_tail"))
sm[, `:=`(bulk_ratio = round(ess_bulk / ndraws, 3),
          tail_ratio = round(ess_tail / ndraws, 3))]
cat("  worst 15 by bulk ESS:\n")
print(head(sm[order(ess_bulk), .(variable, rhat = round(rhat, 4),
                                 ess_bulk = round(ess_bulk), bulk_ratio,
                                 ess_tail = round(ess_tail), tail_ratio)], 15))

# ---------------------------------------------------------------------
sec("SUMMARY")
cat(sprintf("  whole-fit gate (as m2_refit.R computes it): Rhat %.4f %s | ESS %.3f %s\n",
            max(rh, na.rm = TRUE), ifelse(max(rh, na.rm = TRUE) < THRESH_RHAT, "PASS", "FAIL"),
            min(nr, na.rm = TRUE), ifelse(min(nr, na.rm = TRUE) > THRESH_ESS, "PASS", "FAIL")))
cat(sprintf("  CORE only (fixed effects, sigma model, group sds, cor, nu): Rhat %.4f %s | ESS %.3f %s\n",
            core_rhat, ifelse(core_rhat < THRESH_RHAT, "PASS", "FAIL"),
            core_ess,  ifelse(core_ess  > THRESH_ESS,  "PASS", "FAIL")))
cat("\n  If CORE passes and the failures sit in r_tile at small tiles, the reported\n")
cat("  quantities are sound and the whole-fit gate is measuring a nuisance level that\n")
cat("  did not exist when the 1.01 / 0.10 thresholds were written.\n")
cat("  If CORE fails, the fit is not usable as primary and a rerun is required.\n")
cat("  This script does not make that call.\n")
cat("\ndone: ", format(Sys.time()), "\n", sep = "")
