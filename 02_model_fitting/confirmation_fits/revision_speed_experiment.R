#!/usr/bin/env Rscript
# =====================================================================
# revision_speed_experiment.R
#
# Sampling-speed experiment. The 16-site CHM fit took 58 h and the 18-site DTM
# fit 161 h on the slowest chain, which rules out cross-validation at that cost.
# This measures where the time goes and whether a better-specified model samples
# faster.
#
# PART 1 (5 min, read-only): sampler diagnostics on the EXISTING fits.
#   Tells us if the cost is bad posterior geometry (treedepth saturation)
#   or just a tiny step size.
#
# PART 2 (overnight): fit 5 specification variants on an identical small
#   balanced subsample and time them. Relative timings are what matter.
#
# RUN:  cd $PROJECT_ROOT/scripts/reviewed
#       nohup Rscript revision_speed_experiment.R > speed_exp.out 2>&1 &
#       tail -f speed_exp.out
#
# Set SKIP_PART2=TRUE to run only the fast diagnostics first.
# =====================================================================

suppressPackageStartupMessages({
  library(brms); library(data.table); library(posterior)
})

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
CHKPT  <- file.path(PROJECT_ROOT, "checkpoints")
TBL    <- file.path(PROJECT_ROOT, "manuscript_tables")
EXPDIR <- file.path(PROJECT_ROOT, "models", "speed_exp")
dir.create(EXPDIR, recursive = TRUE, showWarnings = FALSE)

SKIP_PART2 <- as.logical(Sys.getenv("SKIP_PART2", "FALSE"))

LOG <- file.path(TBL, sprintf("speed_experiment_%s.log", format(Sys.Date(), "%Y%m%d")))
con <- file(LOG, open = "wt"); sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

hr <- function() cat(strrep("=", 96), "\n")
sec <- function(x) { cat("\n"); hr(); cat("## ", x, "\n"); hr() }
unwrap <- function(p) if (is.list(p) && !is.null(p$data)) p$data else p

cat("revision_speed_experiment.R\n"); cat("run at: ", format(Sys.time()), "\n")


# =====================================================================
sec("PART 1  --  WHY ARE THE EXISTING FITS SLOW?")
# =====================================================================

diagnose_fit <- function(fit, label) {
  cat("\n---------- ", label, " ----------\n")
  np <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)
  if (is.null(np)) { cat("  nuts_params unavailable\n"); return(invisible(NULL)) }
  np <- as.data.table(np)

  td <- np[Parameter == "treedepth__", Value]
  dv <- np[Parameter == "divergent__", Value]
  ss <- np[Parameter == "stepsize__",  Value]

  cat(sprintf("  post-warmup iterations logged : %d\n", length(td)))
  cat(sprintf("  DIVERGENT transitions         : %d (%.2f%%)\n",
              sum(dv), 100 * mean(dv)))
  cat(sprintf("  treedepth  mean / median / max: %.2f / %.0f / %.0f\n",
              mean(td), median(td), max(td)))
  cat(sprintf("  %% iterations AT max_treedepth=15: %.1f%%   <-- KEY NUMBER\n",
              100 * mean(td >= 15)))
  cat(sprintf("  %% iterations at treedepth >= 12 : %.1f%%\n",
              100 * mean(td >= 12)))
  cat(sprintf("  adapted step size per chain    : %s\n",
              paste(signif(unique(ss), 3), collapse = ", ")))
  cat(sprintf("  implied leapfrog steps/iter    : ~%.0f (2^mean treedepth)\n",
              2^mean(td)))
  cat("\n  INTERPRETATION:\n")
  if (mean(td >= 15) > 0.10) {
    cat("    >10% of iterations hit the treedepth ceiling. The cost is BAD\n")
    cat("    POSTERIOR GEOMETRY (weakly identified parameters), not raw data\n")
    cat("    size. Simplifying the model will help enormously.\n")
  } else if (mean(td) > 9) {
    cat("    High but not saturating treedepth. Partly geometry, partly the\n")
    cat("    hardcoded step_size=0.0005. Both are fixable.\n")
  } else {
    cat("    Treedepth is moderate. The cost is per-gradient (data size and\n")
    cat("    parameter count), so THREADING is the main lever.\n")
  }
  invisible(np)
}

for (f in list(
  list(p = file.path(CHKPT, "10b_chm_sensitivity.rds"), obj = "fit_chm_16", lab = "CHM 16-site (58 h)"),
  list(p = file.path(CHKPT, "10_models_stage2.rds"),    obj = "fit_dtm_s2", lab = "DTM 18-site (161 h)")
)) {
  if (!file.exists(f$p)) { cat("missing: ", f$p, "\n"); next }
  pay <- unwrap(readRDS(f$p))
  fit <- pay[[f$obj]]
  if (!is.null(fit)) diagnose_fit(fit, f$lab)
  rm(pay, fit); gc()
}


# =====================================================================
sec("PART 1b  --  WHERE ARE THE WASTED PARAMETERS?")
# =====================================================================

ing <- unwrap(readRDS(file.path(CHKPT, "01_data_ingest.rds")))
chm <- as.data.table(ing$chm_df)
rm(ing); gc()

SITES_FLAGGED <- c("1", "2", "3")
chm16 <- chm[!as.character(site) %in% SITES_FLAGGED & is.finite(chm_error_mean) & !is.na(lc_l1_code)]

cat("\n-- land-cover class frequencies in the 16-site CHM pool\n")
lc_tab <- chm16[, .N, by = lc_l1_code][order(-N)]
lc_tab[, pct := round(100 * N / sum(N), 3)]
print(lc_tab)
cat("\n  Each of these levels currently gets:\n")
cat("    - a slope_mean_z interaction coefficient\n")
cat("    - a wsci_z interaction coefficient\n")
cat("    - a random intercept and a random wsci_z slope (with LKJ correlation)\n")
cat("    - a sigma random intercept\n")
cat(sprintf("  => roughly %d parameters tied to land cover, for %d classes,\n",
            2 * nrow(lc_tab) + 3 * nrow(lc_tab), nrow(lc_tab)))
cat(sprintf("     of which %d classes hold under 0.5%% of the data.\n",
            sum(lc_tab$pct < 0.5)))

cat("\n-- site imbalance in the 16-site CHM pool\n")
st <- chm16[, .N, by = site][order(-N)]
st[, pct := round(100 * N / sum(N), 1)]
print(st)
cat(sprintf("\n  LARGEST SITE = %.1f%% of all training rows.  <-- KEY NUMBER\n",
            max(st$pct)))
cat("  If this is over 50%, the 'multi-site' fixed effects are largely one site.\n")

cat("\n-- redundancy check: s(slope,wsci) spline vs the slope:wsci linear term\n")
cat("   correlation between slope_mean_z and wsci_z: ",
    round(cor(chm16$slope_mean_z, chm16$wsci_z, use = "complete.obs"), 3), "\n")
cat("   The 2D thin-plate spline already spans the linear interaction surface.\n")
cat("   Keeping both creates a ridge in the posterior.\n")

if (SKIP_PART2) { cat("\nSKIP_PART2=TRUE, stopping here.\n"); hr(); quit(save = "no") }


# =====================================================================
sec("PART 2  --  SPECIFICATION SPEED EXPERIMENT")
# =====================================================================

# Balanced subsample: every site represented, no site allowed to dominate.
set.seed(2026)
CAP <- 400L
samp <- chm16[, .SD[sample(.N, min(.N, CAP))], by = site]
samp[, `:=`(site = factor(site), ecoregion = factor(ecoregion),
            lc_l1_code = factor(lc_l1_code))]

# Collapsed land cover: the four forest types the paper is actually about
FOREST <- c("BDF", "DNF", "EBF", "ENF")
samp[, lc_collapsed := factor(fifelse(as.character(lc_l1_code) %in% FOREST,
                                      as.character(lc_l1_code), "OTHER"))]

cat(sprintf("\nExperiment sample: %d rows, %d sites, %d LC classes (%d collapsed)\n",
            nrow(samp), uniqueN(samp$site),
            uniqueN(samp$lc_l1_code), uniqueN(samp$lc_collapsed)))
print(samp[, .N, by = site][order(site)])

priors_base <- c(
  prior(normal(0, 2),        class = "b"),
  prior(normal(0, 10),       class = "Intercept"),
  prior(gamma(2, 0.1),       class = "nu"),
  prior(student_t(3, 0, 2),  class = "sd"),
  prior(normal(0, 1.5),      class = "sds")
)

base_terms <- paste(
  "slope_mean_z + slope_sd_z + wsci_z + rh_98_z + cover_z +",
  "aspect_sin_z + aspect_cos_z + meta_offnad_z + meta_sunel_z +",
  "meta_az_conc_z + meta_stereo_z + meta_leafon_z + meta_fwdrev_z +",
  "meta_absgeo_z + meta_relgeo_z + view_az_sin_z + view_az_cos_z"
)
re_terms <- "(1 + slope_mean_z | ecoregion) + (1 | site) + (1 + wsci_z | lc_l1_code)"
re_terms_c <- "(1 + slope_mean_z | ecoregion) + (1 | site) + (1 + wsci_z | lc_collapsed)"

variants <- list(
  V0_current = list(
    f = bf(as.formula(paste("chm_error_mean ~", base_terms,
        "+ s(slope_mean_z, wsci_z, k = 20) + slope_mean_z:lc_l1_code +",
        "wsci_z:lc_l1_code + slope_mean_z:wsci_z +", re_terms)),
        sigma ~ 1 + slope_mean_z + (1 | site) + (1 | lc_l1_code)),
    ctrl = list(adapt_delta = 0.98, max_treedepth = 15, step_size = 0.0005),
    thr = 2, note = "exactly what was published"),

  V1_no_redundant = list(
    f = bf(as.formula(paste("chm_error_mean ~", base_terms,
        "+ s(slope_mean_z, wsci_z, k = 20) + slope_mean_z:lc_l1_code +",
        "wsci_z:lc_l1_code +", re_terms)),
        sigma ~ 1 + slope_mean_z + (1 | site) + (1 | lc_l1_code)),
    ctrl = list(adapt_delta = 0.98, max_treedepth = 15, step_size = 0.0005),
    thr = 2, note = "drop the linear slope:wsci (spline already spans it)"),

  V2_collapsed_lc = list(
    f = bf(as.formula(paste("chm_error_mean ~", base_terms,
        "+ s(slope_mean_z, wsci_z, k = 20) + slope_mean_z:lc_collapsed +",
        "wsci_z:lc_collapsed +", re_terms_c)),
        sigma ~ 1 + slope_mean_z + (1 | site) + (1 | lc_collapsed)),
    ctrl = list(adapt_delta = 0.98, max_treedepth = 15, step_size = 0.0005),
    thr = 2, note = "16 LC classes -> 4 forest types + OTHER"),

  V3_sampler = list(
    f = bf(as.formula(paste("chm_error_mean ~", base_terms,
        "+ s(slope_mean_z, wsci_z, k = 20) + slope_mean_z:lc_collapsed +",
        "wsci_z:lc_collapsed +", re_terms_c)),
        sigma ~ 1 + slope_mean_z + (1 | site) + (1 | lc_collapsed)),
    ctrl = list(adapt_delta = 0.90, max_treedepth = 12),
    thr = 2, note = "let Stan adapt its own step size; adapt_delta 0.90"),

  V4_threaded = list(
    f = bf(as.formula(paste("chm_error_mean ~", base_terms,
        "+ s(slope_mean_z, wsci_z, k = 20) + slope_mean_z:lc_collapsed +",
        "wsci_z:lc_collapsed +", re_terms_c)),
        sigma ~ 1 + slope_mean_z + (1 | site) + (1 | lc_collapsed)),
    ctrl = list(adapt_delta = 0.90, max_treedepth = 12),
    thr = 16, note = "same as V3 but 16 threads/chain (node has 192 cores)")
)

results <- data.table()
fits <- list()

for (nm in names(variants)) {
  v <- variants[[nm]]
  cat("\n"); hr()
  cat("FITTING ", nm, "  --  ", v$note, "\n"); hr()
  t0 <- Sys.time()
  fit <- tryCatch(
    brm(formula = v$f, data = samp, family = student(), prior = priors_base,
        chains = 2, iter = 1000, warmup = 500,
        cores = 2, threads = threading(v$thr),
        control = v$ctrl, backend = "cmdstanr",
        file = file.path(EXPDIR, nm), seed = 2026, refresh = 200),
    error = function(e) { cat("  FAILED: ", conditionMessage(e), "\n"); NULL })
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (is.null(fit)) next
  fits[[nm]] <- fit

  np <- as.data.table(brms::nuts_params(fit))
  td <- np[Parameter == "treedepth__", Value]
  dv <- np[Parameter == "divergent__", Value]

  results <- rbind(results, data.table(
    variant     = nm,
    note        = v$note,
    minutes     = round(el, 1),
    n_params    = length(brms::variables(fit)),
    divergences = sum(dv),
    pct_maxtd   = round(100 * mean(td >= v$ctrl$max_treedepth), 1),
    mean_td     = round(mean(td), 2),
    max_rhat    = round(max(brms::rhat(fit), na.rm = TRUE), 4),
    min_ess     = round(min(brms::neff_ratio(fit), na.rm = TRUE), 3)
  ))
  cat(sprintf("\n  >>> %s finished in %.1f minutes\n", nm, el))
  print(results)
}

sec("RESULTS")
print(results)
fwrite(results, file.path(TBL, "speed_experiment_results.csv"))

if (nrow(results) > 1) {
  spd <- results[1, minutes] / results[, minutes]
  cat("\nSPEEDUP vs the published specification:\n")
  for (i in seq_len(nrow(results)))
    cat(sprintf("  %-16s %6.1f min   %5.1fx faster\n",
                results$variant[i], results$minutes[i], spd[i]))
  cat(sprintf("\nEXTRAPOLATION to the full 16-site fit (124,966 rows):\n"))
  cat(sprintf("  published fit took 58 h. Best variant implies roughly %.1f h.\n",
              58 / max(spd)))
  cat("  (crude: geometry gains scale better than linearly, so this is pessimistic)\n")
}

sec("DOES SIMPLIFYING COST US ANYTHING? (LOO comparison)")
if (length(fits) >= 2) {
  loos <- lapply(fits, function(f) tryCatch(loo(f), error = function(e) NULL))
  loos <- Filter(Negate(is.null), loos)
  if (length(loos) >= 2) print(loo_compare(loos))
  cat("\nIf the simplified variants are within ~2 SE of V0, the dropped\n")
  cat("parameters were not earning their keep and we simplify with a clear\n")
  cat("conscience. Report this comparison in the paper.\n")
}

hr(); cat("Log: ", LOG, "\n"); hr()
