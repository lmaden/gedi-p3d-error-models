#!/usr/bin/env Rscript
# =====================================================================
# tile_pilot_chm.R  --  price the tile-level random intercept (M2) on a
# seven-site subset before committing the full refits.
#
# KEY DECISION THIS IMPLEMENTS (tile_key_probe.R, 16 Aug):
#   - No true tile/strip ID exists in the data; shot_number is
#     footprint-unique and nothing else is ID-like.
#   - The six geometry columns form exactly ONE partition per site
#     (joint == max_single at 19/19 CHM and 18/18 DTM sites), and the
#     ten modeled acquisition predictors define the SAME partition.
#   - So the tile key = unique geometry-value combination within site:
#     the resolution at which acquisition metadata is assigned.
#     CHM 19-site: 19,450 cells, median 3 footprints, max 1,653.
#
# WHAT THIS PILOT DOES
#   Seven CHM sites in full (4 7 9 11 12 13 17; ~24.6k rows, ~1.9k
#   tiles), two fits at 2 chains x 1000 iter (500 warmup):
#     tile:    published CHM spec + (1 | tile)
#     control: published CHM spec unchanged, same subset
#   Prints diagnostics for both, the tile sd posterior, and the
#   widening of the ten acquisition fixed effects (posterior SD ratio
#   tile/control). The pilot PASSES if the tile fit is clean (zero
#   divergences, no treedepth pile-up, sane Rhat) at the published
#   sampler settings and runtime is acceptable.
#
# Also verifies two things recorded in the plan:
#   - GATE: the modeled-ten partition equals the geometry-six
#     partition on the full 16-site data (hard stop if not)
#   - the vestigial _avg/_ct/_ratio columns are all-NA, not
#     site-constant (corrects a premise-check misreading)
#
# USAGE  (expect tens of minutes per fit; over ~2 min => run in tmux)
#   export PROJECT_ROOT=/gpfs/data1/vclgp/lmaden/chpt1
#   export TMPDIR=$PROJECT_ROOT/tmp
#   export R_LIBS=/gpfs/data1/vclgp/lmaden/Rlib
#   export CMDSTAN=/gpfs/data1/vclgp/lmaden/cmdstan/cmdstan-2.37.0
#   Rscript scripts/reviewed/tile_pilot_chm.R
#   (2 chains x 8 threads = 16 cores)
# =====================================================================

suppressPackageStartupMessages({ library(data.table); library(brms) })

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
OUTMOD <- file.path(PROJECT_ROOT, "models", "m2_pilot")
dir.create(OUTMOD, showWarnings = FALSE, recursive = TRUE)

SITES_PILOT   <- c("4", "7", "9", "11", "12", "13", "17")
SITES_FLAGGED <- c("1", "2", "3")          # CHM exclusions, as in loo_fold_chm.R
ITER <- 1000; WARMUP <- 500; CHAINS <- 2; THREADS <- 8
ADAPT_DELTA <- 0.95; MAX_TREEDEPTH <- 12   # the CHM confirmation-fit settings

hr  <- function() cat(strrep("=", 92), "\n")
sec <- function(x) { cat("\n"); hr(); cat("## ", x, "\n"); hr() }

sec("STEP 0 -- WRITE-LOCATION GUARD")
ALLOWED <- "/gpfs/data1/vclgp/lmaden"
cat("tempdir(): ", tempdir(), "\n")
if (!startsWith(normalizePath(tempdir(), mustWork = FALSE), ALLOWED))
  stop("tempdir() outside ", ALLOWED, " -- export TMPDIR=$PROJECT_ROOT/tmp first")
cat("OK\n")

sec("STEP 1 -- DATA AND THE TILE KEY")
unwrap <- function(p) if (is.list(p) && !is.null(p$data)) p$data else p
mp  <- unwrap(readRDS(file.path(PROJECT_ROOT, "checkpoints", "08_model_prep.rds")))
dat <- as.data.table(mp$mod_chm_s2)
rm(mp); gc()

cat("\nNA share of the 'site-constant' columns (expect 1.000 = empty, not constant):\n")
for (cn in c("meta_abs_geoacc_avg", "meta_az_concentration", "meta_tot_ct", "rh_98"))
  if (cn %in% names(dat))
    cat(sprintf("  %-24s %.3f\n", cn, mean(is.na(dat[[cn]]))))

dat <- dat[!as.character(site) %in% SITES_FLAGGED]
KEEP <- c("ENF", "IWL", "BDF", "RCP", "GRS", "UNK", "DNF", "WTR", "IMP", "EBF")
dat[, lc_grp := factor(fifelse(as.character(lc_l1_code) %in% KEEP,
                               as.character(lc_l1_code), "OTHER"))]
dat[, lc_grp := relevel(lc_grp, ref = "ENF")]
dat[, `:=`(site = factor(site), ecoregion = factor(ecoregion))]

MODEL_VARS <- c("chm_error_mean", "slope_mean_z", "slope_sd_z", "wsci_z",
                "rh_98_z", "cover_z", "aspect_sin_z", "aspect_cos_z",
                "meta_offnad_z", "meta_sunel_z", "meta_az_conc_z",
                "meta_stereo_z", "meta_leafon_z", "meta_fwdrev_z",
                "meta_absgeo_z", "meta_relgeo_z", "view_az_sin_z",
                "view_az_cos_z", "site", "ecoregion", "lc_grp")
dat <- dat[complete.cases(dat[, ..MODEL_VARS])]
cat("\n16-site CHM analysis rows after NA drop:", format(nrow(dat), big.mark = ","),
    "(published fit used 124,966)\n")

GEO   <- c("meta_offnad_z", "meta_sunel_z", "meta_az_conc_z", "meta_absgeo_z",
           "view_az_sin_z", "view_az_cos_z")
ACQ10 <- c(GEO, "meta_relgeo_z", "meta_stereo_z", "meta_fwdrev_z", "meta_leafon_z")

dat[, tile   := .GRP, by = c("site", GEO)]
dat[, tile10 := .GRP, by = c("site", ACQ10)]
same <- all(dat[, uniqueN(tile10), by = tile]$V1 == 1) &&
        all(dat[, uniqueN(tile),   by = tile10]$V1 == 1)
cat("\nGATE: geometry-six partition == modeled-ten partition on 16-site data:",
    ifelse(same, "PASS", "FAIL"), "\n")
if (!same) stop("Partitions differ -- re-examine the key before piloting.")
dat[, tile10 := NULL]
dat[, tile := factor(tile)]
cat("tile levels, 16-site CHM:", format(nlevels(dat$tile), big.mark = ","), "\n")

pil <- dat[as.character(site) %in% SITES_PILOT]
pil[, `:=`(site = droplevels(site), ecoregion = droplevels(ecoregion),
           lc_grp = droplevels(lc_grp), tile = droplevels(tile))]
ntile <- pil[, .N, by = tile]$N
cat(sprintf("\nPILOT subset: %s rows | %d sites | %d ecoregions | %d land-cover | %s tiles\n",
            format(nrow(pil), big.mark = ","), nlevels(pil$site),
            nlevels(pil$ecoregion), nlevels(pil$lc_grp),
            format(nlevels(pil$tile), big.mark = ",")))
cat(sprintf("footprints per tile: median %.0f | p90 %.0f | max %.0f\n",
            median(ntile), quantile(ntile, .9), max(ntile)))

main_effects <- paste(
  "slope_mean_z + slope_sd_z + wsci_z + rh_98_z + cover_z +",
  "aspect_sin_z + aspect_cos_z + meta_offnad_z + meta_sunel_z +",
  "meta_az_conc_z + meta_stereo_z + meta_leafon_z + meta_fwdrev_z +",
  "meta_absgeo_z + meta_relgeo_z + view_az_sin_z + view_az_cos_z"
)
priors <- c(
  prior(normal(0, 2),       class = "b"),
  prior(normal(0, 10),      class = "Intercept"),
  prior(gamma(2, 0.1),      class = "nu"),
  prior(student_t(3, 0, 2), class = "sd"),
  prior(normal(0, 1.5),     class = "sds"),
  prior(normal(0, 1),       class = "b",         dpar = "sigma"),
  prior(student_t(3, 0, 5), class = "Intercept", dpar = "sigma")
)

run_fit <- function(label, extra_re) {
  sec(paste("FIT --", label))
  f <- bf(
    as.formula(paste(
      "chm_error_mean ~", main_effects,
      "+ s(slope_mean_z, wsci_z, k = 20)",
      "+ slope_mean_z:lc_grp + wsci_z:lc_grp",
      "+ (1 + slope_mean_z | ecoregion) + (1 | site) + (1 + wsci_z | lc_grp)",
      extra_re
    )),
    sigma ~ 1 + slope_mean_z + (1 | site) + (1 | lc_grp)
  )
  t0 <- Sys.time()
  fit <- brm(formula = f, data = pil, family = student(), prior = priors,
             chains = CHAINS, iter = ITER, warmup = WARMUP,
             cores = CHAINS, threads = threading(THREADS),
             control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
             backend = "cmdstanr", seed = 2026, refresh = 50,
             file = file.path(OUTMOD, paste0("fit_pilot_", label)))
  cat(sprintf("\n%s fit: %.2f hours\n", label,
              as.numeric(difftime(Sys.time(), t0, units = "hours"))))
  np <- as.data.table(brms::nuts_params(fit))
  td <- np[Parameter == "treedepth__", Value]
  cat(sprintf("  divergences %d | max Rhat %.4f | min ESS ratio %.3f\n",
              sum(np[Parameter == "divergent__", Value]),
              max(brms::rhat(fit), na.rm = TRUE),
              min(brms::neff_ratio(fit), na.rm = TRUE)))
  cat(sprintf("  mean treedepth %.2f | %% at ceiling(%d) %.2f%%\n",
              mean(td), MAX_TREEDEPTH, 100 * mean(td >= MAX_TREEDEPTH)))
  fit
}

fit_tile <- run_fit("tile", "+ (1 | tile)")
fit_ctrl <- run_fit("control", "")

sec("COMPARISON -- what the tile level does")
vc <- VarCorr(fit_tile)$tile$sd
cat("tile sd posterior (Intercept):\n"); print(round(vc, 4))
fx_t <- fixef(fit_tile); fx_c <- fixef(fit_ctrl)
acq <- intersect(ACQ10, rownames(fx_t))
cat("\nacquisition fixed effects: estimate (control -> tile) | posterior SD ratio tile/control\n")
for (v in acq)
  cat(sprintf("  %-16s %+.3f -> %+.3f   | SD x%.2f\n",
              v, fx_c[v, "Estimate"], fx_t[v, "Estimate"],
              fx_t[v, "Est.Error"] / fx_c[v, "Est.Error"]))
cat("\nDone. Interpret against plan section 3 (Path A vs B) before launching full refits.\n")
