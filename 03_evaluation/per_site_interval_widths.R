# =============================================================================
# per_site_interval_widths.R
#
# Part of the 16-site CHM reconciliation.
# Recomputes 95% posterior prediction-interval widths + empirical coverage
# on a HOLDOUT defined as the exact complement of fit_chm_16's training rows,
# using the PRIMARY 16-site fit (fit_chm_16 in 10b_chm_sensitivity), replacing
# the deprecated 19-site fit_chm_s2 numbers that currently feed manuscript
# Section 3.5 and the Section-2A "94% coverage" anchor.
#
# METHOD MATCHES the deprecated run (compute_prediction_interval_widths.R):
#   - posterior_predict() full predictive (Student-t, identity link)
#   - 200 draws, allow_new_levels = TRUE, chunked
#   - width = q0.975 - q0.025 ; coverage = mean(obs in [q0.025, q0.975])
# DIFFERS only in: model = fit_chm_16 (16-site), and holdout = anti-join
#   complement of fit_chm_16$data (NOT a re-run seeded split, which misses by 6).
#
# HARD GATES (abort before any predict if violated):
#   G1  model is fit_chm_16, N == 124,966, 16 sites (no 1/2/3)
#   G2  anti-join matched(train) == 124,966   (holdout provably disjoint)
#   G3  holdout N == 253,700
#   G4  no new site/ecoregion levels in holdout (LMS in lc_l1_code is 2 rows, OK)
#
# Output: tables/prediction_interval_widths_chm16.csv  (overall + by site + by forest)
#         console: overall median/IQR/range + NEW 16-site coverage (CHM)
#         (DTM is unchanged elsewhere and intentionally NOT recomputed here)
# =============================================================================

Sys.setenv(DISPLAY = "")
options(device = pdf)

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tidyr)
  library(readr)
})

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1")
CK           <- file.path(PROJECT_ROOT, "checkpoints")
OUT          <- file.path(PROJECT_ROOT, "tables")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

N_DRAWS    <- 200
CHUNK_SIZE <- 15000

FOREST_CLASSES <- c("BDF", "DNF", "EBF", "ENF")
FOREST_LABELS  <- c(BDF = "Broadleaf Deciduous", DNF = "Deciduous Needleleaf",
                    EBF = "Evergreen Broadleaf", ENF = "Evergreen Needleleaf")

cat("=== 16-site CHM Prediction-Interval Widths + Coverage (PRIMARY fit) ===\n\n")

# ---------------------------------------------------------------------------
# LOAD primary 16-site fit + raw data
# ---------------------------------------------------------------------------
cat("Loading fit_chm_16 (10b_chm_sensitivity) ...\n")
m   <- readRDS(file.path(CK, "10b_chm_sensitivity.rds"))
fit <- m$data$fit_chm_16
stopifnot(inherits(fit, "brmsfit"))

cp  <- readRDS(file.path(CK, "01_data_ingest.rds"))
chm <- cp$data$chm_df

mvars <- all.vars(m$data$chm_formula_s2$formula)   # response + 17 _z + 3 factors
SITES_16 <- as.character(c(4:15, 17:20))

# ---------------------------------------------------------------------------
# GATE G1: model identity
# ---------------------------------------------------------------------------
fit_sites <- sort(unique(as.character(fit$data$site)))
cat(sprintf("G1 model: N=%s  n_sites=%s\n", nrow(fit$data), length(fit_sites)))
stopifnot(nrow(fit$data) == 124966L)
stopifnot(setequal(fit_sites, SITES_16))
stopifnot(!any(c("1","2","3") %in% fit_sites))
cat("   G1 PASS: fit_chm_16, N=124,966, 16 sites (no 1/2/3)\n")

# ---------------------------------------------------------------------------
# Build the 16-site predictor-complete universe and the anti-join HOLDOUT
# ---------------------------------------------------------------------------
U <- chm %>%
  filter(if_all(all_of(intersect(mvars, colnames(chm))), ~ !is.na(.))) %>%
  mutate(site = as.character(site)) %>%
  filter(site %in% SITES_16)

key <- function(d) paste(round(d$chm_error_mean, 6), round(d$rh_98_z, 6),
                         round(d$slope_mean_z, 6), as.character(d$site), sep = "|")

train_keys <- key(fit$data)
is_train   <- key(U) %in% train_keys
H          <- U[!is_train, ]

# ---------------------------------------------------------------------------
# GATE G2/G3: holdout disjointness + size
# ---------------------------------------------------------------------------
matched <- sum(is_train)
cat(sprintf("G2 anti-join matched(train) = %s  (target 124,966)\n", matched))
stopifnot(matched == 124966L)
cat(sprintf("G3 holdout N = %s  (target 253,700)\n", nrow(H)))
stopifnot(nrow(H) == 253700L)
cat("   G2/G3 PASS: holdout is the provable complement of training\n")

# ---------------------------------------------------------------------------
# GATE G4: no unseen site/ecoregion levels (LMS in lc_l1_code is 2 rows, OK)
# ---------------------------------------------------------------------------
new_site <- setdiff(unique(as.character(H$site)),      as.character(levels(fit$data$site)))
new_eco  <- setdiff(unique(as.character(H$ecoregion)), as.character(levels(fit$data$ecoregion)))
new_lc   <- setdiff(unique(as.character(H$lc_l1_code)),as.character(levels(fit$data$lc_l1_code)))
stopifnot(length(new_site) == 0, length(new_eco) == 0)
cat(sprintf("G4 PASS: no new site/eco levels; new lc levels = {%s} (n=%s rows, population-level)\n",
            paste(new_lc, collapse = ","),
            sum(as.character(H$lc_l1_code) %in% new_lc)))

# brms validate_newdata() rejects ANY lc_l1_code level absent from training,
# even under allow_new_levels (that flag covers grouping factors like site/
# ecoregion, not population-level categorical predictors). LMS is 2 rows
# (0.001%); drop them so the holdout uses exactly the trained lc strata.
n_drop_lc <- sum(as.character(H$lc_l1_code) %in% new_lc)
if (n_drop_lc > 0) {
  cat(sprintf("Dropping %s holdout row(s) with untrained lc_l1_code level {%s} (%.4f%%).\n",
              n_drop_lc, paste(new_lc, collapse = ","), 100 * n_drop_lc / nrow(H)))
  H <- H[!(as.character(H$lc_l1_code) %in% new_lc), ]
}

# Match factor coding to the fitted model so posterior_predict aligns terms
H <- H %>%
  mutate(
    site       = factor(as.character(site),       levels = levels(fit$data$site)),
    ecoregion  = factor(as.character(ecoregion),  levels = levels(fit$data$ecoregion)),
    lc_l1_code = factor(as.character(lc_l1_code), levels = levels(fit$data$lc_l1_code))
  )
stopifnot(!any(is.na(H$lc_l1_code)), !any(is.na(H$site)), !any(is.na(H$ecoregion)))

cat(sprintf("\nAll gates passed. Holdout N = %s. Proceeding to posterior_predict.\n\n", nrow(H)))

# ---------------------------------------------------------------------------
# CHUNKED posterior_predict -> per-obs 95% PI
# ---------------------------------------------------------------------------
predict_intervals_chunked <- function(model, newdata, ndraws = N_DRAWS,
                                       chunk_size = CHUNK_SIZE) {
  n      <- nrow(newdata)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  cat(sprintf("  Predicting %d obs in %d chunks (%d draws)...\n",
              n, length(chunks), ndraws))
  res <- vector("list", length(chunks))
  for (j in seq_along(chunks)) {
    idx <- chunks[[j]]
    pp  <- posterior_predict(model, newdata = newdata[idx, ],
                             ndraws = ndraws, allow_new_levels = TRUE)
    q025 <- apply(pp, 2, quantile, probs = 0.025, na.rm = TRUE)
    q975 <- apply(pp, 2, quantile, probs = 0.975, na.rm = TRUE)
    res[[j]] <- data.frame(pi_lower = q025, pi_upper = q975,
                           pi_width = q975 - q025, stringsAsFactors = FALSE)
    cat(sprintf("    Chunk %d/%d done\n", j, length(chunks)))
  }
  bind_rows(res)
}

cat("--- CHM (16-site) prediction intervals ---\n")
chm_pi <- predict_intervals_chunked(fit, H)
chm_pi$site       <- as.character(H$site)
chm_pi$lc_l1_code <- as.character(H$lc_l1_code)
chm_pi$observed   <- H$chm_error_mean

# ---------------------------------------------------------------------------
# SUMMARIES (identical definitions to the deprecated run)
# ---------------------------------------------------------------------------
summarize_pi <- function(d, product, stratum) {
  data.frame(
    product = product, stratum = stratum, n = nrow(d),
    median_width = median(d$pi_width, na.rm = TRUE),
    mean_width   = mean(d$pi_width,   na.rm = TRUE),
    q25_width    = quantile(d$pi_width, 0.25, na.rm = TRUE),
    q75_width    = quantile(d$pi_width, 0.75, na.rm = TRUE),
    min_width    = min(d$pi_width, na.rm = TRUE),
    max_width    = max(d$pi_width, na.rm = TRUE),
    coverage_95  = mean(d$observed >= d$pi_lower & d$observed <= d$pi_upper, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

overall  <- summarize_pi(chm_pi, "CHM", "Overall")

by_site  <- do.call(rbind, lapply(split(chm_pi, chm_pi$site), function(d)
  summarize_pi(d, "CHM", paste0("site=", d$site[1]))))
by_site  <- by_site[order(as.numeric(sub("site=", "", by_site$stratum))), ]

chm_forest <- chm_pi %>% filter(lc_l1_code %in% FOREST_CLASSES)
by_forest  <- do.call(rbind, lapply(split(chm_forest, chm_forest$lc_l1_code), function(d)
  summarize_pi(d, "CHM", paste0("lc_l1_code=", FOREST_LABELS[d$lc_l1_code[1]]))))

# ---------------------------------------------------------------------------
# REPORT (these are the numbers that replace the §3.5 / §2A CHM values)
# ---------------------------------------------------------------------------
cat("\n=====================================================================\n")
cat("  16-site CHM (PRIMARY) — 95% PI widths + coverage on holdout\n")
cat("=====================================================================\n")
cat(sprintf("  Overall: median = %.2f m  IQR = [%.2f, %.2f] m  range = [%.2f, %.2f] m\n",
            overall$median_width, overall$q25_width, overall$q75_width,
            overall$min_width, overall$max_width))
cat(sprintf("  Overall 95%% coverage = %.4f  (= %.1f%%)\n",
            overall$coverage_95, 100 * overall$coverage_95))
cat(sprintf("  Median-of-site-medians = %.2f m ; site-median range = %.2f - %.2f m\n",
            median(by_site$median_width),
            min(by_site$median_width), max(by_site$median_width)))

cat("\n  By site (tracker numbering; map to manuscript later):\n")
for (i in seq_len(nrow(by_site))) {
  r <- by_site[i, ]
  cat(sprintf("    %-9s median=%6.2f  IQR=[%5.2f,%6.2f]  cov=%.3f  n=%s\n",
              r$stratum, r$median_width, r$q25_width, r$q75_width,
              r$coverage_95, format(r$n, big.mark = ",")))
}

cat("\n  Narrowest / widest site median:\n")
o <- by_site[order(by_site$median_width), ]
cat(sprintf("    narrowest: %s  %.2f m (n=%s, cov=%.3f)\n",
            o$stratum[1], o$median_width[1], format(o$n[1], big.mark=","), o$coverage_95[1]))
cat(sprintf("    widest:    %s  %.2f m (n=%s, cov=%.3f)\n",
            o$stratum[nrow(o)], o$median_width[nrow(o)],
            format(o$n[nrow(o)], big.mark=","), o$coverage_95[nrow(o)]))
cat(sprintf("    fold-range (max/min of site medians) = %.2fx\n",
            max(by_site$median_width) / min(by_site$median_width)))

cat("\n  By forest type:\n")
for (i in seq_len(nrow(by_forest))) {
  r <- by_forest[i, ]
  cat(sprintf("    %-40s median=%6.2f  cov=%.3f  n=%s\n",
              r$stratum, r$median_width, r$coverage_95, format(r$n, big.mark = ",")))
}

# ---------------------------------------------------------------------------
# SAVE (16-site CSV, parallel schema to the deprecated file; CHM only)
# ---------------------------------------------------------------------------
all_results <- bind_rows(
  overall %>% mutate(level = "overall"),
  by_site %>% mutate(level = "site"),
  by_forest %>% mutate(level = "forest_type")
)
out_csv <- file.path(OUT, "prediction_interval_widths_chm16.csv")
write_csv(all_results, out_csv)
cat(sprintf("\n[OK] wrote %s  (%d rows)\n", out_csv, nrow(all_results)))

# raw per-obs PI for any follow-up (e.g., the J4a figure)
saveRDS(list(chm16_pi = chm_pi),
        file.path(CK, "holdout_prediction_intervals_chm16.rds"), compress = "xz")
cat(sprintf("[OK] wrote %s\n", file.path(CK, "holdout_prediction_intervals_chm16.rds")))
cat("\nDone.\n")
