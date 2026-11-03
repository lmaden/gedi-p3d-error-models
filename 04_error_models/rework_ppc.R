# rework_ppc.R  (supplement: posterior predictive check of the three full fits)
# Paper: Fig. S1 and Text S2.
# For each full fit, draws 200 replicate error sets from the posterior predictive distribution
# at the training footprints (fit$data, 68,786 rows) and compares them with the observed errors:
#   1. share of footprints per 0.5 m bin from -30 to 30 m (the Fig. 3a bins), observed and
#      the median and 5th to 95th percentile over the replicates;
#   2. per replicate: mean, SD, quantiles, share beyond the observed mean +/- 3 observed SD
#      (low and high tail), share with |error| > 5 m and > 10 m; the same for the observed errors.
# Reads  $CH1/models/rework/fit_<resp>_15.rds (read only; nothing is refit)
# Writes $CH1/manuscript_tables/rework/ppc/{bins,stats}_<resp>.csv and ppc.log
# One core, about 4-8 GB RAM (as rework_divcheck.R), a few minutes per model.
# RUN: Rscript rework_ppc.R dsm dtm chm
suppressPackageStartupMessages({ library(brms); library(data.table) })

RESPS <- commandArgs(trailingOnly = TRUE)
if (!length(RESPS)) RESPS <- c("dsm", "dtm", "chm")
CH1 <- Sys.getenv("PROJECT_ROOT", "<CLUSTER_ROOT>/chpt1")
OUT <- file.path(CH1, "manuscript_tables", "rework", "ppc")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
NDRAWS <- 200L; CHUNK_ROWS <- 5000L
BW <- 0.5; EDGES <- seq(-30, 30, by = BW)

shares <- function(e) {
  k <- findInterval(e, EDGES, rightmost.closed = TRUE)
  100 * tabulate(k[k >= 1 & k < length(EDGES)], nbins = length(EDGES) - 1) / length(e)
}
stats <- function(e, m0, s0) {
  q <- quantile(e, c(0.01, 0.05, 0.5, 0.95, 0.99), names = FALSE)
  data.table(mean = mean(e), sd = sd(e), q01 = q[1], q05 = q[2], q50 = q[3], q95 = q[4], q99 = q[5],
             pct_low3sd = 100 * mean(e < m0 - 3 * s0), pct_high3sd = 100 * mean(e > m0 + 3 * s0),
             pct_over5 = 100 * mean(abs(e) > 5), pct_over10 = 100 * mean(abs(e) > 10))
}

for (r in RESPS) {
  f <- file.path(CH1, "models", "rework", sprintf("fit_%s_15.rds", r))
  cat(sprintf("\n==== %s  (%s)  %s ====\n", r, f, format(Sys.time())))
  fit <- readRDS(f)
  y <- fit$data$y
  n <- length(y)
  cat(sprintf("training footprints %d | observed mean %+.3f  SD %.3f\n", n, mean(y), sd(y)))
  if (n != 68786L) cat("** NOTE: expected 68,786 training footprints **\n")
  set.seed(2026)
  ids <- sort(sample.int(ndraws(fit), NDRAWS))
  starts <- seq(1L, n, by = CHUNK_ROWS)
  yrep <- do.call(cbind, lapply(starts, function(a) posterior_predict(
    fit, newdata = fit$data[a:min(a + CHUNK_ROWS - 1L, n), ], draw_ids = ids)))
  stopifnot(nrow(yrep) == NDRAWS, ncol(yrep) == n)
  m0 <- mean(y); s0 <- sd(y)

  sh <- t(apply(yrep, 1, shares))
  bins <- data.table(response = r, x = head(EDGES, -1) + BW / 2, obs = shares(y),
                     rep_median = apply(sh, 2, median),
                     rep_q05 = apply(sh, 2, quantile, 0.05), rep_q95 = apply(sh, 2, quantile, 0.95))
  fwrite(bins, file.path(OUT, sprintf("bins_%s.csv", r)))

  st <- rbind(cbind(data.table(response = r, set = "observed", draw = NA_integer_), stats(y, m0, s0)),
              rbindlist(lapply(seq_len(NDRAWS), function(i)
                cbind(data.table(response = r, set = "replicate", draw = ids[i]), stats(yrep[i, ], m0, s0)))))
  fwrite(st, file.path(OUT, sprintf("stats_%s.csv", r)))

  sm <- st[set == "replicate", lapply(.SD, function(v) sprintf("%.2f [%.2f, %.2f]", median(v),
           quantile(v, 0.05), quantile(v, 0.95))), .SDcols = names(st)[-(1:3)]]
  cat("observed:\n"); print(st[set == "observed", -(1:3)], digits = 3)
  cat("replicates, median [5th, 95th percentile]:\n"); print(t(sm))
  rm(fit, yrep, sh); gc()
}
cat("\nfinished:", format(Sys.time()), "\n")
