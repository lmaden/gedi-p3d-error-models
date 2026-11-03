# rework_divcheck.R  (supplement: are the reported numbers moved by the divergent draws?)
# Paper: Text S2 (divergent draws).
# For the full fits with divergences (DTM 5, DSM 2 of 6,000 post-warmup draws) this
#   1. computes R2 total and R2 fixed-only for EVERY draw (in chunks of draws, so memory
#      stays small) and summarises them with all draws and without the divergent ones;
#   2. does the same for the group SDs, RH98 slope and nu;
#   3. places each divergent draw within the posterior of every group SD (percentile),
#      which shows whether the divergences sit in a funnel (group SD near zero).
# Reads  $CH1/models/rework/fit_<resp>_15.rds (read only; nothing is refit)
# Writes $CH1/manuscript_tables/rework/divcheck/{summary,divergent_draws}_<resp>.csv
# One core, about 4-8 GB RAM, about 10-30 minutes per model.
# RUN: Rscript rework_divcheck.R dtm dsm > divcheck.out 2>&1
suppressPackageStartupMessages({ library(brms); library(data.table); library(posterior) })

RESPS <- commandArgs(trailingOnly = TRUE)
if (!length(RESPS)) RESPS <- c("dtm", "dsm")
CH1 <- Sys.getenv("PROJECT_ROOT", "<CLUSTER_ROOT>/chpt1")
OUT <- file.path(CH1, "manuscript_tables", "rework", "divcheck")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CHUNK <- 250L

for (r in RESPS) {
  f <- file.path(CH1, "models", "rework", sprintf("fit_%s_15.rds", r))
  cat(sprintf("\n==== %s  (%s) ====\n", r, f))
  fit <- readRDS(f)
  dr  <- as_draws_df(fit)
  ids <- data.table(Chain = dr$.chain, Iteration = dr$.iteration, draw = dr$.draw)
  np  <- as.data.table(nuts_params(fit))[Parameter == "divergent__"]
  np  <- merge(np, ids, by = c("Chain", "Iteration"))
  div <- np[Value == 1, draw]
  N   <- nrow(ids)
  cat(sprintf("draws %d | divergent %d (%.3f%%) | draw ids: %s\n", N, length(div),
              100 * length(div) / N, paste(sort(div), collapse = " ")))

  # R2 for every draw, in chunks of draws
  r2t <- r2f <- numeric(N)
  for (a in seq(1L, N, by = CHUNK)) {
    d_ <- a:min(a + CHUNK - 1L, N)
    r2t[d_] <- bayes_R2(fit, draw_ids = d_, summary = FALSE)[, 1]
    r2f[d_] <- bayes_R2(fit, draw_ids = d_, summary = FALSE, re_formula = NA)[, 1]
    if (a %% 1500L == 1L) cat(sprintf("  R2 draws %d-%d done\n", a, max(d_)))
  }
  pars <- grep("^sd_|^b_rh_98_z$|^nu$", variables(dr), value = TRUE)
  vals <- c(list(R2_total = r2t, R2_fixed_only = r2f),
            lapply(setNames(pars, pars), function(p) as.numeric(dr[[p]])))
  keep <- !(seq_len(N) %in% div)
  sm <- rbindlist(lapply(names(vals), function(q) {
    v <- vals[[q]]
    data.table(response = r, quantity = q,
               median_all = median(v), median_without = median(v[keep]),
               q2.5_all = quantile(v, 0.025), q2.5_without = quantile(v[keep], 0.025),
               q97.5_all = quantile(v, 0.975), q97.5_without = quantile(v[keep], 0.975))
  }))
  sm[, change_in_median := median_without - median_all]
  fwrite(sm, file.path(OUT, sprintf("summary_%s.csv", r)))
  cat("\nmedian with all draws vs without the divergent draws:\n")
  print(sm[, .(quantity, all = signif(median_all, 4), without = signif(median_without, 4),
               change = signif(change_in_median, 2))])

  if (!length(div)) { cat("no divergent draws\n"); rm(fit, dr); gc(); next }
  # where the divergent draws sit (percentile within each posterior)
  sdp <- grep("^sd_", names(vals), value = TRUE)
  dd <- rbindlist(lapply(div, function(k) {
    data.table(response = r, draw = k, parameter = c("R2_total", sdp),
               value = sapply(c("R2_total", sdp), function(q) vals[[q]][k]),
               percentile = sapply(c("R2_total", sdp), function(q) round(100 * mean(vals[[q]] <= vals[[q]][k]))))
  }))
  fwrite(dd, file.path(OUT, sprintf("divergent_draws_%s.csv", r)))
  cat("\npercentile of each divergent draw within the posterior (low = funnel neck):\n")
  print(dcast(dd, parameter ~ draw, value.var = "percentile"))
  rm(fit, dr); gc()
}
cat("\ndone\n")
quit(save = "no")
