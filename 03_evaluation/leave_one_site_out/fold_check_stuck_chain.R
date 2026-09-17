suppressPackageStartupMessages({library(brms); library(posterior); library(data.table)})
M <- file.path(Sys.getenv("PROJECT_ROOT"), "models", "loo")
for (S in c(11, 20)) {
  f <- file.path(M, sprintf("fit_chm_site%d.rds", S))
  if (!file.exists(f)) { cat("missing:", f, "\n"); next }
  cat("\n################ site", S, "################\n")
  fit <- readRDS(f)
  np  <- as.data.table(brms::nuts_params(fit))
  tdc <- np[Parameter == "treedepth__",
            .(mean_td = round(mean(Value), 2),
              pct_ceiling = round(100 * mean(Value >= 14), 1)), by = Chain][order(Chain)]
  print(tdc)
  stuck   <- tdc[pct_ceiling > 50, Chain]
  healthy <- setdiff(tdc$Chain, stuck)
  cat("stuck chain(s):", paste(stuck, collapse = " "),
      "| healthy:", paste(healthy, collapse = " "), "\n")
  if (!length(stuck) || !length(healthy)) { cat("no clean split\n"); next }

  dr <- posterior::as_draws_array(fit)
  b  <- grep("^b_", posterior::variables(dr), value = TRUE)
  sa <- as.data.table(posterior::summarise_draws(
          posterior::subset_draws(dr, variable = b), "rhat", "ess_bulk", "mean"))
  sh <- as.data.table(posterior::summarise_draws(
          posterior::subset_draws(dr, variable = b, chain = healthy),
          "rhat", "ess_bulk", "mean"))
  cat(sprintf("\nmax Rhat, all 4 chains ....... %.4f\n", max(sa$rhat, na.rm = TRUE)))
  cat(sprintf("max Rhat, healthy chains only  %.4f   <-- under 1.01 means unimodal\n",
              max(sh$rhat, na.rm = TRUE)))
  cat(sprintf("min ESS ratio, healthy only .. %.3f\n",
              min(sh$ess_bulk / (length(healthy) * 1500), na.rm = TRUE)))
  m <- merge(sa[, .(variable, mean_all = mean)],
             sh[, .(variable, mean_healthy = mean)], by = "variable")
  m[, shift := round(mean_healthy - mean_all, 4)]
  cat("\nlargest coefficient shifts when the stuck chain is dropped:\n")
  print(head(m[order(-abs(shift)), .(variable, mean_all = round(mean_all, 3),
                                     mean_healthy = round(mean_healthy, 3), shift)], 8))
  rm(fit, dr); gc(verbose = FALSE)
}
