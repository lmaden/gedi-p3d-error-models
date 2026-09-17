library(data.table)
L <- file.path(Sys.getenv("PROJECT_ROOT"), "manuscript_tables", "loo")
rd <- function(pat) {
  f <- list.files(L, pat, full.names = TRUE)
  if (!length(f)) return(NULL)
  d <- rbindlist(lapply(f, fread), fill = TRUE); d[, site := as.integer(site)]; d
}
new <- rd("^loo_metrics_dtm_site[0-9]+_ad099\\.csv$")
old <- rd("^loo_metrics_dtm_site[0-9]+\\.csv$")
cat("rerun folds:", uniqueN(new$site), " original folds:", uniqueN(old$site), "\n\n")

cat("=== DID 0.99 FIX THE DIVERGENCES? ===\n")
a <- unique(old[site %in% unique(new$site),
                .(site, div_098 = divergences, rhat_098 = round(max_rhat, 4),
                  hours_098 = round(fit_hours, 2))])
b <- unique(new[, .(site, div_099 = divergences, rhat_099 = round(max_rhat, 4),
                    hours_099 = round(fit_hours, 2), ok = fold_diagnostics_ok)])
cmp <- merge(a, b, by = "site")
print(cmp[order(-div_098)])
cat("\ndivergences before:", sum(cmp$div_098), " after:", sum(cmp$div_099), "\n")
cat("folds now passing:", sum(cmp$ok), "of", nrow(cmp), "\n")
cat("extra compute:", round(sum(cmp$hours_099) - sum(cmp$hours_098), 1), "hours\n\n")

cat("=== DID THE ESTIMATES MOVE? (RMSE per site per rung) ===\n")
m <- merge(old[site %in% unique(new$site), .(site, rung, rmse_098 = rmse, cov90_098 = cov90)],
           new[, .(site, rung, rmse_099 = rmse, cov90_099 = cov90)],
           by = c("site", "rung"))
m[, d_rmse := round(rmse_099 - rmse_098, 4)]
m[, d_pct  := round(100 * (rmse_099 - rmse_098) / rmse_098, 2)]
print(m[order(-abs(d_pct)), .(site, rung, rmse_098 = round(rmse_098, 3),
                              rmse_099 = round(rmse_099, 3), d_rmse, d_pct)])
cat("\nlargest RMSE change:", round(max(abs(m$d_pct)), 2), "%\n")
cat("median absolute change:", round(median(abs(m$d_pct)), 2), "%\n")
