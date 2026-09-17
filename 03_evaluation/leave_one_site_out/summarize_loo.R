library(data.table)
L <- file.path(Sys.getenv("PROJECT_ROOT"), "manuscript_tables", "loo")
d <- rbindlist(lapply(list.files(L, "^loo_metrics_chm_site[0-9]+\\.csv$",
                                 full.names = TRUE), fread), fill = TRUE)
d[, site := as.integer(site)]
in0 <- list.files(L, "^loo_metrics_chm_site[0-9]+_init0\\.csv$", full.names = TRUE)
if (length(in0)) {                         # prefer the clean init = 0 reruns
  r0 <- rbindlist(lapply(in0, fread), fill = TRUE)
  r0[, site := as.integer(site)]
  d <- rbind(d[!site %in% unique(r0$site)], r0, fill = TRUE)
  cat("(CHM uses the init = 0 reruns for", uniqueN(r0$site), "folds)\n")
}
cat("folds:", uniqueN(d$site), " rows:", nrow(d), " n:", sum(d[rung == unique(rung)[1], n]), "\n\n")

f <- unique(d[, .(site, n_held, fit_hours = round(fit_hours, 2), divergences,
                  max_rhat = round(max_rhat, 4), min_ess = round(min_ess_ratio, 3),
                  pct_ceiling = round(pct_at_ceiling, 2), ok = fold_diagnostics_ok)])
cat("=== DIAGNOSTICS ===\n"); print(f[order(-divergences, site)])
cat("\nfolds with divergences:", sum(f$divergences > 0),
    " total:", sum(f$divergences),
    " | folds failing:", sum(!f$ok),
    " | compute:", round(sum(f$fit_hours), 1), "h\n\n")

cat("=== CHM RESULTS BY RUNG (n-weighted) ===\n")
print(d[, .(folds = .N, n = sum(n),
            rmse  = round(sqrt(sum(rmse^2 * n) / sum(n)), 3),
            bias  = round(sum(bias  * n) / sum(n), 3),
            crps  = round(sum(crps  * n) / sum(n), 3),
            cov50 = round(sum(cov50 * n) / sum(n), 3),
            cov90 = round(sum(cov90 * n) / sum(n), 3),
            cov95 = round(sum(cov95 * n) / sum(n), 3)), by = rung][order(rung)])

cat("\n=== ORACLE GAP PER SITE, largest sites first ===\n")
w <- dcast(d, site + n_held ~ rung, value.var = "rmse", fun.aggregate = mean)
setnames(w, 3:5, c("oracle", "rung2", "rung3"))
w[, gap_m := round(rung2 - oracle, 3)][, gap_pct := round(100 * (rung2 - oracle) / oracle, 1)]
print(w[order(-n_held), .(site, n_held, oracle = round(oracle, 3),
                          rung2 = round(rung2, 3), rung3 = round(rung3, 3), gap_m, gap_pct)])
cat("\nmedian gap:", round(median(w$gap_m), 3), "m |", round(median(w$gap_pct), 1), "%\n")
cat("negative gaps:", sum(w$gap_m < 0), "of", nrow(w), "\n")

cat("\n\n=== THE COMPARISON: CHM vs DTM ===\n")
dt <- rbindlist(lapply(list.files(L, "^loo_metrics_dtm_site[0-9]+\\.csv$",
                                  full.names = TRUE), fread), fill = TRUE)
ad <- list.files(L, "^loo_metrics_dtm_site[0-9]+_ad099\\.csv$", full.names = TRUE)
if (length(ad)) {                       # prefer the clean adapt_delta 0.99 reruns
  r <- rbindlist(lapply(ad, fread), fill = TRUE)
  dt <- rbind(dt[!site %in% unique(r$site)], r, fill = TRUE)
  cat("(DTM uses the adapt_delta 0.99 reruns for", uniqueN(r$site), "folds)\n")
}
pool <- function(x, lab) x[, .(product = lab, folds = .N, n = sum(n),
    rmse = round(sqrt(sum(rmse^2 * n) / sum(n)), 3),
    crps = round(sum(crps * n) / sum(n), 3),
    cov90 = round(sum(cov90 * n) / sum(n), 3)), by = rung]
print(rbind(pool(d, "CHM"), pool(dt, "DTM"))[order(rung, product)])

g <- function(x) { v <- dcast(x, site ~ rung, value.var = "rmse", fun.aggregate = mean)
                   100 * (v[[3]] - v[[2]]) / v[[2]] }
cat("\nmedian oracle gap:  CHM", round(median(g(d)), 1), "%   DTM", round(median(g(dt)), 1), "%\n")
