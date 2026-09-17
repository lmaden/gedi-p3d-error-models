library(data.table)
L <- file.path(Sys.getenv("PROJECT_ROOT"), "manuscript_tables", "loo")
g <- function(sfx) { f <- list.files(L, sprintf("^loo_metrics_dtm_site[0-9]+%s\\.csv$", sfx),
                                     full.names = TRUE)
  if (length(f)) rbindlist(lapply(f, fread), fill = TRUE)[, site := as.integer(site)][] }
d <- g(""); r <- g("_ad099")
d[, src := "orig"]; if (!is.null(r)) { r[, src := "ad099"]
  d <- rbind(d[!site %in% unique(r$site)], r, fill = TRUE) }
f <- unique(d[, .(site, src, n_held, fit_hours = round(fit_hours, 2), divergences,
                  max_rhat = round(max_rhat, 4), min_ess = round(min_ess_ratio, 3),
                  pct_ceiling = round(pct_at_ceiling, 2), ok = fold_diagnostics_ok)])
print(f[order(-divergences, -max_rhat)])
cat("\nfolds:", nrow(f), " divergences:", sum(f$divergences),
    " failing:", sum(!f$ok), " ad099 used:", sum(f$src == "ad099"),
    " compute:", round(sum(f$fit_hours), 1), "h\n")
