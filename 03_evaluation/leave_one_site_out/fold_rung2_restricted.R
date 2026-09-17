library(data.table)
L <- file.path(Sys.getenv("PROJECT_ROOT"), "manuscript_tables", "loo")

grab <- function(prod, suffix) {
  f <- list.files(L, sprintf("^loo_metrics_%s_site[0-9]+%s\\.csv$", prod, suffix),
                  full.names = TRUE)
  if (!length(f)) return(NULL)
  x <- rbindlist(lapply(f, fread), fill = TRUE)
  x[, site := as.integer(site)][]
}

pool <- function(x) x[, .(folds = .N, n = sum(n),
    rmse    = round(sqrt(sum(rmse^2 * n) / sum(n)), 3),
    bias    = round(sum(bias    * n) / sum(n), 3),
    crps    = round(sum(crps    * n) / sum(n), 3),
    width90 = round(sum(width90 * n) / sum(n), 3),
    cov50   = round(sum(cov50   * n) / sum(n), 3),
    cov90   = round(sum(cov90   * n) / sum(n), 3),
    cov95   = round(sum(cov95   * n) / sum(n), 3)), by = rung][order(rung)]

for (prod in c("chm", "dtm")) {
  d  <- grab(prod, "")
  rr <- grab(prod, if (prod == "chm") "_init0" else "_ad099")
  if (!is.null(rr)) d <- rbind(d[!site %in% unique(rr$site)], rr, fill = TRUE)

  cat("\n\n##########", toupper(prod), "##########\n")
  cat("\n-- which rung carries the eco_known flag --\n")
  print(dcast(d, rung ~ eco_known, fun.aggregate = length, value.var = "site"))

  eco <- d[!is.na(eco_known), .(eco_ok = any(eco_known)), by = site]
  eco <- merge(eco, d[, .(n_held = n_held[1]), by = site], by = "site")
  tot <- eco[, sum(n_held)]

  cat("\n-- folds, each site counted once --\n")
  cat("valid      :", eco[eco_ok == TRUE,  .N], "folds,",
      eco[eco_ok == TRUE,  sum(n_held)], "footprints\n")
  cat("degenerate :", eco[eco_ok == FALSE, .N], "folds,",
      eco[eco_ok == FALSE, sum(n_held)], "footprints",
      sprintf("(%.1f%%) -> sites %s\n", 100 * eco[eco_ok == FALSE, sum(n_held)] / tot,
              paste(sort(eco[eco_ok == FALSE, site]), collapse = " ")))

  cat("\n-- ALL FOLDS (what is reported now) --\n"); print(pool(d))
  cat("\n-- VALID FOLDS ONLY (rung 2 is a real measurement) --\n")
  print(pool(d[site %in% eco[eco_ok == TRUE, site]]))
}
