library(data.table)
L <- file.path(Sys.getenv("PROJECT_ROOT"), "manuscript_tables", "loo")

grab <- function(prod, suffix) {
  f <- list.files(L, sprintf("^loo_metrics_%s_site[0-9]+%s\\.csv$", prod, suffix),
                  full.names = TRUE)
  if (!length(f)) return(NULL)
  x <- rbindlist(lapply(f, fread), fill = TRUE)
  x[, site := as.integer(site)][]
}

for (prod in c("chm", "dtm")) {
  base <- grab(prod, "")
  rr   <- grab(prod, if (prod == "chm") "_init0" else "_ad099")
  if (!is.null(rr)) base <- rbind(base[!site %in% unique(rr$site)], rr, fill = TRUE)

  cat("\n=====", toupper(prod), "=====\n")
  cat("eco_known constant within fold:",
      all(base[, uniqueN(eco_known), by = site]$V1 == 1), "\n\n")

  s <- unique(base[, .(site, n_held, eco_known)])[order(-n_held)]
  print(s)

  bad <- s[eco_known == FALSE]
  cat("\nfolds with eco_known FALSE:", nrow(bad), "of", nrow(s), "-> sites",
      paste(sort(bad$site), collapse = " "), "\n")
  cat("held-out footprints in those folds:", sum(bad$n_held), "of", sum(s$n_held),
      sprintf("(%.1f%%)\n", 100 * sum(bad$n_held) / sum(s$n_held)))
}
