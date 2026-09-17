library(data.table)
PR <- Sys.getenv("PROJECT_ROOT")
pick <- function(o, nm) if (!is.null(o[[nm]])) o[[nm]] else o$data[[nm]]

cat("=== ML packages ===\n")
for (pk in c("ranger","randomForest","xgboost","gbm")) 
  cat(sprintf("  %-14s %s\n", pk, requireNamespace(pk, quietly = TRUE)))

held <- function(prod, sfx) {
  L <- file.path(PR, "manuscript_tables", "loo")
  g <- function(x) { f <- list.files(L, sprintf("^loo_metrics_%s_site[0-9]+%s\\.csv$", prod, x),
                                     full.names = TRUE)
    if (length(f)) rbindlist(lapply(f, fread), fill = TRUE)[, site := as.integer(site)][] }
  b <- g(""); r <- g(sfx)
  if (!is.null(r)) b <- rbind(b[!site %in% unique(r$site)], r, fill = TRUE)
  unique(b[, .(site, n_held)])
}

check <- function(lab, fd, prod, sfx) {
  d <- as.data.table(fd)
  cat("\n===", lab, "=== rows:", nrow(d), " cols:", ncol(d), "\n")
  print(names(d))
  a <- d[, .(fit_rows = .N), by = .(site = as.integer(as.character(site)))]
  m <- merge(a, held(prod, sfx), by = "site", all = TRUE)[order(-fit_rows)]
  m[, ok := fit_rows == n_held]
  print(m)
  cat("total fit rows:", sum(m$fit_rows), "  total n_held:", sum(m$n_held, na.rm = TRUE),
      "  ALL MATCH:", isTRUE(all(m$ok)), "\n")
}

cat("\n\n######## CHM ########\n")
o <- readRDS(file.path(PR, "checkpoints", "10b_chm_sensitivity.rds"))
check("fit_chm_16$data", pick(o, "fit_chm_16")$data, "chm", "_init0"); rm(o); gc()

cat("\n\n######## DTM ########\n")
o <- readRDS(file.path(PR, "checkpoints", "10_models_stage2.rds"))
check("fit_dtm_s2$data", pick(o, "fit_dtm_s2")$data, "dtm", "_ad099"); rm(o); gc()
