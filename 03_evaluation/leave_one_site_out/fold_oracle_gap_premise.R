library(data.table)
p <- readRDS(file.path(Sys.getenv("PROJECT_ROOT"), "checkpoints", "08_model_prep.rds"))$data

chk <- function(lab, dat, resp, drop_sites) {
  d <- as.data.table(dat)
  sc <- grep("^site", names(d), value = TRUE, ignore.case = TRUE)[1]
  cat("\n==", lab, "==\n")
  cat("site column     :", sc, "\n")
  d[, .sid := as.character(get(sc))]
  cat("site labels     :", paste(head(sort(unique(d$.sid)), 25), collapse = " "), "\n")
  if (length(drop_sites)) d <- d[!.sid %in% drop_sites]
  y <- d[[resp]]; y <- y[is.finite(y)]
  cat("rows used       :", length(y), "\n")
  cat("SD of error     :", round(sd(y), 3), "m\n")
  cat("TOTAL error var :", round(var(y), 3), "m^2\n")
  invisible(var(y))
}

vc <- chk("CHM (16 sites)", p$mod_chm_s2, "chm_error_mean", c("1","2","3"))
vd <- chk("DTM (18 sites)", p$mod_dtm_s2, "dtm_error_mean", character(0))

cat("\n\n== IS THE RECONCILIATION TRUE? ==\n")
cat("(RE variances established: CHM 2.38 m^2, site 67.2%; DTM 1.00 m^2, site 58.3%)\n\n")
for (r in list(list("CHM", 2.38, 0.672, vc, 4.445, 4.831),
               list("DTM", 1.00, 0.583, vd, 3.003, 3.362))) {
  sv <- r[[2]] * r[[3]]
  cat(sprintf("%s  total err var %.3f m^2\n", r[[1]], r[[4]]))
  cat(sprintf("     all RE  = %5.1f%% of total error variance\n", 100 * r[[2]] / r[[4]]))
  cat(sprintf("     site RE = %5.1f%% of total error variance\n", 100 * sv / r[[4]]))
  cat(sprintf("     residual= %5.1f%% of total error variance\n", 100 * (1 - r[[2]] / r[[4]])))
  cat(sprintf("     predicted RMSE cost of losing site: %.1f%%   observed gap: %.1f%%\n\n",
              100 * (1 / sqrt(1 - sv / r[[4]]) - 1), 100 * (r[[6]] - r[[5]]) / r[[5]]))
}
