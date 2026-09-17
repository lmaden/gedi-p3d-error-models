library(data.table)
L <- file.path(Sys.getenv("PROJECT_ROOT"), "manuscript_tables", "loo")
g <- function(p, x) { f <- list.files(L, sprintf("^loo_metrics_%s_site[0-9]+%s\\.csv$", p, x),
                                      full.names = TRUE)
  if (length(f)) rbindlist(lapply(f, fread), fill = TRUE)[, site := as.integer(site)][] }

for (p in c("chm", "dtm")) {
  b <- fread(file.path(L, sprintf("m7_baselines_%s_persite.csv", p)))
  m <- g(p, ""); r <- g(p, if (p == "chm") "_init0" else "_ad099")
  if (!is.null(r)) m <- rbind(m[!site %in% unique(r$site)], r, fill = TRUE)
  ok <- m[!is.na(eco_known), .(eco_ok = any(eco_known)), by = site]
  b <- merge(b, ok, by = "site")
  cols <- setdiff(names(b), c("site", "n", "eco_ok"))
  po <- function(x) sapply(cols, function(cc) round(sqrt(sum(x[[cc]]^2*x$n)/sum(x$n)), 3))
  cat("\n########", toupper(p), "########\n")
  print(data.table(baseline = cols, all_folds = po(b), valid_only = po(b[eco_ok == TRUE])))
  cat("valid folds:", b[eco_ok == TRUE, .N], "n =", b[eco_ok == TRUE, sum(n)], "\n")
  cat("\n-- per-site, largest first --\n")
  print(b[order(-n), .(site, n, site_offset = round(site_offset,2), global = round(global,2),
                       lin = round(lin_blind,2), rf = round(rf_blind,2), eco_ok)])
}
