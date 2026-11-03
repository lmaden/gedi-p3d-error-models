# tableS5_coefficients.R  --  Table S5, effects of the tested variables in the three error models
# Rows: the 17 variables of Table S1 in the same order, then the group offset SDs.
# Cells: posterior median (95% interval), metres of error per SD of the variable (group rows: SD in m).
# Inputs  (P3D_DERIVED) full_fits/{dsm,dtm,chm}/fixef.csv (rework_fullfit.R), zconst.csv
#         (P3D_RESULTS) tableS1/tableS1_predictors.csv (tableS1_predictors.R)
# Outputs (P3D_RESULTS) tableS5/tableS5_coefficients.{csv,md}
# Checks  every Table S1 variable found once per product; the numbers Results 3.2 cites
#         (RH98 effects +0.79 DTM, +0.42 DSM, -0.24 CHM; no other median beyond about 0.4 m).
suppressPackageStartupMessages(library(data.table))
DERIVED  <- Sys.getenv("P3D_DERIVED", "data_derived")
OUT_ROOT <- Sys.getenv("P3D_RESULTS", "results")
OUT  <- file.path(OUT_ROOT, "tableS5"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PRODUCTS <- c("DSM", "DTM", "CHM")

s1 <- fread(file.path(OUT_ROOT, "tableS1", "tableS1_predictors.csv"))[variable != ""]
zc <- fread(file.path(DERIVED, "zconst.csv"))
vars <- merge(s1[, .(variable, group, name, ord = .I)], zc[, .(variable, z)], by = "variable")
setorder(vars, ord); stopifnot(nrow(vars) == 17L)
vars[, par := paste0("b_", z)]
GROUPS <- data.table(par = c("sd_site__Intercept", "sd_ecoregion__Intercept", "sd_lc_grp__Intercept"),
                     group = "Group offsets (SD, m)", name = c("Site", "Level II ecoregion", "Land-cover class"))

fx <- rbindlist(lapply(PRODUCTS, function(p) {
  x <- fread(file.path(DERIVED, "full_fits", tolower(p), "fixef.csv"))
  setnames(x, c("q2.5", "q97.5"), c("lo", "hi"))
  x[parameter %in% c(vars$par, GROUPS$par), .(product = p, par = parameter, median, lo, hi)]
}))
stopifnot(fx[, .N, by = product]$N == 20L)

# checks (Results 3.2)
rh <- fx[par == "b_rh_98_z"]; setkey(rh, product)
cat(sprintf("RH98 per SD: DSM %+.2f, DTM %+.2f, CHM %+.2f (text +0.42, +0.79, -0.24)\n",
            rh["DSM"]$median, rh["DTM"]$median, rh["CHM"]$median))
oth <- fx[par %in% vars$par & par != "b_rh_98_z"][order(-abs(median))][1]
cat(sprintf("largest other median: %s %s %+.3f (text: no more than about 0.4 m)  %s\n",
            oth$product, oth$par, oth$median, ifelse(abs(oth$median) < 0.45, "PASS", "FAIL")))
cat(sprintf("RH98 SD = %.2f m (text: 8 m)\n", zc[variable == "rh_98"]$sd))

# assemble
rows <- rbind(vars[, .(par, group, name)], GROUPS)
rows[, ord := .I]
W <- dcast(fx, par ~ product, value.var = c("median", "lo", "hi"))
T5 <- merge(rows, W, by = "par")[order(ord)]
fwrite(T5, file.path(OUT, "tableS5_coefficients.csv"))

f2 <- function(x) { s <- sprintf("%.2f", x); s[s == "-0.00"] <- "0.00"; sub("-", "\u2212", s, fixed = TRUE) }
cell <- function(m, l, h) sprintf("%s (%s to %s)", f2(m), f2(l), f2(h))
md <- c("| Group | Variable | DSM | DTM | CHM |", "|---|---|---|---|---|")
prev <- ""
for (i in seq_len(nrow(T5))) {
  r <- T5[i]; g <- ifelse(r$group == prev, "", r$group); prev <- r$group
  md <- c(md, sprintf("| %s | %s | %s | %s | %s |", g, r$name,
                      cell(r$median_DSM, r$lo_DSM, r$hi_DSM), cell(r$median_DTM, r$lo_DTM, r$hi_DTM),
                      cell(r$median_CHM, r$lo_CHM, r$hi_CHM)))
}
writeLines(md, file.path(OUT, "tableS5_coefficients.md"))
cat(md, sep = "\n")
