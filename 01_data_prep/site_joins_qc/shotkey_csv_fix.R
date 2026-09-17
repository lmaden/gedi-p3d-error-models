library(data.table)

DIR <- "/gpfs/data1/vclgp/lmaden/chpt1/data/enriched_by_site"
fs <- list.files(DIR, pattern="^site_[0-9]{2}_enriched\\.csv(\\.gz)?$", full.names = TRUE)

for (f in fs) {
  dt <- fread(f, colClasses = c(shot_key="character", shot_number="character"))
  if ("shot_key" %in% names(dt)) {
    # If shot_number is all NA or missing, backfill from shot_key
    if (!"shot_number" %in% names(dt) || all(is.na(dt$shot_number))) {
      dt[, shot_number := shot_key]
    }
    dt[, shot_key := NULL]  # drop the helper
    # rewrite (same name)
    fwrite(dt, f, sep = ",", quote = FALSE, na = "", bom = FALSE)
    cat("fixed:", f, "rows=", nrow(dt), "\n")
  }
}
