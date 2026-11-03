# dsm_direct_rowpreds.R  --  held-out DSM rows with the directly sampled DSM error
# Sourced by table3_transfer.R and tableS2_transfer_groups.R (Table 3b, Table S2).
# The DSM error models were fit to the DSM error from Eq. 4 (CHM error plus DTM error). The paper samples the
# delivered DSM directly (p3d_dsm_direct.R; merge_dsm_direct.R writes dsm_direct/dsm_direct_frame.csv.gz), so
# every held-out DSM quantity that uses the observed error is recomputed with the direct error (observed + delta):
# no correction, the single average correction (mean error of the training footprints at the other 15 sites, which
# are the held-out rows of the other groups), the model's correction (direct error minus the unchanged prediction),
# the site's own mean, and interval coverage. Predictions and interval limits are unchanged.
# Input  (P3D_DATA_ROOT) dsm_direct/dsm_direct_frame.csv.gz (columns site, shot_number, delta)
# Check  the constant corrections recomputed from the Eq. 4 errors reproduce the stored ones (1e-9) before the swap.
use_direct_dsm <- function(R, foldcol, root) {
  dd <- fread(file.path(root, "dsm_direct", "dsm_direct_frame.csv.gz"),
              colClasses = list(character = "shot_number"))[, .(site = as.character(site), shot_number, delta)]
  ds <- copy(R[product == "DSM"]); setnames(ds, foldcol, "g_")
  site_m <- ds[, .(m = mean(observed)), by = .(g_, site)]
  e_site <- merge(ds[, .(g_, site, pred_const_site)], site_m, by = c("g_", "site"))[, max(abs(pred_const_site - m))]
  tot <- ds[, .(s = sum(observed), n = .N), by = g_]
  tot[, m_other := (sum(s) - s) / (sum(n) - n)]
  e_glob <- merge(ds[, .(g_, pred_const_global)], tot[, .(g_, m_other)], by = "g_")[, max(abs(pred_const_global - m_other))]
  cat(sprintf("GATE DSM constant corrections recomputed from the Eq. 4 errors: site mean %.1e, other-sites mean %.1e  %s\n",
              e_site, e_glob, ifelse(max(e_site, e_glob) < 1e-9, "PASS", "FAIL")))
  stopifnot(max(e_site, e_glob) < 1e-9)
  n0 <- nrow(ds)
  ds <- merge(ds, dd, by = c("site", "shot_number"), all.x = TRUE)
  stopifnot(nrow(ds) == n0, !anyNA(ds$delta))
  ds[, observed := observed + delta][, delta := NULL]
  ds[, pred_const_site := mean(observed), by = .(g_, site)]
  tot <- ds[, .(s = sum(observed), n = .N), by = g_]
  tot[, m_other := (sum(s) - s) / (sum(n) - n)]
  ds <- merge(ds, tot[, .(g_, m_other)], by = "g_")
  ds[, pred_const_global := m_other][, m_other := NULL]
  setnames(ds, "g_", foldcol)
  cat(sprintf("DSM held-out rows switched to the direct DSM error: %s rows, mean change %+.4f m\n",
              format(n0, big.mark = ","), mean(ds$observed) - mean(R[product == "DSM"]$observed)))
  rbind(R[product != "DSM"], ds[, names(R), with = FALSE], use.names = TRUE)
}
