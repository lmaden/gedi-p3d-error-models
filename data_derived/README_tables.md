# Derived tables

Errors are P3D minus lidar, in meters. Site identifiers are internal ids. `site_names.csv` gives the
matching site number in the paper (`paper_site`) and the site name for each internal id (`tracker`).
Columns named `site`, `site_id` or `tracker` hold the internal id. None of these files holds a
per-footprint P3D value.

## Inputs read by the scripts

| File | Contents | Made by |
|---|---|---|
| `site_names.csv` | internal id, paper site number, short and formal names, state, lidar flight | lookup |
| `ecoregion_by_site.csv` | footprints per site and CEC Level II ecoregion | `05_tables_figures/table1_facts.R` |
| `lidar_dates.csv` | acquisition dates from the GPS times of the points and the file headers, ten files per flight | `05_tables_figures/table1_facts.R` |
| `p3d_year_summary.csv` | first and last image year per site and P3D product, with post spacing and datums | image metadata of the footprint sampling (`01_data_prep/p3d_metadata/`) |
| `zconst.csv` | mean and SD used to standardize each tested variable | `03_analysis_frame/rework_model_prep.R` |
| `full_fits/{dsm,dtm,chm}/fixef.csv`, `params.csv` | posterior summaries (median, mean, SD, 95% interval, R-hat, bulk ESS) of the population-level effects, group SDs and correlations and the other distributional parameters (`fixef.csv`), and of all parameters including the group offsets (`params.csv`) | `04_error_models/rework_fullfit.R` |
| `full_fits/{dsm,dtm,chm}/r2.csv` | Bayesian R2 with and without the group offsets | `04_error_models/rework_fullfit.R` |
| `full_fits/{dsm,dtm,chm}/priors.csv` | priors of each model | `04_error_models/rework_fullfit.R` |
| `folds/metrics_{dsm,dtm,chm}_fold{1..6}.csv` | held-out scores per fold and site (RMSE, bias, interval coverage and width), with sampler settings and diagnostics; the DSM rows use the DSM error of Eq. 4 | `04_error_models/rework_fold.R` |
| `ppc/bins_{dsm,dtm,chm}.csv`, `ppc/stats_{dsm,dtm,chm}.csv` | binned error densities and per-draw summaries of the observed and replicated errors | `04_error_models/rework_ppc.R` |

## Outputs as published (`results/`)

These are the outputs of `05_tables_figures/` used in the paper, at the paths the scripts write to
under `P3D_RESULTS`. The exception is `fig5/fig5_transects_selected.csv`, which was assembled from the
outputs of `fig5_select_transects.R` and `fig5_profiles.R`.

| File | Contents |
|---|---|
| `table1/table1_sites.{csv,md}` | Table 1 |
| `table2/table2_site_weighted.csv` | Table 2: bias, SD, RMSE, 95th percentile of the absolute error and shares over 5 m and 10 m, with sites weighted equally, with all footprints weighted equally, and for the median site (not printed) |
| `table2/fig2_panel_stats.csv` | bias, RMSE and correlation of the two panels of Fig. 2 |
| `table2/error_tails.csv` | shares of errors beyond three SD and share of the error variance in the 2% of footprints with the largest errors |
| `table3/table3_transfer.{csv,md}`, `table3/table3_panels.md` | Table 3: variance shares and corrections at held-out sites |
| `tableS1/tableS1_predictors.{csv,md}` | Table S1 |
| `tableS2/tableS2_transfer_groups.{csv,md}` | Table S2: corrections at each held-out site |
| `tableS3/*` | Table S3 and its two parts |
| `tableS4/*` | Table S4: site rows, forest-type rows (sites weighted equally) and pooled forest-type rows |
| `tableS5/*` | Table S5: posterior medians and 95% intervals of the effects |
| `persite/persite.csv` | bias, SD and RMSE of each product per site, share of CHM errors over 5 m, median lidar canopy height, DSM and DTM error correlation |
| `persite/height_class.csv`, `persite/forest_class.csv`, `persite/accuracy_summary.csv` | pooled error statistics by lidar height class, by forest class and overall |
| `persite/lp_checks.csv` | per-site checks of the lidar references: ground surface minus 3DEP, rebuilt against delivered lidar CHM, lidar against GEDI canopy height |
| `fig2/fig2_scatter_stats.csv` | Fig. 2 statistics and the footprints outside the plotted range |
| `fig3/fig3b_class_means.csv` | Fig. 3b: mean error and 10th to 90th percentiles of the site means by height class |
| `fig5/fig5_transects_selected.csv` | the three transects of Fig. 5, with the first and last GEDI shot number, the number of footprints and the mean errors printed in the panel subtitles (DSM errors from the delivered DSM sampled at each footprint) |
| `figS1/figS1_numbers.csv` | observed and replicated shares of large errors and tail quantiles (Fig. S1, Text S2) |
| `figS2/variogram_persite.csv`, `figS2/variogram_summary.csv` | normalized semivariance of the errors by distance class, per site and as medians over sites (Fig. S2, Text S3) |
