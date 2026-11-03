# Tables and figures: inputs, run order, outputs

The scripts in this folder build every table and figure of the paper from the outputs of the cluster
pipeline (folders 02 to 04 and two cluster scripts here) and from the derived files shipped in
`data_derived/`. Inside every script and file, sites are internal ids (1 to 9, 11 to 15, 17 to 20).
`data_derived/site_names.csv` maps them (column `tracker`) to the paper's site numbers 1 to 18.

## Roots

Every script reads and writes through three roots, set by environment variables:

| Variable | Default | Holds |
|---|---|---|
| `P3D_DATA_ROOT` | `<LOCAL_DATA_ROOT>` (must be set) | footprint-level inputs, not distributed: a local copy of the cluster folder `$CH1/manuscript_tables` plus `external/` |
| `P3D_DERIVED` | `data_derived` | inputs distributed with this release |
| `P3D_RESULTS` | `results` | outputs, one folder per paper item |

Relative defaults are resolved from the working directory, so run the scripts from the repository root.
Example (bash):

```
export P3D_DATA_ROOT=/path/to/manuscript_tables P3D_DERIVED=data_derived P3D_RESULTS=results
Rscript 05_tables_figures/persite_class_tables.R
python 05_tables_figures/site_weighted_accuracy.py
```

Checked with R 4.4.0 (data.table 1.15.4, ggplot2 3.5.2, patchwork 1.2.0, sf 1.0.16, ggrepel 0.9.5,
maps 3.4.3, hexbin 1.28.3) and Python 3.12 (pandas 3.0.5, numpy 2.5.2, scipy 1.18.0), which reproduce the
published tables and figures exactly. Set `PYTHONIOENCODING=utf-8` on Windows consoles.

### `P3D_DATA_ROOT` (not distributed)

| File | Produced by |
|---|---|
| `rework/frame_footprints.csv.gz` | `03_analysis_frame/rework_model_prep.R`, then `03_analysis_frame/merge_dsm_direct.R` (directly sampled DSM error in `dsm_error_mean`, Eq. 4 value in `dsm_error_eq4`) |
| `rework/folds/rowpreds_{dsm,dtm,chm}_fold{1..6}.csv.gz` | `04_error_models/rework_fold.R` |
| `lidar_products/lp_site_NN.csv.gz`, `lidar_products/lp_summary.csv` | `02_lidar_references/lidar_products_rebuild.R` |
| `dsm_direct/dsm_direct_frame.csv.gz` | `03_analysis_frame/p3d_dsm_direct.R`, stacked by `merge_dsm_direct.R` |
| `lidar_ground/lg_site_06_neon_jerc2019.csv.gz`, `lg_site_04.csv.gz`, `lg_site_14.csv.gz` | an earlier ground-only lidar rebuild (not in this release); used only for the footprint coordinates of the Fig. 5 window selection |
| `fig5/fig5_profile_site{04,06,14}.csv`, `fig5/fig5_fpcheck_site{04,06,14}.csv` | `fig5_profiles.R` on the cluster (`$CH1/manuscript_tables/fig5/`) |
| `external/section2_1_site_summary_full_v2.csv` | site centroids from `01_data_prep/r_workflow/extract_section2_1_study_sites.R` (Fig. 1 markers) |
| `external/na_cec_eco_l2/NA_CEC_Eco_Level2.shp` (with .dbf, .shx, .prj) | CEC North American Level II ecoregions (Commission for Environmental Cooperation) |

### `P3D_DERIVED` (`data_derived/`)

| File | Produced by |
|---|---|
| `site_names.csv` | tracker id, paper site number, names, lidar flight |
| `ecoregion_by_site.csv`, `lidar_dates.csv` | `table1_facts.R` (cluster) |
| `p3d_year_summary.csv` | first and last image year per site and P3D product (image metadata of the sampling stage) |
| `zconst.csv` | `03_analysis_frame/rework_model_prep.R` (standardization constants) |
| `full_fits/{dsm,dtm,chm}/{fixef,params,r2,priors}.csv` | `04_error_models/rework_fullfit.R` |
| `folds/metrics_{dsm,dtm,chm}_fold{1..6}.csv` | `04_error_models/rework_fold.R` |
| `ppc/{bins,stats}_{dsm,dtm,chm}.csv` | `04_error_models/rework_ppc.R` |

## Upstream (cluster) order

The cluster scripts read the project root, written `$CH1` in comments and in this note, from an
environment variable. Most scripts use `PROJECT_ROOT`, `lidar_products_rebuild.R` and `p3d_dsm_direct.R`
use `CH1`, and `fig5_profiles.R` uses `CH1_ROOT`. `p3d_dsm_direct.R` and `fig5_profiles.R` also read a
second root, where the P3D rasters and the lidar tiles were kept (`CH2` and `CH2_ROOT`). Paths specific
to our cluster are written as placeholders: `<CLUSTER_ROOT>` for the project storage, `<R_LIBRARY>` for the
R package library and `<LIDAR_ARCHIVE>` for the folder of lidar point clouds. Before a run, set `IMPORTED`
in `lidar_products_rebuild.R` and `table1_facts.R` to the lidar folder, and `ALLOWED` in `rework_fullfit.R`
and `rework_fold.R` to the storage root. `ALLOWED` is a guard that stops a fit whose temporary files would
land outside that root. The module names in `run_rework_queue.sh` are those of our cluster.

1. `02_lidar_references/lidar_products_rebuild.R` (lidar ground and pit-free CHM at every footprint)
2. `03_analysis_frame/rework_model_prep.R` (frame, datum steps, standardization, training sample and folds). It needs the
   NOAA GEOID18 grid `us_noaa_g2018u0.tif`, from the PROJ data files, in `$CH1/data/geoid/`.
3. `04_error_models/run_rework_queue.sh` (3 full fits with `rework_fullfit.R`, 18 held-out fits with `rework_fold.R`)
4. `04_error_models/rework_divcheck.R`, `04_error_models/rework_ppc.R` (read the full fits)
5. `03_analysis_frame/p3d_dsm_direct.R` (the delivered DSM at every frame footprint, read with the frame of step 2)
6. `05_tables_figures/table1_facts.R`
7. Copy `$CH1/manuscript_tables` to `P3D_DATA_ROOT`, then run `03_analysis_frame/merge_dsm_direct.R` once on that copy.

## Run order in this folder

| Step | Script | Paper item | Reads | Writes (under `P3D_RESULTS`) |
|---|---|---|---|---|
| 1 | `persite_class_tables.R` | per-site and class tables behind Table 1, Table S4 and Fig. 4 | frame, `lidar_products/lp_summary.csv` | `persite/persite.csv`, `height_class.csv`, `forest_class.csv`, `accuracy_summary.csv`, `lp_checks.csv`, `error_densities.{png,pdf}`, log |
| 2 | `table1_sites.R` | Table 1 | `site_names.csv`, `ecoregion_by_site.csv`, `p3d_year_summary.csv`, `persite/persite.csv`, frame | `table1/table1_sites.{csv,md}` |
| 3 | `site_weighted_accuracy.py` | Table 2, Fig. 2 panel statistics, tail statistics | frame, `lidar_products/lp_site_NN.csv.gz` | `table2/table2_site_weighted.csv`, `fig2_panel_stats.csv`, `error_tails.csv` |
| 4 | `table3_transfer.R` (sources `dsm_direct_rowpreds.R`) | Table 3 | `rework/folds/rowpreds_*`, `dsm_direct/dsm_direct_frame.csv.gz`, `folds/metrics_*`, `full_fits/*/r2.csv` | `table3/table3_transfer.{csv,md}` |
| 5 | `table3_panels.R` | Table 3 as printed | `table3/table3_transfer.csv` | `table3/table3_panels.md` |
| 6 | `tableS1_predictors.R` | Table S1 | `zconst.csv` | `tableS1/tableS1_predictors.{csv,md}` |
| 7 | `tableS2_transfer_groups.R` (sources `dsm_direct_rowpreds.R`) | Table S2 | `rework/folds/rowpreds_*`, `dsm_direct/dsm_direct_frame.csv.gz`, `site_names.csv`, `table3/table3_transfer.csv` | `tableS2/tableS2_transfer_groups.{csv,md}` |
| 8 | `fig3_error_distributions.R` | Fig. 3 | frame, `table2/error_tails.csv` | `fig3/fig3_error_distributions.{png,pdf}`, `fig3b_class_means.csv` |
| 9 | `tableS3_height_classes.R` | Table S3 class means | frame, `fig3/fig3b_class_means.csv` | `tableS3/tableS3_height_classes.{csv,md}` |
| 10 | `tableS3_rmse_over5.R` | Table S3a RMSE and share > 5 m | frame, `persite/height_class.csv`, `tableS3/tableS3_height_classes.csv` | `tableS3/tableS3_rmse_over5.csv` |
| 11 | `tableS3_assemble.py` | Table S3 as printed | `tableS3/tableS3_height_classes.md`, `tableS3/tableS3_rmse_over5.csv` | `tableS3/tableS3.md` |
| 12 | `tableS4_persite.R` | Table S4 site rows | frame, `site_names.csv`, `persite/persite.csv` | `tableS4/tableS4_persite.{csv,md}` |
| 13 | `tableS4_forest_rows.py` | Table S4 forest-type rows | frame, `persite/forest_class.csv` | `tableS4/tableS4_forest_rows.{csv,md}`, `forest_rows_pooled.md` |
| 14 | `tableS4_assemble.py` | Table S4 as printed | `tableS4/tableS4_persite.md`, `tableS4/tableS4_forest_rows.md`, `table1/table1_sites.csv`, frame | `tableS4/tableS4.md` |
| 15 | `tableS5_coefficients.R` | Table S5 | `full_fits/*/fixef.csv`, `zconst.csv`, `tableS1/tableS1_predictors.csv` | `tableS5/tableS5_coefficients.{csv,md}` |
| 16 | `tableS5_assemble.py` | Table S5 as printed | `tableS5/tableS5_coefficients.md`, `tableS1/tableS1_predictors.md` | `tableS5/tableS5.md` |
| 17 | `fig1_site_map.R` | Fig. 1 | `site_names.csv`, `external/*`, `table1/table1_sites.csv` | `fig1/fig1_site_map.{png,pdf}` |
| 18 | `fig2_height_scatter.R` | Fig. 2 | frame, `lidar_products/lp_site_NN.csv.gz` | `fig2/fig2_height_scatter.{png,pdf}`, `fig2_scatter_stats.csv` |
| 19 | `fig4_site_accuracy.R` | Fig. 4 | `persite/persite.csv`, `site_names.csv` | `fig4/fig4_site_accuracy.{png,pdf}` |
| 20 | `fig5_select_transects.R` | Fig. 5, transect selection | frame, `lidar_ground/lg_site_*` | `fig5/fig5_windows.csv`, `fig5_footprints_siteNN.csv`, `fig5_footprints.csv` |
| 21 | `fig5_profiles.R` (cluster) | Fig. 5, profiles | `fig5_footprints.csv` in the script folder, lidar and P3D rasters | `$CH1/manuscript_tables/fig5/` (copy to `P3D_DATA_ROOT/fig5/`) |
| 22 | `fig5_transects_plot.R` | Fig. 5 | `fig5/fig5_profile_*`, `fig5/fig5_fpcheck_*` | `fig5/fig5_transects.{png,pdf}` |
| 23 | `figS1_ppc.R` | Fig. S1 | `ppc/*` | `figS1/figS1_ppc.{png,pdf}`, `figS1_numbers.csv` |
| 24 | `within_site_variogram.py` | Fig. S2 and Text S3 numbers | frame, `lidar_products/lp_site_NN.csv.gz` | `figS2/variogram_persite.csv`, `variogram_summary.csv` |
| 25 | `figS2_variogram.R` | Fig. S2 | `figS2/variogram_summary.csv` | `figS2/figS2_variogram.{png,pdf}` |

"frame" is `P3D_DATA_ROOT/rework/frame_footprints.csv.gz`. Steps 1 to 19 and 23 to 25 need nothing from
Fig. 5. Step 22 needs the cluster output of step 21 (the published profiles were cut from the windows that
step 20 reproduces).

## Notes

- Table 2 is `table2/table2_site_weighted.csv`, rows "Sites weighted equally" and "All footprints" (printed as
  "Footprints weighted equally"), rounded to two decimals (percentages to one). The "Median site" rows are not printed.
- The DSM error models were fit to the DSM error of Eq. 4 (CHM error plus DTM error), and Table 3a, Table S5 and
  Fig. S1 use it. Table 2, Table 3b, Tables S2 to S4, Figs. 2 to 4 and Fig. S2 use the directly sampled
  DSM error (Eq. 1). The Fig. 5 transects were selected with the Eq. 4 DSM error (see `fig5_select_transects.R`). The
  errors printed in Fig. 5 are recomputed on the cluster from the rasters.
- Table S1, Table S5 and Fig. S1 can be rebuilt from `data_derived/` alone. The Table 3a shares are the medians in
  `data_derived/full_fits/*/r2.csv`. For panel (b), `table3_transfer.R` also needs the held-out predictions.
- `merge_dsm_direct.R` reproduces the frame and `dsm_direct_frame.csv.gz` used for the paper to within 1e-16 m in every
  value. The files used for the paper were written by a merge that read the per-site files with a float parser that is
  not exactly round-trip, so a few trailing digits differ. No reported number changes.
