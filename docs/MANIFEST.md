# Manifest: paper items, scripts and outputs

Script paths are relative to the repository root. Output paths are relative to the results root
(`P3D_RESULTS`). A copy of each published table output is in `data_derived/results/` under the same
relative path, so a rerun can be compared file by file. `05_tables_figures/RUN_ORDER.md` gives the run
order and the inputs of every script.

## Tables

| Paper item | Script | Output |
|---|---|---|
| Table 1 | `05_tables_figures/table1_sites.R` (site facts from `05_tables_figures/table1_facts.R`) | `table1/table1_sites.{csv,md}` |
| Table 2 | `05_tables_figures/site_weighted_accuracy.py` | `table2/table2_site_weighted.csv` |
| Table 3 | `05_tables_figures/table3_transfer.R`, then `table3_panels.R` | `table3/table3_transfer.{csv,md}`, `table3/table3_panels.md` |
| Table S1 | `05_tables_figures/tableS1_predictors.R` | `tableS1/tableS1_predictors.{csv,md}` |
| Table S2 | `05_tables_figures/tableS2_transfer_groups.R` | `tableS2/tableS2_transfer_groups.{csv,md}` |
| Table S3 | `05_tables_figures/tableS3_height_classes.R`, `tableS3_rmse_over5.R`, then `tableS3_assemble.py` | `tableS3/tableS3.md` (parts in `tableS3/`) |
| Table S4 | `05_tables_figures/tableS4_persite.R`, `tableS4_forest_rows.py`, then `tableS4_assemble.py` | `tableS4/tableS4.md` (parts in `tableS4/`) |
| Table S5 | `05_tables_figures/tableS5_coefficients.R`, then `tableS5_assemble.py` | `tableS5/tableS5.md` |

## Figures

| Paper item | Script | Output |
|---|---|---|
| Fig. 1 | `05_tables_figures/fig1_site_map.R` | `fig1/fig1_site_map.{png,pdf}` |
| Fig. 2 | `05_tables_figures/fig2_height_scatter.R` (panel statistics from `site_weighted_accuracy.py`) | `fig2/fig2_height_scatter.{png,pdf}`, `fig2/fig2_scatter_stats.csv` |
| Fig. 3 | `05_tables_figures/fig3_error_distributions.R` | `fig3/fig3_error_distributions.{png,pdf}`, `fig3/fig3b_class_means.csv` |
| Fig. 4 | `05_tables_figures/fig4_site_accuracy.R` (per-site values from `persite_class_tables.R`) | `fig4/fig4_site_accuracy.{png,pdf}` |
| Fig. 5 | `05_tables_figures/fig5_select_transects.R`, `fig5_profiles.R` (cluster), `fig5_transects_plot.R` | `fig5/fig5_transects.{png,pdf}`; the three transects are listed in `data_derived/results/fig5/fig5_transects_selected.csv` |
| Fig. S1 | `04_error_models/rework_ppc.R` (cluster), then `05_tables_figures/figS1_ppc.R` | `figS1/figS1_ppc.{png,pdf}`, `figS1/figS1_numbers.csv` |
| Fig. S2 | `05_tables_figures/within_site_variogram.py`, then `figS2_variogram.R` | `figS2/variogram_{persite,summary}.csv`, `figS2/figS2_variogram.{png,pdf}` |

## Numbers in the text

| Number | Where | Source |
|---|---|---|
| 630,390 footprints at 18 sites; 68,786 training footprints | Abstract, 2.3, 2.4, Text S2 | `03_analysis_frame/rework_model_prep.R`; fit logs in `docs/run_logs/` |
| Bias, RMSE and shares of errors over 5 m and 10 m, sites and footprints weighted equally | Abstract, 3.1, 4.1 | `table2/table2_site_weighted.csv` |
| RMSE 5.37 m to 3.63 m and r 0.80 to 0.90 with the lidar ground surface | Abstract, 3.1 | `table2/fig2_panel_stats.csv`, `fig2/fig2_scatter_stats.csv` |
| Errors beyond three SD 7 to 9 times the normal share; 2% of footprints hold 42%, 30% and 29% of the error variance | 3.1, 4.2 | `table2/error_tails.csv` |
| DTM error by canopy height class; share of CHM errors over 5 m by forest type | 3.1 | `tableS3/`, `tableS4/` |
| Error variance explained by the tested variables (5%, 10%, 4%) and by the group offsets | 3.2 | `data_derived/full_fits/*/r2.csv`, `table3/table3_transfer.csv` |
| Corrections at held-out sites | 3.3 | `table3/table3_transfer.csv`, `tableS2/tableS2_transfer_groups.csv` |
| Transect means at three sites | 3.1 | `data_derived/results/fig5/fig5_transects_selected.csv` |
| Within-site correlation 0.57 at 60 m, 0.25 at 300 m, near 0 beyond 1 km | 3.1 | `figS2/variogram_summary.csv` |
| Sampler settings and diagnostics | Text S2 | `docs/run_logs/` |
| Posterior predictive shares | Text S2, Fig. S1 | `figS1/figS1_numbers.csv` |

A few counts and shares in the text were computed directly from the footprint-level analysis frame,
which is not distributed. They are the share of the DSM errors below -10 m that fall at Coastal South
Carolina (77%), the median lidar canopy height of non-forest footprints (16 m), the share of footprints
where the two DSM errors agree within 0.1 m (97%), the footprints removed by the screening against 3DEP
(Text S1), and the footprint counts in the notes of Fig. 3 and Table S4.

## Upstream (cluster)

| Step | Script | Output used downstream |
|---|---|---|
| Footprint sampling of the P3D DTM and height model, 3DEP, GEDI variables, image metadata, land cover and ecoregions | `01_data_prep/` | per-site footprint tables |
| Lidar ground surface and pit-free CHM at every footprint | `02_lidar_references/lidar_products_rebuild.R` | `lidar_products/lp_site_NN.csv.gz`, `lp_summary.csv` |
| Analysis frame, datum steps, standardization, training sample and folds | `03_analysis_frame/rework_model_prep.R` | `rework/frame_footprints.csv.gz`, `data_derived/zconst.csv` |
| Full fits and held-out fits | `04_error_models/run_rework_queue.sh`, `rework_fullfit.R`, `rework_fold.R` | `data_derived/full_fits/`, `data_derived/folds/`, held-out predictions |
| Divergent-draw check and posterior predictive check | `04_error_models/rework_divcheck.R`, `rework_ppc.R` | `data_derived/ppc/` |
| Delivered DSM at every footprint | `03_analysis_frame/p3d_dsm_direct.R`, then `merge_dsm_direct.R` | `dsm_direct/dsm_direct_frame.csv.gz`, updated frame |
| Site facts for Table 1 | `05_tables_figures/table1_facts.R` | `data_derived/ecoregion_by_site.csv`, `data_derived/lidar_dates.csv` |
