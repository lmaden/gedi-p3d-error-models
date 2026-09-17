# Manifest: which script produces which artifact

Paths are repository-relative. Scripts marked as running on the University of Maryland cluster
read their project root from `PROJECT_ROOT`. Site
numbers inside scripts and per-site CSVs are TRACKER numbers; the manuscript numbering is given by
`data_derived/site_id_lookup.csv` (tracker 1 to 15 = manuscript 1 to 15; tracker 17 to 20 =
manuscript 16 to 19; tracker 16 and 21 are not in the study). The helper `tracker_to_manuscript()`
in `05_figures_tables/fig_common.R` applies the mapping.

The published CHM fit is the 16-site fit stored as `fit_chm_16` inside `checkpoints/10b_chm_sensitivity.rds`
(produced by `02_model_fitting/section_10b_chm_sensitivity_refit.R`). The published DTM fit is
`fit_dtm_s2` inside `checkpoints/10_models_stage2.rds` (produced by `02_model_fitting/section_10_models_stage2.R`).

## Main text

| Artifact | Script(s) | Notes |
|---|---|---|
| Table 1 | `01_data_prep/r_workflow/extract_section2_1_study_sites.R` | outputs are tracker-indexed |
| Table 2 | `01_data_prep/r_workflow/extract_table3_predictor_stats.R` | the script name carries the table's old number; the published run used the pre-filter footprint set |
| Table 3 | not on disk (see Known gaps in README); same basis as Table S8; `05_figures_tables/figS15_chm_dtm_diagnostics.R` reproduces the n and the descriptive statistics printed in its console log | |
| Table 4 | `02_model_fitting/section_10b_full_comparison.R` (CHM R², ν, holdout metrics), `03_evaluation/compute_definitive_r2.R` (DTM row only), `03_evaluation/verify_predictive_accuracy.R`, `03_evaluation/holdout_02_scoring.R`, `03_evaluation/holdout_03_cluster_values.R`, `03_evaluation/pipeline_preflight.R` | composite table; cells transcribed from the listed outputs |
| Table 5 | `02_model_fitting/section_10b_chm_sensitivity_refit.R` (CHM column, `est_16`), `05_figures_tables/fig3_coefficients.R` (DTM column) | |
| Table 6 | `02_model_fitting/section_10b_full_comparison.R` (CHM), `03_evaluation/forest_type_accuracy_holdout_full.R` (DTM), `03_evaluation/pipeline_preflight.R` (CRPS) | |
| Table 7 | `03_evaluation/leave_one_site_out/loo_table_build.R` pooling `loo_fold_chm.R` and `loo_fold_dtm.R` outputs | run log in `docs/run_logs/loo_table_build.txt` |
| Fig. 1 | `05_figures_tables/fig1_sitemap.R` via `fig1_sitemap_run_local.R` | |
| Fig. 2 | `05_figures_tables/fig2_workflow.svg` (hand-drawn, rendered with rsvg-convert) | |
| Fig. 3 | `05_figures_tables/fig3_coefficients.R` | |
| Fig. 4 | `05_figures_tables/fig4_coupling_15site.R` (panels a to c from `fig7_within_site_coupling.R`, panel d from `figure6_panel_d.R`) | the Section 3.4 cross-site correlation and its jackknife are reproduced by `03_evaluation/coupling_geoaccuracy_jackknife.R` |
| Fig. 5 | `05_figures_tables/fig5_pred_vs_obs.R` | |
| Fig. 6 | `05_figures_tables/fig_persite_prediction_intervals.R` | upstream `03_evaluation/per_site_interval_widths.R` |
| Fig. 7 | `05_figures_tables/fig5_ppc_refined_p3d.R` | |

## Supplement

| Artifact | Script(s) |
|---|---|
| Figs. S1, S2 | `05_figures_tables/publication_figures.R` (conditional effects) and `supp_figs_plotting.R` |
| Fig. S3 | `05_figures_tables/s3_vif_corr_probdir_16site.R` (exports) and `render_figS3_from_csv.R` |
| Figs. S4, S6, S12 | `05_figures_tables/supp_figs_rerender.R`, `supp_figs_rerender_fixes.R` |
| Fig. S5 | `05_figures_tables/publication_figures.R` |
| Figs. S7, S8 | `03_evaluation/temporal_mismatch_analysis.R`, `figS13_variogram_values.R` |
| Fig. S9 | `05_figures_tables/figS15_chm_dtm_diagnostics.R` |
| Fig. S10 | `01_data_prep/r_workflow/section_03_eda.R` |
| Fig. S11 | `05_figures_tables/figS11_site_random_effects.R` |
| Fig. S13 | `05_figures_tables/figS20_dtm_spherical.R`, data from `03_evaluation/holdout_02_scoring.R` and `holdout_04_dtm_panel_data.R` |
| Fig. S14 | `05_figures_tables/supp_figs_plotting.R` |
| Fig. S15 | `05_figures_tables/fig6_spatial_residuals.R` |
| Fig. S16 | `05_figures_tables/suppfig_S16_15panel_hexbin.R` |
| Figs. S17, S18 | `05_figures_tables/suppfig_25_site5_cover_decile.R`, `suppfig_26_site5_lc_strata.R`; data from `04_coupling_analysis/site5_cover/site5_03_cover_decile_detail.R` |
| Fig. S19 | `05_figures_tables/fig3_nu_posteriors.R` |
| Tables S1, S4, S5 | `02_model_fitting/section_10b_chm_sensitivity_refit.R` (CHM), `fixef()` of the DTM fit (DTM) |
| Table S2 | `03_evaluation/holdout_01_posterior_values.R`, `holdout_04_dtm_panel_data.R`; CHM column patch `03_evaluation/site_random_effect_correlations.R` |
| Table S3 | prior block of `02_model_fitting/section_10_models_stage2.R`; verified by `03_evaluation/verify_priors.R` |
| Tables S6a, S6b | `03_evaluation/figS9_16site_group_values.R`, `holdout_01_posterior_values.R` |
| Table S7 | `03_evaluation/section_12_spatial.R` |
| Table S8 | same as Table 3 |
| Table S9 | `03_evaluation/holdout_02_scoring.R` (CHM), `section_11b_variance_decomposition.R` (DTM) |
| Table S10 | `01_data_prep/site_joins_qc/extract_als_dates.py`; point-cloud inventory `04_coupling_analysis/reference_quality/als_pointcloud_inventory.sh`, `als_lasinfo_verify.sh` |
| Table S11 | `03_evaluation/baselines/baseline_01_ladder.R`, `baseline_03_conformal_intervals.R`, `baseline_02_restricted.R`; run log in `docs/run_logs/baseline_conformal_intervals.log` |
| Text S2 (reference quality) | `04_coupling_analysis/reference_quality/gedi_als_offset_distributions.R`, `gedi_als_offset_spatial_structure.R`; offsets recorded in the header of `section_10b_chm_sensitivity_refit.R` |
| Text S2.1 | `04_coupling_analysis/reference_quality/als_vs_als_noise_floor.R` |
| Text S2.2 | `02_model_fitting/counterfactual_chm_refit.R`; audit `03_evaluation/leave_one_site_out/loo_decomposition_audit.R` |
| Text S2.3 | `04_coupling_analysis/site5_cover/*.R` |
| Text S2.4 | prior predictive block of `02_model_fitting/section_08_model_prep.R` |
| Text S2.5 | `03_evaluation/leave_one_site_out/*`, `02_model_fitting/confirmation_fits/*`, `02_model_fitting/tile_sensitivity_pilot/tile_pilot_chm.R`; run logs in `docs/run_logs/` |

## Pipeline order

1. `01_data_prep/3dep_and_reference_rasters/` (download, derivatives, resampling), then `01_data_prep/footprint_extraction/gedi_sample_chunked.py` with `helper_speed.py` (footprint means within a 12.5 m radius, `als_chm_p90`, ecoregion join).
2. `01_data_prep/p3d_metadata/` (acquisition-geometry predictors from tile metadata) and `01_data_prep/land_cover/` (GLC-FCS30 class per footprint), joined per site in `01_data_prep/site_joins_qc/`.
3. `01_data_prep/r_workflow/` sections 00 to 07 (ingest, filters, z-scaling, exploratory analysis, collinearity, summaries), configured by `analysis_config.R` and `analysis_utils.R`.
4. `02_model_fitting/section_08_model_prep.R` (seed 2025; 20% Stage 1 and 33% Stage 2 site-stratified samples; holdout is the complement), `section_09_models_stage1.R`, `section_10_models_stage2.R` (DTM fit), `section_10b_chm_sensitivity_refit.R` (16-site CHM fit).
5. `03_evaluation/` (holdout scoring, PPCs, convergence, spatial residuals, leave-one-site-out folds and baselines), `04_coupling_analysis/`, then `05_figures_tables/`.
