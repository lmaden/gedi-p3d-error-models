# Analysis code and derived data for "Bayesian hierarchical characterization of satellite photogrammetric elevation errors in forests using GEDI"

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22831703.svg)](https://doi.org/10.5281/zenodo.22831703)

Levi Madenberg, Paul B. May, Hao Tang, Adrian Pascual, Ralph Dubayah.
Submitted to *ISPRS Journal of Photogrammetry and Remote Sensing*, September 2026.
Corresponding author: Levi Madenberg (lmadenbe@umd.edu), Department of Geographical Sciences, University of Maryland.

This repository holds the R and Python code behind the tables and figures in the paper and its
supplement, together with the derived site-level and forest-type-level tables the paper reports.
It is organised so that a reader can see what was done and could redo it given the input data
described below. It is not a one-command pipeline: the model fits took days on a university
cluster and several inputs are proprietary or must be downloaded from their providers.
`docs/MANIFEST.md` maps every
table and figure to the script that produced it.

## What is here

| Folder | Contents |
|---|---|
| `01_data_prep/` | GEDI footprint quality checks, 3DEP download and resampling, footprint-mean extraction of the P3D and reference rasters, P3D tile-metadata aggregation, land-cover and ecoregion assignment, per-site joins and quality control, and the R ingest workflow (sections 00 to 07) |
| `02_model_fitting/` | Stage 1 and Stage 2 `brms` fits for the CHM (16 sites) and DTM (18 sites) error models, the counterfactual reference-substitution refit, the tuned confirmation fits, the tile-level sensitivity pilot, and the cluster environment setup |
| `03_evaluation/` | Holdout scoring, posterior predictive checks, convergence diagnostics, spatial residual analysis, the temporal-mismatch analysis, the leave-one-site-out folds and their pooled table |
| `04_coupling_analysis/` | Within-site CHM-DTM error coupling and its cross-site correlates, the Site 5 cover analysis, the ALS-versus-ALS reference noise floor and the reference-quality screening |
| `05_figures_tables/` | Scripts that render the main and supplementary figures, the Fig. 2 SVG source, and shared plotting helpers |
| `data_derived/` | Machine-readable copies of Tables 1 to 7 and S1 to S11 exactly as printed, the site-level coupling results, and the tracker-to-manuscript site-ID lookup |
| `docs/` | `MANIFEST.md` mapping each table and figure to its script, and run logs of the leave-one-site-out table build, the baseline ladder, the tile pilot and the fold gates |

## What is not here

- **Vantor Precision3D imagery and the delivered DSM, DTM and CHM rasters.** They are proprietary,
  were obtained under a data sharing agreement, and cannot be redistributed. No per-footprint
  P3D value is included in this repository.
- **Fitted model objects.** Stan fits are large and version-fragile. The posterior summaries the
  paper reports are in `data_derived/`.
- **Public input data**, to be obtained from the providers: GEDI L2A, L2B, L4A and L4C
  (NASA LP DAAC), 3DEP lidar-derived elevation (USGS), NEON airborne lidar canopy height models
  (NEON Data Portal), GLC-FCS30 land cover (Zhang et al., 2021), and the CEC Level II ecoregions.
- **The per-site GEDI footprint files** (`GEDI_site<N>_hq_ALL.gpkg`) were built upstream in the
  UMD GEDI lab archive with the quality filters stated in Section 2.2.3 of the paper; the builder is
  not part of this repository.

## Software environment

R 4.5.0 with `brms` 2.23.0, `posterior` 1.6.1, `cmdstanr` on CmdStan 2.37.0, `scoringRules` 1.1.3,
`ranger` 0.18.0 and the tidyverse; Python 3 with `rasterio`, `numpy`, `pandas` and `h5py` for the
raster and GEDI extraction steps. Sampler settings for the published fits: four chains of 3,000
iterations (1,500 warmup), no thinning, `adapt_delta = 0.98`, `max_treedepth = 15`; the CHM fit
took about 58 h and the DTM fit about 161 h on the cluster; the 34 leave-one-site-out folds ran as
a work queue over several days.

Cluster scripts read the project root from the environment variable `PROJECT_ROOT` (or
`CHPT1_ROOT`) and fall back to the original cluster path. The Python extraction scripts in
`01_data_prep/` carry the original local Windows and external-drive paths at the top of each file
and must be edited before use. No script contains credentials; the one API token reference is a
placeholder.

## Reproducibility notes and known gaps

The data-preparation scripts under `01_data_prep/` name their input directories with the placeholder
`<LOCAL_DATA_ROOT>`; substitute the directory that holds the inputs described above.

- Seeds: `set.seed(2025)` in `02_model_fitting/section_08_model_prep.R` draws the 20% Stage 1 and
  33% Stage 2 site-stratified samples; the holdout is the complement. The two samples were drawn
  independently and overlap by 19.9%. Downstream scripts regenerate the holdout from the seed; one
  re-render in August 2026 drew a slightly different holdout under a newer dplyr, so the published
  Fig. 5 image is the original render.
- Site identifiers inside scripts and per-site files are tracker numbers; see
  `data_derived/site_id_lookup.csv`.
- The 16-site CHM fit is the only CHM fit the paper reports and is stored under the file name
  `10b_chm_sensitivity.rds` (element `fit_chm_16`); the 18-site DTM fit is `fit_dtm_s2` in
  `10_models_stage2.rds`.
- One artifact has no located producing script: the descriptive statistics of Table 3 and Table S8,
  a May 2026 recompute on the forest subset. `05_figures_tables/figS15_chm_dtm_diagnostics.R`
  applies the same filters and prints the same statistics.

## License and citation

Code is released under the MIT License (`LICENSE`). Derived tables in `data_derived/` are released
under CC BY 4.0 (`data_derived/LICENSE`). Please cite the paper (see `CITATION.cff`). Releases are
archived at Zenodo under the concept DOI https://doi.org/10.5281/zenodo.22831703.

## Funding

NASA GEDI mission, contract NNL15AA03C to the University of Maryland.
