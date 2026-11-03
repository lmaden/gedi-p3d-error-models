# Analysis code and derived data for "Assessing Precision3D satellite photogrammetric elevation models for forest canopy height retrieval"

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22831703.svg)](https://doi.org/10.5281/zenodo.22831703)

Levi Madenberg, Paul B. May, Hao Tang, Adrian Pascual and Ralph Dubayah. Submitted to the ISPRS Open
Journal of Photogrammetry and Remote Sensing, November 2026. The corresponding author is Levi
Madenberg (lmadenbe@umd.edu), Department of Geographical Sciences, University of Maryland.

The paper compares the Precision3D (P3D) digital surface model (DSM), digital terrain model (DTM) and
canopy height model (CHM) with airborne lidar at 630,390 GEDI footprints at 18 sites in the
conterminous United States. Bayesian hierarchical models relate the errors to environmental,
structural, terrain and image variables, and sites held out of the fits test whether the models
correct the errors elsewhere.

The code in this repository produced every table and figure in the paper and its supplement.
`data_derived/` holds the site-level, class-level and model-summary tables the paper reports, and
`docs/MANIFEST.md` gives the script behind each table, figure and reported number. Each model fit ran
for about a day on a university cluster, and some inputs are proprietary or must be downloaded from
their providers, so a rerun needs the inputs listed below and access to a cluster.

## Contents

| Folder | Contents |
|---|---|
| `01_data_prep/` | Footprint sampling from the earlier analysis (GEDI quality checks, 3DEP download and resampling, footprint means of the P3D DTM and height model and of the reference rasters, P3D image metadata, land cover and ecoregions), with its R workflow and join checks. This paper uses the per-site tables it produced and the Fig. 1 site centroids from `r_workflow/extract_section2_1_study_sites.R`. |
| `02_lidar_references/` | Lidar ground surface and pit-free canopy height model built from each site's point cloud and averaged over every footprint |
| `03_analysis_frame/` | The analysis frame (datum steps, screening against 3DEP, standardization, training sample and held-out groups), the delivered DSM sampled at every footprint, and the step that merges the two |
| `04_error_models/` | The Bayesian hierarchical error models: full fits, held-out fits, the divergent-draw check, the posterior predictive check and the cluster job queue |
| `05_tables_figures/` | The scripts that build every table and figure, with the run order in `RUN_ORDER.md` |
| `data_derived/` | The inputs these scripts read and their outputs as published, described in `data_derived/README_tables.md` |
| `docs/` | `MANIFEST.md` and the run logs of the model fits |

## Not included

- Vantor Precision3D imagery and the delivered DSM, DTM and height-model rasters. They are
  proprietary and were provided under a data sharing agreement, so they cannot be redistributed.
  NASA-funded researchers can request Precision3D through NASA's Commercial Satellite Data
  Acquisition program.
- Footprint-level tables (the analysis frame, the lidar values at each footprint and the
  per-footprint model predictions). Together with the public lidar they would give back the P3D
  values.
- Fitted model objects, which are large and contain the footprint-level data. Their posterior
  summaries are in `data_derived/`.
- Public inputs. GEDI L2A and L2B come from the NASA LP DAAC and GEDI L4A and L4C from the ORNL DAAC.
  3DEP elevation comes from the USGS, NEON lidar from the NEON Data Portal, and the USDA Forest
  Service and NASA lidar from the acquisitions listed in Table 1 of the paper. Land cover is
  the 2022 layer of GLC_FCS30D, the ecoregions are the CEC Level II map, and the geoid grid is
  NOAA's GEOID18.
- The per-site GEDI footprint files (`GEDI_site<N>_hq_ALL.gpkg`). They were built in the University
  of Maryland GEDI lab archive with the filters given in Section 2.2.3 of the paper. The script that
  built them is not in this repository.

## Software

The cluster steps ran in R 4.5.0 with `brms` and `cmdstanr` on CmdStan 2.37.0. The tables and figures were
rebuilt from the released scripts with R 4.4.0 and Python 3.12, and they reproduce the published
values and figures exactly. `05_tables_figures/RUN_ORDER.md` lists the package versions used for the tables and figures. Each full
fit ran four chains of 3,000 iterations, half of them warmup, with `adapt_delta` 0.98,
`max_treedepth` 14 and seed 2026, and took 18 to 24 hours on 32 cores. The 18 held-out fits used
`adapt_delta` 0.99 and ran as a job queue over several days. Their logs are in `docs/run_logs/`.

## Notes on provenance

- Two DSM errors appear in the code. The models were fit to the DSM error of Eq. 4 in the paper
  (the CHM error plus the DTM error), which Table 3a, Table S5 and Fig. S1 use. Table 2, Table 3b,
  Tables S2 to S4, Figs. 2 to 4 and Fig. S2 use the error of the DSM sampled directly from the
  delivered raster, as Text S2 explains.
- The scripts in `01_data_prep/` are unchanged from version 1. They keep the original cluster paths,
  and their local input folders are written as `<LOCAL_DATA_ROOT>`, so both need to be set before use.
- The analysis frame reuses the footprint sampling of the earlier analysis (`01_data_prep/`) for the
  P3D DTM and height model, 3DEP, the GEDI variables, the image metadata, land cover and ecoregions.
  The lidar references and the P3D DSM were sampled again for this paper.
- The direct DSM sampling ran with the current `p3d_dsm_direct.R` at Sites 1, 5 and 17 (internal ids
  1, 5 and 19). The other 15 sites ran with an earlier revision of the same script that was not kept.
  The current revision adds a check on which footprint positions to read. At those 15 sites the first
  positions it tries reproduce the P3D DTM values of the frame, so it should read the same positions,
  but the 15 sites were not rerun with it.
- The markers in Fig. 1 are site centroids from the footprint set of the earlier analysis. They lie
  within about 20 km of the centroids of the footprints analyzed here.
- The Fig. 5 transects were selected with the Eq. 4 DSM error, before the DSM was sampled directly.
  With the directly sampled error the same windows rank first at Sites 6 and 13 and second at Site 4.
  The lidar surfaces in Fig. 5 come from 1 km tiles at 1 m resolution, built with the same ground and CHM steps as
  `02_lidar_references/lidar_products_rebuild.R`, run over whole tiles. The tile script is not
  included.
- Site numbers inside the scripts and per-site files are internal site ids. The `tracker` column of
  `data_derived/site_names.csv` holds them, next to the site numbers and names used in the paper.
- Names such as `rework_*`, `chpt1` and `chpt2` in scripts and paths come from the project's own
  folder layout. They are kept because the run logs print them.
- The outputs of the divergent-draw check and some cluster logs (the job queue, the lidar rebuild and
  the DSM sampling) are not included.

## Versions

Version 2.0.0 goes with the submitted manuscript named above. Versions 1.0.0 to 1.0.2 (September
2026) archived the analysis of an earlier version of the manuscript, which used different models and
site numbers. They remain available from the Zenodo record and as tagged releases here.

## License and citation

The code is under the MIT License (`LICENSE`) and the tables in `data_derived/` under CC BY 4.0
(`data_derived/LICENSE`). If you use either, please cite the paper (see `CITATION.cff`). Releases are
archived at Zenodo under the concept DOI https://doi.org/10.5281/zenodo.22831703.

## Funding

This work was carried out within the NASA GEDI mission (contract NNL15AA03C to the University of
Maryland).
