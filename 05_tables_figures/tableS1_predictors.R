# tableS1_predictors.R  --  Table S1, the 17 variables and 3 groupings of the error models
# Rows follow the model formula in rework_fullfit.R (main effects; groups site, ecoregion, lc_grp);
# descriptions and sources as printed in Table S1. Mean and SD = the constants that standardize
# each variable across all 630,390 frame footprints (rework_model_prep.R STEP 4), so a coefficient
# per SD converts to units.
# Input   (P3D_DERIVED) zconst.csv (rework_model_prep.R, cluster $CH1/manuscript_tables/rework/zconst.csv)
# Outputs (P3D_RESULTS) tableS1/tableS1_predictors.{csv,md}
# Check   the 17 variables = the formula's.
suppressPackageStartupMessages(library(data.table))
DERIVED  <- Sys.getenv("P3D_DERIVED", "data_derived")
OUT_ROOT <- Sys.getenv("P3D_RESULTS", "results")
OUT  <- file.path(OUT_ROOT, "tableS1"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
z  <- fread(file.path(DERIVED, "zconst.csv"))

D <- fread(text = "variable|group|name|description|units|source
slope_mean|Terrain|Slope|Mean terrain slope over the footprint|degrees|3DEP
slope_sd|Terrain|Slope variability|Standard deviation of slope over the footprint|degrees|3DEP
aspect_sin_mean|Terrain|Aspect (sine)|East-west orientation|unitless|3DEP
aspect_cos_mean|Terrain|Aspect (cosine)|North-south orientation|unitless|3DEP
rh_98|Canopy|Canopy height (RH98)|GEDI 98th percentile relative height|m|GEDI L2A (Hofton and Blair, 2019)
cover|Canopy|Canopy cover|GEDI canopy cover fraction|proportion|GEDI L2B (Tang and Armston, 2019)
wsci|Canopy|Structural complexity (WSCI)|Waveform structural complexity index|unitless|GEDI L4C (de Conto et al., 2024)
meta_off_nadir_avg|Image|Off-nadir angle|Mean off-nadir viewing angle of the images covering the footprint|degrees|P3D metadata
meta_sun_elev_avg|Image|Sun elevation|Mean sun elevation of the images covering the footprint|degrees|P3D metadata
meta_az_concentration|Image|Azimuth concentration|Mean resultant length of the images' viewing azimuths (0 when evenly spread, 1 when all from one direction)|unitless|P3D metadata
view_az_sin|Image|Viewing azimuth (sine)|East-west component of the mean viewing azimuth|unitless|P3D metadata
view_az_cos|Image|Viewing azimuth (cosine)|North-south component of the mean viewing azimuth|unitless|P3D metadata
meta_stereo_ratio|Image|Stereo ratio|Share of images acquired as in-track stereo|proportion|P3D metadata
fwdrev|Image|Forward-reverse balance|Share of forward-looking minus share of reverse-looking images|unitless|P3D metadata
meta_leaf_on_ratio|Image|Leaf-on fraction|Share of images acquired between 9 March and 13 October (leaf-on season)|proportion|P3D metadata
meta_abs_geoacc_avg|Image|Absolute geoaccuracy|Mean absolute geolocation accuracy of the images covering the footprint|m|P3D metadata
meta_rel_geoacc_avg|Image|Relative geoaccuracy|Mean relative geolocation accuracy of the images covering the footprint|m|P3D metadata", sep = "|")

stopifnot(nrow(D) == 17L, setequal(D$variable, z$variable), nrow(z) == 17L)
m <- merge(D, z[, .(variable, mean, sd)], by = "variable", sort = FALSE)
m <- m[match(D$variable, variable)]

G <- data.table(group = "Grouping", name = c("Site", "Level II ecoregion", "Land-cover class"),
  description = c("The 18 sites (Table 1)",
                  "Ecoregion of the footprint",
                  "Class of the footprint, with the ten most common classes kept separate and all other classes combined into one group"),
  units = "category", source = c("This study", "CEC Level II (Commission for Environmental Cooperation, 1997)",
                                 "GLC_FCS30 (Zhang et al., 2021)"))
fwrite(rbind(m, G, fill = TRUE), file.path(OUT, "tableS1_predictors.csv"))

fm <- function(x, v) {
  s <- sprintf(if (v %in% c("cover", "meta_az_concentration", "meta_stereo_ratio", "meta_leaf_on_ratio", "fwdrev",
                            "aspect_sin_mean", "aspect_cos_mean", "view_az_sin", "view_az_cos")) "%.2f" else "%.1f", x)
  s <- sub("^-(0\\.0+)$", "\\1", s)                  # no "-0.00"
  sub("-", "\u2212", s, fixed = TRUE)
}
md <- c("| Group | Variable | Description | Units | Source | Mean | SD |", "|---|---|---|---|---|---|---|",
        m[, sprintf("| %s | %s | %s | %s | %s | %s | %s |", group, name, description, units, source,
                    mapply(fm, mean, variable), mapply(fm, sd, variable))],
        G[, sprintf("| %s | %s | %s | %s | %s |  |  |", group, name, description, units, source)])
writeLines(md, file.path(OUT, "tableS1_predictors.md"))
cat(md, sep = "\n")
