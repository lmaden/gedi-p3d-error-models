#!/usr/bin/env python3
"""
gedi_sample_chunked.py
==========================

Full GEDI sampling pipeline with:
- Complete sampling from v09 (terrain, DTM, CHM, errors)
- CHUNKED processing for aligned pair metrics (from resample_3dep script)
- Aggressive memory management for large sites (1, 5, etc.)
- Fine-grained per-step + per-chunk checkpointing

Samples:
- GEDI attributes (rh_98, cover, pft, wsci, elev_low, delta_time, shot_number)
- Ecoregion join
- Slope (mean, sd, valid_frac)
- Aspect (sin_mean, cos_mean, valid_frac)
- 3DEP DTM (mean, median, percentiles, valid_frac, center point)
- P3D DTM (mean, median, percentiles, valid_frac, center point)
- ALS CHM (mean, median, sd, percentiles, valid_frac, center point)
- P3D CHM (mean, median, sd, percentiles, valid_frac, center point)
- CHM error metrics (RMSE, error_mean, error_median, valid_frac) [CHUNKED]
- DTM error metrics (RMSE, error_mean, error_median, valid_frac) [CHUNKED]

Usage:
  python gedi_sample_chunked.py --site 5
  python gedi_sample_chunked.py --site 5 --cleanup
"""

from __future__ import annotations

import gc
import json
import os
import pickle
from pathlib import Path
from typing import Dict, Any, Tuple, List, Optional, Set

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio

# Import from helper_speed
from helper_speed import (
    build_vrt,
    zonal_stats,
    zonal_aligned_pair_metrics,
    sample_points,
)

# ---------------------------------------------------------------------
# CONFIG - UPDATED FOR YOUR DRIVE SETUP
# ---------------------------------------------------------------------
os.environ.setdefault("GDAL_CACHEMAX", "2048")  # MB

base_gedi_dir = Path(r"C:\Users\levim\OneDrive\Desktop\work\gedi\dissertation\chpt1\data\gedi\sites")
eco_path      = Path(r"C:\Users\levim\OneDrive\Desktop\work\gedi\dissertation\chpt1\data\general\na_cec_eco_l2\NA_CEC_Eco_Level2.shp")

# Drive paths - CORRECTED for your setup
slope_base    = Path(r"E:\slope")
aspect_base   = Path(r"E:\aspect")
dep_base      = Path(r"E:\3dep")
als_chm_base  = Path(r"E:\chms")
mhrsi_base    = Path(r"F:\p3d")

OUTPUT_CSV = Path(r"C:\Users\levim\OneDrive\Desktop\work\gedi\dissertation\chpt1\data\sampled_data_FINAL.csv")
META_JSON  = OUTPUT_CSV.with_suffix(".json")

# Footprint radius & controls
FOOTPRINT_RADIUS = 12.5  # meters
MIN_VALID_FRAC: float | None = None

# Diagnostics toggles
BYPASS_CHM_MASK = False

# Strict alignment for error metrics
STRICT_ALIGN_ERROR_MEAN   = True
STRICT_ALIGN_ERROR_MEDIAN = True

# Percentiles
CHM_PERCENTILES = (0.90, 0.75, 0.25)
DTM_PERCENTILES = (0.90, 0.75, 0.25)

# Chunk size for aligned metrics (100K footprints per chunk)
ALIGNED_CHUNK_SIZE = 100_000

# Sites to process
SITES = list(range(1, 21))

# ---------------------------------------------------------------------
# Ecoregion data (loaded once at import)
# ---------------------------------------------------------------------
print("Reading ecoregion shapefile...")
eco_l2_src = gpd.read_file(eco_path)
eco_cache: dict[str, gpd.GeoDataFrame] = {}
print("Ecoregion data loaded.")

# ---------------------------------------------------------------------
# Checkpoint helpers
# ---------------------------------------------------------------------
def get_checkpoint_dir(site: int) -> Path:
    """Get checkpoint directory for a site."""
    return Path(f"./checkpoints_site_{site:02d}")


def save_checkpoint(site: int, name: str, data: Dict):
    """Save checkpoint data for a specific sampling step."""
    ckpt_dir = get_checkpoint_dir(site)
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    
    ckpt_file = ckpt_dir / f"{name}.pkl"
    with open(ckpt_file, 'wb') as f:
        pickle.dump(data, f)
    
    print(f"   💾 Checkpoint saved: {ckpt_file.name}")


def load_checkpoint(site: int, name: str) -> Optional[Dict]:
    """Load checkpoint data if it exists."""
    ckpt_file = get_checkpoint_dir(site) / f"{name}.pkl"
    if ckpt_file.exists():
        with open(ckpt_file, 'rb') as f:
            return pickle.load(f)
    return None


def list_checkpoints(site: int) -> List[str]:
    """List available checkpoints for a site."""
    ckpt_dir = get_checkpoint_dir(site)
    if not ckpt_dir.exists():
        return []
    return [p.stem for p in ckpt_dir.glob("*.pkl")]


def get_existing_chunks(site: int, prefix: str) -> Set[int]:
    """Get set of existing chunk checkpoints."""
    ckpt_dir = get_checkpoint_dir(site)
    existing = set()
    if ckpt_dir.exists():
        for p in ckpt_dir.glob(f"{prefix}_chunk_*.pkl"):
            try:
                chunk_num = int(p.stem.split('_')[-1])
                existing.add(chunk_num)
            except:
                pass
    return existing


def cleanup_chunk_checkpoints(site: int, prefix: str, n_chunks: int):
    """Remove individual chunk checkpoints after combining."""
    ckpt_dir = get_checkpoint_dir(site)
    for i in range(n_chunks):
        chunk_file = ckpt_dir / f'{prefix}_chunk_{i}.pkl'
        if chunk_file.exists():
            try:
                chunk_file.unlink()
            except:
                pass


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------
def _rgp(root: Path, *sub: str) -> List[str]:
    """Recursively glob for .tif files."""
    base = root.joinpath(*sub)
    if not base.exists():
        return []
    out = [str(p) for p in base.rglob("*.tif")]
    out += [str(p) for p in base.rglob("*.tiff")]
    return out


def _quick_band_stats(ds) -> str:
    """Quick statistics for a raster dataset."""
    from rasterio.enums import Resampling
    h = max(1, ds.height // 512)
    w = max(1, ds.width // 512)
    arr = ds.read(1, out_shape=(1, h, w), resampling=Resampling.bilinear).astype(np.float32)
    msk = ds.read_masks(1, out_shape=(1, h, w), resampling=Resampling.nearest)
    arr[msk == 0] = np.nan
    finite = np.isfinite(arr)
    if finite.any():
        vmin, vmax = float(np.nanmin(arr)), float(np.nanmax(arr))
    else:
        vmin, vmax = float("nan"), float("nan")
    invalid_frac = float(np.isnan(arr).sum()) / arr.size
    return f"dtype={ds.dtypes[0]}, nodata={ds.nodata}, min≈{vmin:.3f}, max≈{vmax:.3f}, invalid≈{invalid_frac:.1%}"


def _vfrac_summary(name: str, vf: np.ndarray) -> str:
    """Summary of valid fraction statistics."""
    any_valid = float(np.isfinite(vf).mean())
    ge10 = float((vf >= 0.10).mean())
    ge50 = float((vf >= 0.50).mean())
    return f"{name}: valid-frac — any={any_valid:.1%}, ≥10%={ge10:.1%}, ≥50%={ge50:.1%}"


# ---------------------------------------------------------------------
# Chunked aligned metrics computation
# ---------------------------------------------------------------------
def compute_aligned_metrics_chunked(
    site: int,
    ds1: rasterio.DatasetReader,
    ds2: rasterio.DatasetReader,
    xs: np.ndarray,
    ys: np.ndarray,
    prefix: str,  # 'chm' or 'dtm'
    use_dataset_mask: bool = True,
    _log = None
) -> Dict[str, np.ndarray]:
    """
    Compute aligned pair metrics in chunks with per-chunk checkpointing.
    
    This prevents memory crashes on large sites by:
    1. Processing footprints in 100K chunks
    2. Saving each chunk immediately after computation
    3. Resuming from last completed chunk on restart
    
    Returns dict with: rmse, vfrac, error_mean, error_median
    """
    if _log is None:
        _log = print
    
    n_footprints = len(xs)
    n_chunks = (n_footprints + ALIGNED_CHUNK_SIZE - 1) // ALIGNED_CHUNK_SIZE
    
    _log(f"   Processing {n_footprints:,} footprints in {n_chunks} chunks")
    
    # Check for existing combined checkpoint
    combined_data = load_checkpoint(site, f'{prefix}_aligned_combined')
    if combined_data:
        _log(f"   ✓ Loaded combined {prefix.upper()} aligned metrics from checkpoint")
        return combined_data
    
    # Check which chunks are already done
    existing_chunks = get_existing_chunks(site, f'{prefix}_aligned')
    if existing_chunks:
        _log(f"   Found {len(existing_chunks)} existing chunk checkpoint(s)")
        if existing_chunks != set(range(n_chunks)):
            missing = set(range(n_chunks)) - existing_chunks
            _log(f"   Will resume from chunk {min(missing) + 1}")
    
    # Pre-allocate output arrays
    rmse_all = np.full(n_footprints, np.nan, dtype=np.float32)
    vfrac_all = np.full(n_footprints, np.nan, dtype=np.float32)
    mu1_all = np.full(n_footprints, np.nan, dtype=np.float32)
    mu2_all = np.full(n_footprints, np.nan, dtype=np.float32)
    med1_all = np.full(n_footprints, np.nan, dtype=np.float32)
    med2_all = np.full(n_footprints, np.nan, dtype=np.float32)
    
    # Track failed chunks
    failed_chunks = []
    
    # Process each chunk
    for i in range(n_chunks):
        start_idx = i * ALIGNED_CHUNK_SIZE
        end_idx = min((i + 1) * ALIGNED_CHUNK_SIZE, n_footprints)
        chunk_size = end_idx - start_idx
        
        # Check if this chunk is already done
        if i in existing_chunks:
            print(f"   Chunk {i+1}/{n_chunks}: ✓ Loading from checkpoint... ", end='', flush=True)
            chunk_data = load_checkpoint(site, f'{prefix}_aligned_chunk_{i}')
            rmse_all[start_idx:end_idx] = chunk_data['rmse']
            vfrac_all[start_idx:end_idx] = chunk_data['vfrac']
            mu1_all[start_idx:end_idx] = chunk_data['mu1']
            mu2_all[start_idx:end_idx] = chunk_data['mu2']
            med1_all[start_idx:end_idx] = chunk_data['med1']
            med2_all[start_idx:end_idx] = chunk_data['med2']
            print(f"({(i+1)/n_chunks*100:.1f}% complete)")
            continue
        
        print(f"   Chunk {i+1}/{n_chunks}: footprints {start_idx:,} to {end_idx:,}...", end='', flush=True)
        
        try:
            # Get coordinates for this chunk
            xs_chunk = xs[start_idx:end_idx]
            ys_chunk = ys[start_idx:end_idx]
            
            # Compute aligned metrics for this chunk
            rmse_chunk, mu1, mu2, vfrac_chunk, med1, med2 = zonal_aligned_pair_metrics(
                ds1, ds2, xs_chunk, ys_chunk, FOOTPRINT_RADIUS,
                want_median=True,
                return_valid_fraction=True,
                min_valid_frac=MIN_VALID_FRAC,
                use_dataset_mask=use_dataset_mask
            )
            
            # Store results in main arrays
            rmse_all[start_idx:end_idx] = rmse_chunk.astype(np.float32)
            vfrac_all[start_idx:end_idx] = vfrac_chunk.astype(np.float32)
            mu1_all[start_idx:end_idx] = mu1.astype(np.float32)
            mu2_all[start_idx:end_idx] = mu2.astype(np.float32)
            med1_all[start_idx:end_idx] = med1.astype(np.float32)
            med2_all[start_idx:end_idx] = med2.astype(np.float32)
            
            # Save chunk checkpoint IMMEDIATELY
            save_checkpoint(site, f'{prefix}_aligned_chunk_{i}', {
                'rmse': rmse_chunk.astype(np.float32),
                'vfrac': vfrac_chunk.astype(np.float32),
                'mu1': mu1.astype(np.float32),
                'mu2': mu2.astype(np.float32),
                'med1': med1.astype(np.float32),
                'med2': med2.astype(np.float32)
            })
            
            print(f" ✓ ({(i+1)/n_chunks*100:.1f}% complete)")
            
            # Explicit cleanup after each chunk
            del rmse_chunk, mu1, mu2, vfrac_chunk, med1, med2, xs_chunk, ys_chunk
            gc.collect()
            
        except Exception as chunk_error:
            # Handle corrupted tiles - log and continue with NaN
            print(f" ⚠️  FAILED (corrupted data)")
            print(f"      Error: {str(chunk_error)[:100]}")
            
            # Arrays already initialized with NaN, so chunk will remain NaN
            failed_chunks.append(i + 1)
            
            # Save checkpoint with NaN values so we don't retry this chunk
            save_checkpoint(site, f'{prefix}_aligned_chunk_{i}', {
                'rmse': np.full(chunk_size, np.nan, dtype=np.float32),
                'vfrac': np.full(chunk_size, np.nan, dtype=np.float32),
                'mu1': np.full(chunk_size, np.nan, dtype=np.float32),
                'mu2': np.full(chunk_size, np.nan, dtype=np.float32),
                'med1': np.full(chunk_size, np.nan, dtype=np.float32),
                'med2': np.full(chunk_size, np.nan, dtype=np.float32)
            })
            
            print(f"      Filled with NaN, continuing... ({(i+1)/n_chunks*100:.1f}% complete)")
            
            # Cleanup
            try:
                del xs_chunk, ys_chunk
                gc.collect()
            except:
                pass
    
    # Compute error metrics from aligned means
    error_mean = (mu2_all - mu1_all).astype(np.float32)
    error_median = (med2_all - med1_all).astype(np.float32)
    
    # Clean up intermediate arrays
    del mu1_all, mu2_all, med1_all, med2_all
    gc.collect()
    
    # Report any failed chunks
    if failed_chunks:
        _log(f"   ⚠️  WARNING: {len(failed_chunks)} chunk(s) failed due to corrupted data:")
        _log(f"      Chunks: {failed_chunks}")
        _log(f"      These chunks are filled with NaN and will be excluded from analysis")
    
    # Report summary
    valid_errors = np.isfinite(error_mean)
    if valid_errors.sum() > 0:
        mean_err = error_mean[valid_errors].mean()
        median_err = np.median(error_mean[valid_errors])
        rmse_val = np.sqrt((error_mean[valid_errors]**2).mean())
        _log(f"   ✓ All chunks complete! {prefix.upper()} errors computed:")
        _log(f"      Mean:   {mean_err:.3f} m")
        _log(f"      Median: {median_err:.3f} m")
        _log(f"      RMSE:   {rmse_val:.3f} m")
        _log(f"      N valid: {valid_errors.sum():,} / {len(error_mean):,} ({100*valid_errors.sum()/len(error_mean):.1f}%)")
    else:
        _log(f"   ⚠️  WARNING: No valid {prefix.upper()} errors computed (all data corrupted)")
    
    # Save combined checkpoint
    combined_results = {
        'rmse': rmse_all,
        'vfrac': vfrac_all,
        'error_mean': error_mean,
        'error_median': error_median
    }
    save_checkpoint(site, f'{prefix}_aligned_combined', combined_results)
    
    # Clean up individual chunk checkpoints
    _log(f"   Cleaning up {n_chunks} chunk checkpoints...")
    cleanup_chunk_checkpoints(site, f'{prefix}_aligned', n_chunks)
    
    return combined_results


# ---------------------------------------------------------------------
# Main worker function
# ---------------------------------------------------------------------
def process_site(site: int) -> Tuple[int, pd.DataFrame | None, str]:
    """
    Process a single site with fine-grained checkpointing.
    
    Returns:
        (site_id, dataframe, log_string)
    """
    log_lines: list[str] = []
    def _log(msg: str):
        log_lines.append(f"[{pd.Timestamp.now()}] {msg}")
        print(f"[{pd.Timestamp.now()}] {msg}")
    
    _log(f"="*80)
    _log(f"STARTING SITE {site}")
    _log(f"="*80)
    
    # Check for existing checkpoints
    existing_ckpts = list_checkpoints(site)
    if existing_ckpts:
        _log(f"📂 Found {len(existing_ckpts)} existing checkpoints: {', '.join(existing_ckpts[:5])}...")
        _log(f"   Will resume from last checkpoint")
    else:
        _log(f"📂 No existing checkpoints - starting fresh")
    
    # ---------------------------------------------------------------------
    # STEP 1: Load GEDI data
    # ---------------------------------------------------------------------
    gedi_data = load_checkpoint(site, 'gedi_loaded')
    if gedi_data:
        _log("1. ✓ Loaded GEDI data from checkpoint")
        gdf_orig = gedi_data['gdf']
    else:
        _log("1. Loading GEDI data...")
        gedi_fp = base_gedi_dir / str(site) / f"GEDI_site{site}_hq_ALL.gpkg"
        if not gedi_fp.exists():
            return site, None, f"GEDI geopackage not found: {gedi_fp}"
        
        wanted_cols = [
            "geometry", "rh_opt_098", "cover_z_000", "land_cover_data/pft_class",
            "elev_lowestmode", "wsci", "delta_time", "shot_number"
        ]
        try:
            gdf_orig = gpd.read_file(gedi_fp, columns=wanted_cols)
        except TypeError:
            tmp = gpd.read_file(gedi_fp)
            keep = [c for c in wanted_cols if c in tmp.columns or c == "geometry"]
            gdf_orig = tmp[keep].copy()
            del tmp
            gc.collect()
        
        if gdf_orig.empty:
            return site, None, "GEDI file is empty."
        
        # Normalize column names
        gdf_orig.columns = [c.replace("/", "_") for c in gdf_orig.columns]
        gdf_orig = gdf_orig.rename(columns={
            "cover_z_000": "cover",
            "land_cover_data_pft_class": "pft",
            "elev_lowestmode": "elev_low"
        })
        if "rh_opt_098" in gdf_orig.columns:
            gdf_orig["rh_98"] = gdf_orig["rh_opt_098"] / 100.0
        
        _log(f"   Loaded {len(gdf_orig):,} GEDI observations")
        save_checkpoint(site, 'gedi_loaded', {'gdf': gdf_orig})
    
    # ---------------------------------------------------------------------
    # STEP 2: Ecoregion join
    # ---------------------------------------------------------------------
    eco_data = load_checkpoint(site, 'ecoregion_joined')
    if eco_data:
        _log("2. ✓ Loaded ecoregion join from checkpoint")
        gdf_orig = eco_data['gdf']
    else:
        _log("2. Joining ecoregions...")
        try:
            eco_key = gdf_orig.crs.to_string()
            if eco_key not in eco_cache:
                eco_cache[eco_key] = eco_l2_src.to_crs(gdf_orig.crs)
            gdf_orig = (
                gpd.sjoin(
                    gdf_orig,
                    eco_cache[eco_key][["geometry", "NA_L2NAME"]],
                    how="left",
                    predicate="intersects",
                )
                .rename(columns={"NA_L2NAME": "ecoregion"})
                .drop(columns="index_right", errors='ignore')
            )
            _log(f"   ✓ Ecoregions joined")
        except Exception as e:
            _log(f"   ⚠ Ecoregion join failed: {e}")
            gdf_orig["ecoregion"] = None
        
        save_checkpoint(site, 'ecoregion_joined', {'gdf': gdf_orig})
    
    # Get shot numbers for later use
    n_shots = len(gdf_orig)
    _log(f"   Processing {n_shots:,} shots")
    
    # ---------------------------------------------------------------------
    # STEP 3: Discover and open datasets
    # ---------------------------------------------------------------------
    datasets_data = load_checkpoint(site, 'datasets_opened')
    if datasets_data:
        _log("3. ✓ Reopening datasets from checkpoint paths")
        data_sources = datasets_data['data_sources']
        datasets_to_sample = {}
        for name, cfg in data_sources.items():
            if cfg["paths"]:
                try:
                    ds = build_vrt(cfg["paths"])
                    datasets_to_sample[name] = ds
                except Exception as e:
                    _log(f"   ⚠ Failed to reopen '{name}': {e}")
    else:
        _log("3. Discovering and opening datasets...")
        data_sources = {
            "slope":   {"paths": _rgp(slope_base, f"polygon_{site}")},
            "aspect":  {"paths": _rgp(aspect_base, f"polygon_{site}")},
            "dep_dtm": {"paths": _rgp(dep_base, f"polygon_{site}")},
            "p3d_dtm": {"paths": _rgp(mhrsi_base, str(site), "vricon_raster_50cm", "dtm")},
            "p3d_chm": {"paths": _rgp(mhrsi_base, str(site), "vricon_raster_50cm", "dhm")},
            "als_chm": {"paths": _rgp(als_chm_base, str(site), "nad83") or _rgp(als_chm_base, str(site), "chm")},
        }
        
        datasets_to_sample = {}
        for name, cfg in data_sources.items():
            if not cfg["paths"]:
                _log(f"   ⚠ No files found for '{name}', skipping")
                continue
            try:
                ds = build_vrt(cfg["paths"])
                _log(f"   ✓ Opened '{name}': {len(cfg['paths'])} tiles, CRS={ds.crs.to_string()}")
                datasets_to_sample[name] = ds
            except Exception as e:
                _log(f"   ⚠ Failed to open '{name}': {e}")
        
        if not datasets_to_sample:
            return site, None, "No datasets to sample."
        
        save_checkpoint(site, 'datasets_opened', {'data_sources': data_sources})
    
    # Coordinate cache (per-CRS)
    coords_cache: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    def _coords_for(crs) -> tuple[np.ndarray, np.ndarray]:
        key = crs.to_string()
        if key not in coords_cache:
            gdf_r = gdf_orig.to_crs(crs)
            coords_cache[key] = (
                gdf_r.geometry.x.values.astype(np.float64),
                gdf_r.geometry.y.values.astype(np.float64),
            )
            del gdf_r
            gc.collect()
        return coords_cache[key]
    
    # Results dictionary
    results: Dict[str, Any] = {}
    
    # ---------------------------------------------------------------------
    # STEP 4: Sample SLOPE
    # ---------------------------------------------------------------------
    if "slope" in datasets_to_sample:
        slope_data = load_checkpoint(site, 'slope_stats')
        if slope_data:
            _log("4. ✓ Loaded slope stats from checkpoint")
            results.update(slope_data['stats'])
        else:
            _log("4. Sampling slope...")
            ds = datasets_to_sample["slope"]
            xs, ys = _coords_for(ds.crs)
            
            ignore_nodata = bool(ds.nodata is not None and float(ds.nodata) == 0.0)
            if ignore_nodata:
                _log("   Slope nodata==0 detected → ignoring nodata value")
            
            st = zonal_stats(
                ds, xs, ys, FOOTPRINT_RADIUS,
                want_std=True, return_valid_fraction=True,
                min_valid_frac=MIN_VALID_FRAC,
                use_dataset_mask=True, ignore_nodata_value=ignore_nodata
            )
            
            slope_results = {
                "slope_mean": st["mean"].astype(np.float32),
                "slope_sd": st["std"].astype(np.float32),
                "slope_valid_frac": st["valid_frac"].astype(np.float32)
            }
            results.update(slope_results)
            _log(f"   {_vfrac_summary('Slope', st['valid_frac'])}")
            save_checkpoint(site, 'slope_stats', {'stats': slope_results})
            del st
            gc.collect()
    
    # ---------------------------------------------------------------------
    # STEP 5: Sample ASPECT
    # ---------------------------------------------------------------------
    if "aspect" in datasets_to_sample:
        aspect_data = load_checkpoint(site, 'aspect_stats')
        if aspect_data:
            _log("5. ✓ Loaded aspect stats from checkpoint")
            results.update(aspect_data['stats'])
        else:
            _log("5. Sampling aspect...")
            ds = datasets_to_sample["aspect"]
            xs, ys = _coords_for(ds.crs)
            
            sin_stats = zonal_stats(
                ds, xs, ys, FOOTPRINT_RADIUS,
                transform_fn=lambda v: np.sin(np.deg2rad(v)),
                return_valid_fraction=True,
                min_valid_frac=MIN_VALID_FRAC, use_dataset_mask=True
            )
            cos_stats = zonal_stats(
                ds, xs, ys, FOOTPRINT_RADIUS,
                transform_fn=lambda v: np.cos(np.deg2rad(v)),
                return_valid_fraction=False,
                min_valid_frac=MIN_VALID_FRAC, use_dataset_mask=True
            )
            
            aspect_results = {
                "aspect_sin_mean": sin_stats["mean"].astype(np.float32),
                "aspect_cos_mean": cos_stats["mean"].astype(np.float32),
                "aspect_valid_frac": sin_stats["valid_frac"].astype(np.float32)
            }
            results.update(aspect_results)
            _log(f"   {_vfrac_summary('Aspect', sin_stats['valid_frac'])}")
            save_checkpoint(site, 'aspect_stats', {'stats': aspect_results})
            del sin_stats, cos_stats
            gc.collect()
    
    # ---------------------------------------------------------------------
    # STEP 6: Sample ALS CHM
    # ---------------------------------------------------------------------
    if "als_chm" in datasets_to_sample:
        als_chm_data = load_checkpoint(site, 'als_chm_stats')
        if als_chm_data:
            _log("6. ✓ Loaded ALS CHM stats from checkpoint")
            results.update(als_chm_data['stats'])
        else:
            _log("6. Sampling ALS CHM...")
            ds = datasets_to_sample["als_chm"]
            xs, ys = _coords_for(ds.crs)
            
            use_mask = not BYPASS_CHM_MASK
            if not use_mask:
                _log("   CHM mask bypass is ON")
            
            st = zonal_stats(
                ds, xs, ys, FOOTPRINT_RADIUS,
                want_median=True,
                want_std=True,
                percentiles=CHM_PERCENTILES,
                return_valid_fraction=True,
                min_valid_frac=MIN_VALID_FRAC, use_dataset_mask=use_mask
            )
            
            als_chm_results = {
                "als_chm_mean": st["mean"].astype(np.float32),
                "als_chm_median": st["median"].astype(np.float32),
                "als_chm_sd": st["std"].astype(np.float32),
                "als_chm_valid_frac": st["valid_frac"].astype(np.float32),
                "als_chm": sample_points(ds, xs, ys, use_dataset_mask=use_mask)
            }
            # Add percentiles
            for pct in CHM_PERCENTILES:
                pct_key = f"p{int(pct*100):02d}"
                if pct_key in st:
                    als_chm_results[f"als_chm_{pct_key}"] = st[pct_key].astype(np.float32)
            
            results.update(als_chm_results)
            _log(f"   {_vfrac_summary('ALS CHM', st['valid_frac'])}")
            save_checkpoint(site, 'als_chm_stats', {'stats': als_chm_results})
            del st
            gc.collect()
    
    # ---------------------------------------------------------------------
    # STEP 7: Sample P3D CHM
    # ---------------------------------------------------------------------
    if "p3d_chm" in datasets_to_sample:
        p3d_chm_data = load_checkpoint(site, 'p3d_chm_stats')
        if p3d_chm_data:
            _log("7. ✓ Loaded P3D CHM stats from checkpoint")
            results.update(p3d_chm_data['stats'])
        else:
            _log("7. Sampling P3D CHM...")
            ds = datasets_to_sample["p3d_chm"]
            xs, ys = _coords_for(ds.crs)
            
            use_mask = not BYPASS_CHM_MASK
            
            st = zonal_stats(
                ds, xs, ys, FOOTPRINT_RADIUS,
                want_median=True,
                want_std=True,
                percentiles=CHM_PERCENTILES,
                return_valid_fraction=True,
                min_valid_frac=MIN_VALID_FRAC, use_dataset_mask=use_mask
            )
            
            p3d_chm_results = {
                "p3d_chm_mean": st["mean"].astype(np.float32),
                "p3d_chm_median": st["median"].astype(np.float32),
                "p3d_chm_sd": st["std"].astype(np.float32),
                "p3d_chm_valid_frac": st["valid_frac"].astype(np.float32),
                "p3d_chm": sample_points(ds, xs, ys, use_dataset_mask=use_mask)
            }
            # Add percentiles
            for pct in CHM_PERCENTILES:
                pct_key = f"p{int(pct*100):02d}"
                if pct_key in st:
                    p3d_chm_results[f"p3d_chm_{pct_key}"] = st[pct_key].astype(np.float32)
            
            results.update(p3d_chm_results)
            _log(f"   {_vfrac_summary('P3D CHM', st['valid_frac'])}")
            save_checkpoint(site, 'p3d_chm_stats', {'stats': p3d_chm_results})
            del st
            gc.collect()
    
    # ---------------------------------------------------------------------
    # STEP 8: Sample 3DEP DTM
    # ---------------------------------------------------------------------
    if "dep_dtm" in datasets_to_sample:
        dep_dtm_data = load_checkpoint(site, 'dep_dtm_stats')
        if dep_dtm_data:
            _log("8. ✓ Loaded 3DEP DTM stats from checkpoint")
            results.update(dep_dtm_data['stats'])
        else:
            _log("8. Sampling 3DEP DTM...")
            ds = datasets_to_sample["dep_dtm"]
            xs, ys = _coords_for(ds.crs)
            
            st = zonal_stats(
                ds, xs, ys, FOOTPRINT_RADIUS,
                want_median=True,
                want_std=False,
                percentiles=DTM_PERCENTILES,
                return_valid_fraction=True,
                min_valid_frac=MIN_VALID_FRAC, use_dataset_mask=True
            )
            
            dep_dtm_results = {
                "dep_dtm_mean": st["mean"].astype(np.float32),
                "dep_dtm_median": st["median"].astype(np.float32),
                "dep_dtm_valid_frac": st["valid_frac"].astype(np.float32),
                "dep_dtm": sample_points(ds, xs, ys, use_dataset_mask=True)
            }
            # Add percentiles
            for pct in DTM_PERCENTILES:
                pct_key = f"p{int(pct*100):02d}"
                if pct_key in st:
                    dep_dtm_results[f"dep_dtm_{pct_key}"] = st[pct_key].astype(np.float32)
            
            results.update(dep_dtm_results)
            _log(f"   {_vfrac_summary('3DEP DTM', st['valid_frac'])}")
            save_checkpoint(site, 'dep_dtm_stats', {'stats': dep_dtm_results})
            del st
            gc.collect()
    
    # ---------------------------------------------------------------------
    # STEP 9: Sample P3D DTM
    # ---------------------------------------------------------------------
    if "p3d_dtm" in datasets_to_sample:
        p3d_dtm_data = load_checkpoint(site, 'p3d_dtm_stats')
        if p3d_dtm_data:
            _log("9. ✓ Loaded P3D DTM stats from checkpoint")
            results.update(p3d_dtm_data['stats'])
        else:
            _log("9. Sampling P3D DTM...")
            ds = datasets_to_sample["p3d_dtm"]
            xs, ys = _coords_for(ds.crs)
            
            st = zonal_stats(
                ds, xs, ys, FOOTPRINT_RADIUS,
                want_median=True,
                want_std=False,
                percentiles=DTM_PERCENTILES,
                return_valid_fraction=True,
                min_valid_frac=MIN_VALID_FRAC, use_dataset_mask=True
            )
            
            p3d_dtm_results = {
                "p3d_dtm_mean": st["mean"].astype(np.float32),
                "p3d_dtm_median": st["median"].astype(np.float32),
                "p3d_dtm_valid_frac": st["valid_frac"].astype(np.float32),
                "mhrsi_dtm": sample_points(ds, xs, ys, use_dataset_mask=True)
            }
            # Add percentiles
            for pct in DTM_PERCENTILES:
                pct_key = f"p{int(pct*100):02d}"
                if pct_key in st:
                    p3d_dtm_results[f"p3d_dtm_{pct_key}"] = st[pct_key].astype(np.float32)
            
            results.update(p3d_dtm_results)
            _log(f"   {_vfrac_summary('P3D DTM', st['valid_frac'])}")
            save_checkpoint(site, 'p3d_dtm_stats', {'stats': p3d_dtm_results})
            del st
            gc.collect()
    
    # ---------------------------------------------------------------------
    # STEP 10: CHM aligned metrics (CHUNKED for memory safety)
    # ---------------------------------------------------------------------
    if "als_chm" in datasets_to_sample and "p3d_chm" in datasets_to_sample:
        chm_error_data = load_checkpoint(site, 'chm_aligned_combined')
        if chm_error_data:
            _log("10. ✓ Loaded CHM aligned metrics from checkpoint")
            chm_metrics = chm_error_data
        else:
            _log("10. Computing CHM aligned metrics (CHUNKED)...")
            try:
                ds1 = datasets_to_sample["als_chm"]
                ds2 = datasets_to_sample["p3d_chm"]
                xs1, ys1 = _coords_for(ds1.crs)
                
                chm_metrics = compute_aligned_metrics_chunked(
                    site, ds1, ds2, xs1, ys1,
                    prefix='chm',
                    use_dataset_mask=not BYPASS_CHM_MASK,
                    _log=_log
                )
            except Exception as e:
                _log(f"   ⚠ CHM aligned metrics failed: {e}")
                chm_metrics = None
        
        if chm_metrics:
            results["chm_rmse"] = chm_metrics['rmse']
            results["chm_rmse_valid_frac"] = chm_metrics['vfrac']
            if STRICT_ALIGN_ERROR_MEAN:
                results["error_mean"] = chm_metrics['error_mean']
            if STRICT_ALIGN_ERROR_MEDIAN:
                results["error_median"] = chm_metrics['error_median']
    
    # ---------------------------------------------------------------------
    # STEP 11: DTM aligned metrics (CHUNKED for memory safety)
    # ---------------------------------------------------------------------
    if "dep_dtm" in datasets_to_sample and "p3d_dtm" in datasets_to_sample:
        dtm_error_data = load_checkpoint(site, 'dtm_aligned_combined')
        if dtm_error_data:
            _log("11. ✓ Loaded DTM aligned metrics from checkpoint")
            dtm_metrics = dtm_error_data
        else:
            _log("11. Computing DTM aligned metrics (CHUNKED)...")
            try:
                ds1 = datasets_to_sample["dep_dtm"]
                ds2 = datasets_to_sample["p3d_dtm"]
                xs1, ys1 = _coords_for(ds1.crs)
                
                dtm_metrics = compute_aligned_metrics_chunked(
                    site, ds1, ds2, xs1, ys1,
                    prefix='dtm',
                    use_dataset_mask=True,
                    _log=_log
                )
            except Exception as e:
                _log(f"   ⚠ DTM aligned metrics failed: {e}")
                dtm_metrics = None
        
        if dtm_metrics:
            results["dtm_rmse"] = dtm_metrics['rmse']
            results["dtm_rmse_valid_frac"] = dtm_metrics['vfrac']
            if STRICT_ALIGN_ERROR_MEAN:
                results["dtm_error_mean"] = dtm_metrics['error_mean']
            if STRICT_ALIGN_ERROR_MEDIAN:
                results["dtm_error_median"] = dtm_metrics['error_median']
    
    # Close all datasets to free memory
    _log("Closing raster datasets...")
    for ds in datasets_to_sample.values():
        try:
            ds.close()
        except:
            pass
    del datasets_to_sample
    gc.collect()
    
    # ---------------------------------------------------------------------
    # STEP 12: Assemble final dataframe
    # ---------------------------------------------------------------------
    final_df_data = load_checkpoint(site, 'final_dataframe')
    if final_df_data:
        _log("12. ✓ Loaded final dataframe from checkpoint")
        df_final = final_df_data['df']
    else:
        _log("12. Assembling final dataframe...")
        
        # Start with GEDI attributes
        base_cols = ["shot_number", "cover", "pft", "wsci", "rh_98", "delta_time", 
                     "elev_low", "ecoregion"]
        available_base_cols = [c for c in base_cols if c in gdf_orig.columns]
        df_final = gdf_orig[available_base_cols].copy()
        
        # Add site column
        df_final.insert(0, "site", site)
        
        # Add all sampled results
        for key, arr in results.items():
            if isinstance(arr, np.ndarray):
                if len(arr) != len(df_final):
                    _log(f"   ⚠ Length mismatch for {key}: {len(arr)} != {len(df_final)}")
                    continue
                df_final[key] = arr.astype(np.float32) if arr.dtype == np.float64 else arr
        
        # Final column order
        final_cols_ordered = [
            # Context
            "site","cover","pft","wsci","rh_98","delta_time","shot_number","elev_low","ecoregion",
            # Terrain
            "slope_mean","slope_sd","slope_valid_frac",
            "aspect_sin_mean","aspect_cos_mean","aspect_valid_frac",
            # CHM summary
            "als_chm_mean","p3d_chm_mean","error_mean",
            "als_chm_median","p3d_chm_median","error_median",
            "als_chm_sd","p3d_chm_sd",
            "als_chm_p25","als_chm_p75","als_chm_p90",
            "p3d_chm_p25","p3d_chm_p75","p3d_chm_p90",
            "chm_rmse","als_chm_valid_frac","p3d_chm_valid_frac","chm_rmse_valid_frac",
            # DTM summary
            "dep_dtm_mean","p3d_dtm_mean","dtm_error_mean",
            "dep_dtm_median","p3d_dtm_median","dtm_error_median",
            "dep_dtm_p25","dep_dtm_p75","dep_dtm_p90",
            "p3d_dtm_p25","p3d_dtm_p75","p3d_dtm_p90",
            "dtm_rmse","dep_dtm_valid_frac","p3d_dtm_valid_frac","dtm_rmse_valid_frac",
            # Point samples
            "dep_dtm","mhrsi_dtm","als_chm","p3d_chm",
        ]
        
        # Fill missing columns with NaN, reorder
        for col in final_cols_ordered:
            if col not in df_final.columns:
                df_final[col] = np.nan
        
        # Keep only ordered columns that exist
        df_final = df_final[[c for c in final_cols_ordered if c in df_final.columns]]
        
        _log(f"   ✓ Final dataframe: {len(df_final):,} rows × {len(df_final.columns)} columns")
        _log(f"   Memory usage: {df_final.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
        save_checkpoint(site, 'final_dataframe', {'df': df_final})
    
    _log(f"="*80)
    _log(f"SITE {site} COMPLETE")
    _log(f"="*80)
    
    return site, df_final, "\n".join(log_lines)


# ---------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Sample GEDI data with chunked processing for large sites"
    )
    parser.add_argument("--site", type=int, required=True,
                       help="Site number to process")
    parser.add_argument("--cleanup", action="store_true",
                       help="Delete checkpoint files after successful completion")
    parser.add_argument("--output-dir", type=str, default=None,
                       help="Output directory for CSV (default: current directory)")
    
    args = parser.parse_args()
    
    try:
        site_id, df, log = process_site(args.site)
        
        if df is not None:
            if args.output_dir:
                output_dir = Path(args.output_dir)
                output_dir.mkdir(parents=True, exist_ok=True)
            else:
                output_dir = Path(".")
            
            output_file = output_dir / f"site_{site_id:02d}_sampled.csv"
            df.to_csv(output_file, index=False, float_format="%.4f")
            print(f"\n✓ Saved {len(df):,} rows to {output_file}")
            
            # Cleanup if requested
            if args.cleanup:
                import shutil
                ckpt_dir = get_checkpoint_dir(args.site)
                if ckpt_dir.exists():
                    shutil.rmtree(ckpt_dir)
                    print(f"🗑️  Cleaned up checkpoint directory: {ckpt_dir}")
        else:
            print("\n⚠ No data returned")
            
    except Exception as e:
        print(f"\n{'!'*80}")
        print(f"ERROR: {e}")
        print(f"{'!'*80}")
        import traceback
        traceback.print_exc()
