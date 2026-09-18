#!/usr/bin/env python3
"""
resample_dtm_slope_aspect.py
============================

Efficiently re-sample 3DEP-derived data (DTM, slope, aspect) while preserving CHM data.

This script:
1. Loads existing checkpoint CSV (with buggy DTM/slope/aspect but good CHM data)
2. Re-samples 3DEP DTM, P3D DTM, slope, and aspect using the FIXED helper_speed.py
3. Replaces DTM/slope/aspect-related columns in the dataframe
4. Recalculates ALL filter/coverage statistics properly
5. Saves corrected checkpoint

Time estimate: 5-8 hours for site 1 (vs 48 hours for full resample)

Why slope/aspect need resampling:
- Slope and aspect rasters are derived from 3DEP DTM
- They have the same CRS, resolution, and extent
- Therefore they have the same alignment bug that affected 3DEP DTM
"""

from __future__ import annotations
import argparse
import pickle
from pathlib import Path
from typing import List, Optional, Dict
import pandas as pd
import geopandas as gpd
import numpy as np
import rasterio

# Import the FIXED helper functions
from helper_speed import (
    build_vrt,
    zonal_stats,
    zonal_aligned_pair_metrics,
    sample_points,
)

# ---------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------
CHECKPOINT_DIR = Path(r"<LOCAL_DATA_ROOT>\sampled_data_FINAL_by_site")
GEDI_BASE = Path(r"<LOCAL_DATA_ROOT>\gedi\sites")
SLOPE_BASE = Path(r"D:\slope")
ASPECT_BASE = Path(r"D:\aspect")
DEP_BASE = Path(r"D:\3dep")
MHRSI_BASE = Path(r"E:\p3d")

FOOTPRINT_RADIUS = 12.5  # meters
MIN_VALID_FRAC = None
DTM_PERCENTILES = (0.90, 0.75, 0.25)

# ---------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------
def _rgp(root: Path, *sub: str) -> List[str]:
    """Recursively glob for .tif files."""
    base = root.joinpath(*sub)
    if not base.exists():
        return []
    return [str(p) for p in base.rglob("*.tif")]


# ---------------------------------------------------------------------
# CHECKPOINT HELPERS
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


def calculate_filter_statistics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Recalculate ALL filter/coverage statistics properly.
    
    This is what the previous patching script was missing!
    """
    total = len(df)
    print(f"\nCalculating filter statistics for {total:,} observations...")
    
    # Define filter columns and their criteria
    filters = {
        # P3D CHM filters
        'p3d_chm_valid_frac': lambda df: df['p3d_chm_valid_frac'].notna(),
        'p3d_chm_valid_frac_gt05': lambda df: (df['p3d_chm_valid_frac'] >= 0.5),
        
        # P3D DTM filters
        'p3d_dtm_valid_frac': lambda df: df['p3d_dtm_valid_frac'].notna(),
        'p3d_dtm_valid_frac_gt05': lambda df: (df['p3d_dtm_valid_frac'] >= 0.5),
        
        # ALS CHM filters
        'als_chm_valid_frac': lambda df: df['als_chm_valid_frac'].notna(),
        'als_chm_valid_frac_gt05': lambda df: (df['als_chm_valid_frac'] >= 0.5),
        
        # 3DEP DTM filters
        'dep_dtm_valid_frac': lambda df: df['dep_dtm_valid_frac'].notna(),
        'dep_dtm_valid_frac_gt05': lambda df: (df['dep_dtm_valid_frac'] >= 0.5),
        
        # Slope filter
        'slope_valid_frac': lambda df: df['slope_valid_frac'].notna(),
        'slope_valid_frac_gt05': lambda df: (df['slope_valid_frac'] >= 0.5),
        
        # ALS CHM p90 filters
        'als_chm_p90_gt2m': lambda df: (df['als_chm_p90'] >= 2),
    }
    
    print("\nFilter survival rates:")
    print("=" * 80)
    for name, filter_fn in filters.items():
        mask = filter_fn(df)
        count = mask.sum()
        pct = 100 * count / total if total > 0 else 0
        print(f"{name:.<40} {count:>10,} ({pct:>5.1f}%)")
    
    # Combined filters
    print("\n" + "=" * 80)
    
    # All CHM filters combined (for CHM error analysis)
    chm_mask = (
        (df['p3d_chm_valid_frac'] >= 0.5) &
        (df['als_chm_valid_frac'] >= 0.5) &
        (df['slope_valid_frac'] >= 0.5) &
        (df['als_chm_p90'] >= 2)
    )
    chm_count = chm_mask.sum()
    print(f"{'All CHM filters combined':.<40} {chm_count:>10,} ({100*chm_count/total:>5.1f}%)")
    
    # All DTM filters combined (for DTM error analysis)
    dtm_mask = (
        (df['p3d_dtm_valid_frac'] >= 0.5) &
        (df['dep_dtm_valid_frac'] >= 0.5) &
        (df['slope_valid_frac'] >= 0.5)
    )
    dtm_count = dtm_mask.sum()
    print(f"{'All DTM filters combined':.<40} {dtm_count:>10,} ({100*dtm_count/total:>5.1f}%)")
    print("=" * 80)
    
    return df


def resample_3dep_derived_data(site: int, checkpoint_csv: Path, output_csv: Path):
    """
    Re-sample 3DEP-derived data (DTM, slope, aspect) for a site.
    
    Args:
        site: Site number
        checkpoint_csv: Path to existing checkpoint CSV with all data
        output_csv: Path to save corrected CSV
    """
    print(f"\n{'='*80}")
    print(f"RE-SAMPLING 3DEP-DERIVED DATA FOR SITE {site}")
    print(f"(DTM, SLOPE, ASPECT)")
    print(f"{'='*80}\n")
    
    # Check for existing checkpoints
    existing_ckpts = list_checkpoints(site)
    if existing_ckpts:
        print(f"📂 Found existing checkpoints: {', '.join(existing_ckpts)}")
        print(f"   Will resume from last checkpoint\n")
    else:
        print(f"📂 No existing checkpoints found - starting fresh\n")
    
    # 1. Load existing checkpoint
    print(f"1. Loading existing checkpoint: {checkpoint_csv.name}")
    if not checkpoint_csv.exists():
        raise FileNotFoundError(f"Checkpoint not found: {checkpoint_csv}")
    
    df = pd.read_csv(checkpoint_csv)
    n_orig = len(df)
    print(f"   Loaded {n_orig:,} observations")
    
    # Optimize memory: convert float64 columns to float32 where possible
    print(f"   Optimizing dtypes to reduce memory...")
    for col in df.select_dtypes(include=['float64']).columns:
        df[col] = df[col].astype(np.float32)
    
    # Check what we have
    cols_to_replace = [c for c in df.columns if any(x in c.lower() for x in ['dtm', 'slope', 'aspect'])]
    print(f"   Found {len(cols_to_replace)} columns to replace: {cols_to_replace[:8]}...")
    print(f"   Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
    
    # 2. Load GEDI data to get coordinates
    print(f"\n2. Loading GEDI geometries...")
    gedi_fp = GEDI_BASE / str(site) / f"GEDI_site{site}_hq_ALL.gpkg"
    if not gedi_fp.exists():
        raise FileNotFoundError(f"GEDI file not found: {gedi_fp}")
    
    # Read just geometry and shot_number to match with df
    gedi_gdf = gpd.read_file(gedi_fp)
    print(f"   Loaded {len(gedi_gdf):,} GEDI footprints")
    
    # Ensure shot_number column exists for matching
    if 'shot_number' not in df.columns:
        raise ValueError("shot_number column not found in checkpoint CSV!")
    
    # Store shot_numbers for later, then we can delete gedi_gdf after coordinate extraction
    shot_numbers = gedi_gdf['shot_number'].values
    
    # 3. Find rasters
    print(f"\n3. Finding rasters...")
    slope_paths = _rgp(SLOPE_BASE, f"polygon_{site}")
    aspect_paths = _rgp(ASPECT_BASE, f"polygon_{site}")
    dep_dtm_paths = _rgp(DEP_BASE, f"polygon_{site}")
    p3d_dtm_paths = _rgp(MHRSI_BASE, str(site), "vricon_raster_50cm", "dtm")
    
    if not slope_paths:
        raise RuntimeError(f"No slope rasters found for site {site}")
    if not aspect_paths:
        raise RuntimeError(f"No aspect rasters found for site {site}")
    if not dep_dtm_paths:
        raise RuntimeError(f"No 3DEP DTM rasters found for site {site}")
    if not p3d_dtm_paths:
        raise RuntimeError(f"No P3D DTM rasters found for site {site}")
    
    print(f"   Slope:    {len(slope_paths)} tiles")
    print(f"   Aspect:   {len(aspect_paths)} tiles")
    print(f"   3DEP DTM: {len(dep_dtm_paths)} tiles")
    print(f"   P3D DTM:  {len(p3d_dtm_paths)} tiles")
    
    # 4. Open datasets
    print(f"\n4. Opening datasets...")
    try:
        slope_ds = build_vrt(slope_paths)
        print(f"   Slope:    CRS={slope_ds.crs.to_string()}")
        aspect_ds = build_vrt(aspect_paths)
        print(f"   Aspect:   CRS={aspect_ds.crs.to_string()}")
        dep_ds = build_vrt(dep_dtm_paths)
        print(f"   3DEP DTM: CRS={dep_ds.crs.to_string()}")
        p3d_ds = build_vrt(p3d_dtm_paths)
        print(f"   P3D DTM:  CRS={p3d_ds.crs.to_string()}")
    except Exception as e:
        raise RuntimeError(f"Failed to open datasets: {e}")
    
    # 5. Reproject coordinates to each dataset's CRS
    print(f"\n5. Reprojecting coordinates...")
    # Slope and aspect have same CRS as 3DEP (they're derived from it)
    gedi_gdf_3dep = gedi_gdf.to_crs(dep_ds.crs)
    xs_3dep = gedi_gdf_3dep.geometry.x.values.astype(np.float64)
    ys_3dep = gedi_gdf_3dep.geometry.y.values.astype(np.float64)
    
    gedi_gdf_p3d = gedi_gdf.to_crs(p3d_ds.crs)
    xs_p3d = gedi_gdf_p3d.geometry.x.values.astype(np.float64)
    ys_p3d = gedi_gdf_p3d.geometry.y.values.astype(np.float64)
    
    print(f"   Coordinates ready for {len(xs_3dep):,} footprints")
    
    # CRITICAL: Delete gedi GeoDataFrames to free memory before sampling
    del gedi_gdf, gedi_gdf_3dep, gedi_gdf_p3d
    import gc
    gc.collect()
    
    # 6. Sample slope
    print(f"\n6. Re-sampling slope...")
    slope_data = load_checkpoint(site, 'slope')
    if slope_data:
        print(f"   ✓ Loaded from checkpoint (skipping sampling)")
        slope_stats = slope_data['stats']
    else:
        print(f"   Starting fresh (this may take a while)...")
        try:
            # Check if slope nodata is 0 (common issue)
            ignore_nodata = bool(slope_ds.nodata is not None and float(slope_ds.nodata) == 0.0)
            if ignore_nodata:
                print(f"   Note: Slope nodata==0 detected → ignoring nodata value")
            
            slope_stats = zonal_stats(
                slope_ds, xs_3dep, ys_3dep, FOOTPRINT_RADIUS,
                want_std=True,  # slope gets SD
                return_valid_fraction=True,
                min_valid_frac=MIN_VALID_FRAC,
                use_dataset_mask=True,
                ignore_nodata_value=ignore_nodata
            )
            
            valid_pct = (slope_stats['valid_frac'] >= 0.5).mean() * 100
            print(f"   ✓ Slope sampled: {valid_pct:.1f}% footprints with >=50% coverage")
            
            # Save checkpoint
            save_checkpoint(site, 'slope', {'stats': slope_stats})
            
        except Exception as e:
            raise RuntimeError(f"Failed to sample slope: {e}")
    
    # 7. Sample aspect  
    print(f"\n7. Re-sampling aspect...")
    aspect_data = load_checkpoint(site, 'aspect')
    if aspect_data:
        print(f"   ✓ Loaded from checkpoint (skipping sampling)")
        sin_stats = aspect_data['sin_stats']
        cos_stats = aspect_data['cos_stats']
    else:
        print(f"   Starting fresh...")
        try:
            # Aspect needs sin/cos transformation
            sin_stats = zonal_stats(
                aspect_ds, xs_3dep, ys_3dep, FOOTPRINT_RADIUS,
                transform_fn=lambda v: np.sin(np.deg2rad(v)),
                return_valid_fraction=True,
                min_valid_frac=MIN_VALID_FRAC,
                use_dataset_mask=True
            )
            cos_stats = zonal_stats(
                aspect_ds, xs_3dep, ys_3dep, FOOTPRINT_RADIUS,
                transform_fn=lambda v: np.cos(np.deg2rad(v)),
                return_valid_fraction=False,
                min_valid_frac=MIN_VALID_FRAC,
                use_dataset_mask=True
            )
            
            valid_pct = (sin_stats['valid_frac'] >= 0.5).mean() * 100
            print(f"   ✓ Aspect sampled: {valid_pct:.1f}% footprints with >=50% coverage")
            
            # Save checkpoint
            save_checkpoint(site, 'aspect', {
                'sin_stats': sin_stats,
                'cos_stats': cos_stats
            })
            
        except Exception as e:
            raise RuntimeError(f"Failed to sample aspect: {e}")
    
    # 8. Sample 3DEP DTM
    print(f"\n8. Re-sampling 3DEP DTM...")
    dep_data = load_checkpoint(site, 'dep_dtm')
    if dep_data:
        print(f"   ✓ Loaded from checkpoint (skipping sampling)")
        dep_stats = dep_data['stats']
        dep_dtm_center = dep_data['center']
    else:
        print(f"   Starting fresh...")
        try:
            dep_stats = zonal_stats(
                dep_ds, xs_3dep, ys_3dep, FOOTPRINT_RADIUS,
                want_median=True,
                want_std=False,
                percentiles=DTM_PERCENTILES,
                return_valid_fraction=True,
                min_valid_frac=MIN_VALID_FRAC,
                use_dataset_mask=True
            )
            
            # Sample center point too
            dep_dtm_center = sample_points(dep_ds, xs_3dep, ys_3dep, use_dataset_mask=True)
            
            valid_pct = np.isfinite(dep_stats['valid_frac']).mean() * 100
            print(f"   ✓ 3DEP DTM sampled: {valid_pct:.1f}% footprints with valid data")
            
            # Save checkpoint
            save_checkpoint(site, 'dep_dtm', {
                'stats': dep_stats,
                'center': dep_dtm_center
            })
            
        except Exception as e:
            raise RuntimeError(f"Failed to sample 3DEP DTM: {e}")
    
    # 9. Sample P3D DTM
    print(f"\n9. Re-sampling P3D DTM...")
    p3d_data = load_checkpoint(site, 'p3d_dtm')
    if p3d_data:
        print(f"   ✓ Loaded from checkpoint (skipping sampling)")
        p3d_stats = p3d_data['stats']
        p3d_dtm_center = p3d_data['center']
    else:
        print(f"   Starting fresh...")
        try:
            p3d_stats = zonal_stats(
                p3d_ds, xs_p3d, ys_p3d, FOOTPRINT_RADIUS,
                want_median=True,
                want_std=False,
                percentiles=DTM_PERCENTILES,
                return_valid_fraction=True,
                min_valid_frac=MIN_VALID_FRAC,
                use_dataset_mask=True
            )
            
            # Sample center point too
            p3d_dtm_center = sample_points(p3d_ds, xs_p3d, ys_p3d, use_dataset_mask=True)
            
            valid_pct = np.isfinite(p3d_stats['valid_frac']).mean() * 100
            print(f"   ✓ P3D DTM sampled: {valid_pct:.1f}% footprints with valid data")
            
            # Save checkpoint
            save_checkpoint(site, 'p3d_dtm', {
                'stats': p3d_stats,
                'center': p3d_dtm_center
            })
            
        except Exception as e:
            raise RuntimeError(f"Failed to sample P3D DTM: {e}")
    
    # 10. Compute aligned pair metrics (CHUNKED to avoid memory crash!)
    print(f"\n10. Computing aligned DTM error metrics (with FIXED code)...")
    aligned_data = load_checkpoint(site, 'aligned_dtm')
    if aligned_data:
        print(f"   ✓ Loaded from checkpoint (skipping computation)")
        rmse = aligned_data['rmse']
        vfrac = aligned_data['vfrac']
        dtm_error_mean = aligned_data['error_mean']
        dtm_error_median = aligned_data['error_median']
    else:
        print(f"   Processing in chunks to avoid memory crash...")
        try:
            # Process in chunks of 100K footprints
            CHUNK_SIZE = 100000
            n_footprints = len(xs_3dep)
            n_chunks = (n_footprints + CHUNK_SIZE - 1) // CHUNK_SIZE
            
            print(f"   Total footprints: {n_footprints:,}")
            print(f"   Chunk size: {CHUNK_SIZE:,}")
            print(f"   Number of chunks: {n_chunks}")
            
            # Check which chunks are already done
            ckpt_dir = get_checkpoint_dir(site)
            existing_chunks = set()
            if ckpt_dir.exists():
                for p in ckpt_dir.glob("aligned_chunk_*.pkl"):
                    try:
                        chunk_num = int(p.stem.split('_')[-1])
                        existing_chunks.add(chunk_num)
                    except:
                        pass
            
            if existing_chunks:
                print(f"   Found {len(existing_chunks)} existing chunk checkpoint(s)")
                if existing_chunks != set(range(n_chunks)):
                    missing = set(range(n_chunks)) - existing_chunks
                    print(f"   Will resume from chunk {min(missing) + 1}")
            
            # Pre-allocate output arrays
            rmse = np.full(n_footprints, np.nan, dtype=np.float32)
            vfrac = np.full(n_footprints, np.nan, dtype=np.float32)
            mu1_all = np.full(n_footprints, np.nan, dtype=np.float32)
            mu2_all = np.full(n_footprints, np.nan, dtype=np.float32)
            med1_all = np.full(n_footprints, np.nan, dtype=np.float32)
            med2_all = np.full(n_footprints, np.nan, dtype=np.float32)
            
            # Track failed chunks
            failed_chunks = []
            
            # Process each chunk
            for i in range(n_chunks):
                start_idx = i * CHUNK_SIZE
                end_idx = min((i + 1) * CHUNK_SIZE, n_footprints)
                
                # Check if this chunk is already done
                if i in existing_chunks:
                    print(f"   Chunk {i+1}/{n_chunks}: ✓ Loading from checkpoint... ", end='', flush=True)
                    chunk_data = load_checkpoint(site, f'aligned_chunk_{i}')
                    rmse[start_idx:end_idx] = chunk_data['rmse']
                    vfrac[start_idx:end_idx] = chunk_data['vfrac']
                    mu1_all[start_idx:end_idx] = chunk_data['mu1']
                    mu2_all[start_idx:end_idx] = chunk_data['mu2']
                    med1_all[start_idx:end_idx] = chunk_data['med1']
                    med2_all[start_idx:end_idx] = chunk_data['med2']
                    print(f"({(i+1)/n_chunks*100:.1f}% complete)")
                    continue
                
                print(f"   Chunk {i+1}/{n_chunks}: footprints {start_idx:,} to {end_idx:,}...", end='', flush=True)
                
                try:
                    # Get coordinates for this chunk
                    xs_chunk = xs_3dep[start_idx:end_idx]
                    ys_chunk = ys_3dep[start_idx:end_idx]
                    
                    # Compute aligned metrics for this chunk
                    rmse_chunk, mu1, mu2, vfrac_chunk, med1, med2 = zonal_aligned_pair_metrics(
                        dep_ds, p3d_ds, xs_chunk, ys_chunk, FOOTPRINT_RADIUS,
                        want_median=True,
                        return_valid_fraction=True,
                        min_valid_frac=MIN_VALID_FRAC,
                        use_dataset_mask=True
                    )
                    
                    # Store results in main arrays
                    rmse[start_idx:end_idx] = rmse_chunk.astype(np.float32)
                    vfrac[start_idx:end_idx] = vfrac_chunk.astype(np.float32)
                    mu1_all[start_idx:end_idx] = mu1.astype(np.float32)
                    mu2_all[start_idx:end_idx] = mu2.astype(np.float32)
                    med1_all[start_idx:end_idx] = med1.astype(np.float32)
                    med2_all[start_idx:end_idx] = med2.astype(np.float32)
                    
                    # Save chunk checkpoint IMMEDIATELY
                    save_checkpoint(site, f'aligned_chunk_{i}', {
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
                    import gc
                    gc.collect()
                    
                except Exception as chunk_error:
                    # Handle corrupted tiles - log and continue with NaN
                    print(f" ⚠️  FAILED (corrupted data)")
                    print(f"      Error: {str(chunk_error)[:100]}")
                    
                    # Arrays already initialized with NaN, so chunk will remain NaN
                    # Log this failure
                    failed_chunks.append(i + 1)
                    
                    # Save checkpoint with NaN values so we don't retry this chunk
                    chunk_size = end_idx - start_idx
                    save_checkpoint(site, f'aligned_chunk_{i}', {
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
                        import gc
                        gc.collect()
                    except:
                        pass
            
            # Compute error metrics from aligned means
            dtm_error_mean = (mu2_all - mu1_all).astype(np.float32)
            dtm_error_median = (med2_all - med1_all).astype(np.float32)
            
            # Clean up intermediate arrays
            del mu1_all, mu2_all, med1_all, med2_all
            gc.collect()
            
            # Report any failed chunks
            if failed_chunks:
                print(f"\n   ⚠️  WARNING: {len(failed_chunks)} chunk(s) failed due to corrupted data:")
                print(f"      Chunks: {failed_chunks}")
                print(f"      These chunks are filled with NaN and will be excluded from analysis")
            
            valid_errors = np.isfinite(dtm_error_mean)
            if valid_errors.sum() > 0:
                mean_err = dtm_error_mean[valid_errors].mean()
                median_err = np.median(dtm_error_mean[valid_errors])
                rmse_val = np.sqrt((dtm_error_mean[valid_errors]**2).mean())
                print(f"\n   ✓ All chunks complete! DTM errors computed:")
                print(f"      Mean:   {mean_err:.3f} m")
                print(f"      Median: {median_err:.3f} m")
                print(f"      RMSE:   {rmse_val:.3f} m")
                print(f"      N valid: {valid_errors.sum():,} / {len(dtm_error_mean):,} ({100*valid_errors.sum()/len(dtm_error_mean):.1f}%)")
            else:
                print(f"\n   ⚠️  WARNING: No valid DTM errors computed (all data corrupted)")

            
            # Save final combined checkpoint
            print(f"   Saving final combined checkpoint...")
            save_checkpoint(site, 'aligned_dtm', {
                'rmse': rmse,
                'vfrac': vfrac,
                'error_mean': dtm_error_mean,
                'error_median': dtm_error_median
            })
            
            # Clean up individual chunk checkpoints
            print(f"   Cleaning up {n_chunks} chunk checkpoints...")
            for i in range(n_chunks):
                chunk_file = get_checkpoint_dir(site) / f'aligned_chunk_{i}.pkl'
                if chunk_file.exists():
                    try:
                        chunk_file.unlink()
                    except:
                        pass  # Ignore errors during cleanup
            
        except Exception as e:
            import traceback
            print(f"\n   ERROR: {e}")
            traceback.print_exc()
            raise RuntimeError(f"Failed to compute aligned DTM metrics: {e}")
    
    # CRITICAL: Close all datasets to free memory before merge
    print(f"\n   Closing raster datasets to free memory...")
    for ds in [slope_ds, aspect_ds, dep_ds, p3d_ds]:
        try:
            ds.close()
        except:
            pass
    del slope_ds, aspect_ds, dep_ds, p3d_ds
    import gc
    gc.collect()
    
    # 11. Replace columns in dataframe (MEMORY-OPTIMIZED)
    print(f"\n11. Replacing DTM/slope/aspect columns in dataframe...")
    print(f"   Using memory-efficient chunked processing...")
    
    # First, check if we have a merge checkpoint
    merge_data = load_checkpoint(site, 'merged_dataframe')
    if merge_data:
        print(f"   ✓ Loaded merged dataframe from checkpoint (skipping merge)")
        df_corrected = merge_data['df']
    else:
        print(f"   Starting merge (this is memory-intensive)...")
        
        # CRITICAL: Convert all arrays to float32 to save memory
        # float64 uses 8 bytes per value, float32 uses 4 bytes (50% savings!)
        print(f"   Optimizing data types...")
        slope_stats_f32 = {k: v.astype(np.float32) if isinstance(v, np.ndarray) else v 
                           for k, v in slope_stats.items()}
        sin_stats_f32 = {k: v.astype(np.float32) if isinstance(v, np.ndarray) else v 
                         for k, v in sin_stats.items()}
        cos_stats_f32 = {k: v.astype(np.float32) if isinstance(v, np.ndarray) else v 
                         for k, v in cos_stats.items()}
        dep_stats_f32 = {k: v.astype(np.float32) if isinstance(v, np.ndarray) else v 
                         for k, v in dep_stats.items()}
        p3d_stats_f32 = {k: v.astype(np.float32) if isinstance(v, np.ndarray) else v 
                         for k, v in p3d_stats.items()}
        
        # Create new data dictionary (still in memory but optimized)
        new_data_dict = {
            'shot_number': shot_numbers,
            # Slope stats
            'slope_mean': slope_stats_f32['mean'],
            'slope_sd': slope_stats_f32['std'],
            'slope_valid_frac': slope_stats_f32['valid_frac'],
            # Aspect stats
            'aspect_sin_mean': sin_stats_f32['mean'],
            'aspect_cos_mean': cos_stats_f32['mean'],
            'aspect_valid_frac': sin_stats_f32['valid_frac'],
            # 3DEP DTM stats
            'dep_dtm_mean': dep_stats_f32['mean'],
            'dep_dtm_median': dep_stats_f32['median'],
            'dep_dtm_valid_frac': dep_stats_f32['valid_frac'],
            'dep_dtm_p25': dep_stats_f32.get('p25', np.full(len(shot_numbers), np.nan, dtype=np.float32)),
            'dep_dtm_p75': dep_stats_f32.get('p75', np.full(len(shot_numbers), np.nan, dtype=np.float32)),
            'dep_dtm_p90': dep_stats_f32.get('p90', np.full(len(shot_numbers), np.nan, dtype=np.float32)),
            'dep_dtm': dep_dtm_center.astype(np.float32),
            # P3D DTM stats
            'p3d_dtm_mean': p3d_stats_f32['mean'],
            'p3d_dtm_median': p3d_stats_f32['median'],
            'p3d_dtm_valid_frac': p3d_stats_f32['valid_frac'],
            'p3d_dtm_p25': p3d_stats_f32.get('p25', np.full(len(shot_numbers), np.nan, dtype=np.float32)),
            'p3d_dtm_p75': p3d_stats_f32.get('p75', np.full(len(shot_numbers), np.nan, dtype=np.float32)),
            'p3d_dtm_p90': p3d_stats_f32.get('p90', np.full(len(shot_numbers), np.nan, dtype=np.float32)),
            'mhrsi_dtm': p3d_dtm_center.astype(np.float32),
            # Aligned error metrics
            'dtm_rmse': rmse.astype(np.float32),
            'dtm_rmse_valid_frac': vfrac.astype(np.float32),
            'dtm_error_mean': dtm_error_mean,
            'dtm_error_median': dtm_error_median,
        }
        
        # Drop old columns from original df
        old_cols = [c for c in df.columns if any(x in c for x in 
                    ['slope_', 'aspect_', 'dep_dtm', 'p3d_dtm', 'mhrsi_dtm', 'dtm_rmse', 'dtm_error'])]
        print(f"   Dropping {len(old_cols)} old columns...")
        df_reduced = df.drop(columns=old_cols, errors='ignore')
        
        # CRITICAL: Set shot_number as index on BOTH dataframes for memory-efficient assignment
        print(f"   Setting up index-based assignment (avoiding merge)...")
        if 'shot_number' in df_reduced.columns:
            df_reduced = df_reduced.set_index('shot_number')
        
        # Assign new columns directly by index (much more memory-efficient than merge!)
        print(f"   Assigning new columns...")
        for col_name, col_data in new_data_dict.items():
            if col_name != 'shot_number':
                df_reduced[col_name] = col_data
        
        # Reset index to get shot_number back as a column
        df_corrected = df_reduced.reset_index()
        
        # Explicit cleanup
        del df_reduced, new_data_dict
        del slope_stats_f32, sin_stats_f32, cos_stats_f32, dep_stats_f32, p3d_stats_f32
        import gc
        gc.collect()
        
        print(f"   ✓ Columns replaced: {len(df_corrected):,} rows")
        print(f"   Final memory usage: {df_corrected.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
        
        # Save merge checkpoint BEFORE filter statistics
        print(f"   Saving merge checkpoint (in case filter stats fail)...")
        save_checkpoint(site, 'merged_dataframe', {'df': df_corrected})
    
    if len(df_corrected) != n_orig:
        raise RuntimeError(f"Row count mismatch after merge: {len(df_corrected)} != {n_orig}")
    
    print(f"   ✓ Columns replaced: {len(df_corrected)} rows")
    
    # 12. Recalculate filter statistics
    print(f"\n12. Recalculating filter statistics...")
    df_corrected = calculate_filter_statistics(df_corrected)
    
    # 13. Save corrected CSV
    print(f"\n13. Saving corrected checkpoint...")
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df_corrected.to_csv(output_csv, index=False, float_format="%.4f")
    
    print(f"   ✓ Saved {len(df_corrected):,} observations to:")
    print(f"      {output_csv}")
    
    # 14. Final summary
    print(f"\n{'='*80}")
    print(f"SITE {site} 3DEP-DERIVED DATA RE-SAMPLING COMPLETE")
    print(f"{'='*80}\n")
    
    print("Next steps:")
    print("  1. Verify the filter statistics look correct")
    print("     - Slope/aspect/DTM filters should all be >0%")
    print("  2. Check DTM error statistics are reasonable (~1-2m RMSE)")
    print("  3. If good, replace the old checkpoint file")
    print(f"  4. Run resume_driver.py --resume to continue with other sites")
    print(f"\n💡 Checkpoint files saved in: {get_checkpoint_dir(site)}")
    print(f"   (You can delete these after verifying the output)\n")
    
    return df_corrected


def main():
    parser = argparse.ArgumentParser(
        description="Re-sample 3DEP-derived data (DTM, slope, aspect) while keeping CHM data"
    )
    parser.add_argument("--site", type=int, required=True,
                       help="Site number to resample (e.g., 1)")
    parser.add_argument("--input", type=str, default=None,
                       help="Input checkpoint CSV (default: auto-detect from CHECKPOINT_DIR)")
    parser.add_argument("--output", type=str, default=None,
                       help="Output CSV path (default: input_path with _3DEP_FIXED suffix)")
    parser.add_argument("--cleanup", action="store_true",
                       help="Delete checkpoint files after successful completion")
    
    args = parser.parse_args()
    
    # Determine input path
    if args.input:
        input_csv = Path(args.input)
    else:
        input_csv = CHECKPOINT_DIR / f"site_{args.site:02d}.csv"
    
    if not input_csv.exists():
        print(f"ERROR: Input checkpoint not found: {input_csv}")
        print(f"\nPlease ensure site {args.site} has been processed previously.")
        print(f"Looking in: {CHECKPOINT_DIR}")
        return 1
    
    # Determine output path
    if args.output:
        output_csv = Path(args.output)
    else:
        output_csv = input_csv.with_name(input_csv.stem + "_3DEP_FIXED.csv")
    
    print(f"\nConfiguration:")
    print(f"  Site:   {args.site}")
    print(f"  Input:  {input_csv}")
    print(f"  Output: {output_csv}")
    print()
    
    try:
        resample_3dep_derived_data(args.site, input_csv, output_csv)
        
        # Clean up checkpoints if requested
        if args.cleanup:
            import shutil
            ckpt_dir = get_checkpoint_dir(args.site)
            if ckpt_dir.exists():
                shutil.rmtree(ckpt_dir)
                print(f"\n🗑️  Cleaned up checkpoint directory: {ckpt_dir}")
        
        return 0
    except Exception as e:
        print(f"\n{'!'*80}")
        print(f"ERROR: {e}")
        print(f"{'!'*80}\n")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
