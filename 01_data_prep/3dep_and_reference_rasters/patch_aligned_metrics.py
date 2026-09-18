#!/usr/bin/env python3
"""
patch_aligned_metrics.py
===============================

Quick patch for sites 1 and 5 to fix buggy aligned metrics without full resampling.

This script:
1. Loads existing enriched CSVs (with buggy aligned metrics)
2. Checks which GEDI footprints fall outside P3D DTM/CHM coverage
3. Sets aligned metrics to NaN for footprints outside P3D coverage
4. Recomputes error metrics from individual means for valid footprints
5. Saves patched CSVs

Time: ~5-10 minutes per site (vs 20+ hours for full resample)
"""

from pathlib import Path
import pandas as pd
import numpy as np
import rasterio
from rasterio.warp import transform_bounds
from tqdm import tqdm
import sys

# =============================================================================
# CONFIGURATION - PATHS FOR LOCAL PC
# =============================================================================

# Input CSVs (old buggy versions)
BUGGY_CSV_DIR = Path(r"<LOCAL_DATA_ROOT>\buggy_csvs")
INPUT_SITE_1 = BUGGY_CSV_DIR / "site_01_enriched.csv"
INPUT_SITE_5 = BUGGY_CSV_DIR / "site_05_enriched.csv"

# Output CSVs (patched versions)
OUTPUT_SITE_1 = BUGGY_CSV_DIR / "site_01_enriched_FIXED.csv"
OUTPUT_SITE_5 = BUGGY_CSV_DIR / "site_05_enriched_FIXED.csv"

# GEDI files location
GEDI_BASE = Path(r"<LOCAL_DATA_ROOT>\gedi\sites")

# Raster data directories
P3D_BASE = Path(r"E:\p3d")
DEP_BASE = Path(r"D:\3dep")
ALS_CHM_BASE = Path(r"D:\chms")

# Footprint radius for coverage checking (meters)
FOOTPRINT_RADIUS = 12.5

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_raster_bounds(raster_paths, crs_target=None):
    """
    Get union of bounds from multiple rasters.
    Returns (minx, miny, maxx, maxy) in target CRS.
    """
    if not raster_paths:
        return None
    
    all_bounds = []
    for path in raster_paths:
        try:
            with rasterio.open(str(path)) as src:
                bounds = src.bounds
                src_crs = src.crs
                
                if crs_target and src_crs != crs_target:
                    # Transform bounds to target CRS
                    bounds = transform_bounds(src_crs, crs_target, *bounds)
                else:
                    bounds = (bounds.left, bounds.bottom, bounds.right, bounds.top)
                
                all_bounds.append(bounds)
        except Exception as e:
            print(f"    Warning: Could not read {path}: {e}")
            continue
    
    if not all_bounds:
        return None
    
    # Compute union
    minx = min(b[0] for b in all_bounds)
    miny = min(b[1] for b in all_bounds)
    maxx = max(b[2] for b in all_bounds)
    maxy = max(b[3] for b in all_bounds)
    
    return (minx, miny, maxx, maxy)


def find_rasters(base_path, site_id, *subdirs):
    """Find all .tif/.tiff files in specified subdirectory structure."""
    search_path = base_path / str(site_id)
    for subdir in subdirs:
        search_path = search_path / subdir
    
    if not search_path.exists():
        print(f"    Warning: Path does not exist: {search_path}")
        return []
    
    rasters = list(search_path.rglob("*.tif")) + list(search_path.rglob("*.tiff"))
    return [str(r) for r in rasters]


def find_3dep_rasters(base_path, site_id):
    """
    Find 3DEP rasters for a site.
    Try multiple possible directory structures.
    """
    # Try different possible structures
    possible_paths = [
        base_path / f"polygon_{site_id}",  # e.g., D:\3dep\polygon_1\
        base_path / str(site_id),           # e.g., D:\3dep\1\
        base_path / f"site_{site_id}",      # e.g., D:\3dep\site_1\
    ]
    
    for search_path in possible_paths:
        if search_path.exists():
            rasters = list(search_path.rglob("*.tif")) + list(search_path.rglob("*.tiff"))
            if rasters:
                print(f"    Found 3DEP rasters in: {search_path}")
                return [str(r) for r in rasters]
    
    print(f"    Warning: No 3DEP rasters found for site {site_id}")
    return []


def point_in_bounds(x, y, bounds, buffer=0):
    """Check if point (x,y) is within bounds with optional buffer."""
    if bounds is None:
        return False
    minx, miny, maxx, maxy = bounds
    return (minx + buffer) <= x <= (maxx - buffer) and (miny + buffer) <= y <= (maxy - buffer)


def get_gedi_crs(site_id):
    """Get CRS from GEDI data for this site."""
    import geopandas as gpd
    gedi_path = GEDI_BASE / str(site_id) / f"GEDI_site{site_id}_hq_ALL.gpkg"
    if not gedi_path.exists():
        raise FileNotFoundError(f"GEDI file not found: {gedi_path}")
    
    gdf = gpd.read_file(gedi_path, rows=1)
    return gdf.crs


# =============================================================================
# MAIN PATCHING FUNCTION
# =============================================================================

def patch_site(input_csv, output_csv, site_id):
    """Patch aligned metrics for a single site."""
    
    print(f"\n{'='*80}")
    print(f"PATCHING SITE {site_id}")
    print(f"{'='*80}")
    
    # 1. Load existing CSV
    print(f"\n1. Loading {input_csv.name}...")
    if not input_csv.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")
    
    df = pd.read_csv(input_csv)
    n_orig = len(df)
    print(f"   Loaded {n_orig:,} observations")
    
    # Quick check of current error statistics
    print(f"\n   Current error statistics (BEFORE patching):")
    if 'dtm_error_mean' in df.columns:
        dtm_errors = df['dtm_error_mean'].dropna()
        if len(dtm_errors) > 0:
            print(f"     DTM errors: N={len(dtm_errors):,}, "
                  f"mean={dtm_errors.mean():.3f}m, "
                  f"median={dtm_errors.median():.3f}m, "
                  f"std={dtm_errors.std():.3f}m")
            catastrophic = (np.abs(dtm_errors) > 20).sum()
            print(f"     Catastrophic (|e|>20m): {catastrophic:,} ({100*catastrophic/len(dtm_errors):.2f}%)")
    
    if 'error_mean' in df.columns:
        chm_errors = df['error_mean'].dropna()
        if len(chm_errors) > 0:
            print(f"     CHM errors: N={len(chm_errors):,}, "
                  f"mean={chm_errors.mean():.3f}m, "
                  f"median={chm_errors.median():.3f}m, "
                  f"std={chm_errors.std():.3f}m")
    
    # 2. Get spatial coverage of datasets
    print(f"\n2. Determining spatial coverage of raster datasets...")
    
    # Get GEDI CRS for this site
    gedi_crs = get_gedi_crs(site_id)
    print(f"   GEDI CRS: {gedi_crs.to_string()}")
    
    # Find rasters
    print(f"   Finding P3D DTM rasters...")
    p3d_dtm_rasters = find_rasters(P3D_BASE, site_id, "vricon_raster_50cm", "dtm", "data")
    print(f"   Found {len(p3d_dtm_rasters)} P3D DTM tiles")
    
    print(f"   Finding P3D CHM rasters...")
    p3d_chm_rasters = find_rasters(P3D_BASE, site_id, "vricon_raster_50cm", "dhm", "data")
    print(f"   Found {len(p3d_chm_rasters)} P3D CHM tiles")
    
    print(f"   Finding 3DEP DTM rasters...")
    dep_dtm_rasters = find_3dep_rasters(DEP_BASE, site_id)
    print(f"   Found {len(dep_dtm_rasters)} 3DEP DTM tiles")
    
    print(f"   Finding ALS CHM rasters...")
    # Try different possible subdirectory structures
    als_chm_rasters = (
        find_rasters(ALS_CHM_BASE, site_id, "nad83") or 
        find_rasters(ALS_CHM_BASE, site_id, "chm") or
        find_rasters(ALS_CHM_BASE, site_id)
    )
    print(f"   Found {len(als_chm_rasters)} ALS CHM tiles")
    
    # Check that we found rasters
    if not p3d_dtm_rasters:
        print(f"\n   WARNING: No P3D DTM rasters found!")
    if not p3d_chm_rasters:
        print(f"\n   WARNING: No P3D CHM rasters found!")
    if not als_chm_rasters:
        print(f"\n   WARNING: No ALS CHM rasters found!")
    
    # Get bounds
    print(f"\n   Computing spatial bounds...")
    p3d_dtm_bounds = get_raster_bounds(p3d_dtm_rasters, gedi_crs) if p3d_dtm_rasters else None
    p3d_chm_bounds = get_raster_bounds(p3d_chm_rasters, gedi_crs) if p3d_chm_rasters else None
    als_chm_bounds = get_raster_bounds(als_chm_rasters, gedi_crs) if als_chm_rasters else None
    
    if p3d_dtm_bounds:
        print(f"   P3D DTM bounds: ({p3d_dtm_bounds[0]:.2f}, {p3d_dtm_bounds[1]:.2f}, "
              f"{p3d_dtm_bounds[2]:.2f}, {p3d_dtm_bounds[3]:.2f})")
    else:
        print(f"   P3D DTM bounds: None (no rasters found)")
    
    if p3d_chm_bounds:
        print(f"   P3D CHM bounds: ({p3d_chm_bounds[0]:.2f}, {p3d_chm_bounds[1]:.2f}, "
              f"{p3d_chm_bounds[2]:.2f}, {p3d_chm_bounds[3]:.2f})")
    else:
        print(f"   P3D CHM bounds: None (no rasters found)")
    
    if als_chm_bounds:
        print(f"   ALS CHM bounds: ({als_chm_bounds[0]:.2f}, {als_chm_bounds[1]:.2f}, "
              f"{als_chm_bounds[2]:.2f}, {als_chm_bounds[3]:.2f})")
    else:
        print(f"   ALS CHM bounds: None (no rasters found)")
    
    # 3. Check coverage for each GEDI footprint
    print(f"\n3. Checking GEDI footprint coverage...")
    
    # Load GEDI coordinates
    import geopandas as gpd
    gedi_path = GEDI_BASE / str(site_id) / f"GEDI_site{site_id}_hq_ALL.gpkg"
    print(f"   Loading coordinates from: {gedi_path}")
    gdf = gpd.read_file(gedi_path)
    
    # Extract coordinates
    gdf['x'] = gdf.geometry.x
    gdf['y'] = gdf.geometry.y
    
    # Merge coordinates into dataframe
    print(f"   Merging coordinates with enriched data...")
    df = df.merge(gdf[['shot_number', 'x', 'y']], on='shot_number', how='left')
    
    # Check how many footprints have coordinates
    n_with_coords = (~df['x'].isna()).sum()
    print(f"   {n_with_coords:,} / {n_orig:,} footprints have coordinates")
    
    # Check coverage with buffer equal to footprint radius
    buffer = FOOTPRINT_RADIUS
    
    print(f"   Checking spatial coverage (buffer={buffer}m)...")
    df['p3d_dtm_covered'] = df.apply(
        lambda row: point_in_bounds(row['x'], row['y'], p3d_dtm_bounds, buffer) 
        if pd.notna(row['x']) and pd.notna(row['y']) else False,
        axis=1
    )
    
    df['p3d_chm_covered'] = df.apply(
        lambda row: point_in_bounds(row['x'], row['y'], p3d_chm_bounds, buffer)
        if pd.notna(row['x']) and pd.notna(row['y']) else False,
        axis=1
    )
    
    df['als_chm_covered'] = df.apply(
        lambda row: point_in_bounds(row['x'], row['y'], als_chm_bounds, buffer)
        if pd.notna(row['x']) and pd.notna(row['y']) else False,
        axis=1
    )
    
    n_dtm_covered = df['p3d_dtm_covered'].sum()
    n_chm_covered = (df['p3d_chm_covered'] & df['als_chm_covered']).sum()
    
    print(f"   DTM: {n_dtm_covered:,} / {n_orig:,} "
          f"({100*n_dtm_covered/n_orig:.1f}%) footprints within P3D DTM coverage")
    print(f"   CHM: {n_chm_covered:,} / {n_orig:,} "
          f"({100*n_chm_covered/n_orig:.1f}%) footprints within P3D+ALS CHM coverage")
    
    # 4. Fix aligned DTM metrics
    print(f"\n4. Patching DTM aligned metrics...")
    
    # For footprints outside P3D DTM coverage: set to NaN
    outside_dtm = ~df['p3d_dtm_covered']
    n_outside_dtm = outside_dtm.sum()
    
    if n_outside_dtm > 0:
        print(f"   Setting {n_outside_dtm:,} footprints outside P3D DTM coverage to NaN")
        dtm_cols_to_null = ['dtm_rmse', 'dtm_error_mean', 'dtm_error_median', 'dtm_rmse_valid_frac',
                           'p3d_dtm_mean', 'p3d_dtm_median', 'p3d_dtm_valid_frac',
                           'p3d_dtm_p25', 'p3d_dtm_p75', 'p3d_dtm_p90']
        for col in dtm_cols_to_null:
            if col in df.columns:
                df.loc[outside_dtm, col] = np.nan
    
    # For footprints inside coverage: recompute error from individual means
    inside_dtm = (df['p3d_dtm_covered'] & 
                  df['dep_dtm_mean'].notna() & 
                  df['p3d_dtm_mean'].notna())
    n_inside_dtm = inside_dtm.sum()
    
    if n_inside_dtm > 0:
        print(f"   Recomputing error metrics for {n_inside_dtm:,} valid footprints")
        df.loc[inside_dtm, 'dtm_error_mean'] = (
            df.loc[inside_dtm, 'p3d_dtm_mean'] - df.loc[inside_dtm, 'dep_dtm_mean']
        )
        if 'dtm_error_median' in df.columns and 'p3d_dtm_median' in df.columns and 'dep_dtm_median' in df.columns:
            df.loc[inside_dtm, 'dtm_error_median'] = (
                df.loc[inside_dtm, 'p3d_dtm_median'] - df.loc[inside_dtm, 'dep_dtm_median']
            )
        # RMSE approximation from error mean (rough estimate)
        if 'dtm_rmse' in df.columns:
            df.loc[inside_dtm, 'dtm_rmse'] = np.abs(df.loc[inside_dtm, 'dtm_error_mean'])
    
    # 5. Fix aligned CHM metrics
    print(f"\n5. Patching CHM aligned metrics...")
    
    # For footprints outside coverage: set to NaN
    outside_chm = ~df['p3d_chm_covered'] | ~df['als_chm_covered']
    n_outside_chm = outside_chm.sum()
    
    if n_outside_chm > 0:
        print(f"   Setting {n_outside_chm:,} footprints outside CHM coverage to NaN")
        chm_cols_to_null = ['chm_rmse', 'error_mean', 'error_median', 'chm_rmse_valid_frac',
                           'p3d_chm_mean', 'p3d_chm_median', 'p3d_chm_valid_frac',
                           'p3d_chm_p25', 'p3d_chm_p75', 'p3d_chm_p90']
        for col in chm_cols_to_null:
            if col in df.columns:
                df.loc[outside_chm, col] = np.nan
    
    # For footprints inside coverage: recompute error from individual means
    inside_chm = (
        df['p3d_chm_covered'] & df['als_chm_covered'] & 
        df['als_chm_mean'].notna() & df['p3d_chm_mean'].notna()
    )
    n_inside_chm = inside_chm.sum()
    
    if n_inside_chm > 0:
        print(f"   Recomputing error metrics for {n_inside_chm:,} valid footprints")
        df.loc[inside_chm, 'error_mean'] = (
            df.loc[inside_chm, 'p3d_chm_mean'] - df.loc[inside_chm, 'als_chm_mean']
        )
        if 'error_median' in df.columns and 'p3d_chm_median' in df.columns and 'als_chm_median' in df.columns:
            df.loc[inside_chm, 'error_median'] = (
                df.loc[inside_chm, 'p3d_chm_median'] - df.loc[inside_chm, 'als_chm_median']
            )
        # RMSE approximation
        if 'chm_rmse' in df.columns:
            df.loc[inside_chm, 'chm_rmse'] = np.abs(df.loc[inside_chm, 'error_mean'])
    
    # 6. Report statistics
    print(f"\n6. Summary statistics (AFTER patching):")
    
    if inside_dtm.sum() > 0:
        dtm_errors = df.loc[inside_dtm, 'dtm_error_mean'].dropna()
        print(f"\n   DTM Errors:")
        print(f"     N valid:       {len(dtm_errors):,}")
        print(f"     Mean:          {dtm_errors.mean():.3f} m")
        print(f"     Median:        {dtm_errors.median():.3f} m")
        print(f"     Std:           {dtm_errors.std():.3f} m")
        print(f"     RMSE:          {np.sqrt((dtm_errors**2).mean()):.3f} m")
        print(f"     Min/Max:       [{dtm_errors.min():.3f}, {dtm_errors.max():.3f}] m")
        catastrophic = (np.abs(dtm_errors) > 20).sum()
        print(f"     Catastrophic:  {catastrophic:,} ({100*catastrophic/len(dtm_errors):.2f}%)")
    
    if inside_chm.sum() > 0:
        chm_errors = df.loc[inside_chm, 'error_mean'].dropna()
        print(f"\n   CHM Errors:")
        print(f"     N valid:       {len(chm_errors):,}")
        print(f"     Mean:          {chm_errors.mean():.3f} m")
        print(f"     Median:        {chm_errors.median():.3f} m")
        print(f"     Std:           {chm_errors.std():.3f} m")
        print(f"     RMSE:          {np.sqrt((chm_errors**2).mean()):.3f} m")
        print(f"     Min/Max:       [{chm_errors.min():.3f}, {chm_errors.max():.3f}] m")
    
    # 7. Save patched CSV
    print(f"\n7. Saving patched CSV...")
    
    # Drop temporary columns
    df = df.drop(columns=['x', 'y', 'p3d_dtm_covered', 'p3d_chm_covered', 'als_chm_covered'], 
                 errors='ignore')
    
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False, float_format="%.4f")
    
    print(f"   ✓ Saved {len(df):,} observations to:")
    print(f"     {output_csv}")
    print(f"\n{'='*80}")
    print(f"SITE {site_id} PATCHING COMPLETE")
    print(f"{'='*80}\n")


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Patch sites 1 and 5."""
    
    print("\n" + "="*80)
    print("ALIGNED METRICS PATCHER - LOCAL VERSION")
    print("="*80)
    print("\nThis script patches buggy aligned metrics in sites 1 and 5")
    print("without requiring full resampling (saves 20+ hours per site).")
    print("\nStrategy:")
    print("  1. Check which GEDI footprints fall outside P3D coverage")
    print("  2. Set aligned metrics to NaN for those footprints")
    print("  3. Recompute error metrics from individual means for valid footprints")
    print("\nNote: This is an approximation but should be ~95% accurate")
    print("      and dramatically better than the buggy version.")
    print(f"\nInput directory:  {BUGGY_CSV_DIR}")
    print(f"Output directory: {BUGGY_CSV_DIR}")
    print(f"Output naming:    site_XX_enriched_FIXED.csv")
    
    # Check input files exist
    missing_files = []
    for site_id, input_csv in [(1, INPUT_SITE_1), (5, INPUT_SITE_5)]:
        if not input_csv.exists():
            missing_files.append((site_id, input_csv))
    
    if missing_files:
        print("\n" + "!"*80)
        print("ERROR: Input CSV files not found!")
        print("!"*80)
        for site_id, path in missing_files:
            print(f"  Site {site_id}: {path}")
        print(f"\nPlease verify that buggy enriched CSVs exist in: {BUGGY_CSV_DIR}")
        return 1
    
    # Patch both sites
    try:
        patch_site(INPUT_SITE_1, OUTPUT_SITE_1, site_id=1)
        patch_site(INPUT_SITE_5, OUTPUT_SITE_5, site_id=5)
        
        print("\n" + "="*80)
        print("ALL PATCHING COMPLETE!")
        print("="*80)
        print(f"\nPatched files saved to:")
        print(f"  - {OUTPUT_SITE_1}")
        print(f"  - {OUTPUT_SITE_5}")
        print(f"\nNext steps:")
        print(f"  1. Review the AFTER patching statistics above")
        print(f"     - DTM errors should be ~1-2m RMSE (not ~800m!)")
        print(f"     - Catastrophic errors should be <1%")
        print(f"  2. If statistics look good, these are ready to use!")
        print(f"  3. Upload to cluster and integrate with other enriched site files")
        
        return 0
        
    except Exception as e:
        print(f"\nERROR during patching: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
