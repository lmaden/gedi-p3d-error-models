# find_overlaps_fast.py
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import rasterio
from rasterio.warp import transform_bounds
from rasterio.errors import RasterioIOError

def bboxes_intersect(a, b):
    """True if axis-aligned boxes a and b have area overlap (> 0)."""
    l1, b1, r1, t1 = a
    l2, b2, r2, t2 = b
    return not (r1 <= l2 or r2 <= l1 or t1 <= b2 or t2 <= b1)

def read_meta(path):
    """Open once to read CRS and bounds; return (path, crs_str, bounds) or (path, None, None) on failure."""
    try:
        with rasterio.open(path) as ds:
            crs = ds.crs
            if crs is None:
                return (path, None, None)
            b = ds.bounds
            return (path, crs.to_string(), (b.left, b.bottom, b.right, b.top))
    except RasterioIOError:
        return (path, None, None)
    except Exception:
        return (path, None, None)

def xform_bounds(src_crs, dst_crs, bounds, densify_pts=0):
    """Reproject bounding box with optional densification (0 = corners only, very fast)."""
    l, b, r, t = bounds
    try:
        return transform_bounds(src_crs, dst_crs, l, b, r, t, densify_pts=densify_pts, always_xy=True)
    except TypeError:
        # Older rasterio without always_xy
        return transform_bounds(src_crs, dst_crs, l, b, r, t, densify_pts=densify_pts)

def find_overlaps_fast(
    dir_path,
    target_path,
    recursive=False,
    exts=(".tif", ".tiff", ".vrt"),
    workers=8,
    densify_pts=0,
    use_gdal_speedups=True,
):
    """
    Return file NAMES in dir_path that spatially overlap target_path.
    Speed tricks:
      * One target->CRS transform per unique candidate CRS
      * Threaded metadata reads
      * Fast bbox reprojection (densify_pts=0 by default)
    """
    # Optional GDAL speedups (safe for typical GeoTIFFs on local disk)
    env_kwargs = {}
    if use_gdal_speedups:
        env_kwargs.update({
            "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
        })

    with rasterio.Env(**env_kwargs):
        # Load target metadata
        with rasterio.open(target_path) as tgt:
            if tgt.crs is None:
                raise ValueError("Target raster has no CRS.")
            tgt_crs = tgt.crs
            tb = tgt.bounds
            tgt_bounds = (tb.left, tb.bottom, tb.right, tb.top)

        # Gather candidate files
        root = Path(dir_path)
        iterator = root.rglob("*") if recursive else root.iterdir()
        candidates = [p for p in iterator if p.is_file() and p.suffix.lower() in exts]

        # Read CRS & bounds in parallel
        metas = []
        max_workers = max(1, workers)
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            futures = {ex.submit(read_meta, p): p for p in candidates}
            for fut in as_completed(futures):
                metas.append(fut.result())

        # Group by CRS string (avoid transform per file)
        by_crs = {}
        for path, crs_str, bounds in metas:
            if crs_str and bounds:
                by_crs.setdefault(crs_str, []).append((path, bounds))

        # Precompute target bounds in each group's CRS once
        overlaps = []
        # Make a small cache of CRS objects so we don't reparse strings repeatedly
        crs_cache = {}

        def get_crs_obj(crs_str):
            obj = crs_cache.get(crs_str)
            if obj is None:
                obj = rasterio.crs.CRS.from_string(crs_str)
                crs_cache[crs_str] = obj
            return obj

        for crs_str, items in by_crs.items():
            group_crs = get_crs_obj(crs_str)
            if group_crs == tgt_crs:
                tgt_in_group = tgt_bounds
            else:
                tgt_in_group = xform_bounds(tgt_crs, group_crs, tgt_bounds, densify_pts=densify_pts)

            # Now cheap AABB tests for all in this CRS
            for path, cand_bounds in items:
                if bboxes_intersect(cand_bounds, tgt_in_group):
                    overlaps.append(Path(path).name)  # return just the name

        return sorted(overlaps)

if __name__ == "__main__":
    directory = r"D:\chms\5\nad83"
    target = r"E:\p3d\5\vricon_raster_50cm\dhm\data\0802733w_330747n_20250401T082952Z_dhm.tif"

    hits = find_overlaps_fast(
        directory,
        target,
        recursive=False,
        workers=8,          # try 8–12; it's I/O-bound
        densify_pts=0,      # corners only; use 5–8 if you need more accuracy
        use_gdal_speedups=True
    )
    for name in hits:
        print(name)

