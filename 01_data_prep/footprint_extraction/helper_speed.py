#!/usr/bin/env python3
"""
helper_speed.py (FIXED)
=======================

High‑performance raster utilities with CRITICAL BUG FIX for spatial alignment.

CRITICAL FIX in v09:
- zonal_aligned_pair_metrics now checks if footprints fall within ds2's valid bounds
- Prevents sampling from outside P3D tile extents when comparing with 3DEP
- Fixes the bug causing ~800m DTM errors in enriched CSVs

This revision:
- Adds bounds checking to prevent sampling outside ds2's valid extent
- Adds `percentiles` support to `zonal_stats` (e.g., p90/p75/p25).
- Adds `zonal_aligned_means` to compute mean/median on *aligned* windows
  (ds2 -> ds1) using identical valid-pixel masks (mirrors RMSE alignment).
- Adds `zonal_aligned_pair_metrics` — single-pass computation of
  RMSE + aligned means (+ aligned medians) using one WarpedVRT + one window loop.
- Keeps mask‑aware zonal reducers, metric‑true kernels, and VRT helpers.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import tempfile
import warnings
from functools import lru_cache
from pathlib import Path
from subprocess import CalledProcessError, run
from typing import List, Sequence, Tuple, Union, Optional, Dict

import numba as nb
import numpy as np
import rasterio
from filelock import FileLock
from rasterio.crs import CRS
from rasterio.enums import Resampling, ColorInterp
from rasterio.errors import NotGeoreferencedWarning
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window
from rasterio.windows import transform as window_transform

# ---------------------------------------------------------------------
# 0. Global sanity checks & warning filters
# ---------------------------------------------------------------------
warnings.filterwarnings("ignore", category=NotGeoreferencedWarning)

_REQUIRED_GDAL_CMDS: Sequence[str] = ("gdalbuildvrt", "gdalwarp")
_missing = [cmd for cmd in _REQUIRED_GDAL_CMDS if shutil.which(cmd) is None]
if _missing:
    raise RuntimeError(
        f"Required GDAL utilities not found on PATH: {', '.join(_missing)}\n"
        "→ Install GDAL (conda `gdal` or OSGeo4W) and ensure its `bin/` is on PATH."
    )

# ---------------------------------------------------------------------
# CRS safety patch – auto‑assign EPSG:4326 if .crs is missing
# ---------------------------------------------------------------------
import geopandas as gpd
if not getattr(gpd.GeoDataFrame, "_gedi_crs_patch", False):
    _orig_gdf_to_crs = gpd.GeoDataFrame.to_crs
    _orig_gs_to_crs  = gpd.GeoSeries.to_crs
    def _safe_to_crs_gdf(self, *args, **kwargs):
        if self.crs is None:
            self = self.set_crs("EPSG:4326", inplace=False)
        return _orig_gdf_to_crs(self, *args, **kwargs)
    def _safe_to_crs_gs(self, *args, **kwargs):
        if self.crs is None:
            self = self.set_crs("EPSG:4326", inplace=False)
        return _orig_gs_to_crs(self, *args, **kwargs)
    gpd.GeoDataFrame.to_crs = _safe_to_crs_gdf
    gpd.GeoSeries.to_crs    = _safe_to_crs_gs
    gpd.GeoDataFrame._gedi_crs_patch = True

# ---------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------
def _run(cmd: List[str]) -> None:
    """Run *cmd* hiding stdout/stderr; raise if exit-code ≠ 0."""
    run(cmd, check=True, capture_output=True)

def _safe(p: str | Path) -> str:
    """Return absolute path, adding \\?\\ on Windows for very long paths."""
    s = str(Path(p).resolve())
    if os.name == "nt" and len(s) >= 250 and not s.startswith("\\\\?\\"):
        s = "\\\\?\\" + s
    return s

def _has_alpha_band(ds: rasterio.io.DatasetReader) -> bool:
    try:
        return any(ci == ColorInterp.alpha for ci in ds.colorinterp)
    except Exception:
        return False

# ---------------------------------------------------------------------
# 1. Build / open a VRT (thread‑safe) — with unique filename hashing
# ---------------------------------------------------------------------
def build_vrt(
    paths: List[str],
    *,
    target_crs: CRS | None = None,
    resampling: str = "nearest",
) -> rasterio.DatasetReader:
    """
    Return an open Rasterio dataset that is a mosaic of *paths*.

    * Thread‑safe via file‑locking.
    * Optional on‑the‑fly reprojection.
    * Unique VRT filename per (paths, target_crs, resampling).
    """
    if not paths:
        raise FileNotFoundError("Empty path list given to build_vrt()")

    # Fast path: single file & no reprojection requested
    if target_crs is None and len(paths) == 1:
        return rasterio.open(_safe(paths[0]))

    tag = target_crs.to_string() if target_crs else "src"
    digest = hashlib.sha1(
        ("|".join(sorted(map(str, paths))) + f"|{tag}|{resampling}").encode("utf-8")
    ).hexdigest()[:16]
    vrt_fp = Path(paths[0]).parent / f"_tmp_{digest}.vrt"
    lock_fp = vrt_fp.with_suffix(".lock")

    with FileLock(str(lock_fp)):
        if not vrt_fp.exists():
            with tempfile.NamedTemporaryFile("w", delete=False, suffix=".txt") as flist:
                for p in paths:
                    flist.write(f"{_safe(p)}\n")
                list_path = flist.name
            try:
                if target_crs is None:
                    _run([
                        "gdalbuildvrt", "-q",
                        "-overwrite",
                        "-input_file_list", list_path, str(vrt_fp)
                    ])
                else:
                    cmd = [
                        "gdalbuildvrt", "-q",
                        "-overwrite",
                        "-t_srs", target_crs.to_string(),
                        "-allow_projection_difference",
                        "-r", resampling,
                        "-input_file_list", list_path, str(vrt_fp),
                    ]
                    try:
                        _run(cmd)
                    except CalledProcessError as e:
                        if b"Unknown argument: -t_srs" not in e.stderr:
                            raise
                        # Fallback: build then warp
                        src_vrt = Path(vrt_fp).with_suffix(".src.vrt")
                        _run([
                            "gdalbuildvrt", "-q",
                            "-allow_projection_difference",
                            "-input_file_list", list_path, str(src_vrt)
                        ])
                        try:
                            _run([
                                "gdalwarp", "-q", "-overwrite", "-of", "VRT",
                                "-t_srs", target_crs.to_string(),
                                "-r", resampling,
                                str(src_vrt), str(vrt_fp),
                            ])
                        finally:
                            src_vrt.unlink(missing_ok=True)
            finally:
                Path(list_path).unlink(missing_ok=True)

    try:
        ds = rasterio.open(_safe(vrt_fp))
        if ds.count == 0:
            raise rasterio.errors.RasterioIOError("empty VRT")
        return ds
    except rasterio.errors.RasterioIOError:
        base = rasterio.open(_safe(paths[0]))
        return WarpedVRT(
            base,
            crs=target_crs,
            resampling=Resampling[resampling],
            src_nodata=base.nodata,
            dst_nodata=base.nodata,
            add_alpha=not _has_alpha_band(base),
            init_dest_nodata=True,
            vrt_path="/vsimem/fallback.vrt",
            warp_mem_limit=256,
        )

# ---------------------------------------------------------------------
# 2. Metric‑true kernel cache (supports non‑square pixels)
# ---------------------------------------------------------------------
@lru_cache(maxsize=8)
def _circular_kernel(px_radius: int) -> np.ndarray:
    y, x = np.ogrid[-px_radius:px_radius + 1, -px_radius:px_radius + 1]
    return (x * x + y * y) <= px_radius * px_radius

@lru_cache(maxsize=64)
def _metric_kernel(rx: int, ry: int, res_x: float, res_y: float) -> np.ndarray:
    """
    Boolean mask approximating a metric circle of R meters when dx != dy:
        (x*res_x)^2 + (y*res_y)^2 <= (rx*res_x)^2
    """
    y, x = np.ogrid[-ry:ry + 1, -rx:rx + 1]
    return ((x * res_x) ** 2 + (y * res_y) ** 2) <= (rx * res_x) ** 2

def kernel_for(ds: rasterio.DatasetReader, radius_m: float) -> np.ndarray:
    """Metric‑accurate circular kernel for *ds* at *radius_m* (meters)."""
    res_x = abs(ds.transform.a)
    res_y = abs(ds.transform.e)
    rx = max(1, int(np.ceil(radius_m / res_x)))
    ry = max(1, int(np.ceil(radius_m / res_y)))
    if abs(res_x - res_y) < 1e-9:
        return _circular_kernel(rx)
    return _metric_kernel(rx, ry, res_x, res_y)

def _window_and_kernel_slices(
    col: int,
    row: int,
    ds: rasterio.DatasetReader,
    kernel: np.ndarray,
) -> Tuple[Window, slice, slice]:
    """Return (window, kernel_row_slice, kernel_col_slice) for footprint at (col, row)."""
    kh, kw = kernel.shape
    ry, rx = kh // 2, kw // 2
    left = max(0, col - rx)
    top = max(0, row - ry)
    right = min(ds.width - 1, col + rx)
    bottom = min(ds.height - 1, row + ry)
    win = Window(left, top, right - left + 1, bottom - top + 1)
    ksr = slice(max(0, ry - row), kh - max(0, row + ry + 1 - ds.height))
    ksc = slice(max(0, rx - col), kw - max(0, col + rx + 1 - ds.width))
    return win, ksr, ksc

# ---------------------------------------------------------------------
# 3. Fast reducers (numba‑accelerated)
# ---------------------------------------------------------------------
@nb.njit
def _nanmean(arr: np.ndarray) -> float:
    s, c = 0.0, 0
    for v in arr.flat:
        if np.isfinite(v):
            s += v
            c += 1
    return (s / c) if c > 0 else np.nan

@nb.njit
def _nanmedian(arr: np.ndarray) -> float:
    finite = arr[np.isfinite(arr)]
    if len(finite) == 0:
        return np.nan
    return float(np.median(finite))

@nb.njit
def _nanstd(arr: np.ndarray) -> float:
    s, s2, c = 0.0, 0.0, 0
    for v in arr.flat:
        if np.isfinite(v):
            s += v
            s2 += v * v
            c += 1
    if c < 2:
        return np.nan
    mean = s / c
    var = (s2 / c) - mean * mean
    return np.sqrt(max(0.0, var))

@nb.njit
def _rmse(a1: np.ndarray, a2: np.ndarray) -> float:
    s, c = 0.0, 0
    for i in range(len(a1)):
        v1, v2 = a1[i], a2[i]
        if np.isfinite(v1) and np.isfinite(v2):
            s += (v1 - v2) ** 2
            c += 1
    return np.sqrt(s / c) if c > 0 else np.nan

# ---------------------------------------------------------------------
# 4. Zonal stats (mean, median, std, percentiles, valid_frac)
# ---------------------------------------------------------------------
def zonal_stats(
    ds: rasterio.DatasetReader,
    xs: np.ndarray,
    ys: np.ndarray,
    radius_m: float,
    *,
    want_std: bool = False,
    want_median: bool = False,
    percentiles: Sequence[float] = (),
    transform_fn=None,
    return_valid_fraction: bool = False,
    min_valid_frac: float | None = None,
    use_dataset_mask: bool = True,
    ignore_nodata_value: bool = False,
) -> Dict[str, np.ndarray]:
    """
    Extract zonal statistics within circular footprints.
    
    Returns dict with keys: mean, [median], [std], [p90], [p75], [p25], [valid_frac]
    """
    kernel = kernel_for(ds, radius_m)
    n = len(xs)
    
    mu = np.full(n, np.nan, dtype=np.float32)
    med = np.full(n, np.nan, dtype=np.float32) if want_median else None
    sd = np.full(n, np.nan, dtype=np.float32) if want_std else None
    vfr = np.full(n, np.nan, dtype=np.float32) if return_valid_fraction else None
    
    pct_arrays = {}
    for p in percentiles:
        pct_arrays[p] = np.full(n, np.nan, dtype=np.float32)
    
    Tinv = ~ds.transform
    nodata = ds.nodata if not ignore_nodata_value else None
    
    cols_f, rows_f = Tinv * (xs, ys)
    inside = (cols_f >= 0) & (cols_f < ds.width) & (rows_f >= 0) & (rows_f < ds.height)
    
    for i in np.nonzero(inside)[0]:
        c, r = int(cols_f[i]), int(rows_f[i])
        win, ks_row, ks_col = _window_and_kernel_slices(c, r, ds, kernel)
        arr = ds.read(1, window=win).astype(np.float32)
        
        if use_dataset_mask:
            m = ds.read_masks(1, window=win)
            arr[m == 0] = np.nan
        
        if nodata is not None:
            arr[arr == nodata] = np.nan
        
        if transform_fn is not None:
            arr = transform_fn(arr)
        
        k = kernel[ks_row, ks_col]
        arr_masked = arr.copy()
        arr_masked[~k] = np.nan
        
        valid = np.isfinite(arr_masked)
        k_count = int(k.sum())
        valid_count = int(valid.sum())
        vfrac = (valid_count / k_count) if k_count > 0 else 0.0
        
        if (min_valid_frac is not None) and (vfrac < min_valid_frac):
            if vfr is not None:
                vfr[i] = vfrac
            continue
        
        vals = arr_masked[valid]
        if len(vals) > 0:
            mu[i] = _nanmean(vals)
            if want_median:
                med[i] = _nanmedian(vals)
            if want_std:
                sd[i] = _nanstd(vals)
            for p in percentiles:
                pct_arrays[p][i] = np.percentile(vals, p * 100)
        
        if vfr is not None:
            vfr[i] = vfrac
    
    result = {"mean": mu}
    if want_median:
        result["median"] = med
    if want_std:
        result["std"] = sd
    for p in percentiles:
        pname = f"p{int(p*100):02d}"
        result[pname] = pct_arrays[p]
    if return_valid_fraction:
        result["valid_frac"] = vfr
    
    return result

# ---------------------------------------------------------------------
# 5. NEW: Bounds checking helper
# ---------------------------------------------------------------------
def _check_point_in_bounds(
    ds: rasterio.DatasetReader,
    x: float,
    y: float,
    crs: CRS,
    buffer_m: float = 0.0
) -> bool:
    """
    Check if point (x, y) in given CRS falls within dataset bounds.
    Optionally add buffer_m to shrink the valid region.
    """
    # Transform point to dataset's CRS if needed
    if crs != ds.crs:
        from rasterio.warp import transform as warp_transform
        xs, ys = warp_transform(crs, ds.crs, [x], [y])
        x_ds, y_ds = xs[0], ys[0]
    else:
        x_ds, y_ds = x, y
    
    # Check bounds with optional buffer
    bounds = ds.bounds
    return (
        (bounds.left + buffer_m) <= x_ds <= (bounds.right - buffer_m) and
        (bounds.bottom + buffer_m) <= y_ds <= (bounds.top - buffer_m)
    )

# ---------------------------------------------------------------------
# 6. FIXED: Single‑pass aligned RMSE + means (+ medians) with bounds checking
# ---------------------------------------------------------------------
def zonal_aligned_pair_metrics(
    ds1: rasterio.DatasetReader,
    ds2: rasterio.DatasetReader,
    xs: np.ndarray,
    ys: np.ndarray,
    radius_m: float,
    *,
    want_median: bool = True,
    return_valid_fraction: bool = True,
    min_valid_frac: float | None = None,
    use_dataset_mask: bool = True,
    ignore_nodata1: bool = False,
    ignore_nodata2: bool = False,
    override_nodata1: float | None = None,
    override_nodata2: float | None = None,
) -> Union[
    Tuple[np.ndarray, np.ndarray, np.ndarray],
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
]:
    """
    Compute, in a SINGLE aligned pass (ds2 → ds1):
      - RMSE(ds1, ds2) on the common valid footprint pixels,
      - mean(ds1), mean(ds2) on the same pixels,
      - [optional] median(ds1), median(ds2) on the same pixels.

    CRITICAL FIX: Now checks if footprints fall within ds2's valid bounds
    before computing metrics. This prevents sampling outside P3D tile extents
    when comparing with 3DEP.

    Returns (by default):
        rmse, mean1, mean2, valid_frac [, median1, median2]
    If `return_valid_fraction=False` and/or `want_median=False`, the tuple shortens accordingly.
    """
    kernel = kernel_for(ds1, radius_m)
    warp_mem_limit = int(os.environ.get("GEODATA_ALIGN_WARP_MEM_LIMIT", "512"))
    
    # CRITICAL FIX: Pre-check which footprints fall within ds2's bounds
    # This prevents sampling from outside ds2's valid extent
    n = len(xs)
    ds2_coverage = np.zeros(n, dtype=bool)
    
    # Add buffer equal to footprint radius to ensure full footprint is within bounds
    for i in range(n):
        ds2_coverage[i] = _check_point_in_bounds(ds2, xs[i], ys[i], ds1.crs, buffer_m=radius_m)
    
    with WarpedVRT(
        ds2,
        crs=ds1.crs, transform=ds1.transform, width=ds1.width, height=ds1.height,
        resampling=Resampling.bilinear,
        src_nodata=ds2.nodata, dst_nodata=np.nan, dtype="float32",
        add_alpha=False, init_dest_nodata=True,
        vrt_path="/vsimem/pair_align.vrt", warp_mem_limit=warp_mem_limit,
    ) as ds2w:
        rmse = np.full(n, np.nan, dtype=np.float32)
        mu1  = np.full(n, np.nan, dtype=np.float32)
        mu2  = np.full(n, np.nan, dtype=np.float32)
        vfr  = np.full(n, np.nan, dtype=np.float32) if return_valid_fraction else None
        md1  = np.full(n, np.nan, dtype=np.float32) if want_median else None
        md2  = np.full(n, np.nan, dtype=np.float32) if want_median else None

        Tinv = ~ds1.transform

        if ignore_nodata1: nod1 = None
        elif override_nodata1 is not None: nod1 = override_nodata1
        else: nod1 = ds1.nodata

        if ignore_nodata2: nod2 = None
        elif override_nodata2 is not None: nod2 = override_nodata2
        else: nod2 = ds2w.nodata

        cols_f, rows_f = Tinv * (xs, ys)
        inside = (cols_f >= 0) & (cols_f < ds1.width) & (rows_f >= 0) & (rows_f < ds1.height)
        
        # CRITICAL FIX: Only process footprints that are inside BOTH datasets
        valid_footprints = inside & ds2_coverage

        for i in np.nonzero(valid_footprints)[0]:
            c, r = int(cols_f[i]), int(rows_f[i])
            win, ks_row, ks_col = _window_and_kernel_slices(c, r, ds1, kernel)
            a1 = ds1.read(1, window=win).astype(np.float32)
            a2 = ds2w.read(1, window=win).astype(np.float32)
            if use_dataset_mask:
                m1 = ds1.read_masks(1, window=win)
                m2 = ds2w.read_masks(1, window=win)
                a1[m1 == 0] = np.nan
                a2[m2 == 0] = np.nan
            if nod1 is not None: a1[a1 == nod1] = np.nan
            if nod2 is not None: a2[a2 == nod2] = np.nan
            k = kernel[ks_row, ks_col]
            vmask = np.isfinite(a1) & np.isfinite(a2) & k
            k_count = int(k.sum()); valid_count = int(vmask.sum())
            vfrac = (valid_count / k_count) if k_count > 0 else 0.0
            if (min_valid_frac is not None) and (vfrac < min_valid_frac):
                if vfr is not None: vfr[i] = vfrac
                continue
            v1 = a1[vmask]; v2 = a2[vmask]
            mu1_i = _nanmean(v1); mu2_i = _nanmean(v2)
            mu1[i] = mu1_i; mu2[i] = mu2_i
            rmse[i] = _rmse(v1, v2)
            if want_median:
                md1[i] = _nanmedian(v1)
                md2[i] = _nanmedian(v2)
            if vfr is not None:
                vfr[i] = vfrac

    outs: Tuple = (rmse, mu1, mu2)
    if return_valid_fraction:
        outs = (*outs, vfr)
    if want_median:
        outs = (*outs, md1, md2)
    return outs  # type: ignore[return-value]

# ---------------------------------------------------------------------
# 7. Point sampler
# ---------------------------------------------------------------------
def sample_points(
    ds: rasterio.DatasetReader,
    xs: np.ndarray,
    ys: np.ndarray,
    *,
    chunk_size: int = 10_000,
    use_dataset_mask: bool = True,
) -> np.ndarray:
    import rasterio
    from rasterio.windows import Window
    n = len(xs)
    out = np.full(n, np.nan, dtype=np.float32)
    if n == 0:
        return out

    Tinv = ~ds.transform
    nodata = ds.nodata
    cols_f, rows_f = Tinv * (xs, ys)
    inside = (cols_f >= 0) & (cols_f < ds.width) & (rows_f >= 0) & (rows_f < ds.height)
    if not inside.any():
        return out

    idx_inside = np.nonzero(inside)[0]
    coords = np.column_stack((xs[idx_inside], ys[idx_inside]))

    for i0 in range(0, len(coords), chunk_size):
        sub_coords = coords[i0:i0 + chunk_size]
        try:
            for j, v in enumerate(ds.sample(map(tuple, sub_coords))):
                val = v[0]; idx = idx_inside[i0 + j]
                if use_dataset_mask:
                    c, r = int(cols_f[idx]), int(rows_f[idx])
                    m = ds.read_masks(1, window=Window(c, r, 1, 1))
                    if m[0, 0] == 0:
                        val = np.nan
                if nodata is not None and val == nodata:
                    val = np.nan
                out[idx] = val
        except Exception:
            # Fallback: read 1x1 windows individually; mark pixel NaN if read fails
            for j, (x, y) in enumerate(sub_coords):
                idx = idx_inside[i0 + j]
                c, r = int(cols_f[idx]), int(rows_f[idx])
                try:
                    arr = ds.read(1, window=Window(c, r, 1, 1)).astype(np.float32)
                    val = float(arr[0, 0])
                    if use_dataset_mask:
                        m = ds.read_masks(1, window=Window(c, r, 1, 1))
                        if m[0, 0] == 0:
                            val = np.nan
                    if nodata is not None and val == nodata:
                        val = np.nan
                    out[idx] = val
                except Exception:
                    out[idx] = np.nan
    return out


# ---------------------------------------------------------------------
# 8. Debug: save a 1 km chip around a point
# ---------------------------------------------------------------------
def save_debug_chip(
    ds: rasterio.DatasetReader,
    x: float,
    y: float,
    out_tif: str | Path,
    side_m: float = 1000.0,
    *,
    use_dataset_mask: bool = True,
) -> Path:
    """Save a square chip of size `side_m` (meters) centered at (x,y) to GeoTIFF."""
    out_tif = Path(out_tif)
    out_tif.parent.mkdir(parents=True, exist_ok=True)

    res_x = abs(ds.transform.a)
    res_y = abs(ds.transform.e)
    half_w = int(np.ceil((side_m / 2) / res_x))
    half_h = int(np.ceil((side_m / 2) / res_y))

    c_f, r_f = (~ds.transform) * (x, y)
    c, r = int(round(c_f)), int(round(r_f))

    left = max(0, c - half_w); top = max(0, r - half_h)
    right = min(ds.width - 1, c + half_w); bottom = min(ds.height - 1, r + half_h)

    win = Window(left, top, right - left + 1, bottom - top + 1)
    arr = ds.read(1, window=win).astype(np.float32)

    if use_dataset_mask:
        m = ds.read_masks(1, window=win)
        arr[m == 0] = np.nan

    prof = ds.profile.copy()
    prof.update(
        driver="GTiff", height=arr.shape[0], width=arr.shape[1],
        count=1, dtype="float32", compress="lzw", tiled=True, nodata=np.nan,
        transform=window_transform(win, ds.transform),
    )
    with rasterio.open(str(out_tif), "w", **prof) as dst:
        dst.write(arr, 1)
    return out_tif
