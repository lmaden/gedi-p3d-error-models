#!/usr/bin/env python
"""
reproject_chms_to_nad83.py
==========================

Batch‑reproject every ALS CHM to the *NAD 83 / UTM* zone used by its
corresponding P3D DTM.

Path conventions
----------------
CHMs:  D:\\chms\\<site>\\…\\*.tif      (may be in a "chm" sub‑folder)
DTMs:  E:\\p3d\\<site>\\vricon_raster_50cm\\dtm\\*.tif

For each site the script:

1. Reads the CRS from the first DTM it finds (this is the NAD 83 / UTM
   target CRS).
2. Finds every CHM GeoTIFF beneath the site folder **except** those that
   are already under a "nad83" directory.
3. Reprojects each CHM with *gdalwarp*, then converts the result into a
   Cloud Optimised GeoTIFF (COG) with:
      - LZW compression
      - 512×512 internal tiles
      - pyramid levels 2‑4‑8‑16
4. Writes outputs to:
      D:\\chms\\<site>\\nad83\\<original>_nad83.tif

Only files that do **not** yet exist are processed, so the script is
re‑entrant.

Requirements
------------
* GDAL command‑line utilities (`gdalwarp`, `gdaladdo`, `gdal_translate`)
  available in the current environment (conda‑forge, OSGeo4W, or QGIS).

Usage
-----
    (gedi_env) > python reproject_chms_to_nad83.py
    # optional arguments
    (gedi_env) > python reproject_chms_to_nad83.py --chm-root D:\\chms --dtm-root E:\\p3d
"""

import argparse
import subprocess
import sys
from pathlib import Path
from typing import List

import rasterio
from rasterio.crs import CRS

# ────────────────────────────────────────────────────────────────────────────
# Configuration constants
# ────────────────────────────────────────────────────────────────────────────
DEFAULT_CHM_ROOT = Path(r"D:\chms")
DEFAULT_DTM_ROOT = Path(r"E:\p3d")
OVERVIEW_LEVELS: List[str] = ["2", "4", "8", "16"]  # internal pyramids
BLOCKSIZE = "512"


# ────────────────────────────────────────────────────────────────────────────
# Utility functions
# ────────────────────────────────────────────────────────────────────────────
def get_target_crs(site: str, dtm_root: Path) -> CRS:
    """
    Read CRS from the first DTM the function finds for *site*.
    Raises RuntimeError if none found or the CRS is missing.
    """
    dtm_dir = dtm_root / site / "vricon_raster_50cm" / "dtm"
    candidates = list(dtm_dir.rglob("*.tif"))
    if not candidates:
        raise RuntimeError(f"{dtm_dir}: no DTM found")
    with rasterio.open(candidates[0]) as ds:
        crs = ds.crs
    if crs is None:
        raise RuntimeError(f"{candidates[0]} has no CRS")
    return crs


def warp_to_cog(src: Path, dst: Path, target_crs: CRS) -> None:
    """
    Reproject *src* to *target_crs* with gdalwarp, add internal overviews
    with gdaladdo, then convert to a Cloud Optimised GeoTIFF (COG) using
    gdal_translate (GDAL ≥ 3.8).  Intermediate files are cleaned up.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix(".tmp.tif")

    # -- step 1: gdalwarp ---------------------------------------------------
    subprocess.run(
        [
            "gdalwarp",
            "-t_srs",
            target_crs.to_wkt(),
            "-r",
            "bilinear",
            "-multi",
            "-wo",
            "NUM_THREADS=ALL_CPUS",
            "-co",
            f"COMPRESS=LZW",
            "-co",
            "TILED=YES",
            "-co",
            f"BLOCKXSIZE={BLOCKSIZE}",
            "-co",
            f"BLOCKYSIZE={BLOCKSIZE}",
            "-co",
            "BIGTIFF=IF_SAFER",
            str(src),
            str(tmp),
        ],
        check=True,
    )

    # -- step 2: build internal overviews ----------------------------------
    subprocess.run(
        ["gdaladdo", "-r", "average", str(tmp), *OVERVIEW_LEVELS], check=True
    )

    # -- step 3: convert to official COG driver ----------------------------
    subprocess.run(
        [
            "gdal_translate",
            "-of",
            "COG",
            "-co",
            "COMPRESS=LZW",
            "-co",
            f"BLOCKSIZE={BLOCKSIZE}",
            str(tmp),
            str(dst),
        ],
        check=True,
    )

    tmp.unlink()  # remove intermediate
    print(f"[WRITE] {dst.relative_to(dst.parents[2])}")


def find_site_chms(site_dir: Path) -> List[Path]:
    """
    Return every *.tif* under *site_dir* that is **not** already inside a
    'nad83' sub‑folder.
    """
    return [
        p
        for p in site_dir.rglob("*.tif")
        if "nad83" not in {part.lower() for part in p.parts}
    ]


# ────────────────────────────────────────────────────────────────────────────
# Main script
# ────────────────────────────────────────────────────────────────────────────
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reproject ALS CHMs to NAD83 / UTM per site"
    )
    parser.add_argument(
        "--chm-root",
        type=Path,
        default=DEFAULT_CHM_ROOT,
        help=f"Root folder that contains per‑site CHM directories "
        f"(default: {DEFAULT_CHM_ROOT})",
    )
    parser.add_argument(
        "--dtm-root",
        type=Path,
        default=DEFAULT_DTM_ROOT,
        help=f"Root folder that contains per‑site DTM directories "
        f"(default: {DEFAULT_DTM_ROOT})",
    )
    args = parser.parse_args()

    if not args.chm_root.exists():
        sys.exit(f"{args.chm_root} does not exist")
    if not args.dtm_root.exists():
        sys.exit(f"{args.dtm_root} does not exist")

    # discover sites by folder name under CHM_ROOT -------------------------
    site_dirs = sorted([p for p in args.chm_root.iterdir() if p.is_dir()])
    if not site_dirs:
        sys.exit(f"No site folders found in {args.chm_root}")

    print(f"🔄  Reprojecting CHMs from WGS84 → NAD83 using P3D DTMs\n")

    for site_dir in site_dirs:
        site = site_dir.name
        try:
            target_crs = get_target_crs(site, args.dtm_root)
        except Exception as err:
            print(f"[SKIP] {site}: {err}")
            continue

        chms = find_site_chms(site_dir)
        if not chms:
            print(f"[SKIP] {site}: no CHMs found")
            continue

        out_dir = site_dir / "nad83"
        processed = 0

        for chm in chms:
            out_name = chm.stem + "_nad83.tif"
            dst = out_dir / out_name
            if dst.exists():
                continue
            warp_to_cog(chm, dst, target_crs)
            processed += 1

        if processed == 0:
            print(f"[SKIP] {site}: all targets already exist")
        else:
            print(f"[DONE] {site}: {processed} file(s) written\n")

    print("\n✔  Script completed.")


if __name__ == "__main__":
    main()
