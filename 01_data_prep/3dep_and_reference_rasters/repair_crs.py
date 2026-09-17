#!/usr/bin/env python
"""
repair_site5_crs.py
===================

Assign EPSG:32617 (WGS 84 / UTM zone 17 N) to every CHM GeoTIFF in
D:\\chms\\5\\chm\\ that lacks a valid CRS or has the wrong one.

Key points
----------
* Uses the GDAL Python API only – no external subprocess calls.
* Opens COGs with `IGNORE_COG_LAYOUT_BREAK=YES`, so the files can be
  updated in‑place even though the operation technically breaks the
  strict COG layout.  The resulting file is a **valid GeoTIFF**; you may
  run `gdal_translate -of COG` later to restore full cloud‑optimisation.
* Prints a clear report of how many files were modified / skipped.

Usage
-----
    (gedi_env) > python repair_site5_crs.py
    # optional
    (gedi_env) > python repair_site5_crs.py --site-dir "D:\\chms\\5\\chm" --epsg 32617
"""

from pathlib import Path
import sys
import argparse

from osgeo import gdal, osr

# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------
def build_target_srs(epsg: int) -> osr.SpatialReference:
    """Return an osr.SpatialReference instance for *epsg*."""
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    return srs


def dataset_has_same_srs(ds: gdal.Dataset, target_srs: osr.SpatialReference) -> bool:
    """Return True if *ds* already matches *target_srs* (axis order ignored)."""
    if not ds:
        return False
    wkt = ds.GetProjectionRef()
    if not wkt.strip():
        return False
    src_srs = osr.SpatialReference()
    src_srs.ImportFromWkt(wkt)
    return bool(src_srs.IsSame(target_srs))


def assign_srs_inplace(tif_path: Path, target_srs: osr.SpatialReference) -> None:
    """
    Open *tif_path* with update permission (GDAL OF_UPDATE) and write
    *target_srs* to the file.  Uses IGNORE_COG_LAYOUT_BREAK=YES so COGs
    can be modified.
    """
    ds = gdal.OpenEx(
        str(tif_path),
        gdal.OF_UPDATE,
        open_options=["IGNORE_COG_LAYOUT_BREAK=YES"],
    )
    if ds is None:
        raise RuntimeError(f"Cannot open {tif_path} for update")

    ds.SetProjection(target_srs.ExportToWkt())
    ds.FlushCache()  # write metadata
    ds = None        # close


# ----------------------------------------------------------------------
# Main routine
# ----------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Repair missing/incorrect CRS headers for Site 5 CHMs"
    )
    parser.add_argument(
        "--site-dir",
        type=Path,
        default=Path(r"D:\chms\5\chm"),
        help="Directory containing the Site 5 CHM GeoTIFFs",
    )
    parser.add_argument(
        "--epsg",
        type=int,
        default=32617,
        help="EPSG code to assign (default: 32617 = WGS84 / UTM 17N)",
    )
    args = parser.parse_args()

    if not args.site_dir.exists():
        sys.exit(f"{args.site_dir} does not exist.")

    tif_files = sorted(args.site_dir.glob("*.tif"))
    if not tif_files:
        sys.exit(f"No .tif files found in {args.site_dir}")

    target_srs = build_target_srs(args.epsg)
    fixed, skipped = 0, 0

    print(f"🛠  Assigning EPSG:{args.epsg} to CHMs in {args.site_dir} …\n")

    for tif in tif_files:
        # First check if file already OK --------------------------------
        ds_ro = gdal.Open(str(tif), gdal.GA_ReadOnly)
        already_ok = dataset_has_same_srs(ds_ro, target_srs)
        ds_ro = None  # close read‑only handle

        if already_ok:
            print(f"[SKIP] {tif.name:35s} – already EPSG:{args.epsg}")
            skipped += 1
            continue

        try:
            assign_srs_inplace(tif, target_srs)
            print(f"[FIX]  {tif.name:35s} -> EPSG:{args.epsg}")
            fixed += 1
        except Exception as exc:
            print(f"[ERR ] {tif.name:35s} – {exc}")

    print(
        f"\n✔  Completed. {fixed} file(s) updated, {skipped} already correct "
        f"in {args.site_dir}."
    )


if __name__ == "__main__":
    main()
