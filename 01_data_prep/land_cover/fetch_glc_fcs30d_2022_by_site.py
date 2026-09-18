#!/usr/bin/env python3
"""
fetch_glc_fcs30d_2022_by_site.py

Streams 2022 GLC_FCS30D (30 m) land-cover COG, clips it by each site's
bounding-box KML, and writes a per-site GeoTIFF:

Output layout:
  C:\...\chpt1\data\general\lc\site_{ID}\glc_fcs30d_2022_site{ID}.tif

Requires: rasterio, geopandas, shapely, requests
"""

import os
from pathlib import Path
import json
import requests
import rasterio
from rasterio.mask import mask
from rasterio.enums import Resampling
import geopandas as gpd
from shapely.geometry import mapping

# ---------- Paths (edit if needed) ----------
SITES_DIR = Path(r"<LOCAL_DATA_ROOT>\gedi\sites")
OUT_ROOT  = Path(r"<LOCAL_DATA_ROOT>\general\lc")

# ---------- GLC_FCS30D 2022 asset (STAC + fallback URL) ----------
STAC_ITEM_2022 = "https://stac.openlandmap.org/lc_glc.fcs30d/lc_glc.fcs30d_20220101_20221231/lc_glc.fcs30d_20220101_20221231.json"
COG_FALLBACK   = "https://s3.openlandmap.org/arco/lc_glc.fcs30d_c_30m_s_20220101_20221231_go_epsg.4326_v20231026.tif"

def glc_2022_href() -> str:
    """Resolve the 2022 COG href from STAC, fallback to known URL."""
    try:
        r = requests.get(STAC_ITEM_2022, timeout=30)
        r.raise_for_status()
        item = r.json()
        assets = item.get("assets", {})
        # Prefer a "c" (classification) asset, else first .tif
        if "c" in assets and assets["c"].get("href", "").lower().endswith(".tif"):
            return assets["c"]["href"]
        for a in assets.values():
            href = a.get("href", "")
            if href.lower().endswith(".tif"):
                return href
    except Exception:
        pass
    return COG_FALLBACK  # robust fallback

GLC_HREF = glc_2022_href()

# Optimize remote COG reads
os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", "tif,TIF,ovr,msk,xml,json,geojson")

def read_kml_polygon(kml_path: Path):
    """Read KML and return a single polygon in EPSG:4326."""
    try:
        gdf = gpd.read_file(kml_path, driver="LIBKML")
    except Exception:
        gdf = gpd.read_file(kml_path)  # let GDAL auto-detect
    if gdf.empty:
        raise ValueError(f"No geometry in {kml_path}")
    gdf = gdf.to_crs(4326)
    geom = gdf.unary_union  # dissolve (in case of multi-features)
    return geom

def clip_site(site_id: int):
    kml = SITES_DIR / f"{site_id}" / f"site_{site_id}_bb.kml"
    if not kml.exists():
        print(f"↪ Site {site_id}: KML not found → {kml}")
        return
    out_dir = OUT_ROOT / f"site_{site_id}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tif = out_dir / f"glc_fcs30d_2022_site{site_id}.tif"

    geom = read_kml_polygon(kml)
    with rasterio.Env():
        with rasterio.open(GLC_HREF) as src:
            # Clip in source CRS (EPSG:4326)
            out_img, out_transform = mask(
                src, [mapping(geom)], crop=True, nodata=src.nodata
            )
            meta = src.meta.copy()
            meta.update(
                driver="GTiff",
                height=out_img.shape[1],
                width=out_img.shape[2],
                transform=out_transform,
                compress="lzw",
                tiled=True,
                bigtiff="IF_SAFER"
            )
            with rasterio.open(out_tif, "w", **meta) as dst:
                dst.write(out_img)
    print(f"✔ Site {site_id}: wrote {out_tif}")

def main():
    # Iterate numeric site folders
    site_ids = sorted(
        int(p.name) for p in SITES_DIR.iterdir() if p.is_dir() and p.name.isdigit()
    )
    print(f"Resolved GLC_FCS30D 2022 href: {GLC_HREF}")
    for sid in site_ids:
        try:
            clip_site(sid)
        except Exception as e:
            print(f"✖ Site {sid}: {e}")

if __name__ == "__main__":
    main()
