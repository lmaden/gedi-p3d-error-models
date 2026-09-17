#!/usr/bin/env python3
"""
sample_glc_fcs30d_2022_to_csv.py  —  memory‑safe (streaming) updater

Adds 2022 GLC_FCS30D land cover to your existing CSV(s) WITHOUT
loading the whole merged file into memory.

Two phases in one command:
  (A) For each site, sample land cover at GEDI shots and write a small
      per‑site lookup file:  site_<id>\lc2022_lookup_site<id>.(parquet|csv)
      Columns: shot_key, lc2022_id_center, lc2022_id_mode, lc2022_mode_frac,
               lc2022_l1_code, lc2022_l1_name
  (B) Stream through the big merged CSV in CHUNKS and left‑join those
      lookups chunk‑by‑chunk → write an updated merged CSV.

Sampling notes:
- Center‑pixel sampling (always).
- Optional 3×3 majority via --majority_3x3 (default OFF).
  (Good edge‑case guard; fast. Your 25 m circle is smaller than 1 pixel at 30 m,
   so center‑pixel is already a consistent proxy. Majority_3x3 is a conservative
   fallback if you want a local mode near pixel edges.)

Dependencies: pandas, geopandas, rasterio, numpy, (optional) pyarrow (for parquet)
"""

from __future__ import annotations
import argparse
from pathlib import Path
import sys
from typing import Dict, List

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.transform import rowcol

# ---------- Level‑1 mapping (from your legend snapshot) ----------
ID_TO_L1_CODE: Dict[int, str] = {
    # Cropland
    10: "RCP", 11: "RCP", 12: "RCP", 20: "ICP",
    # Forests
    51: "EBF", 52: "EBF", 61: "BDF", 62: "BDF",
    71: "ENF", 72: "ENF", 81: "DNF", 82: "DNF",
    91: "MFT", 92: "MFT",
    # Shrubland / Grassland
    120: "SHR", 121: "SHR", 122: "SHR",
    130: "GRS",
    # Tundra (lichens & mosses)
    140: "LMS",
    # Wetlands (inland / coastal)
    180: "IWL", 181: "IWL", 182: "IWL", 183: "IWL", 184: "IWL",
    185: "CWL", 186: "CWL", 187: "CWL", 189: "CWL",
    # Impervious / Sparse vegetation / Bare
    190: "IMP",
    150: "SVG", 151: "SVG", 152: "SVG", 153: "SVG",
    200: "BAL", 201: "BAL", 202: "BAL",
    # Water / Snow & Ice
    210: "WTR",
    220: "PSI",
}

L1_CODE_TO_NAME: Dict[str, str] = {
    "RCP": "Rainfed cropland",
    "ICP": "Irrigated cropland",
    "EBF": "Evergreen broadleaved forest",
    "BDF": "Deciduous broadleaved forest",
    "ENF": "Evergreen needleleaved forest",
    "DNF": "Deciduous needleleaved forest",
    "MFT": "Mixed-leaf forest",
    "SHR": "Shrubland",
    "GRS": "Grassland",
    "LMS": "Lichens and mosses (tundra)",
    "IWL": "Inland wetland",
    "CWL": "Coastal wetland",
    "IMP": "Impervious surface",
    "SVG": "Sparse vegetation",
    "BAL": "Bare areas",
    "WTR": "Water body",
    "PSI": "Permanent snow and ice",
}

# ---------- I/O helpers ----------

def site_csv_path(data_dir: Path, site_id: int) -> Path:
    p = data_dir / f"site_{site_id:02d}.csv"
    if not p.exists():
        p = data_dir / f"site_{site_id}.csv"
    return p

def site_lc_tif(lc_root: Path, site_id: int) -> Path:
    d = lc_root / f"site_{site_id:02d}"
    if not d.exists():
        d = lc_root / f"site_{site_id}"
    # expect exactly one tif produced by your fetch script
    cands = list(d.glob("*.tif"))
    if not cands:
        raise FileNotFoundError(f"No land‑cover tif under {d}")
    return cands[0]

def site_gpkg(gedi_root: Path, site_id: int) -> Path:
    return gedi_root / f"{site_id}" / f"GEDI_site{site_id}_hq_ALL.gpkg"

def site_lookup_path(lc_root: Path, site_id: int, as_parquet: bool) -> Path:
    d = lc_root / f"site_{site_id:02d}"
    if not d.exists():
        d = lc_root / f"site_{site_id}"
    d.mkdir(parents=True, exist_ok=True)
    ext = "parquet" if as_parquet else "csv"
    return d / f"lc2022_lookup_site{site_id}.{ext}"

def have_pyarrow() -> bool:
    try:
        import pyarrow  # noqa: F401
        return True
    except Exception:
        return False

# ---------- sampling ----------

def sample_center(ds: rasterio.io.DatasetReader, xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    """Center‑pixel class for each (x, y) in ds CRS."""
    vals = []
    for x, y in zip(xs, ys):
        for v in ds.sample([(x, y)]):
            vals.append(v[0])
    out = np.array(vals, dtype="float32")
    nodata = ds.nodata
    if nodata is not None:
        out[out == nodata] = np.nan
    return out

def sample_majority_3x3(ds: rasterio.io.DatasetReader, xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    """3×3 mode around center pixel (good edge guard; fast)."""
    vals = np.full(len(xs), np.nan, dtype="float32")
    Tinv = ~ds.transform
    nodata = ds.nodata
    for i, (x, y) in enumerate(zip(xs, ys)):
        c, r = rowcol(ds.transform, x, y)
        rr = slice(max(0, r-1), min(ds.height, r+2))
        cc = slice(max(0, c-1), min(ds.width,  c+2))
        arr = ds.read(1, window=((rr.start, rr.stop), (cc.start, cc.stop)))
        if nodata is not None:
            arr = np.where(arr == nodata, np.nan, arr)
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            continue
        uu, ccnts = np.unique(arr.astype(np.int64), return_counts=True)
        vals[i] = uu[np.argmax(ccnts)]
    return vals

def build_site_lookup(
    site_id: int,
    gedi_root: Path,
    lc_root: Path,
    majority_3x3: bool = False
) -> pd.DataFrame:
    """Return per‑site lookup table with LC classes for each shot."""
    tpath = site_lc_tif(lc_root, site_id)
    gpath = site_gpkg(gedi_root, site_id)

    # Read GEDI points
    try:
        g = gpd.read_file(gpath, columns=["geometry", "shot_number"])
    except TypeError:
        g = gpd.read_file(gpath)[["geometry", "shot_number"]]
    if g.crs is None:
        g = g.set_crs(4326)

    with rasterio.open(tpath) as ds:
        g = g.to_crs(ds.crs)
        xs = g.geometry.x.values.astype("float64")
        ys = g.geometry.y.values.astype("float64")

        center = sample_center(ds, xs, ys)
        if majority_3x3:
            mode3 = sample_majority_3x3(ds, xs, ys)
            mode_id = mode3
            mode_frac = np.full_like(mode3, np.nan, dtype="float32")  # not computed here
        else:
            # no majority window requested: treat center as mode proxy
            mode_id = center.copy()
            mode_frac = np.full_like(center, np.nan, dtype="float32")

    df = pd.DataFrame({
        "shot_key": g["shot_number"].astype("string"),
        "lc2022_id_center": pd.Series(center).astype("Int32"),
        "lc2022_id_mode":   pd.Series(mode_id).astype("Int32"),
        "lc2022_mode_frac": pd.Series(mode_frac).astype("Float32"),
    })
    # Map Level‑1 code/name from modal ID
    df["lc2022_l1_code"] = df["lc2022_id_mode"].map(lambda v: ID_TO_L1_CODE.get(int(v)) if pd.notna(v) else pd.NA)
    df["lc2022_l1_name"] = df["lc2022_l1_code"].map(lambda c: L1_CODE_TO_NAME.get(c) if pd.notna(c) else pd.NA)
    return df

# ---------- streaming merge ----------

def stream_update_merged(
    merged_csv: Path,
    out_csv: Path,
    site_ids: List[int],
    lc_root: Path,
    chunksize: int = 1_000_000,
):
    """
    Read the big merged CSV in chunks, attach per‑site LC lookups, and append.
    """
    as_parquet = have_pyarrow()
    out_csv.unlink(missing_ok=True)
    header_written = False

    # pre-resolve lookup paths
    lk_paths: Dict[int, Path] = {
        s: site_lookup_path(lc_root, s, as_parquet) for s in site_ids
    }

    # sanity: ensure all lookups exist
    missing = [s for s, p in lk_paths.items() if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing per‑site lookup files for sites: "
            + ", ".join(map(str, missing))
            + ". Run this script once without --update-only to build them."
        )

    # streaming
    chunk_iter = pd.read_csv(
        merged_csv,
        chunksize=chunksize,
        low_memory=False,
        dtype={"shot_number": "string"},
        memory_map=True,
    )

    for ich, chunk in enumerate(chunk_iter, start=1):
        # prepare join key
        if "site" not in chunk.columns or "shot_number" not in chunk.columns:
            raise ValueError("Merged CSV must have columns 'site' and 'shot_number'.")
        chunk["shot_key"] = chunk["shot_number"].astype("string")

        # attach LC per site present in this chunk
        out_pieces = []
        for s in chunk["site"].astype("Int64").dropna().unique():
            s = int(s)
            sub = chunk.loc[chunk["site"].astype("Int64") == s].copy()
            lk = (pd.read_parquet(lk_paths[s]) if as_parquet
                  else pd.read_csv(lk_paths[s], dtype={"shot_key": "string"}))
            sub = sub.merge(lk, on="shot_key", how="left")
            out_pieces.append(sub.drop(columns=["shot_key"]))
        if out_pieces:
            out_chunk = pd.concat(out_pieces, ignore_index=True)
        else:
            # no site match; just drop the helper column
            out_chunk = chunk.drop(columns=["shot_key"])

        # append to output
        out_chunk.to_csv(out_csv, mode="a", index=False, header=not header_written)
        header_written = True
        print(f"  wrote chunk {ich:,}")

# ---------- CLI ----------

def parse_sites(s: str) -> List[int]:
    s = s.strip()
    if "-" in s:
        a, b = s.split("-")
        return list(range(int(a), int(b) + 1))
    return [int(x) for x in s.split(",") if x.strip()]

def main():
    ap = argparse.ArgumentParser(description="Append 2022 GLC_FCS30D to existing CSVs (memory‑safe).")
    ap.add_argument("--merged-csv", required=True, help="Path to the big merged CSV.")
    ap.add_argument("--gedi-sites-root", required=True, help="Root: .../gedi/sites")
    ap.add_argument("--lc-root", required=True, help="Root where per‑site LC GeoTIFFs (and lookups) live.")
    ap.add_argument("--sites", default="1-20", help="Sites to process, e.g. '1-20' or '3,7,12'.")
    ap.add_argument("--radius-m", type=float, default=12.5, help="(kept for compatibility; not used unless majority_3x3).")
    ap.add_argument("--majority_3x3", action="store_true", help="Use local 3×3 majority around center pixel.")
    ap.add_argument("--chunksize", type=int, default=1_000_000, help="Rows per chunk when streaming merged CSV.")
    ap.add_argument("--out-csv", default=None, help="Output merged path (default: <merged>_with_lc.csv).")
    ap.add_argument("--build-only", action="store_true", help="Only build per‑site lookup files; do not update merged.")
    ap.add_argument("--update-only", action="store_true", help="Only update merged (assumes lookups already exist).")
    args = ap.parse_args()

    merged_csv = Path(args.merged_csv)
    gedi_root  = Path(args.gedi_sites_root)
    lc_root    = Path(args.lc_root)
    site_ids   = parse_sites(args.sites)
    out_csv    = Path(args.out_csv) if args.out_csv else merged_csv.with_name(merged_csv.stem + "_with_lc.csv")

    # Step A — build per‑site lookups (unless update‑only)
    if not args.update_only:
        as_parquet = have_pyarrow()
        for s in site_ids:
            print(f"[Site {s}] building lookup…")
            df = build_site_lookup(s, gedi_root, lc_root, majority_3x3=args.majority_3x3)
            lk_path = site_lookup_path(lc_root, s, as_parquet)
            if as_parquet:
                df.to_parquet(lk_path, index=False)
            else:
                df.to_csv(lk_path, index=False)
            print(f"[Site {s}] wrote {lk_path}")

    # Step B — stream‑update the merged CSV (unless build‑only)
    if not args.build_only:
        print(f"\nStreaming update → {out_csv}")
        stream_update_merged(merged_csv, out_csv, site_ids, lc_root, chunksize=args.chunksize)
        print(f"✔ Done: {out_csv}")

if __name__ == "__main__":
    # Safer CSV field size (Windows sometimes hits defaults on huge lines)
    pd.options.mode.copy_on_write = True
    main()

