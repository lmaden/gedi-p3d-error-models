#!/usr/bin/env python
"""
Batch‑convert selected raster folders to Cloud‑Optimized GeoTIFF (COG),
verify the new files, delete the originals **only after** a successful
check, and rename each COG so your downstream code still sees *.tif.

P3D rasters remain untouched: we simply do not include that directory
in the ROOTS list.

──────────────────────────────────────────────────────────────────────
EDIT SECTION 1 (ROOTS) so the list contains *only* directories you want
to convert—e.g. 3DEP DTMs, slope, aspect, GEDI CHMs.  Leave out E:\p3d.
──────────────────────────────────────────────────────────────────────
"""

import os
import random
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

import rasterio
from tqdm import tqdm


# ────────────────────────────────────────────────────────────────────
# 1)  DIRECTORIES TO CONVERT  (leave the P3D drive/path OUT)
# ────────────────────────────────────────────────────────────────────
ROOTS = [
    r"D:\slope",
    r"D:\aspect",
    r"D:\3dep",
    r"D:\chms",
]

# ────────────────────────────────────────────────────────────────────
# 2)  GDAL COG OPTIONS  (tweak if desired)
# ────────────────────────────────────────────────────────────────────
GDAL_OPTS = [
    "-of", "COG",
    "-co", "COMPRESS=LZW",
    "-co", "BLOCKSIZE=512",
    "-co", "OVERVIEWS=IGNORE_EXISTING",
]

# ────────────────────────────────────────────────────────────────────
# 3)  PARALLEL SETTINGS
# ────────────────────────────────────────────────────────────────────
MAX_WORKERS = os.cpu_count() or 8        # use all logical cores


# -------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------
def gdal_translate_to_cog(src: Path, dst: Path) -> None:
    """Run gdal_translate to create a COG copy."""
    cmd = ["gdal_translate", *GDAL_OPTS, str(src), str(dst)]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def verify_cog(cog: Path) -> bool:
    """Open the COG and read a small window to assert basic integrity."""
    try:
        with rasterio.open(cog) as src:
            h_win = min(100, src.height)
            w_win = min(100, src.width)
            src.read(1, window=((0, h_win), (0, w_win)))
        return True
    except Exception:
        return False


def process_one(src: Path) -> None:
    """
    1. Create <src>.cog.tif
    2. Verify it opens and reads
    3. Delete original
    4. Rename COG to original name
    """
    cog_tmp = src.with_suffix(".cog.tif")
    if cog_tmp.exists():                 # conversion already done
        return

    gdal_translate_to_cog(src, cog_tmp)

    if verify_cog(cog_tmp):
        src.unlink()                     # remove original safely
        cog_tmp.rename(src)              # keep same basename *.tif
    else:
        cog_tmp.unlink(missing_ok=True)
        raise RuntimeError(f"Verification failed for {cog_tmp}")


# -------------------------------------------------------------------
# Main
# -------------------------------------------------------------------
if __name__ == "__main__":
    # Find all *.tif that are NOT already COG outputs
    tiff_paths = [
        p for root in ROOTS
        for p in Path(root).rglob("*.tif")
        if not p.name.endswith(".cog.tif")
    ]

    if not tiff_paths:
        print("No TIFFs found under the specified ROOTS.")
        raise SystemExit(0)

    print(f"→ Found {len(tiff_paths):,} TIFFs to convert")
    random.shuffle(tiff_paths)  # better load‑balance among workers

    # Process in parallel with a progress bar
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        for _ in tqdm(pool.map(process_one, tiff_paths),
                      total=len(tiff_paths),
                      desc="Converting to COG",
                      unit="file"):
            pass

    print("\n✔ All conversions completed successfully")
