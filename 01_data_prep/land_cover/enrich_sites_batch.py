#!/usr/bin/env python3
"""
Enrich individual site CSVs with landcover (sites 2-4, 6-20 only)
"""
import subprocess
from pathlib import Path

SITE_CSV_DIR = Path(r"<LOCAL_DATA_ROOT>\sampled_data_FINAL_by_site")
GEDI_ROOT = Path(r"<LOCAL_DATA_ROOT>\gedi\sites")
LC_ROOT = Path(r"<LOCAL_DATA_ROOT>\general\lc")

# Sites to process (exclude 1 and 5)
SITES = [2, 3, 4] + list(range(6, 21))

for site_id in SITES:
    site_csv = SITE_CSV_DIR / f"site_{site_id:02d}.csv"
    out_csv = SITE_CSV_DIR / f"site_{site_id:02d}_with_lc.csv"
    
    if not site_csv.exists():
        print(f"⚠️  Skipping site {site_id}: CSV not found")
        continue
    
    print(f"\n{'='*60}")
    print(f"Processing Site {site_id}")
    print(f"{'='*60}")
    
    cmd = [
        "python", "sample_glc_fcs30d_2022_to_csv.py",
        "--merged-csv", str(site_csv),
        "--gedi-sites-root", str(GEDI_ROOT),
        "--lc-root", str(LC_ROOT),
        "--sites", str(site_id),
        "--out-csv", str(out_csv),
        "--chunksize", "500000"
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"✔ Site {site_id} complete")
    except subprocess.CalledProcessError as e:
        print(f"✖ Site {site_id} failed: {e}")

print("\n" + "="*60)
print("ALL SITES COMPLETE")
print("="*60)