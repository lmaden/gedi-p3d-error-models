#!/usr/bin/env python3
"""
resume_driver.py — Per-site checkpointing + resume for gedi_sample_v09_checkpointed.py

Usage examples:
  python resume_driver.py --resume
  python resume_driver.py --resume --sites 5,6,10
  python resume_driver.py --workers 6
"""

from __future__ import annotations
import argparse, os, re, json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Iterable, List, Set
import pandas as pd
from tqdm import tqdm

# Import from the v09 checkpointed sampler
from gedi_sample_v09_checkpointed import process_site, SITES as DEFAULT_SITES, OUTPUT_CSV, META_JSON

# Default checkpoint directory lives next to your final CSV
CHECKPOINT_DIR = OUTPUT_CSV.with_name(OUTPUT_CSV.stem + "_by_site")

def parse_sites_arg(arg: str, default_sites: Iterable[int]) -> List[int]:
    """
    Parse strings like '1-5,7,9-10' to a sorted unique site list.
    If arg is empty/None, return list(default_sites).
    """
    if not arg:
        return sorted(set(default_sites))
    sites: Set[int] = set()
    for part in arg.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            a = int(a); b = int(b)
            lo, hi = min(a, b), max(a, b)
            sites.update(range(lo, hi + 1))
        else:
            sites.add(int(part))
    return sorted(sites)

def checkpoint_path(site: int) -> Path:
    return CHECKPOINT_DIR / f"site_{site:02d}.csv"

def existing_checkpoint_sites() -> Set[int]:
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
    done: Set[int] = set()
    for p in CHECKPOINT_DIR.glob("site_*.csv"):
        m = re.match(r"site_(\d+)\.csv$", p.name)
        if not m:
            continue
        try:
            site = int(m.group(1))
        except ValueError:
            continue
        # Consider only non-empty files as valid checkpoints
        if p.stat().st_size > 0:
            done.add(site)
    return done

def assemble_final(output_csv: Path, meta_json: Path) -> int:
    files = sorted(CHECKPOINT_DIR.glob("site_*.csv"))
    if not files:
        raise RuntimeError(f"No per-site checkpoints found in {CHECKPOINT_DIR}")
    dfs = [pd.read_csv(fp) for fp in files]
    final_df = pd.concat(dfs, ignore_index=True)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(output_csv, index=False, float_format="%.4f")
    with open(meta_json, "w", encoding="utf8") as f:
        json.dump({"columns": list(final_df.columns)}, f, indent=2)
    return len(final_df)

def main():
    parser = argparse.ArgumentParser(description="Resume-able driver for GEDI sampling (v09)")
    parser.add_argument("--sites", type=str, default="", help="e.g., '1-20' or '6,8,12-15' (default: use SITES from gedi_sample_v09_checkpointed.py)")
    parser.add_argument("--resume", action="store_true", help="Skip sites with an existing checkpoint")
    parser.add_argument("--workers", type=int, default=int(os.environ.get("GEDI_MAX_WORKERS", "4")),
                        help="Number of processes (default: env GEDI_MAX_WORKERS or 4)")
    parser.add_argument("--no-assemble", action="store_true", help="Do not build the combined CSV at the end")
    args = parser.parse_args()

    all_sites = parse_sites_arg(args.sites, DEFAULT_SITES)
    done_sites = existing_checkpoint_sites() if args.resume else set()
    to_run = [s for s in all_sites if s not in done_sites]

    print(f"Checkpoint dir: {CHECKPOINT_DIR}")
    if done_sites:
        print(f"Already have checkpoints for: {sorted(done_sites)}")
    print(f"Will run sites: {to_run if to_run else '[] (nothing to do)'}")

    if not to_run and not args.no_assemble:
        n = assemble_final(OUTPUT_CSV, META_JSON)
        print(f"✔ Wrote {n:,} rows → {OUTPUT_CSV}")
        print(f"ℹ Column list → {META_JSON}")
        return

    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    with ProcessPoolExecutor(max_workers=max(1, min(args.workers, 8))) as ex:
        futs = {ex.submit(process_site, s): s for s in to_run}
        for fut in tqdm(as_completed(futs), total=len(futs), desc="Sites", unit="site"):
            site = futs[fut]
            try:
                site_id, df, log_str = fut.result()
                if log_str:
                    print(f"\n--- Log for Site {site_id} ---\n{log_str}\n--------------------------\n")
                if df is not None and not df.empty:
                    # Atomic write for the checkpoint
                    tmp = checkpoint_path(site_id).with_suffix(".csv.tmp")
                    df.to_csv(tmp, index=False, float_format="%.4f")
                    os.replace(tmp, checkpoint_path(site_id))
                    print(f"✔ Checkpoint saved → {checkpoint_path(site_id)} ({len(df):,} rows)")
                else:
                    print(f"⚠ Site {site_id}: no dataframe returned (skipped/empty).")
            except Exception as e:
                print(f"[Site {site}] ERROR: {e}")
                import traceback; traceback.print_exc()

    if not args.no_assemble:
        n = assemble_final(OUTPUT_CSV, META_JSON)
        print(f"✔ Wrote {n:,} rows → {OUTPUT_CSV}")
        print(f"ℹ Column list → {META_JSON}")

if __name__ == "__main__":
    main()
