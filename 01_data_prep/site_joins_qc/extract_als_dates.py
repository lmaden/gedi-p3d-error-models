#!/usr/bin/env python3
"""
extract_als_dates.py

Pulls per-site ALS acquisition dates from calval_dates.csv and formats
them for Supplementary Table 2.

Usage:
    python extract_als_dates.py path/to/calval_dates.csv

Output: prints one line per manuscript site with:
  MS_id | tracker_id | site_name | n_tiles | n_unique_dates | formatted_date

Formatting rules:
  - Single date: "YYYY-MM-DD"
  - Multi-day within same month: "YYYY-MM (D1-D2)"  e.g. "2017-05 (9-12)"
  - Multi-month: "YYYY-MM-DD to YYYY-MM-DD"  (earliest to latest)
  - Cross-year: full "YYYY-MM-DD to YYYY-MM-DD"

If a site name doesn't match any row in the CSV, prints "[NOT FOUND]"
so you can add an alias to SITE_ALIASES below and rerun.
"""

import sys
import csv
from collections import defaultdict
from datetime import datetime

# ---------------------------------------------------------------
# Site list in manuscript order. Each entry:
#   (ms_id, tracker_id, manuscript_site_name, [aliases to match in CSV])
# The name_als column in calval_dates.csv may use a different string
# than our internal tracker site names; add aliases as needed.
# ---------------------------------------------------------------
SITES = [
    # ms, tracker, manuscript_name,  aliases (match any, case-insensitive, EXACT match on name_als)
    ( 1,  1, "usda_me",           ["usda_me"]),
    ( 2,  2, "nasa_howland",      ["nasa_howland"]),                  # excludes nasa_howland_2009
    ( 3,  3, "neon_sawb",         ["neon_sawb"]),
    ( 4,  4, "neon_harv2019",     ["neon_harv2019"]),                 # excludes older neon_harv
    ( 5,  5, "usda_sc",           ["usda_sc"]),
    ( 6,  6, "neon_jerc2021",     ["neon_jerc2021"]),                 # excludes older neon_jerc
    ( 7,  7, "neon_tall2021",     ["neon_tall2021"]),
    ( 8,  8, "neon_dela2021",     ["neon_dela2021"]),
    ( 9,  9, "neon_leno2021",     ["neon_leno2021"]),                 # excludes older neon_leno
    (10, 10, "neon_clbj2022",     ["neon_clbj2022"]),
    (11, 11, "neon_konz2020",     ["neon_konz2020"]),
    (12, 12, "neon_stei2022",     ["neon_stei2022"]),
    (13, 13, "neon_unde2022",     ["neon_unde2022"]),
    (14, 14, "neon_steicheq2022", ["neon_steicheq2022"]),
    (15, 15, "neon_wood2021",     ["neon_wood2021"]),
    (16, 17, "neon_ster2022",     ["neon_ster2022"]),                 # excludes older neon_ster
    (17, 18, "neon_cper2021",     ["neon_cper2021"]),                 # excludes older neon_cper
    (18, 19, "ltbmu_20180710",    ["ltbmu_20180710", "ltbmu"]),
    (19, 20, "neon_sjer",         ["neon_sjer"]),
]


def parse_date(s):
    """Try several common date formats; return datetime.date or None."""
    s = s.strip()
    for fmt in ("%m/%d/%Y", "%Y-%m-%d", "%m/%d/%y", "%Y/%m/%d", "%d/%m/%Y"):
        try:
            return datetime.strptime(s, fmt).date()
        except ValueError:
            pass
    return None


def format_range(dates):
    """Format a sorted list of unique date objects per the rules above."""
    if not dates:
        return "[no dates]"
    if len(dates) == 1:
        return dates[0].strftime("%Y-%m-%d")
    first, last = dates[0], dates[-1]
    if first.year == last.year and first.month == last.month:
        return f"{first.strftime('%Y-%m')} ({first.day}\u2013{last.day})"
    return f"{first.strftime('%Y-%m-%d')} to {last.strftime('%Y-%m-%d')}"


def match_any(name_als, aliases):
    """Strict case-insensitive equality match against any alias."""
    lo = name_als.lower().strip()
    return any(lo == a.lower() for a in aliases)


def main(csv_path):
    # Read the CSV and group dates by name_als
    by_site = defaultdict(list)
    with open(csv_path, "r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            nm = row.get("name_als", "").strip()
            ds = row.get("als_date", "").strip()
            if not nm or not ds:
                continue
            d = parse_date(ds)
            if d is not None:
                by_site[nm].append(d)

    # Report per manuscript site
    print(f"{'MS':>2}  {'trk':>3}  {'site_name':<20}  "
          f"{'n_tiles':>7}  {'n_dates':>7}  date")
    print("-" * 78)

    unmatched = []
    for ms_id, trk_id, ms_name, aliases in SITES:
        matched_keys = [k for k in by_site.keys() if match_any(k, aliases)]
        if not matched_keys:
            print(f"{ms_id:>2}  {trk_id:>3}  {ms_name:<20}  "
                  f"{'-':>7}  {'-':>7}  [NOT FOUND]")
            unmatched.append((ms_id, ms_name, aliases))
            continue

        all_dates = []
        for k in matched_keys:
            all_dates.extend(by_site[k])
        unique_dates = sorted(set(all_dates))
        n_tiles = len(all_dates)
        n_uniq = len(unique_dates)
        date_str = format_range(unique_dates)

        match_note = ""
        if len(matched_keys) > 1:
            match_note = f"  (matched: {', '.join(matched_keys)})"

        print(f"{ms_id:>2}  {trk_id:>3}  {ms_name:<20}  "
              f"{n_tiles:>7}  {n_uniq:>7}  {date_str}{match_note}")

    if unmatched:
        print("\n---")
        print(f"{len(unmatched)} site(s) did not match. Candidate keys in CSV:")
        keys = sorted(by_site.keys())
        for k in keys:
            print(f"  {k}")
        print("\nAdd the correct alias to SITE_ALIASES in this script and rerun.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python extract_als_dates.py path/to/calval_dates.csv",
              file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1])
