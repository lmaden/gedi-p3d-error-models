#!/usr/bin/env python3
"""
query_3dep_dates.py
Query USGS TNM API to get 3DEP acquisition/publication dates for each site.

Run from the directory containing your site boundary files or with paths specified.

Usage:
  python query_3dep_dates.py --kml_path /path/to/simple_sites.kml
  or
  python query_3dep_dates.py --bounds_csv /path/to/site_bounds.csv

Output: 3dep_date_summary.csv with dates per site
"""

import os
import re
import json
import requests
import argparse
from datetime import datetime
import csv

# Site bounding boxes (approximate) - UPDATE THESE from your data
# Format: site_id: (min_lon, min_lat, max_lon, max_lat)
# You can get these from your site shapefiles/KML
SITE_BOUNDS = {
    1: (-70.62, 45.27, -69.10, 46.82),   # ME usda
    2: (-68.78, 45.17, -68.74, 45.25),   # ME howland  
    3: (-71.26, 42.46, -71.11, 42.58),   # MA neon_sawb
    4: (-72.17, 42.38, -72.10, 42.58),   # MA neon_harv
    5: (-80.88, 32.54, -79.32, 33.98),   # SC usda
    6: (-84.57, 31.16, -84.35, 31.33),   # GA neon_jerc
    7: (-87.49, 32.86, -87.35, 32.99),   # AL neon_tall
    8: (-87.88, 32.49, -87.76, 32.61),   # AL neon_dela
    9: (-88.26, 31.78, -88.13, 31.89),   # AL neon_leno
    10: (-97.68, 33.32, -97.53, 33.42),  # TX neon_clbj
    11: (-96.67, 39.04, -96.48, 39.25),  # KS neon_konz
    12: (-89.62, 45.44, -89.38, 45.56),  # WI neon_stei
    13: (-89.60, 46.14, -89.42, 46.28),  # MN neon_unde
    14: (-90.06, 45.77, -90.03, 45.85),  # MN neon_steicheq
    15: (-99.31, 47.08, -99.03, 47.24),  # ND neon_wood
    # Site 16 excluded (non-forest)
    17: (-103.06, 40.46, -103.02, 40.52),  # CO neon_nogp -> becomes site 16
    18: (-104.75, 40.78, -104.71, 40.87),  # CO neon_ster -> becomes site 17
    19: (-120.28, 38.68, -119.87, 39.34),  # CA ltbmu -> becomes site 18
    20: (-119.80, 37.04, -119.66, 37.14),  # CA neon_sjer -> becomes site 19
}


def query_3dep_tiles(bbox, max_results=50):
    """Query USGS TNM API for 3DEP 1m DEM tiles within a bounding box."""
    url = 'https://tnmaccess.nationalmap.gov/api/v1/products'
    
    params = {
        'datasets': 'Digital Elevation Model (DEM) 1 meter',
        'bbox': ','.join(map(str, bbox)),
        'outputFormat': 'json',
        'max': max_results,
        'offset': 0
    }
    
    all_tiles = []
    
    while True:
        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            tiles = data.get('items', [])
            
            if not tiles:
                break
            
            all_tiles.extend(tiles)
            
            if len(tiles) < max_results:
                break
                
            params['offset'] += max_results
            
        except Exception as e:
            print(f"  Warning: API error - {e}")
            break
    
    return all_tiles


def extract_dates(tiles):
    """Extract publication and source dates from tiles."""
    pub_dates = []
    source_dates = []
    
    for tile in tiles:
        # Publication date
        pub = tile.get('publicationDate')
        if pub:
            try:
                pub_dates.append(datetime.strptime(pub, "%Y-%m-%d"))
            except:
                pass
        
        # Source/collection date (may be in different fields)
        src = tile.get('sourceDate') or tile.get('dateCreated') or tile.get('lastUpdated')
        if src:
            try:
                if 'T' in str(src):
                    source_dates.append(datetime.fromisoformat(src.replace('Z', '+00:00')))
                else:
                    source_dates.append(datetime.strptime(str(src)[:10], "%Y-%m-%d"))
            except:
                pass
    
    return pub_dates, source_dates


def main():
    print("=" * 60)
    print("3DEP Date Query Tool")
    print("=" * 60)
    print()
    
    results = []
    
    for site_id, bbox in sorted(SITE_BOUNDS.items()):
        print(f"Site {site_id}: Querying 3DEP tiles...")
        
        tiles = query_3dep_tiles(bbox)
        pub_dates, source_dates = extract_dates(tiles)
        
        result = {
            'site_id': site_id,
            'n_tiles': len(tiles),
            'pub_earliest': min(pub_dates).strftime('%Y-%m-%d') if pub_dates else 'N/A',
            'pub_latest': max(pub_dates).strftime('%Y-%m-%d') if pub_dates else 'N/A',
            'source_earliest': min(source_dates).strftime('%Y-%m-%d') if source_dates else 'N/A',
            'source_latest': max(source_dates).strftime('%Y-%m-%d') if source_dates else 'N/A',
        }
        results.append(result)
        
        print(f"  Found {len(tiles)} tiles")
        print(f"  Publication: {result['pub_earliest']} to {result['pub_latest']}")
        print()
    
    # Write results
    output_file = '3dep_date_summary.csv'
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['site_id', 'n_tiles', 'pub_earliest', 
                                                'pub_latest', 'source_earliest', 'source_latest'])
        writer.writeheader()
        writer.writerows(results)
    
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    # Overall date range
    all_pub = []
    for r in results:
        if r['pub_earliest'] != 'N/A':
            all_pub.append(datetime.strptime(r['pub_earliest'], '%Y-%m-%d'))
        if r['pub_latest'] != 'N/A':
            all_pub.append(datetime.strptime(r['pub_latest'], '%Y-%m-%d'))
    
    if all_pub:
        print(f"Overall 3DEP date range: {min(all_pub).strftime('%Y-%m-%d')} to {max(all_pub).strftime('%Y-%m-%d')}")
    
    print(f"\nResults saved to: {output_file}")


if __name__ == '__main__':
    main()
