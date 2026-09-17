#!/usr/bin/env python3
"""
aggregate_p3d_metadata.py
Aggregate acquisition dates from P3D metadata.json files across all sites.

Run from local machine with access to P3D data on hard drive.

Usage:
  python aggregate_p3d_metadata.py --p3d_root F:/p3d

Output: p3d_date_summary.csv with acquisition date ranges per site
"""

import os
import json
import glob
import argparse
from datetime import datetime
import csv


def find_metadata_files(root_dir):
    """Find all metadata.json files in the P3D directory structure."""
    patterns = [
        os.path.join(root_dir, '*', 'metadata.json'),
        os.path.join(root_dir, '*', '*', 'metadata.json'),
        os.path.join(root_dir, '*', 'vricon_raster_*', 'metadata.json'),
        os.path.join(root_dir, 'polygon_*', 'metadata.json'),
        os.path.join(root_dir, '*', 'vricon_raster_*', '*', 'metadata.json'),
    ]
    
    found = set()
    for pattern in patterns:
        found.update(glob.glob(pattern))
    
    return sorted(found)


def extract_site_id(path, root_dir):
    """Try to extract site ID from path."""
    rel_path = os.path.relpath(path, root_dir)
    parts = rel_path.replace('\\', '/').split('/')
    
    for part in parts:
        # Try to find numeric site ID
        if part.isdigit():
            return int(part)
        # Try polygon_N format
        if part.startswith('polygon_'):
            try:
                return int(part.replace('polygon_', ''))
            except:
                pass
        # Try site_N format
        if part.startswith('site_'):
            try:
                return int(part.replace('site_', ''))
            except:
                pass
    
    return None


def parse_metadata(filepath):
    """Parse a P3D metadata.json file."""
    try:
        with open(filepath, 'r') as f:
            meta = json.load(f)
        
        acq = meta.get('acquisition', {})
        first_date = acq.get('first-date')
        last_date = acq.get('last-date')
        
        # Parse dates
        first_dt = None
        last_dt = None
        
        if first_date:
            try:
                first_dt = datetime.fromisoformat(first_date.replace('Z', '+00:00'))
            except:
                pass
        
        if last_date:
            try:
                last_dt = datetime.fromisoformat(last_date.replace('Z', '+00:00'))
            except:
                pass
        
        return {
            'product': meta.get('product', 'Unknown'),
            'post_spacing': meta.get('post-spacing', None),
            'create_date': meta.get('create-date', None),
            'first_date': first_dt,
            'last_date': last_dt,
            'h_datum': meta.get('srs', {}).get('horizontal-datum', None),
            'v_datum': meta.get('srs', {}).get('vertical-datum', None),
        }
    
    except Exception as e:
        print(f"  Error parsing {filepath}: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(description='Aggregate P3D metadata')
    parser.add_argument('--p3d_root', required=True, help='Root directory of P3D data (e.g., F:/p3d)')
    args = parser.parse_args()
    
    print("=" * 60)
    print("P3D Metadata Aggregator")
    print("=" * 60)
    print(f"\nSearching in: {args.p3d_root}")
    print()
    
    # Find metadata files
    meta_files = find_metadata_files(args.p3d_root)
    print(f"Found {len(meta_files)} metadata.json files\n")
    
    if not meta_files:
        # Manual fallback - check direct subdirectories
        print("Trying direct subdirectory search...")
        for item in os.listdir(args.p3d_root):
            item_path = os.path.join(args.p3d_root, item)
            if os.path.isdir(item_path):
                for root, dirs, files in os.walk(item_path):
                    if 'metadata.json' in files:
                        meta_files.append(os.path.join(root, 'metadata.json'))
        print(f"Found {len(meta_files)} metadata.json files\n")
    
    results = []
    all_first = []
    all_last = []
    
    for mf in meta_files:
        site_id = extract_site_id(mf, args.p3d_root)
        meta = parse_metadata(mf)
        
        if meta:
            print(f"Site {site_id if site_id else '?':3}: "
                  f"{meta['first_date'].strftime('%Y-%m-%d') if meta['first_date'] else 'N/A'} to "
                  f"{meta['last_date'].strftime('%Y-%m-%d') if meta['last_date'] else 'N/A'} "
                  f"({meta['product']})")
            
            results.append({
                'site_id': site_id,
                'path': mf,
                'product': meta['product'],
                'post_spacing': meta['post_spacing'],
                'first_date': meta['first_date'].strftime('%Y-%m-%d') if meta['first_date'] else None,
                'last_date': meta['last_date'].strftime('%Y-%m-%d') if meta['last_date'] else None,
                'h_datum': meta['h_datum'],
                'v_datum': meta['v_datum'],
            })
            
            if meta['first_date']:
                all_first.append(meta['first_date'])
            if meta['last_date']:
                all_last.append(meta['last_date'])
    
    # Write results
    output_file = 'p3d_date_summary.csv'
    if results:
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['site_id', 'product', 'post_spacing', 
                                                    'first_date', 'last_date', 
                                                    'h_datum', 'v_datum', 'path'])
            writer.writeheader()
            for r in sorted(results, key=lambda x: x['site_id'] if x['site_id'] else 999):
                writer.writerow(r)
    
    print()
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    if all_first and all_last:
        earliest = min(all_first)
        latest = max(all_last)
        print(f"P3D source imagery date range: {earliest.strftime('%Y-%m-%d')} to {latest.strftime('%Y-%m-%d')}")
        print(f"  Earliest acquisition: {earliest.strftime('%B %d, %Y')}")
        print(f"  Latest acquisition: {latest.strftime('%B %d, %Y')}")
    
    print(f"\nResults saved to: {output_file}")


if __name__ == '__main__':
    main()
