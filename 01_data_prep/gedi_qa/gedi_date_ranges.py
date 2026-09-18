import geopandas as gpd
from datetime import datetime, timedelta
import os
from pathlib import Path

# GEDI epoch: January 1, 2018 00:00:00 UTC
GEDI_EPOCH = datetime(2018, 1, 1, 0, 0, 0)

def delta_time_to_datetime(delta_time):
    """Convert GEDI delta_time (seconds since epoch) to datetime"""
    return GEDI_EPOCH + timedelta(seconds=delta_time)

def process_site(site_num, base_path):
    """Process a single site's geopackage"""
    site_folder = os.path.join(base_path, str(site_num))
    gpkg_path = os.path.join(site_folder, f"GEDI_site{site_num}_hq_ALL.gpkg")
    
    if not os.path.exists(gpkg_path):
        print(f"⚠️  Site {site_num}: Geopackage not found at {gpkg_path}")
        return None
    
    try:
        print(f"📂 Processing Site {site_num}...")
        
        # Read the geopackage
        gdf = gpd.read_file(gpkg_path)
        
        # Check if delta_time column exists
        if 'delta_time' not in gdf.columns:
            print(f"⚠️  Site {site_num}: 'delta_time' column not found")
            return None
        
        # Get min and max delta_time
        min_delta = gdf['delta_time'].min()
        max_delta = gdf['delta_time'].max()
        
        # Convert to datetime
        earliest_date = delta_time_to_datetime(min_delta)
        latest_date = delta_time_to_datetime(max_delta)
        
        # Get total number of shots
        num_shots = len(gdf)
        
        result = {
            'site': site_num,
            'earliest': earliest_date,
            'latest': latest_date,
            'num_shots': num_shots,
            'output_folder': site_folder
        }
        
        print(f"✓ Site {site_num}: {earliest_date.strftime('%Y-%m-%d')} to {latest_date.strftime('%Y-%m-%d')} ({num_shots:,} shots)")
        
        return result
        
    except Exception as e:
        print(f"❌ Site {site_num}: Error processing - {str(e)}")
        return None

def main():
    # Base path to sites folder
    base_path = r"<LOCAL_DATA_ROOT>\gedi\sites"
    
    print("=" * 70)
    print("GEDI Date Range Extractor")
    print("=" * 70)
    print()
    
    # Auto-detect site numbers by scanning the folder
    site_numbers = []
    if os.path.exists(base_path):
        for item in os.listdir(base_path):
            item_path = os.path.join(base_path, item)
            if os.path.isdir(item_path) and item.isdigit():
                site_numbers.append(int(item))
    
    if not site_numbers:
        print("No site folders found. Please check the base path.")
        return
    
    site_numbers.sort()
    print(f"Found {len(site_numbers)} site folders: {site_numbers}")
    print()
    
    # Process each site
    results = []
    for site_num in site_numbers:
        result = process_site(site_num, base_path)
        if result:
            results.append(result)
        print()
    
    # Write summary file in each site folder
    print("=" * 70)
    print("Writing output files...")
    print("=" * 70)
    
    for result in results:
        output_file = os.path.join(result['output_folder'], 'date_range_summary.txt')
        
        with open(output_file, 'w') as f:
            f.write(f"GEDI Site {result['site']} - Date Range Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Geopackage: GEDI_site{result['site']}_hq_ALL.gpkg\n\n")
            f.write(f"Earliest Shot: {result['earliest'].strftime('%Y-%m-%d %H:%M:%S')} UTC\n")
            f.write(f"Latest Shot:   {result['latest'].strftime('%Y-%m-%d %H:%M:%S')} UTC\n\n")
            f.write(f"Date Range: {result['earliest'].strftime('%Y-%m-%d')} to {result['latest'].strftime('%Y-%m-%d')}\n")
            f.write(f"Total Shots: {result['num_shots']:,}\n")
            
            # Calculate duration
            duration = result['latest'] - result['earliest']
            f.write(f"Duration: {duration.days} days\n")
        
        print(f"✓ Written: {output_file}")
    
    # Also create a master summary
    master_file = os.path.join(base_path, 'all_sites_date_summary.txt')
    with open(master_file, 'w') as f:
        f.write("GEDI Date Range Summary - All Sites\n")
        f.write("=" * 70 + "\n\n")
        
        for result in results:
            f.write(f"Site {result['site']}:\n")
            f.write(f"  Earliest: {result['earliest'].strftime('%Y-%m-%d %H:%M:%S')} UTC\n")
            f.write(f"  Latest:   {result['latest'].strftime('%Y-%m-%d %H:%M:%S')} UTC\n")
            f.write(f"  Shots:    {result['num_shots']:,}\n")
            f.write("\n")
    
    print(f"\n✓ Master summary written: {master_file}")
    print("\n" + "=" * 70)
    print("✓ Done!")
    print("=" * 70)

if __name__ == "__main__":
    main()
