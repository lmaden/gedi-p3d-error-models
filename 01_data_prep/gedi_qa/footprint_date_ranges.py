import geopandas as gpd
from datetime import datetime
import os
from pathlib import Path
import re

def parse_ts_utc(ts_value):
    """Parse ts_utc value - handles both datetime objects and strings"""
    # If it's already a datetime object, return it directly
    if isinstance(ts_value, datetime):
        return ts_value
    
    # If it's a pandas Timestamp, convert to datetime
    if hasattr(ts_value, 'to_pydatetime'):
        return ts_value.to_pydatetime()
    
    # If it's a string, parse it
    if isinstance(ts_value, str):
        date_part = ts_value.replace('(UTC)', '').strip()
        return datetime.strptime(date_part, '%m/%d/%Y %H:%M:%S')
    
    # Otherwise, try to convert it
    return datetime.fromisoformat(str(ts_value))

def process_footprint(site_num, input_path, output_base_path):
    """Process a single site's footprint geopackage"""
    
    # Construct input file path with zero-padded site number
    site_num_str = str(site_num).zfill(2)  # e.g., 1 -> 01, 5 -> 05
    gpkg_path = os.path.join(input_path, f"site_{site_num_str}_footprints.gpkg")
    
    if not os.path.exists(gpkg_path):
        print(f"⚠️  Site {site_num}: Footprint geopackage not found at {gpkg_path}")
        return None
    
    try:
        print(f"📂 Processing Site {site_num} footprints...")
        
        # Read the geopackage
        gdf = gpd.read_file(gpkg_path)
        
        # Check if ts_utc column exists
        if 'ts_utc' not in gdf.columns:
            print(f"⚠️  Site {site_num}: 'ts_utc' column not found")
            return None
        
        # Parse all dates
        dates = [parse_ts_utc(ts) for ts in gdf['ts_utc']]
        
        # Get min and max dates
        earliest_date = min(dates)
        latest_date = max(dates)
        
        # Get total number of footprints
        num_footprints = len(gdf)
        
        # Output folder (same as GEDI outputs)
        output_folder = os.path.join(output_base_path, str(site_num))
        
        result = {
            'site': site_num,
            'earliest': earliest_date,
            'latest': latest_date,
            'num_footprints': num_footprints,
            'output_folder': output_folder
        }
        
        print(f"✓ Site {site_num}: {earliest_date.strftime('%Y-%m-%d')} to {latest_date.strftime('%Y-%m-%d')} ({num_footprints:,} footprints)")
        
        return result
        
    except Exception as e:
        print(f"❌ Site {site_num}: Error processing - {str(e)}")
        return None

def main():
    # Path to footprint geopackages
    input_path = r"<LOCAL_DATA_ROOT>\sites"
    
    # Output path (same as GEDI outputs)
    output_base_path = r"<LOCAL_DATA_ROOT>\gedi\sites"
    
    print("=" * 70)
    print("Image Footprint Date Range Extractor")
    print("=" * 70)
    print()
    
    # Auto-detect site numbers by scanning for footprint files
    site_numbers = []
    if os.path.exists(input_path):
        for filename in os.listdir(input_path):
            # Match pattern: site_XX_footprints.gpkg
            match = re.match(r'site_(\d+)_footprints\.gpkg', filename)
            if match:
                site_num = int(match.group(1))  # Convert to int (removes leading zeros)
                site_numbers.append(site_num)
    
    if not site_numbers:
        print("No footprint geopackages found. Please check the input path.")
        return
    
    site_numbers.sort()
    print(f"Found {len(site_numbers)} footprint files for sites: {site_numbers}")
    print()
    
    # Process each site
    results = []
    for site_num in site_numbers:
        result = process_footprint(site_num, input_path, output_base_path)
        if result:
            results.append(result)
        print()
    
    # Write summary file in each site folder
    print("=" * 70)
    print("Writing output files...")
    print("=" * 70)
    
    for result in results:
        # Create output folder if it doesn't exist
        os.makedirs(result['output_folder'], exist_ok=True)
        
        output_file = os.path.join(result['output_folder'], 'footprint_date_range_summary.txt')
        
        with open(output_file, 'w') as f:
            f.write(f"Image Footprint Site {result['site']} - Date Range Summary\n")
            f.write("=" * 50 + "\n\n")
            site_num_str = str(result['site']).zfill(2)
            f.write(f"Geopackage: site_{site_num_str}_footprints.gpkg\n\n")
            f.write(f"Earliest Image: {result['earliest'].strftime('%Y-%m-%d %H:%M:%S')} UTC\n")
            f.write(f"Latest Image:   {result['latest'].strftime('%Y-%m-%d %H:%M:%S')} UTC\n\n")
            f.write(f"Date Range: {result['earliest'].strftime('%Y-%m-%d')} to {result['latest'].strftime('%Y-%m-%d')}\n")
            f.write(f"Total Footprints: {result['num_footprints']:,}\n")
            
            # Calculate duration
            duration = result['latest'] - result['earliest']
            f.write(f"Duration: {duration.days} days\n")
        
        print(f"✓ Written: {output_file}")
    
    # Also create a master summary
    master_file = os.path.join(output_base_path, 'all_sites_footprint_summary.txt')
    with open(master_file, 'w') as f:
        f.write("Image Footprint Date Range Summary - All Sites\n")
        f.write("=" * 70 + "\n\n")
        
        for result in results:
            f.write(f"Site {result['site']}:\n")
            f.write(f"  Earliest: {result['earliest'].strftime('%Y-%m-%d %H:%M:%S')} UTC\n")
            f.write(f"  Latest:   {result['latest'].strftime('%Y-%m-%d %H:%M:%S')} UTC\n")
            f.write(f"  Footprints: {result['num_footprints']:,}\n")
            f.write("\n")
    
    print(f"\n✓ Master summary written: {master_file}")
    print("\n" + "=" * 70)
    print("✓ Done!")
    print("=" * 70)

if __name__ == "__main__":
    main()
