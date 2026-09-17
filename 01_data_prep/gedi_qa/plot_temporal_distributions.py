import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import os
import re
import numpy as np
import pandas as pd

# GEDI epoch: January 1, 2018 00:00:00 UTC
GEDI_EPOCH = datetime(2018, 1, 1, 0, 0, 0)

def delta_time_to_datetime(delta_time):
    """Convert GEDI delta_time (seconds since epoch) to datetime"""
    return GEDI_EPOCH + timedelta(seconds=delta_time)

def parse_ts_utc(ts_value):
    """Parse ts_utc value - handles both datetime objects and strings"""
    if isinstance(ts_value, datetime):
        return ts_value
    if hasattr(ts_value, 'to_pydatetime'):
        return ts_value.to_pydatetime()
    if isinstance(ts_value, str):
        date_part = ts_value.replace('(UTC)', '').strip()
        return datetime.strptime(date_part, '%m/%d/%Y %H:%M:%S')
    return datetime.fromisoformat(str(ts_value))

def load_gedi_dates(site_num, base_path):
    """Load GEDI shot dates for a site"""
    site_folder = os.path.join(base_path, str(site_num))
    gpkg_path = os.path.join(site_folder, f"GEDI_site{site_num}_hq_ALL.gpkg")
    
    if not os.path.exists(gpkg_path):
        return None
    
    try:
        gdf = gpd.read_file(gpkg_path)
        if 'delta_time' not in gdf.columns:
            return None
        
        dates = [delta_time_to_datetime(dt) for dt in gdf['delta_time']]
        return dates
    except Exception as e:
        print(f"⚠️  Error loading GEDI data for site {site_num}: {e}")
        return None

def load_footprint_dates(site_num, input_path):
    """Load footprint dates for a site"""
    site_num_str = str(site_num).zfill(2)
    gpkg_path = os.path.join(input_path, f"site_{site_num_str}_footprints.gpkg")
    
    if not os.path.exists(gpkg_path):
        return None
    
    try:
        gdf = gpd.read_file(gpkg_path)
        if 'ts_utc' not in gdf.columns:
            return None
        
        dates = [parse_ts_utc(ts) for ts in gdf['ts_utc']]
        return dates
    except Exception as e:
        print(f"⚠️  Error loading footprint data for site {site_num}: {e}")
        return None

def plot_site_distribution(site_num, gedi_dates, footprint_dates, output_folder):
    """Create distribution plot for a single site"""
    
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    fig.suptitle(f'Site {site_num} - Temporal Distribution', fontsize=16, fontweight='bold')
    
    # GEDI plot
    ax1 = axes[0]
    if gedi_dates:
        ax1.hist(gedi_dates, bins=50, color='#2E7D32', alpha=0.7, edgecolor='black')
        ax1.set_title(f'GEDI Shots (n={len(gedi_dates):,})', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Frequency', fontsize=11)
        ax1.grid(True, alpha=0.3)
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
    else:
        ax1.text(0.5, 0.5, 'No GEDI Data Available', 
                ha='center', va='center', fontsize=14, color='red')
        ax1.set_title('GEDI Shots', fontsize=12, fontweight='bold')
    
    # Footprint plot
    ax2 = axes[1]
    if footprint_dates:
        ax2.hist(footprint_dates, bins=50, color='#1565C0', alpha=0.7, edgecolor='black')
        ax2.set_title(f'Image Footprints (n={len(footprint_dates):,})', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Frequency', fontsize=11)
        ax2.set_xlabel('Date', fontsize=11)
        ax2.grid(True, alpha=0.3)
        ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        ax2.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
    else:
        ax2.text(0.5, 0.5, 'No Footprint Data Available', 
                ha='center', va='center', fontsize=14, color='red')
        ax2.set_title('Image Footprints', fontsize=12, fontweight='bold')
        ax2.set_xlabel('Date', fontsize=11)
    
    plt.tight_layout()
    
    # Save plot
    os.makedirs(output_folder, exist_ok=True)
    output_file = os.path.join(output_folder, f'site_{site_num}_temporal_distribution.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_file

def plot_combined_timeline(site_num, gedi_dates, footprint_dates, output_folder):
    """Create combined timeline plot showing both datasets"""
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    y_pos = 0
    colors = []
    labels = []
    
    if gedi_dates:
        # Plot GEDI as scatter points
        y_gedi = np.ones(len(gedi_dates)) * y_pos
        ax.scatter(gedi_dates, y_gedi, alpha=0.5, s=10, c='#2E7D32', label=f'GEDI Shots (n={len(gedi_dates):,})')
        y_pos += 1
    
    if footprint_dates:
        # Plot footprints as scatter points
        y_footprint = np.ones(len(footprint_dates)) * y_pos
        ax.scatter(footprint_dates, y_footprint, alpha=0.5, s=10, c='#1565C0', label=f'Image Footprints (n={len(footprint_dates):,})')
    
    ax.set_yticks([])
    ax.set_xlabel('Date', fontsize=12)
    ax.set_title(f'Site {site_num} - Combined Timeline', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3, axis='x')
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    plt.tight_layout()
    
    # Save plot
    os.makedirs(output_folder, exist_ok=True)
    output_file = os.path.join(output_folder, f'site_{site_num}_combined_timeline.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_file

def create_overview_plot(all_sites_data, output_path):
    """Create overview plot showing all sites"""
    
    # Count sites with data
    sites_with_data = [(site, data['gedi'], data['footprint']) 
                       for site, data in all_sites_data.items() 
                       if data['gedi'] is not None or data['footprint'] is not None]
    
    if not sites_with_data:
        print("No data to plot in overview")
        return None
    
    # Create subplots
    n_sites = len(sites_with_data)
    n_cols = 3
    n_rows = (n_sites + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 4*n_rows))
    fig.suptitle('Temporal Distribution Overview - All Sites', fontsize=16, fontweight='bold')
    
    if n_sites == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    for idx, (site_num, gedi_dates, footprint_dates) in enumerate(sites_with_data):
        ax = axes[idx]
        
        # Plot both datasets on same axis
        if gedi_dates:
            ax.hist(gedi_dates, bins=30, alpha=0.6, color='#2E7D32', label=f'GEDI (n={len(gedi_dates):,})')
        if footprint_dates:
            ax.hist(footprint_dates, bins=30, alpha=0.6, color='#1565C0', label=f'Footprints (n={len(footprint_dates):,})')
        
        ax.set_title(f'Site {site_num}', fontsize=11, fontweight='bold')
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=8)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # Hide unused subplots
    for idx in range(len(sites_with_data), len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    
    # Save plot
    output_file = os.path.join(output_path, 'all_sites_temporal_overview.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_file

def main():
    # Paths
    gedi_base_path = r"C:\Users\levim\OneDrive\Desktop\work\gedi\dissertation\chpt1\data\gedi\sites"
    footprint_path = r"C:\Users\levim\OneDrive\Desktop\work\gedi\dissertation\chpt1\data\sites"
    
    print("=" * 70)
    print("Temporal Distribution Plot Generator")
    print("=" * 70)
    print()
    
    # Find all sites
    site_numbers = set()
    
    # Find GEDI sites
    if os.path.exists(gedi_base_path):
        for item in os.listdir(gedi_base_path):
            item_path = os.path.join(gedi_base_path, item)
            if os.path.isdir(item_path) and item.isdigit():
                site_numbers.add(int(item))
    
    # Find footprint sites
    if os.path.exists(footprint_path):
        for filename in os.listdir(footprint_path):
            match = re.match(r'site_(\d+)_footprints\.gpkg', filename)
            if match:
                site_numbers.add(int(match.group(1)))
    
    if not site_numbers:
        print("No sites found!")
        return
    
    site_numbers = sorted(list(site_numbers))
    print(f"Found {len(site_numbers)} sites: {site_numbers}")
    print()
    
    # Process each site
    all_sites_data = {}
    
    for site_num in site_numbers:
        print(f"Processing Site {site_num}...")
        
        # Load data
        gedi_dates = load_gedi_dates(site_num, gedi_base_path)
        footprint_dates = load_footprint_dates(site_num, footprint_path)
        
        if gedi_dates is None and footprint_dates is None:
            print(f"  ⚠️  No data found for site {site_num}")
            continue
        
        # Store for overview plot
        all_sites_data[site_num] = {
            'gedi': gedi_dates,
            'footprint': footprint_dates
        }
        
        # Create output folder
        output_folder = os.path.join(gedi_base_path, str(site_num))
        os.makedirs(output_folder, exist_ok=True)
        
        # Generate plots
        dist_plot = plot_site_distribution(site_num, gedi_dates, footprint_dates, output_folder)
        print(f"  ✓ Distribution plot: {dist_plot}")
        
        timeline_plot = plot_combined_timeline(site_num, gedi_dates, footprint_dates, output_folder)
        print(f"  ✓ Timeline plot: {timeline_plot}")
        
        print()
    
    # Create overview plot
    print("Creating overview plot...")
    overview_plot = create_overview_plot(all_sites_data, gedi_base_path)
    if overview_plot:
        print(f"✓ Overview plot: {overview_plot}")
    
    print()
    print("=" * 70)
    print("✓ Done!")
    print("=" * 70)

if __name__ == "__main__":
    main()
