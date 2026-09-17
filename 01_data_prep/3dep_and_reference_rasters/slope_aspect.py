import os
import argparse
from osgeo import gdal

def process_file(in_fp, slope_root, aspect_root, input_root):
    # build the matching output paths
    rel = os.path.relpath(in_fp, input_root)
    slope_fp  = os.path.join(slope_root,  rel)
    aspect_fp = os.path.join(aspect_root, rel)
    os.makedirs(os.path.dirname(slope_fp),  exist_ok=True)
    os.makedirs(os.path.dirname(aspect_fp), exist_ok=True)

    ds = gdal.Open(in_fp)
    # Compute slope in degrees
    gdal.DEMProcessing(
        slope_fp,  ds, 'slope',
        format='GTiff',
        computeEdges=True
    )
    # Compute aspect in degrees
    gdal.DEMProcessing(
        aspect_fp, ds, 'aspect',
        format='GTiff',
        computeEdges=True
    )
    ds = None
    print(f"→ {in_fp} → slope: {slope_fp}, aspect: {aspect_fp}")

def main():
    p = argparse.ArgumentParser(
        description="Batch compute slope & aspect from DEMs via GDAL."
    )
    p.add_argument('--input_dir',  required=True, help='Root DEM folder (e.g. D:\\3dep)')
    p.add_argument('--slope_dir',  required=True, help='Where to save slope rasters (e.g. D:\\slope)')
    p.add_argument('--aspect_dir', required=True, help='Where to save aspect rasters (e.g. D:\\aspect)')
    args = p.parse_args()

    for root, _, files in os.walk(args.input_dir):
        for f in files:
            if f.lower().endswith('.tif'):
                process_file(
                    os.path.join(root, f),
                    args.slope_dir,
                    args.aspect_dir,
                    args.input_dir
                )

if __name__ == '__main__':
    main()