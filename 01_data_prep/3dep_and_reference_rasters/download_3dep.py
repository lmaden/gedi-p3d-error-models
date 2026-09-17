from pykml import parser
from shapely.geometry import Polygon, Point, LineString, box
import requests
import os
import geopandas as gpd
import json
import datetime
import re


def sanitize_name(name):
    """Make a filesystem-safe name from a placemark name."""
    # replace spaces and illegal chars with underscore
    return re.sub(r"[^\w\-]+", "_", name)


def parse_kml(file_path):
    """
    Parse KML and return a list of dicts with 'name' and 'geometry'.
    """
    with open(file_path) as f:
        root = parser.parse(f).getroot()

    features = []
    placemarks = root.findall('.//{http://www.opengis.net/kml/2.2}Placemark')
    for idx, placemark in enumerate(placemarks, start=1):
        # get name or fallback to index
        name_elem = placemark.find('.//{http://www.opengis.net/kml/2.2}name')
        name = name_elem.text.strip() if name_elem is not None else f'feature_{idx}'

        # extract geometry
        geom = None
        if placemark.find('.//{http://www.opengis.net/kml/2.2}Polygon') is not None:
            coords = placemark.find('.//{http://www.opengis.net/kml/2.2}coordinates').text.strip().split()
            coords = [tuple(map(float, coord.split(','))) for coord in coords]
            geom = Polygon(coords)
        elif placemark.find('.//{http://www.opengis.net/kml/2.2}LineString') is not None:
            coords = placemark.find('.//{http://www.opengis.net/kml/2.2}coordinates').text.strip().split()
            coords = [tuple(map(float, coord.split(','))) for coord in coords]
            geom = LineString(coords)
        elif placemark.find('.//{http://www.opengis.net/kml/2.2}Point') is not None:
            coord = placemark.find('.//{http://www.opengis.net/kml/2.2}coordinates').text.strip().split(',')
            geom = Point(tuple(map(float, coord)))

        if geom is not None:
            features.append({'name': sanitize_name(name), 'geometry': geom})
    return features


def split_bbox(bbox, rows, cols):
    minx, miny, maxx, maxy = bbox
    width = (maxx - minx) / cols
    height = (maxy - miny) / rows
    bboxes = []
    for i in range(rows):
        for j in range(cols):
            new_minx = minx + j * width
            new_miny = miny + i * height
            new_maxx = new_minx + width
            new_maxy = new_miny + height
            bboxes.append((new_minx, new_miny, new_maxx, new_maxy))
    return bboxes


def find_dem_tiles(polygon, rows=2, cols=2):
    url = 'https://tnmaccess.nationalmap.gov/api/v1/products'
    bbox_segments = split_bbox(polygon.bounds, rows, cols)
    all_tiles = {}

    for bbox in bbox_segments:
        params = {
            'datasets': 'Digital Elevation Model (DEM) 1 meter',
            'bbox': ','.join(map(str, bbox)),
            'outputFormat': 'json',
            'max': 50,
            'offset': 0
        }

        while True:
            response = requests.get(url, params=params)
            response.raise_for_status()
            data = response.json()
            tiles = data.get('items', [])
            if not tiles:
                break

            for tile in tiles:
                bb = tile.get('boundingBox', {})
                footprint = box(bb['minX'], bb['minY'], bb['maxX'], bb['maxY'])
                if polygon.intersects(footprint):
                    pub_date = datetime.datetime.strptime(tile['publicationDate'], "%Y-%m-%d")
                    key = (bb['minX'], bb['minY'], bb['maxX'], bb['maxY'])
                    if key not in all_tiles or all_tiles[key]['publicationDate'] < pub_date:
                        all_tiles[key] = {'tile': tile, 'publicationDate': pub_date}

            params['offset'] += params['max']

    return [item['tile'] for item in all_tiles.values()]


def download_dem_files(dem_tiles, download_dir):
    os.makedirs(download_dir, exist_ok=True)
    for item in dem_tiles:
        url = item['downloadURL']
        out_file = os.path.join(download_dir, os.path.basename(url))
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(out_file, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        print(f'Downloaded {out_file}')


if __name__ == '____main__':
    # Site boundary KML (see docs/MANIFEST.md for the input inventory).
    kml_path = os.environ.get("SITES_KML", "simple_sites.kml")
    base_download_dir = 'D:/3dep'

    features = parse_kml(kml_path)
    for feat in features:
        name = feat['name']
        geom = feat['geometry']
        site_dir = os.path.join(base_download_dir, name)

        print(f"Processing feature: {name}")
        tiles = find_dem_tiles(geom)
        if tiles:
            download_dem_files(tiles, site_dir)
        else:
            print(f'No DEM tiles for feature: {name}')
