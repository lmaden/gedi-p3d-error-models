import re, json, time, requests, pandas as pd
from pathlib import Path

# --- inputs you edit ---
TOKEN = "YOUR_BEARER_TOKEN"
PATH_TO_IDS = "cats_1_within_aoi.txt"         # your uploaded file
OUT_CSV = "maxar_catalog_metadata.csv"
BATCH_SIZE = 100                              # Discovery/eAPI returns up to 100 per page
BASE_URL = "https://api.maxar.com/discovery/v1/search"  # STAC Discovery search
# -----------------------

# Map collection -> friendly sensor label
SENSOR_LABEL = {
    "ge01": "GeoEye-1",
    "wv01": "WorldView-1",
    "wv02": "WorldView-2",
    "wv03-vnir": "WorldView-3 (VNIR)",
    "wv03-swir": "WorldView-3 (SWIR)",
    "wv04": "WorldView-4",
    "lg01": "WorldView Legion-1",
    "lg02": "WorldView Legion-2",
    "lg03": "WorldView Legion-3",
    "lg04": "WorldView Legion-4",
    "lg05": "WorldView Legion-5",
    "lg06": "WorldView Legion-6",
}

# Read and sanitize IDs from the text file
raw = Path(PATH_TO_IDS).read_text().splitlines()
# Keep only classic 16-hex-char Maxar inventory IDs (skip GUIDs like ...-inv)
valid = [ln.strip() for ln in raw if re.fullmatch(r"[0-9A-F]{16}", ln.strip())]
ids = sorted(set(valid))  # de-duplicate

session = requests.Session()
session.headers.update({
    "Authorization": f"Bearer {TOKEN}",
    "Content-Type": "application/json"
})

rows, not_found = [], []

def fetch_batch(id_batch):
    # POST body with ids filter
    payload = {"ids": id_batch, "limit": 100}  # page size hint
    r = session.post(BASE_URL, data=json.dumps(payload), timeout=60)
    r.raise_for_status()
    return r.json()

# Query in chunks
for i in range(0, len(ids), BATCH_SIZE):
    batch = ids[i:i+BATCH_SIZE]
    data = fetch_batch(batch)
    found_ids = set()
    for feat in data.get("features", []):
        fid = feat.get("id")
        coll = feat.get("collection")
        props = feat.get("properties", {})
        rows.append({
            "id": fid,
            "sensor_collection": coll,
            "sensor": SENSOR_LABEL.get(coll, coll),
            "acquisition_datetime": props.get("datetime"),
            "gsd": props.get("gsd"),
            "cloud_cover": props.get("eo:cloud_cover") or props.get("cloud_cover"),
            "off_nadir": props.get("view:off_nadir") or props.get("offNadirAngle"),
            "sun_elevation": props.get("view:sun_elevation") or props.get("sun_elevation"),
            "bbox": feat.get("bbox"),
        })
        if fid:
            found_ids.add(fid)
    missing = set(batch) - found_ids
    if missing:
        not_found.extend(sorted(missing))
    # polite pacing if needed
    time.sleep(0.1)

# Write results
df = pd.DataFrame(rows).sort_values(["sensor_collection","acquisition_datetime","id"])
df.to_csv(OUT_CSV, index=False)

# Optionally write not-found IDs
if not_found:
    Path("not_found_ids.txt").write_text("\n".join(not_found))

print(f"Wrote {len(df)} rows to {OUT_CSV}")
print(f"Not found: {len(not_found)} (see not_found_ids.txt)" if not_found else "All IDs resolved.")
