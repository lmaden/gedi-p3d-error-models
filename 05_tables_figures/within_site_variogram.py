"""within_site_variogram.py -- Fig. S2 and Text S3 numbers, the within-site error variogram.

Per site and product, error minus the site mean; classical semivariogram gamma(h) = mean of half squared
differences over within-site pairs in each lag bin, to 3 km; divided by the site's within-site variance
(1 = no correlation); median and quartiles over sites (sites with at least 100 pairs in the bin).
Correlation at lag h = 1 - gamma(h) / variance. Sites with more than CAP footprints are thinned at random
(seed SEED).

Inputs  (P3D_DATA_ROOT) rework/frame_footprints.csv.gz (site = tracker id)
                        lidar_products/lp_site_NN.csv.gz (NN = tracker id; x, y in site UTM, m)
Outputs (P3D_RESULTS)   figS2/variogram_persite.csv, figS2/variogram_summary.csv
"""
import os
import pathlib
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

DATA_ROOT = pathlib.Path(os.environ.get("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>"))
OUT = pathlib.Path(os.environ.get("P3D_RESULTS", "results")) / "figS2"
OUT.mkdir(parents=True, exist_ok=True)

SEED, CAP, CUT = 7, 30000, 3000.0
EDGES = np.array([30, 90, 150, 210, 270, 330, 450, 600, 800, 1000, 1500, 2000, 3000], float)

f = pd.read_csv(DATA_ROOT / "rework/frame_footprints.csv.gz", dtype={"shot_number": str},
                usecols=["site", "shot_number", "dsm_error_mean", "dtm_error_mean", "chm_error_mean"])
parts = []
for p in sorted((DATA_ROOT / "lidar_products").glob("lp_site_*.csv.gz")):
    d = pd.read_csv(p, dtype={"shot_number": str}, usecols=["shot_number", "x", "y"])
    d["site"] = int(p.name[8:10])
    parts.append(d)
xy = pd.concat(parts).drop_duplicates(["site", "shot_number"])
f = f.merge(xy, on=["site", "shot_number"], how="left", validate="one_to_one")
assert len(f) == 630390 and f.x.notna().all()
f["track"] = f.shot_number.astype("int64") // 10**11

rng = np.random.default_rng(SEED)
rows = []
for site, g in f.groupby("site"):
    if len(g) > CAP:
        g = g.sample(CAP, random_state=rng.integers(1 << 31))
    pts = g[["x", "y"]].to_numpy()
    pr = cKDTree(pts).query_pairs(CUT, output_type="ndarray")
    d = np.sqrt(((pts[pr[:, 0]] - pts[pr[:, 1]]) ** 2).sum(1))
    k = np.digitize(d, EDGES) - 1
    keep = (k >= 0) & (k < len(EDGES) - 1)
    pr, d, k = pr[keep], d[keep], k[keep]
    tr = g.track.to_numpy()
    same = tr[pr[:, 0]] == tr[pr[:, 1]]
    for prod in ("dsm", "dtm", "chm"):
        a = g[f"{prod}_error_mean"].to_numpy(float)
        a = a - a.mean()
        var = a.var()
        h2 = 0.5 * (a[pr[:, 0]] - a[pr[:, 1]]) ** 2
        n = np.bincount(k, minlength=len(EDGES) - 1)
        gs = np.bincount(k, weights=h2, minlength=len(EDGES) - 1)
        ds = np.bincount(k, weights=d, minlength=len(EDGES) - 1)
        ss = np.bincount(k, weights=same.astype(float), minlength=len(EDGES) - 1)
        for b in range(len(EDGES) - 1):
            if n[b]:
                rows.append(dict(tracker=site, product=prod, n_fp=len(g), bin_lo=EDGES[b], bin_hi=EDGES[b + 1],
                                 n_pairs=int(n[b]), h_mean=ds[b] / n[b], same_track=ss[b] / n[b],
                                 gamma_norm=gs[b] / n[b] / var))
v = pd.DataFrame(rows)
v.to_csv(OUT / "variogram_persite.csv", index=False, float_format="%.6f")
ok = v[v.n_pairs >= 100]
s = (ok.groupby(["product", "bin_lo", "bin_hi"])
     .agg(n_sites=("tracker", "nunique"), h_med=("h_mean", "median"), g_med=("gamma_norm", "median"),
          g_q1=("gamma_norm", lambda x: x.quantile(0.25)), g_q3=("gamma_norm", lambda x: x.quantile(0.75)),
          same_track_med=("same_track", "median"))
     .reset_index())
s["corr_med"] = 1 - s.g_med
s.to_csv(OUT / "variogram_summary.csv", index=False, float_format="%.4f")
pd.set_option("display.width", 200)
print(s.to_string(index=False))
