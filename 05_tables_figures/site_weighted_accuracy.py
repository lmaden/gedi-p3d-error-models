"""site_weighted_accuracy.py -- Table 2, the Fig. 2 panel statistics and the tail statistics of Results 3.1.

Accuracy statistics with every site weighted equally: each footprint gets weight 1/n_site, so each of
the 18 sites carries the same total weight. Computed from the analysis frame and the lidar CHM means:
  - Table 2 rows: bias, SD, RMSE, 95th percentile of |error|, % |error| > 5 m and > 10 m.
    "Sites weighted equally" and "All footprints" are the two halves of Table 2; "Median site" is
    written for reference only.
  - Fig. 2 panel statistics: bias, RMSE and Pearson r, weighted (and unweighted).
  - Tail statistics: share beyond 3 SD (and multiple of the normal 0.27%), share of the error
    variance held by the 2% of weight with the largest squared deviations.
Inputs  (P3D_DATA_ROOT) rework/frame_footprints.csv.gz; lidar_products/lp_site_NN.csv.gz (chm_mean)
Outputs (P3D_RESULTS)   table2/table2_site_weighted.csv, table2/fig2_panel_stats.csv, table2/error_tails.csv
"""
import glob
import os
import numpy as np
import pandas as pd

DATA_ROOT = os.environ.get("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>")
OUT_ROOT = os.environ.get("P3D_RESULTS", "results")
OUT = os.path.join(OUT_ROOT, "table2")
os.makedirs(OUT, exist_ok=True)

f = pd.read_csv(os.path.join(DATA_ROOT, "rework", "frame_footprints.csv.gz"),
                dtype={"shot_number": str})
n_site = f.groupby("site")["site"].transform("size")
f["w"] = 1.0 / n_site
print(f"footprints {len(f):,}, sites {f.site.nunique()}, weight sum {f.w.sum():.6f} (= sites)")


def wq(x, w, q):
    o = np.argsort(x)
    x, w = np.asarray(x)[o], np.asarray(w)[o]
    c = (np.cumsum(w) - 0.5 * w) / w.sum()
    return float(np.interp(q, c, x))


def stats(e, w):
    b = np.average(e, weights=w)
    sd = np.sqrt(np.average((e - b) ** 2, weights=w))
    rmse = np.sqrt(np.average(e ** 2, weights=w))
    a = np.abs(e)
    return dict(bias=b, sd=sd, rmse=rmse, q95=wq(a, w, 0.95),
                over5=100 * np.average(a > 5, weights=w), over10=100 * np.average(a > 10, weights=w))


rows = []
for prod in ("dsm", "dtm", "chm"):
    e = f[f"{prod}_error_mean"].to_numpy()
    rows.append(dict(frame="Sites weighted equally", product=prod.upper(), **stats(e, f.w.to_numpy())))
    s = stats(e, np.ones_like(e))
    s["q95"] = float(np.quantile(np.abs(e), 0.95))  # numpy's default quantile, the same as R's type 7
    rows.append(dict(frame="All footprints", product=prod.upper(), **s))
    per = pd.DataFrame([dict(**stats(g[f"{prod}_error_mean"].to_numpy(), np.ones(len(g))),
                             q95_r=float(np.quantile(np.abs(g[f"{prod}_error_mean"]), 0.95)))
                        for _, g in f.groupby("site")])
    per["q95"] = per["q95_r"]
    rows.append(dict(frame="Median site", product=prod.upper(),
                     **per[["bias", "sd", "rmse", "q95", "over5", "over10"]].median().to_dict()))
t2 = pd.DataFrame(rows)
t2.to_csv(os.path.join(OUT, "table2_site_weighted.csv"), index=False)
pd.set_option("display.width", 200)
print("\nTABLE 2")
print(t2.round(2).to_string(index=False))

# Fig. 2 panel statistics (weighted Pearson r)
lp = []
for p in sorted(glob.glob(os.path.join(DATA_ROOT, "lidar_products", "lp_site_*.csv.gz"))):
    d = pd.read_csv(p, usecols=["shot_number", "chm_mean"], dtype={"shot_number": str})
    lp.append(d)
lp = pd.concat(lp).dropna()
j = f.merge(lp, on="shot_number", how="inner")
print(f"\njoined with lidar CHM means: {len(j):,}")
j_n = j.groupby("site")["site"].transform("size")
jw = (1.0 / j_n).to_numpy()


def wr(x, y, w):
    mx, my = np.average(x, weights=w), np.average(y, weights=w)
    return np.average((x - mx) * (y - my), weights=w) / np.sqrt(
        np.average((x - mx) ** 2, weights=w) * np.average((y - my) ** 2, weights=w))


hl = j.chm_mean.to_numpy()
fig = []
for panel, err in (("a P3D CHM", "chm_error_mean"), ("b DSM minus lidar ground", "dsm_error_mean")):
    e = j[err].to_numpy()
    hp = hl + e
    for lab, w in (("Sites weighted equally", jw), ("All footprints", np.ones_like(jw))):
        fig.append(dict(panel=panel, frame=lab, bias=np.average(e, weights=w),
                        rmse=np.sqrt(np.average(e ** 2, weights=w)), r=wr(hp, hl, w)))
fig = pd.DataFrame(fig)
fig.to_csv(os.path.join(OUT, "fig2_panel_stats.csv"), index=False)
print("\nFIG. 2 PANEL STATISTICS")
print(fig.round(3).to_string(index=False))

# tail statistics (Results 3.1)
tails = []
w = f.w.to_numpy()
for prod in ("dsm", "dtm", "chm"):
    e = f[f"{prod}_error_mean"].to_numpy()
    for lab, ww in (("Sites weighted equally", w), ("All footprints", np.ones_like(w))):
        m = np.average(e, weights=ww)
        sd = np.sqrt(np.average((e - m) ** 2, weights=ww))
        beyond = np.average(np.abs(e - m) > 3 * sd, weights=ww)
        d2 = (e - m) ** 2
        o = np.argsort(-d2)
        cw = np.cumsum(ww[o]) / ww.sum()
        top = o[cw <= 0.02]
        share_var = (ww[top] * d2[top]).sum() / (ww * d2).sum()
        tails.append(dict(product=prod.upper(), frame=lab, beyond3sd_pct=100 * beyond,
                          times_normal=beyond / 0.0026998, top2pct_var_share=100 * share_var))
tails = pd.DataFrame(tails)
tails.to_csv(os.path.join(OUT, "error_tails.csv"), index=False)
print("\nTAILS (Results 3.1)")
print(tails.round(2).to_string(index=False))
