# tableS4_forest_rows.py  --  Table S4 forest-type rows, recomputed from the frame
# Forest types as in Table 1 (table1_sites.R: FOREST = BDF, ENF, DNF, EBF, MFT), so the all-forest row matches the
# 75.6% forest share in Section 2.1 and Table 1. The printed rows weight each site equally, counting only sites with at
# least 30 footprints of the type (rule of Table S3 and Fig. 3b): bias and share > 5 m are means of the site values,
# RMSE is the root of the mean of the site MSEs (Table 2 rule). Mixed forest gets no row of its own (only 2 sites reach
# 30 footprints) and stays inside "All forest types".
# Inputs  (P3D_DATA_ROOT) rework/frame_footprints.csv.gz (630,390 footprints)
#         (P3D_RESULTS)   persite/forest_class.csv (persite_class_tables.R, for the check)
# Outputs (P3D_RESULTS)   tableS4/tableS4_forest_rows.md and tableS4_forest_rows.csv (site-weighted rows, with site
#                         counts and pooled values for comparison); tableS4/forest_rows_pooled.md (pooled rows)
# Checks  the four classes of persite/forest_class.csv are reproduced; n sums to 630,390;
#         all-forest plus non-forest pooled CHM bias and RMSE reproduce Table 2, footprints weighted equally (−2.40, 5.13).
import gzip, csv, os, pathlib, math, sys
sys.stdout.reconfigure(encoding="utf-8")
DATA_ROOT = pathlib.Path(os.environ.get("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>"))
OUT_ROOT = pathlib.Path(os.environ.get("P3D_RESULTS", "results"))
OUT = OUT_ROOT / "tableS4"
OUT.mkdir(parents=True, exist_ok=True)
rows = {}
with gzip.open(DATA_ROOT / "rework" / "frame_footprints.csv.gz", "rt") as fh:
    for r in csv.DictReader(fh):
        rows.setdefault(r["lc_l1_code"], []).append((float(r["dsm_error_mean"]), float(r["dtm_error_mean"]), float(r["chm_error_mean"]), r["site"]))
FOREST = ["ENF", "BDF", "EBF", "DNF", "MFT"]
names = {"ENF": "Evergreen needleleaf forest", "BDF": "Deciduous broadleaf forest", "EBF": "Evergreen broadleaf forest",
         "DNF": "Deciduous needleleaf forest", "MFT": "Mixed forest", "ALLF": "All forest types", "NONF": "Non-forest classes"}
groups = {k: rows[k] for k in FOREST}
groups["ALLF"] = [x for k in FOREST for x in rows[k]]
groups["NONF"] = [x for k, v in rows.items() if k not in FOREST for x in v]
assert len(groups["ALLF"]) + len(groups["NONF"]) == 630390

def stats(v, j):
    e = [x[j] for x in v]; n = len(e)
    return sum(e) / n, math.sqrt(sum(x * x for x in e) / n), 100 * sum(abs(x) > 5 for x in e) / n

# check against persite/forest_class.csv (four classes)
fc = {(r["stratum"], r["product"]): r for r in csv.DictReader(open(OUT_ROOT / "persite" / "forest_class.csv", encoding="utf-8"))}
dmax = 0
for k in ["ENF", "BDF", "EBF", "DNF"]:
    for j, p in enumerate(["DSM", "DTM", "CHM"]):
        b, rm, o5 = stats(groups[k], j); ref = fc[(k, p)]
        dmax = max(dmax, abs(b - float(ref["mean"])), abs(rm - float(ref["rmse"])), abs(o5 - float(ref["pct_over5"])))
allv = groups["ALLF"] + groups["NONF"]; cb, cr, _ = stats(allv, 2)
print(f"GATE forest_class.csv four classes reproduced: max |diff| {dmax:.1e} {'PASS' if dmax < 1e-9 else 'FAIL'}; "
      f"all footprints CHM bias {cb:+.2f} RMSE {cr:.2f} (Table 2: -2.40, 5.13) {'PASS' if abs(cb+2.40)<0.005 and abs(cr-5.13)<0.005 else 'FAIL'}")
print(f"forest share {100*len(groups['ALLF'])/630390:.2f}% (Section 2.1: 75.6%)")
m = lambda x: f"{x:.2f}".replace("-", "−")
out = []
for k in FOREST + ["ALLF", "NONF"]:
    cells = []
    for j in range(3):
        b, rm, o5 = stats(groups[k], j); cells += [m(b), f"{rm:.2f}", f"{o5:.1f}"]
    out.append(f"| {names[k]} | " + " | ".join(cells) + " |")
    print(f"{names[k]:30s} n {len(groups[k]):8,d}")
dest = OUT / "forest_rows_pooled.md"
dest.write_text("\n".join(out) + "\n", encoding="utf-8")
print("\n".join(out))

# ---- each site weighted equally, sites with >= 30 footprints of the type
MIN_N = 30
def sw_stats(v, j):
    by = {}
    for x in v: by.setdefault(x[3], []).append(x[j])
    q = {s_: e for s_, e in by.items() if len(e) >= MIN_N}
    bias = [sum(e) / len(e) for e in q.values()]
    mse = [sum(x * x for x in e) / len(e) for e in q.values()]
    o5 = [100 * sum(abs(x) > 5 for x in e) / len(e) for e in q.values()]
    k = len(q)
    return sum(bias) / k, math.sqrt(sum(mse) / k), sum(o5) / k, k, sum(len(e) for e in q.values()), len(by)
SW_ROWS = ["ENF", "BDF", "EBF", "DNF", "ALLF", "NONF"]          # mixed forest only inside ALLF
csv_rows, out_sw = [], []
for k in SW_ROWS:
    cells = []
    for j, p in enumerate(["DSM", "DTM", "CHM"]):
        b, rm, o5, ns, n30, nany = sw_stats(groups[k], j); pb, prm, po5 = stats(groups[k], j)
        cells += [m(b), f"{rm:.2f}", f"{o5:.1f}"]
        csv_rows.append(dict(cls=names[k], product=p, n=len(groups[k]), n_at_counted_sites=n30, n_sites_any=nany, n_sites_30=ns,
                             sw30_bias=b, sw30_rmse=rm, sw30_o5=o5, pooled_bias=pb, pooled_rmse=prm, pooled_o5=po5))
    out_sw.append(f"| {names[k]} | " + " | ".join(cells) + " |")
dest2 = OUT / "tableS4_forest_rows.md"
dest2.write_text("\n".join(out_sw) + "\n", encoding="utf-8")
with open(OUT / "tableS4_forest_rows.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(csv_rows[0])); w.writeheader(); w.writerows(csv_rows)
print("\n".join(out_sw))
for k in SW_ROWS:
    r = next(x for x in csv_rows if x["cls"] == names[k])
    print(f"{names[k]:30s} sites >= {MIN_N}: {r['n_sites_30']:2d} of {r['n_sites_any']:2d}; footprints {r['n']:,} ({r['n_at_counted_sites']:,} at counted sites)")
mft = {}
for x in groups["MFT"]: mft[x[3]] = mft.get(x[3], 0) + 1
print("Mixed forest footprints by tracker site:", dict(sorted(mft.items(), key=lambda t: -t[1])))
