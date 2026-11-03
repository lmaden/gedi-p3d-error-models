# tableS3_assemble.py  --  Table S3 as printed, by class definition
# (a) lidar classes: mean error, RMSE and |error| > 5 m for DSM, DTM, CHM, then Sites; (b) GEDI classes: mean error, then Sites.
# No "+" signs (as in Tables S4 and S5). Mean errors and site counts come from the Table S3 class means; RMSE cells are
# the root of the mean of the site MSEs (the Table 2 rule, rmse_rms) and |error| > 5 m the mean of the site shares.
# Inputs  (P3D_RESULTS) tableS3/tableS3_height_classes.md (tableS3_height_classes.R),
#                       tableS3/tableS3_rmse_over5.csv (tableS3_rmse_over5.R)
# Output  (P3D_RESULTS) tableS3/tableS3.md
# Check   the site counts of the two inputs agree for all three products.
import csv, os, pathlib, sys
sys.stdout.reconfigure(encoding="utf-8")
OUT = pathlib.Path(os.environ.get("P3D_RESULTS", "results")) / "tableS3"
def rows(p): return [[c.strip() for c in l.strip().strip("|").split("|")] for l in p.read_text(encoding="utf-8").splitlines() if l.startswith("|") and not l.startswith("|---")]
A = rows(OUT / "tableS3_height_classes.md")[1:]
np_ = lambda v: v.lstrip("+")
X = {(r["product"], r["hc"]): r for r in csv.DictReader(open(OUT / "tableS3_rmse_over5.csv", encoding="utf-8"))}
HC = {"2–5": "2-5", "5–10": "5-10", "10–15": "10-15", "15–20": "15-20", "20–30": "20-30", ">30": ">30"}
for a in A:   # check: site counts agree
    for p_ in ("DSM", "DTM", "CHM"):
        r = X[(p_, HC[a[0]])]
        assert r["n_sites"] == a[4], (a[0], p_)
rr = lambda hc, p_: f"{float(X[(p_, HC[hc])]['rmse_rms']):.1f}"
o5 = lambda hc, p_: f"{float(X[(p_, HC[hc])]['over5']):.0f}"
E = r"\|error\| > 5 m (%)"
out = ["*(a) Classes defined by lidar canopy height*", "",
       f"| Lidar canopy height class (m) | DSM bias (m) | DSM RMSE (m) | DSM {E} | DTM bias (m) | DTM RMSE (m) | DTM {E} | CHM bias (m) | CHM RMSE (m) | CHM {E} | Sites |",
       "|" + "---|" * 11]
for a in A:
    out.append("| " + " | ".join([a[0], np_(a[1]), rr(a[0], "DSM"), o5(a[0], "DSM"), np_(a[2]), rr(a[0], "DTM"), o5(a[0], "DTM"), np_(a[3]), rr(a[0], "CHM"), o5(a[0], "CHM"), a[4]]) + " |")
out += ["", "*(b) Classes defined by GEDI canopy height*", "",
        "| GEDI canopy height class (m) | DSM bias (m) | DTM bias (m) | CHM bias (m) | Sites |", "|---|---|---|---|---|"]
for a in A: out.append("| " + " | ".join([a[0], np_(a[5]), np_(a[6]), np_(a[7]), a[8]]) + " |")
(OUT / "tableS3.md").write_text("\n".join(out) + "\n", encoding="utf-8"); print("\n".join(out))
