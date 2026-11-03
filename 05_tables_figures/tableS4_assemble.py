# tableS4_assemble.py  --  Table S4 as printed: a Footprints column and a subhead row before the forest-type rows,
# forest types ordered by footprints. The forest rows weight each site equally (sites with >= 30 footprints of the
# type, RMSE by the Table 2 rule); mixed forest has no row of its own and stays inside "All forest types". The
# Footprints column gives all footprints of the class.
# Inputs  (P3D_RESULTS)   tableS4/tableS4_persite.md (site rows, tableS4_persite.R),
#                         tableS4/tableS4_forest_rows.md (forest rows, tableS4_forest_rows.py),
#                         table1/table1_sites.csv (site footprints, table1_sites.R)
#         (P3D_DATA_ROOT) rework/frame_footprints.csv.gz (class footprints)
# Output  (P3D_RESULTS)   tableS4/tableS4.md
# Checks  site footprints sum to 630,390; forest-type footprints sum to all forest; all forest + non-forest = 630,390.
import gzip, csv, os, pathlib, sys
sys.stdout.reconfigure(encoding="utf-8")
DATA_ROOT = pathlib.Path(os.environ.get("P3D_DATA_ROOT", "<LOCAL_DATA_ROOT>"))
OUT_ROOT = pathlib.Path(os.environ.get("P3D_RESULTS", "results"))
site = [l for l in (OUT_ROOT / "tableS4" / "tableS4_persite.md").read_text(encoding="utf-8").splitlines() if not l.startswith("|---")]
forest = [l for l in (OUT_ROOT / "tableS4" / "tableS4_forest_rows.md").read_text(encoding="utf-8").splitlines() if l.startswith("|")]
n_site = {int(r["paper_site"]): int(r["n"]) for r in csv.DictReader(open(OUT_ROOT / "table1" / "table1_sites.csv", encoding="utf-8"))}
assert sum(n_site.values()) == 630390
cnt = {}
with gzip.open(DATA_ROOT / "rework" / "frame_footprints.csv.gz", "rt") as fh:
    for r in csv.DictReader(fh): cnt[r["lc_l1_code"]] = cnt.get(r["lc_l1_code"], 0) + 1
FOREST = {"Evergreen needleleaf forest": "ENF", "Deciduous broadleaf forest": "BDF", "Evergreen broadleaf forest": "EBF",
          "Deciduous needleleaf forest": "DNF", "Mixed forest": "MFT"}
n_f = {k: cnt[v] for k, v in FOREST.items()}
n_f["All forest types"] = sum(n_f.values()); n_f["Non-forest classes"] = 630390 - n_f["All forest types"]
assert n_f["All forest types"] + n_f["Non-forest classes"] == 630390
fmt = lambda n: f"{n:,}"
def cells(l): return [c.strip() for c in l.strip().strip("|").split("|")]
out = []
h = cells(site[0]); h[0] = "Site or forest type"; h = [x.replace("> 5 m (%)", "\\|error\\| > 5 m (%)") for x in h]; out.append("| " + " | ".join([h[0], "Footprints"] + h[1:]) + " |"); out.append("|" + "---|" * (len(h) + 1))
for l in site[1:]:
    c = cells(l); k = int(c[0].split()[0]); out.append("| " + " | ".join([c[0], fmt(n_site[k])] + c[1:]) + " |")
out.append("| *Forest types, sites weighted equally* |" + " |" * len(h))
fr = [cells(l) for l in forest]
types = sorted([c for c in fr if c[0] in FOREST], key=lambda c: -n_f[c[0]]) + [c for c in fr if c[0] not in FOREST]
assert "Mixed forest" not in [c[0] for c in types] and len(types) == 6, [c[0] for c in types]
for c in types: out.append("| " + " | ".join([c[0], fmt(n_f[c[0]])] + c[1:]) + " |")
(OUT_ROOT / "tableS4" / "tableS4.md").write_text("\n".join(out) + "\n", encoding="utf-8")
print("\n".join(out)); print({k: n_f[k] for k in n_f})
