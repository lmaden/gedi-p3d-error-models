# tableS5_assemble.py  --  Table S5 as printed: Table S5 with a column giving one SD of each variable in its units,
# taken from Table S1 (the SD column as printed), so each effect reads in physical units without turning to Table S1.
# Inputs  (P3D_RESULTS) tableS5/tableS5_coefficients.md (tableS5_coefficients.R), tableS1/tableS1_predictors.md (tableS1_predictors.R)
# Output  (P3D_RESULTS) tableS5/tableS5.md
import os, pathlib, sys
sys.stdout.reconfigure(encoding="utf-8")
OUT_ROOT = pathlib.Path(os.environ.get("P3D_RESULTS", "results"))
def rows(p): return [[c.strip() for c in l.strip().strip("|").split("|")] for l in p.read_text(encoding="utf-8").splitlines() if l.startswith("|") and not l.startswith("|---")]
s1 = rows(OUT_ROOT / "tableS1" / "tableS1_predictors.md"); h1 = s1[0]
iv, iu, isd = h1.index("Variable"), h1.index("Units"), h1.index("SD")
U = {"degrees": "°", "m": " m", "unitless": "", "proportion": ""}
sd = {r[iv]: r[isd] + U[r[iu]] for r in s1[1:] if r[iu] in U and r[isd]}
s5 = rows(OUT_ROOT / "tableS5" / "tableS5_coefficients.md")
out = ["| " + " | ".join(s5[0][:2] + ["SD of variable"] + s5[0][2:]) + " |", "|" + "---|" * (len(s5[0]) + 1)]
n = 0
for r in s5[1:]:
    v = sd.get(r[1], "")
    n += bool(v)
    out.append("| " + " | ".join(r[:2] + [v] + r[2:]) + " |")
assert n == 17, n
(OUT_ROOT / "tableS5" / "tableS5.md").write_text("\n".join(out) + "\n", encoding="utf-8"); print("\n".join(out))
