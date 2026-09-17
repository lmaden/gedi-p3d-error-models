"""build_tableS1_chm16.py (15 Sep 2026)

Rebuild the CHM block of Supplementary Table S1 from the exact posterior draws exported on
the cluster by s3_vif_corr_probdir_16site.R.

  in :  tableS1_chm16_interactions.csv   (term, mean, median, lo95, hi95, p_positive)
  out:  tableS1_chm_16site_rows.csv      (29 non-credible rows, document column order)

Conventions: posterior MEDIAN as the estimate, matching the Table S4 note ("estimates are
posterior medians"). The superseded 19-site block printed the posterior mean (fixef default),
under which WSCI x Coastal Wetland and WSCI x Shrubland print a negative estimate beside a
"positive" probability of direction; the median removes that. U+2212 for minus, rows ordered
by interaction type then land-cover code, exact P(direction) from the draws.
Credible terms (95% CI excludes zero) are dropped here; they belong in Table S4.
"""
import csv, sys, os

MINUS = "−"
NAMES = {
    "BDF": "Broadleaf Deciduous Forest", "CWL": "Coastal Wetland",
    "DNF": "Deciduous Needleleaf Forest", "EBF": "Evergreen Broadleaf Forest",
    "ENF": "Evergreen Needleleaf Forest", "GRS": "Grassland",
    "ICP": "Irrigated Cropland", "IMP": "Impervious Surface",
    "IWL": "Inland Wetland", "MFT": "Mixed Forest",
    "RCP": "Rainfed Cropland", "SHR": "Shrubland",
    "SVG": "Sparse Vegetation", "UNK": "Unknown/Unclassified", "WTR": "Water",
}


def num(x):
    """Two decimals, explicit sign, typographic minus."""
    s = "%+.2f" % x
    return s.replace("-", MINUS)


def pct(p_positive):
    """Exact posterior probability of direction, as the table phrases it."""
    if p_positive >= 0.5:
        return "%d%% positive" % int(100 * p_positive + 0.5)
    return "%d%% negative" % int(100 * (1 - p_positive) + 0.5)


def parse(term):
    """slope_mean_z:lc_l1_codeXXX -> ('Slope', 'XXX'); slope_mean_z:wsci_z -> ('Slope', None)."""
    left, right = term.split(":", 1)
    base = {"slope_mean_z": "Slope", "wsci_z": "WSCI"}[left]
    if right.startswith("lc_l1_code"):
        return base, right[len("lc_l1_code"):]
    return base, None


def main(src, dst):
    rows, dropped = [], []
    with open(src, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            lo, hi = float(r["lo95"]), float(r["hi95"])
            base, code = parse(r["term"])
            if lo > 0 or hi < 0:                      # credible -> Table S4, not S1
                dropped.append((base, code))
                continue
            if code is None:                          # Slope x WSCI, no land-cover code
                itype, name = "%s × WSCI" % base, "–"
                code = "–"
            else:
                itype, name = "%s × Land Cover" % base, NAMES[code]
            rows.append({
                "Product": "CHM", "Interaction Type": itype, "Land Cover": name, "Code": code,
                "Est. (m)": num(float(r["median"])),
                "95% CI": "[%s, %s]" % (num(lo), num(hi)),
                "Prob. Dir.": pct(float(r["p_positive"])),
            })

    order = {"Slope × Land Cover": 0, "Slope × WSCI": 1, "WSCI × Land Cover": 2}
    rows.sort(key=lambda d: (order[d["Interaction Type"]], d["Code"]))

    with open(dst, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print("wrote %d rows to %s" % (len(rows), os.path.basename(dst)))
    print("dropped as credible (Table S4): %s" % ", ".join(
        "%s x %s" % (b, c or "WSCI") for b, c in dropped))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
