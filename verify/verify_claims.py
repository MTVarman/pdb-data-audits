#!/usr/bin/env python3
"""Re-derive every headline number in README.md from the tables in this repository.

Standard library only. No network. No arguments. Run from anywhere:

    python3 verify/verify_claims.py

Each claim is either RE-DERIVED from a shipped table, or reported as NOT
RE-DERIVABLE with the missing input named. A missing input is a stated gap,
not a crash: the script always runs every other check and always prints a
full report. Exit status is 1 only if a number that should re-derive does not.
"""

import csv
import os
import sys
from collections import Counter

ROOT = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))

STRICT = "unrecorded-covalent/strict-cut-491/unrecorded_covalent_tight_cut_491_sites.csv"
DELETED = "unrecorded-covalent/remediation-history/deleted_rows_strict.csv"
NONCYS = "unrecorded-covalent/reactive-residues/non_cys_covalent_targets.csv"
REMOVED = "unrecorded-covalent/reactive-residues/removed_811_classified.csv"

results = []  # (status, claim, expected, got, note)


def record(status, claim, expected, got, note=""):
    results.append((status, claim, expected, got, note))


def load(rel):
    """Read a CSV. Returns None if absent rather than raising."""
    path = os.path.join(ROOT, rel)
    if not os.path.exists(path):
        return None
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def missing(claim, expected, rel, why):
    record("GAP", claim, expected, "not shipped",
           "needs %s - %s" % (rel, why))


# ---------------------------------------------------------------- the four numbers

strict = load(STRICT)

# Row 1 and row 2 of the README table are scan-wide counts. The scan output
# itself is not in this repository; only the strict cut below it is.
missing("scan floor", "16,827 sites / 2,568 entries",
        "unrecorded-covalent/scan-floor/<table>.csv",
        "full scan output of 199,603 entries is not shipped")
missing("chemistry-consistent (triage class 2)", "1,705 sites / 477 entries",
        "unrecorded-covalent/chemistry-consistent/<table>.csv",
        "class-2 table before the X-ray and resolution cuts is not shipped")

if strict is None:
    missing("strict cut", "485 sites / 289 entries", STRICT, "file absent")
    missing("density-supported", "397 sites / 233 entries", STRICT, "file absent")
else:
    n_sites = len(strict)
    n_entries = len(set(r["entry"] for r in strict))
    record("PASS" if (n_sites, n_entries) == (485, 289) else "FAIL",
           "strict cut", "485 sites / 289 entries",
           "%d sites / %d entries" % (n_sites, n_entries),
           "file is still named ...491... from an earlier count")

    # The strict cut's own definition: triage class 2, X-ray, resolution <= 2.5 A.
    classes = set(r["triage_class"] for r in strict)
    methods = set(r["method"] for r in strict)
    res = [num(r["resolution_A"]) for r in strict]
    worst = max(v for v in res if v is not None)
    ok = classes == {"2"} and methods == {"X-ray"} and worst <= 2.5
    record("PASS" if ok else "FAIL",
           "strict cut matches its own definition",
           "class 2 only, X-ray only, resolution <= 2.5 A",
           "classes %s, methods %s, worst resolution %.2f A"
           % (sorted(classes), sorted(methods), worst))

    # Density support is recomputed from the raw wwPDB columns, not read from
    # the precomputed density_class column, so the column is checked too.
    def supported(r):
        v = [num(r["protein_rscc"]), num(r["ligand_rscc"]),
             num(r["protein_rsr"]), num(r["ligand_rsr"]),
             num(r["protein_occupancy"]), num(r["ligand_occupancy"])]
        if any(x is None for x in v):
            return False
        prs, lrs, prr, lrr, po, lo = v
        return prs >= 0.8 and lrs >= 0.8 and prr <= 0.3 and lrr <= 0.3 and po > 0 and lo > 0

    sup = [r for r in strict if supported(r)]
    s_sites = len(sup)
    s_entries = len(set(r["entry"] for r in sup))
    record("PASS" if (s_sites, s_entries) == (397, 233) else "FAIL",
           "density-supported", "397 sites / 233 entries",
           "%d sites / %d entries" % (s_sites, s_entries),
           "recomputed from RSCC/RSR/occupancy, thresholds as stated")

    by_col = [r for r in strict if r["density_class"] == "density_supported"]
    record("PASS" if len(by_col) == s_sites else "FAIL",
           "density_class column agrees with the rule",
           "%d" % s_sites, "%d" % len(by_col))

# ---------------------------------------------------------------- prose claims

if strict is not None:
    # "None of the 485 strict sites is reported as a steric-clash outlier."
    cl = Counter(r["clash_class"] for r in strict)
    pair_outliers = sum(n for k, n in cl.items() if k.startswith("clash_for_pair"))
    record("PASS" if pair_outliers == 0 else "FAIL",
           "no strict site is a clash outlier for its own bonded pair",
           "0", str(pair_outliers),
           "%d sites do carry a clash between the same residue and the "
           "ligand's other atoms" % cl.get(
               "clash_between_same_residue_and_ligand_other_atoms", 0))

    # "More than half of the strict 289 were released since 2020."
    first = {}
    for r in strict:
        first.setdefault(r["entry"], r["release_date"])
    since = sum(1 for d in first.values() if d >= "2020-01-01")
    record("PASS" if since > len(first) / 2 else "FAIL",
           "more than half of the strict entries released since 2020",
           "> %d of %d" % (len(first) // 2, len(first)),
           "%d of %d" % (since, len(first)))

# "135 of the 485 sites (88 entries) did have the bond recorded in an earlier
#  release; 98 of those 135 are density-supported."
deleted = load(DELETED)
if deleted is None or strict is None:
    missing("previously-recorded sites", "135 sites / 88 entries", DELETED, "file absent")
else:
    index = {r["site_id"]: r for r in strict}
    matched = set(r["matched_site_id"] for r in deleted if r["matched_site_id"].strip())
    hit = matched & set(index)
    h_entries = len(set(index[i]["entry"] for i in hit))
    record("PASS" if (len(hit), h_entries) == (135, 88) else "FAIL",
           "sites whose bond record was removed by a later revision",
           "135 sites / 88 entries",
           "%d sites / %d entries" % (len(hit), h_entries))
    hd = sum(1 for i in hit if index[i]["density_class"] == "density_supported")
    record("PASS" if hd == 98 else "FAIL",
           "of those, density-supported", "98", str(hd))

# "Across all 811 protein-ligand covalent sites removed ... 23 density-supported
#  non-cysteine covalent sites on real drug targets."
removed = load(REMOVED)
if removed is None:
    missing("archive-wide removed sites", "811 sites", REMOVED, "file absent")
else:
    record("PASS" if len(removed) == 804 else "FAIL",
           "archive-wide removed sites, located", "804 of 811",
           "%d rows" % len(removed),
           "the 811 total needs archive_removed_811_sites.csv, which is not "
           "shipped; reactive-residues/README.md states 7 atoms were not found")

noncys = load(NONCYS)
if noncys is None:
    missing("non-cysteine targets", "23 sites", NONCYS, "file absent")
else:
    record("PASS" if len(noncys) == 23 else "FAIL",
           "density-supported non-cysteine sites on drug targets",
           "23", "%d rows" % len(noncys))
    kras = [r for r in noncys if r["entry"].upper() == "8JHL"]
    record("PASS" if kras else "FAIL",
           "KRAS G12D 8JHL present among them", "1 row", "%d row(s)" % len(kras))

# ---------------------------------------------------------------- report

W = 62
print("Re-derivation of README.md headline numbers")
print("repository: %s" % ROOT)
print("=" * (W + 26))
for status, claim, expected, got, note in results:
    print("%-5s %s" % (status, claim))
    print("      expected: %s" % expected)
    print("      got:      %s" % got)
    if note:
        print("      note:     %s" % note)
print("=" * (W + 26))

counts = Counter(r[0] for r in results)
print("%d re-derived, %d mismatched, %d not re-derivable from this repository"
      % (counts["PASS"], counts["FAIL"], counts["GAP"]))

if counts["GAP"]:
    print()
    print("Not re-derivable here. Each needs the input named above:")
    for status, claim, expected, got, note in results:
        if status == "GAP":
            print("  - %s (%s)" % (claim, expected))

sys.exit(1 if counts["FAIL"] else 0)
