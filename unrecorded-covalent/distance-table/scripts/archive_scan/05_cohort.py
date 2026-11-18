"""Define the at-risk population and the scan cohorts."""
import csv, json, os, random, re, sys
import pandas as pd

D = "/home/claude/work/data"
REPO = "/home/claude/dl/pdb-data-audits"
meta = pd.read_pickle(f"{D}/meta_all.pkl")
comps = pd.read_pickle(f"{D}/comps_all.pkl")

# --- organic ligand test: >1 heavy atom, contains carbon, not a lone metal/ion/solvent
SOLVENT = {"HOH", "DOD", "H2O", "WAT"}
def heavy_atoms(formula):
    n = 0
    for el, cnt in re.findall(r"([A-Z][a-z]?)(\d*)", formula or ""):
        if not el or el in ("H", "D"):
            continue
        n += int(cnt) if cnt else 1
    return n
def has_c(formula):
    return bool(re.search(r"\bC(\d|\b)", formula or "")) or bool(re.match(r"^C\d*(\s|$)", (formula or "").strip()))

comps["heavy"] = comps["formula"].map(heavy_atoms)
comps["carbon"] = comps["formula"].map(lambda f: any(el == "C" for el, _ in re.findall(r"([A-Z][a-z]?)(\d*)", f or "")))
comps["organic"] = (comps["heavy"] > 1) & comps["carbon"] & ~comps["comp"].isin(SOLVENT)
print("distinct components:", comps["comp"].nunique(), "organic:", comps.loc[comps.organic, "comp"].nunique())

organic_entries = set(comps.loc[comps.organic, "entry"])
meta["has_organic_lig"] = meta["entry"].isin(organic_entries)
meta["protein"] = meta["n_prot"].fillna(0) > 0
atrisk = meta[meta.protein & meta.has_organic_lig]
print("experimental entries:", len(meta))
print("protein entries:", int(meta.protein.sum()))
print("at-risk (protein + organic non-polymer ligand):", len(atrisk))
print(atrisk["site"].value_counts().to_string())

# --- components named in the published archive-wide unrecorded tables (pdb-data-audits)
want = set()
with open(f"{REPO}/unrecorded-covalent/strict-cut-491/unrecorded_covalent_tight_cut_491_sites.csv", newline="") as fh:
    for r in csv.DictReader(fh):
        want.add(r["ligand_comp"])
with open(f"{REPO}/unrecorded-covalent/two-rulers/probe_line_unrecorded_pairs.csv", newline="") as fh:
    for r in csv.DictReader(fh):
        want.add(r["ligand"].rsplit("_", 2)[0])
with open(f"{REPO}/unrecorded-covalent/reactive-residues/removed_811_classified.csv", newline="") as fh:
    for r in csv.DictReader(fh):
        want.add(r["ligand"])
want = {w for w in want if w}
print("components named in his unrecorded/removed tables:", len(want))

cm = comps[comps["comp"].isin(want) & comps["organic"]]
per_comp = cm.groupby("comp")["entry"].nunique().sort_values(ascending=False)
print("entries containing any of them:", cm["entry"].nunique())
print(per_comp.head(25).to_string())
per_comp.to_csv(f"{D}/cohortA_components.csv", header=["entries"])

cohortA = sorted(set(cm["entry"]) & set(atrisk["entry"]))
print("cohort A size:", len(cohortA))

rest = sorted(set(atrisk["entry"]) - set(cohortA))
random.seed(20261115)
cohortB = sorted(random.sample(rest, min(20000, len(rest))))
print("cohort B size:", len(cohortB))

json.dump({"cohortA": cohortA, "cohortB": cohortB}, open(f"{D}/cohorts.json", "w"))
atrisk.to_pickle(f"{D}/atrisk.pkl")
meta.to_pickle(f"{D}/meta_all.pkl")
comps.to_pickle(f"{D}/comps_all.pkl")
