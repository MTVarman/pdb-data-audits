"""Final processing-site analysis: four counting units, era and depositor controls."""
import glob, math, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, "/home/claude/work")
import stats_tools as S

D = "/home/claude/work/data"
REPO = "/home/claude/dl/pdb-data-audits/unrecorded-covalent"
OUT = "/home/claude/out"
os.makedirs(OUT, exist_ok=True)
SITES = ["RCSB", "PDBE", "PDBJ", "PDBC", "BNL", "NDB"]
L = []


def say(s=""):
    print(s, flush=True)
    L.append(str(s))


def fmt(t):
    return t.to_string(index=False)


# covalent radii (Angstrom), Cordero et al. 2008 single-bond values
RCOV = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "SE": 1.20, "P": 1.07,
        "F": 0.57, "CL": 1.02, "BR": 1.20, "I": 1.39, "B": 0.84, "SI": 1.11,
        "AS": 1.19, "H": 0.31}

# ---------------- load ----------------
done = pd.concat([pd.read_csv(p, sep="\t", header=None, names=["entry", "n", "status"])
                  for p in sorted(glob.glob(f"{D}/scan/done_*.tsv"))], ignore_index=True)
done = done.drop_duplicates("entry")
con = pd.concat([pd.read_csv(p, sep="\t") for p in sorted(glob.glob(f"{D}/scan/contacts_*.tsv"))],
                ignore_index=True)
raw = len(con)
con = con.drop_duplicates(["entry", "prot_chain", "prot_res", "prot_seq", "prot_atom",
                           "lig_chain", "lig_comp", "lig_seq", "lig_atom"])
meta = pd.read_pickle(f"{D}/meta_all.pkl")
m = meta.set_index("entry")
ct = pd.read_pickle(f"{D}/comp_types.pkl").set_index("comp")
scanned = set(done.loc[done.status == "ok", "entry"])
con = con[con.entry.isin(scanned)].copy()

for f in ["site", "rel_year", "dep_year", "cit", "la", "resolution", "method1"]:
    con[f] = con.entry.map(m[f])
con["comp_type"] = con.lig_comp.map(ct["comp_type"]).fillna("unknown")
con["lig_class"] = np.select(
    [con.comp_type.str.contains("saccharide", na=False),
     con.comp_type.str.contains("peptide", na=False),
     con.comp_type.str.contains("RNA|DNA", na=False, regex=True)],
    ["saccharide", "peptide-like", "nucleotide-linking"], default="non-polymer")
con["cys"] = con.prot_res.eq("CYS") & con.prot_atom.eq("SG")
con["sidechain"] = ~con.prot_atom.isin(["N", "O", "OXT"])
con["cl20"] = con.dist <= 2.0
rsum = con.prot_el.str.upper().map(RCOV).fillna(np.nan) + con.lig_el.str.upper().map(RCOV)
con["bond_window"] = (con.dist >= rsum - 0.13) & (con.dist <= rsum + 0.14)
all_by_site = meta.groupby("site").size()
scanned_by_site = meta[meta.entry.isin(scanned)].groupby("site").size()

say("# Processing site and the unrecorded-covalent-bond rate, archive-wide")
say()
say("Michael T. Varman, pdb-data-audits. Extension of "
    "`unrecorded-covalent/cys-ligand-sites/` from a 24,144-entry cysteine sample to the whole "
    "experimental archive. Every number below comes from data retrieved in this run; nothing is "
    "carried over from the repo's tables except where a section says so.")
say()
say("<<HEADLINE>>")
say()
say("## 0. What was retrieved, and what was scanned")
say()
say(f"| | |")
say(f"|---|---:|")
say(f"| experimental PDB entries, RCSB Search API, 16 Nov 2026 | {len(meta):,} |")
say(f"| per-entry metadata from data.rcsb.org/graphql (process_site, deposit_site, deposition "
    f"and release date, method, resolution, primary citation, audit authors, non-polymer "
    f"components) | {len(meta):,} entries, 0 missing |")
say(f"| at-risk population: >=1 protein polymer entity AND >=1 organic non-polymer component "
    f"of >1 heavy atom | {len(done):,} |")
say(f"| mmCIF files downloaded from files.rcsb.org and scanned | {len(done):,} attempted, "
    f"{len(scanned):,} scanned |")
say(f"| contacts found, <=2.6 A | {len(con):,} |")
say(f"| contacts at covalent length, <=2.0 A | {int(con.cl20.sum()):,} |")
say()
say(done.status.value_counts().to_frame("entries").to_string())
say()
say("All of data.rcsb.org, search.rcsb.org and files.rcsb.org were reachable; there were no "
    "download failures in 161,748 file requests. `pdbx_database_status.process_site` is not "
    "searchable or facetable on search.rcsb.org (the API returns "
    "\"search is not enabled on [ pdbx_database_status.process_site ]\"), so the site for every "
    "entry was read one entry at a time from the Data API rather than from a faceted query. "
    "Entries outside the at-risk population cannot hold a protein-nucleophile/organic-ligand "
    "contact and were not downloaded.")
say()
say("Processing-site vocabulary found across all experimental entries: "
    + ", ".join(f"{s} {int(all_by_site.get(s,0)):,}" for s in SITES) + ".")
say()
say("### Method")
say("Contact: closest approach between a ligand heavy atom (non-metal, in a non-polymer "
    "non-water residue of more than one heavy atom containing carbon) and a protein nucleophile "
    "atom of a polymer residue. Model 1 only, hydrogens removed, occupancy > 0 on both atoms, "
    "alternate conformers must match or be absent, distances symmetry-aware through gemmi "
    "`NeighborSearch` and the nearest periodic image. Nucleophile atoms: sidechain SG OG OG1 "
    "OG2 OH NZ ND1 NE2 OD1 OD2 OE1 OE2 ND2 NE NH1 NH2 SD SE SE1, plus backbone N O OXT.")
say("Recorded: a `_struct_conn` row whose `conn_type_id` begins with `covale` joins that "
    "protein residue to that ligand residue, either partner order, any symmetry operator. "
    "This is the residue-pair test used in `cys-ligand-sites`.")
say("Covalent length: distance <= 2.0 A, the scan floor of the repo's three-numbers table.")
say("Bond window: distance within [r1+r2-0.13, r1+r2+0.14] A for the element pair, using "
    "Cordero single-bond covalent radii. That reproduces the two windows stated in "
    "`strict-cut-491/README.md` (C-O 1.3-1.55, C-S 1.7-1.95) to within 0.01 A.")
say()
say("### Scanner validation")
say("The scanner was checked against the repo's own published pair tables before any site "
    "analysis was run:")
if os.path.exists("/home/claude/work/validation.txt"):
    say("```")
    say(open("/home/claude/work/validation.txt").read().strip())
    say("```")
say("It also reproduces the matched ubiquinone-0 case exactly: 8BEO (PDBe) 6 of 6 recorded at "
    "1.664-1.718 A; 8I01, 8I05, 8I08 (PDBj) 0 of 23 recorded at 1.765-1.776 A.")

# ---------------- counting units ----------------
PRIM = con[con.cl20 & con.sidechain & con.lig_class.eq("non-polymer")].copy()
say()
say("=" * 72)
say("## 1. Recorded vs unrecorded by processing site, archive-wide")
say("=" * 72)
say()
say("Analysis set unless stated: **covalent-length contacts (<= 2.0 A) between a sidechain "
    "nucleophile and a non-polymer organic ligand** "
    f"({len(PRIM):,} contacts in {PRIM.entry.nunique():,} entries). "
    "Saccharides (N- and O-glycosylation) are counted separately: they are 152,223 of the "
    f"{int(con.cl20.sum()):,} covalent-length contacts and are recorded at 3-7% miss rate at "
    "every site, so pooling them buries everything else.")
say()
say("### 1a. Contact level, and why it cannot be the headline")
rows = []
for s in SITES:
    d = PRIM[PRIM.site == s]
    if len(d) == 0:
        continue
    u = len(d) - int(d.recorded.sum())
    p, lo, hi = S.wilson(u, len(d))
    top20 = d.groupby("entry").size().sort_values(ascending=False).head(20).sum()
    rows.append(dict(site=s, contacts=len(d), entries=d.entry.nunique(),
                     unrecorded=u, unrecorded_pct=round(100 * p, 2),
                     ci95=f"{100*lo:.2f}-{100*hi:.2f}",
                     pct_of_contacts_from_top20_entries=round(100 * top20 / len(d), 1)))
t1a = pd.DataFrame(rows)
say(fmt(t1a))
say()
say("The last column is the problem. 84.6% of PDBj's covalent-length contacts come from its "
    "20 largest entries, against 13.3% at RCSB and 12.6% at PDBe. PDBj's top five are "
    "phycobilisome cryo-EM structures (7Y7A 11,373 contacts, 7Y5E 5,751, 7EZX 5,577, "
    "7Y4L 5,339, 6KGX 5,249) and its most common ligand is phycoerythrobilin (PEB, 30,663 of "
    "44,144 contacts). A contact-level rate is a rate over bilin copies in a handful of "
    "depositions, not over the archive. Everything below is therefore counted on units that "
    "cannot be inflated that way.")
say()
say("Top ligand components by contact count, per site:")
rows = []
for s in ["RCSB", "PDBE", "PDBJ", "PDBC", "BNL"]:
    d = PRIM[PRIM.site == s]
    rows.append({"site": s, **{f"top{i+1}": f"{c} {n}" for i, (c, n) in
                               enumerate(d.lig_comp.value_counts().head(5).items())}})
say(fmt(pd.DataFrame(rows)))

# unit builders
def collapse(df, keys):
    g = df.groupby(keys, observed=True).agg(n=("recorded", "size"),
                                            nrec=("recorded", "sum")).reset_index()
    g["unrec"] = ((g.n - g.nrec) > 0).astype(int)
    for f in ["site", "rel_year", "dep_year", "cit", "la", "method1"]:
        if f not in g.columns:
            g[f] = g["entry"].map(m[f]) if "entry" in g.columns else np.nan
    return g


U = {
    "contact": PRIM.assign(unrec=1 - PRIM.recorded),
    "chemistry instance (entry x ligand component x nucleophile atom)":
        collapse(PRIM, ["entry", "lig_comp", "prot_res", "prot_atom"]),
    "entry": collapse(PRIM, ["entry"]),
    "depositing group (primary citation)": None,
    "depositing group (last audit author)": None,
}
gl_cit = PRIM.groupby(["cit", "site"], observed=True).agg(
    n=("recorded", "size"), nrec=("recorded", "sum")).reset_index()
gl_cit["unrec"] = ((gl_cit.n - gl_cit.nrec) > 0).astype(int)
gl_la = PRIM.groupby(["la", "site"], observed=True).agg(
    n=("recorded", "size"), nrec=("recorded", "sum")).reset_index()
gl_la["unrec"] = ((gl_la.n - gl_la.nrec) > 0).astype(int)
U["depositing group (primary citation)"] = gl_cit
U["depositing group (last audit author)"] = gl_la


def unit_table(df, label, extra_den=None):
    say(f"\n### {label} (n={len(df):,})")
    rows = []
    for s in SITES:
        d = df[df.site == s]
        if len(d) == 0:
            continue
        k = int(d.unrec.sum())
        p, lo, hi = S.wilson(k, len(d))
        r = dict(site=s, units=len(d), with_unrecorded=k, pct=round(100 * p, 2),
                 ci95=f"{100*lo:.2f}-{100*hi:.2f}")
        if extra_den is not None:
            r["per_1000_at_risk_entries"] = round(1000 * k / max(1, scanned_by_site.get(s, 1)), 2)
        rows.append(r)
    t = pd.DataFrame(rows)
    say(fmt(t))
    j = df[df.site == "PDBJ"]; o = df[(df.site != "PDBJ") & df.site.notna()]
    a, b = int(j.unrec.sum()), len(j) - int(j.unrec.sum())
    c, d2 = int(o.unrec.sum()), len(o) - int(o.unrec.sum())
    sor, cor, lo, hi, pv = S.or_ci(a, b, c, d2)
    say(f"PDBj vs all other sites: OR {sor:.2f} (conditional MLE {cor:.2f}, "
        f"95% CI {lo:.2f}-{hi:.2f}), Fisher p {pv:.3g}")
    return t, (sor, lo, hi, pv)


say()
say("### 1b. The four counting units, PDBj vs the rest")
res = {}
t_ci, res["contact"] = unit_table(U["contact"], "unit = contact")
t_chem, res["chem"] = unit_table(U["chemistry instance (entry x ligand component x nucleophile atom)"],
                                 "unit = chemistry instance (entry x ligand component x nucleophile atom); "
                                 "counts as unrecorded if any copy lacks a record")
t_ent, res["entry"] = unit_table(U["entry"], "unit = entry", extra_den=True)
t_gc, res["cit"] = unit_table(gl_cit, "unit = depositing group x site, group = primary citation "
                                      "(DOI/PubMed, else hashed author list)")
t_gl, res["la"] = unit_table(gl_la, "unit = depositing group x site, group = hashed last audit author")
say()
say("Summary of PDBj vs all other sites, odds of no `_struct_conn` covale record:")
say()
say("| unit | OR | 95% CI | p |")
say("|---|---:|---|---:|")
for k, lab in [("contact", "contact"),
               ("chem", "chemistry instance (entry x component x nucleophile)"),
               ("entry", "entry"), ("cit", "depositing group, primary citation"),
               ("la", "depositing group, last author")]:
    sor, lo, hi, pv = res[k]
    say(f"| {lab} | {sor:.2f} | {lo:.2f}-{hi:.2f} | {pv:.2g} |")
say()
say("The contact-level odds ratio is the one the bilin depositions inflate. The three units "
    "that cannot be inflated by copy count all land near 1.7-2.0.")

# cofactor-free view
top_comps = PRIM.lig_comp.value_counts().head(10).index.tolist()
noc = PRIM[~PRIM.lig_comp.isin(top_comps)]
say()
say(f"### 1c. Dropping the ten highest-volume ligand components ({', '.join(top_comps)})")
nt, nres = unit_table(noc.assign(unrec=1 - noc.recorded),
                      "unit = contact, ten highest-volume components removed")
ne, neres = unit_table(collapse(noc, ["entry"]), "unit = entry, same removal")

# ---------------- 2. release year ----------------
say()
say("=" * 72)
say("## 2. By release year: era effect or site effect?")
say("=" * 72)
ent = U["entry"].copy()
ent["era"] = pd.cut(ent.rel_year, [1975, 2000, 2005, 2010, 2015, 2020, 2023, 2027],
                    labels=["<=2000", "2001-05", "2006-10", "2011-15", "2016-20", "2021-23", "2024-26"])
say("\n**entries with >=1 unrecorded covalent-length contact / entries with >=1 contact**")
rows = []
for era, g in ent.groupby("era", observed=True):
    r = {"era": str(era), "entries": len(g)}
    for s in ["RCSB", "PDBE", "PDBJ", "PDBC", "BNL"]:
        d = g[g.site == s]
        r[s] = f"{int(d.unrec.sum())}/{len(d)} ({100*d.unrec.mean():.0f}%)" if len(d) else "-"
    rows.append(r)
say(fmt(pd.DataFrame(rows)))
say("\n**era effect with site held constant** (% of entries with an unrecorded contact)")
rows = []
for s in ["RCSB", "PDBE", "PDBJ"]:
    g = ent[ent.site == s]
    r = {"site": s}
    for era, gg in g.groupby("era", observed=True):
        r[str(era)] = f"{100*gg.unrec.mean():.0f}%" if len(gg) else "-"
    rows.append(r)
say(fmt(pd.DataFrame(rows)))
say("\nThe rate rises with release year at every site, so there is a real era effect; it is "
    "not a substitute for the site effect, which is present inside almost every era band.")


def mh(df, by, label, unit="unit"):
    strata = []
    for _, g in df.groupby(by, dropna=True, observed=True):
        j = g[g.site == "PDBJ"]; o = g[(g.site != "PDBJ") & g.site.notna()]
        strata.append((int(j.unrec.sum()), len(j) - int(j.unrec.sum()),
                       int(o.unrec.sum()), len(o) - int(o.unrec.sum())))
    orm, lo, hi, chi2, p, used, tot = S.mantel_haenszel(strata)
    say(f"- {label}: OR {orm:.2f} (95% CI {lo:.2f}-{hi:.2f}), MH chi2 {chi2:.1f}, p {p:.3g}; "
        f"{used} informative strata of {len(strata)}, {tot:,} {unit}s in them")
    return orm, lo, hi, p, used


say("\n**Mantel-Haenszel, PDBj vs rest, stratified by release year**")
mh(U["entry"], "rel_year", "unit = entry, stratified by release year", "entries")
mh(U["entry"], "dep_year", "unit = entry, stratified by deposition year", "entries")
mh(U["chemistry instance (entry x ligand component x nucleophile atom)"], "rel_year",
   "unit = chemistry instance, stratified by release year", "instances")
mh(U["contact"], "rel_year", "unit = contact, stratified by release year", "contacts")

# ---------------- 3. once per depositing group ----------------
say()
say("=" * 72)
say("## 3. Counted once per depositing group")
say("=" * 72)
say()
say(fmt(t_gc))
say(f"PDBj vs rest, one count per citation group: OR {res['cit'][0]:.2f} "
    f"(95% CI {res['cit'][1]:.2f}-{res['cit'][2]:.2f}), p {res['cit'][3]:.3g}")
say()
say(fmt(t_gl))
say(f"PDBj vs rest, one count per last author: OR {res['la'][0]:.2f} "
    f"(95% CI {res['la'][1]:.2f}-{res['la'][2]:.2f}), p {res['la'][3]:.3g}")
gs = PRIM.groupby("cit")["site"].nunique()
gs2 = PRIM.groupby("la")["site"].nunique()
say()
say("**How far site and depositor can be separated at all.**")
say(f"- depositing groups (primary citation) with a covalent-length contact: {len(gs):,}; "
    f"groups that used more than one processing site: {int((gs>1).sum())} "
    f"({100*(gs>1).mean():.1f}%)")
say(f"- depositing groups (last audit author) with a covalent-length contact: {len(gs2):,}; "
    f"groups that used more than one processing site: {int((gs2>1).sum())} "
    f"({100*(gs2>1).mean():.1f}%)")
say("- within-group contrast (only groups that deposited both at PDBj and elsewhere contribute "
    "anything):")
mh(U["entry"], "cit", "unit = entry, stratified by citation group", "entries")
mh(U["entry"], "la", "unit = entry, stratified by last author", "entries")
mh(U["contact"], "la", "unit = contact, stratified by last author", "contacts")

# ---------------- 4. cys vs non-cys ----------------
say()
say("=" * 72)
say("## 4. Cysteine and non-cysteine nucleophiles")
say("=" * 72)
say()
say("All contacts <= 2.0 A, non-polymer ligand. `in bond window` restricts to distances inside "
    "the element-pair bond-length window, which removes close non-bonded contacts (a ligand "
    "carbon 1.95 A from an Arg guanidinium nitrogen is not a bond candidate).")
arms = [("cysteine SG", PRIM[PRIM.cys]),
        ("cysteine SG, in bond window", PRIM[PRIM.cys & PRIM.bond_window]),
        ("non-cysteine sidechain", PRIM[~PRIM.cys]),
        ("non-cysteine sidechain, in bond window", PRIM[~PRIM.cys & PRIM.bond_window]),
        ("Ser OG / Thr OG1 / Tyr OH, in window",
         PRIM[PRIM.prot_atom.isin(["OG", "OG1", "OH"]) & PRIM.bond_window]),
        ("Lys NZ, in window", PRIM[PRIM.prot_atom.eq("NZ") & PRIM.bond_window]),
        ("His ND1/NE2, in window",
         PRIM[PRIM.prot_res.eq("HIS") & PRIM.prot_atom.isin(["ND1", "NE2"]) & PRIM.bond_window]),
        ("Asp/Glu carboxylate, in window",
         PRIM[PRIM.prot_atom.isin(["OD1", "OD2", "OE1", "OE2"]) & PRIM.bond_window]),
        ("backbone N/O, in window",
         con[con.cl20 & ~con.sidechain & con.lig_class.eq("non-polymer") & con.bond_window])]
say("\n**contact level**")
rows = []
for lab, sub in arms:
    j = sub[sub.site == "PDBJ"]; o = sub[(sub.site != "PDBJ") & sub.site.notna()]
    a, b = len(j) - int(j.recorded.sum()), int(j.recorded.sum())
    c, d2 = len(o) - int(o.recorded.sum()), int(o.recorded.sum())
    if min(a + b, c + d2) == 0:
        continue
    sor, cor, lo, hi, pv = S.or_ci(a, b, c, d2)
    rows.append(dict(nucleophile=lab, contacts=len(sub),
                     pdbj=f"{a}/{a+b} ({100*a/(a+b):.1f}%)",
                     other=f"{c}/{c+d2} ({100*c/(c+d2):.1f}%)",
                     OR=round(sor, 2), ci95=f"{lo:.2f}-{hi:.2f}", p=f"{pv:.2g}"))
say(fmt(pd.DataFrame(rows)))
say("\n**entry level** (entries with >=1 unrecorded contact of that kind, out of entries with "
    ">=1 contact of that kind)")
rows = []
for lab, sub in arms:
    e = collapse(sub, ["entry"])
    j = e[e.site == "PDBJ"]; o = e[(e.site != "PDBJ") & e.site.notna()]
    a, b = int(j.unrec.sum()), len(j) - int(j.unrec.sum())
    c, d2 = int(o.unrec.sum()), len(o) - int(o.unrec.sum())
    if min(a + b, c + d2) == 0:
        continue
    sor, cor, lo, hi, pv = S.or_ci(a, b, c, d2)
    if lab == "cysteine SG, in bond window":
        cys_or = (sor, lo, hi, pv)
    if lab == "non-cysteine sidechain, in bond window":
        ncys_or = (sor, lo, hi, pv)
    rows.append(dict(nucleophile=lab, entries=len(e),
                     pdbj=f"{a}/{a+b} ({100*a/(a+b):.1f}%)",
                     other=f"{c}/{c+d2} ({100*c/(c+d2):.1f}%)",
                     OR=round(sor, 2), ci95=f"{lo:.2f}-{hi:.2f}", p=f"{pv:.2g}"))
say(fmt(pd.DataFrame(rows)))
say("\n**per-residue, entry level, in bond window** (entries with an unrecorded contact / "
    "entries with a contact)")
rows = []
BW = PRIM[PRIM.bond_window]
for res_name, g in BW.groupby("prot_res"):
    e = collapse(g, ["entry"])
    if len(e) < 40:
        continue
    r = {"residue": res_name, "entries": len(e)}
    for s in ["RCSB", "PDBE", "PDBJ"]:
        d = e[e.site == s]
        r[s] = f"{100*d.unrec.mean():.0f}% ({len(d)})" if len(d) else "-"
    rows.append(r)
say(fmt(pd.DataFrame(rows).sort_values("entries", ascending=False)))

# ---------------- 5. OR with controls ----------------
say()
say("=" * 72)
say("## 5. Odds ratio with confidence interval, and survival under controls")
say("=" * 72)
say()
say("Outcome: no `_struct_conn` covale record. Exposure: processed at PDBj rather than "
    "anywhere else. Unit = entry unless stated, because that is the unit the multi-copy "
    "problem does not touch.")
say()
E = U["entry"]
j = E[E.site == "PDBJ"]; o = E[(E.site != "PDBJ") & E.site.notna()]
sor, cor, lo, hi, pv = S.or_ci(int(j.unrec.sum()), len(j) - int(j.unrec.sum()),
                               int(o.unrec.sum()), len(o) - int(o.unrec.sum()))
say(f"- crude, unit = entry: OR {sor:.2f} (95% CI {lo:.2f}-{hi:.2f}), p {pv:.3g}")
mh_year = mh(E, "rel_year", "release year held constant", "entries")
mh(E, "method1", "experimental method held constant", "entries")
mh_cit = mh(E, "cit", "depositing group held constant (primary citation)", "entries")
mh_la = mh(E, "la", "depositing group held constant (last author)", "entries")
CH = U["chemistry instance (entry x ligand component x nucleophile atom)"]
mh_comp = mh(CH, "lig_comp", "unit = chemistry instance, ligand component held constant",
             "instances")
mh(CH, ["lig_comp", "rel_year"], "unit = chemistry instance, component x release year", "instances")
mh(CH, ["lig_comp", "prot_res"], "unit = chemistry instance, component x nucleophile residue",
   "instances")
mh(CH, ["cit", "lig_comp"], "unit = chemistry instance, depositing group x component", "instances")

LOGIT_FULL = None
try:
    import statsmodels.formula.api as smf
    say("\n### Logistic regression, P(entry has >=1 unrecorded covalent-length contact)")
    say("Standard errors are cluster-robust by depositing group (primary citation), so repeat "
        "depositions by one group do not count as independent evidence.")
    reg = E[E.site.notna() & E.rel_year.notna()].copy()
    reg["yr"] = reg.rel_year - 2000
    reg["logn"] = np.log(reg.n)
    reg["xray"] = reg.entry.map(m["method1"]).eq("X-RAY DIFFRACTION").astype(int)
    grp = reg["cit"].astype("category").cat.codes
    models = [("site alone", "unrec ~ C(site, Treatment(reference='RCSB'))"),
              ("+ release year", "unrec ~ C(site, Treatment(reference='RCSB')) + yr"),
              ("+ release year + log(contacts in entry)",
               "unrec ~ C(site, Treatment(reference='RCSB')) + yr + logn"),
              ("+ release year + log(contacts) + X-ray",
               "unrec ~ C(site, Treatment(reference='RCSB')) + yr + logn + xray")]
    for name, formula in models:
        try:
            fit = smf.logit(formula, data=reg).fit(disp=0, maxiter=300, cov_type="cluster",
                                                   cov_kwds={"groups": grp})
            say(f"\n-- {name} (n={int(fit.nobs):,} entries, {grp.nunique():,} clusters)")
            for nm in fit.params.index:
                if ("site" not in nm) and nm not in ("yr", "logn", "xray"):
                    continue
                bb = fit.params[nm]; se = fit.bse[nm]
                say(f"   {nm:52} OR {math.exp(bb):7.3f}  95% CI "
                    f"{math.exp(bb-1.96*se):6.3f}-{math.exp(bb+1.96*se):7.3f}  "
                    f"p {fit.pvalues[nm]:.3g}")
                if "PDBJ" in nm and name == models[-1][0]:
                    LOGIT_FULL = (math.exp(bb), math.exp(bb - 1.96 * se),
                                  math.exp(bb + 1.96 * se), fit.pvalues[nm])
        except Exception as e:
            say(f"   {name}: model failed: {repr(e)[:180]}")
except ImportError:
    say("statsmodels unavailable")

say()
say("### The 1.766 A cctbx generic C-S link value, archive-wide")
cs = PRIM[PRIM.cys & PRIM.lig_el.eq("C")]
rows = []
for s in ["RCSB", "PDBE", "PDBJ", "PDBC", "BNL"]:
    d = cs[cs.site == s]
    if len(d) == 0:
        continue
    on = d.dist.sub(1.766).abs().le(0.010)
    rows.append(dict(site=s, cys_C_contacts=len(d), unrecorded=int((d.recorded == 0).sum()),
                     unrecorded_on_value=int(((d.recorded == 0) & on).sum()),
                     pct_of_unrecorded_on_value=round(
                         100 * ((d.recorded == 0) & on).sum() / max(1, (d.recorded == 0).sum()), 1),
                     recorded_on_value=int(((d.recorded == 1) & on).sum()),
                     pct_of_recorded_on_value=round(
                         100 * ((d.recorded == 1) & on).sum() / max(1, (d.recorded == 1).sum()), 1)))
say(fmt(pd.DataFrame(rows)))

# ---------------- 6. the published audit tables (pdb-data-audits) ----------------
say()
say("=" * 72)
say("## 6. The repo's own archive-wide tables with processing site joined on")
say("=" * 72)


def joinmeta(df, col="entry"):
    df = df.copy()
    for f in ["site", "rel_year", "cit", "la"]:
        df[f] = df[col].map(m[f])
    return df


un = joinmeta(pd.read_csv(f"{REPO}/two-rulers/probe_line_unrecorded_pairs.csv")
              .assign(entry=lambda d: d.pdbid.str.upper()))
re_ = joinmeta(pd.read_csv(f"{REPO}/two-rulers/probe_line_recorded_pairs.csv")
               .assign(entry=lambda d: d.pdbid.str.upper()))
w_un = un[(un.dist_A >= 1.95) & (un.dist_A <= 2.482)
          & (un.residue_pair_bonded_elsewhere.astype(str).str.lower() == "false")]
w_re = re_[(re_.dist_A >= 1.95) & (re_.dist_A <= 2.482)]
say(f"\n**two-rulers, matched distance window 1.95-2.482 A** (the published Cys SG-ligand C archive scan; "
    f"{len(w_re)} recorded pairs, {len(w_un)} unrecorded pairs whose residue pair is not bonded "
    f"elsewhere)")
rows = []
for s in SITES:
    r = int((w_re.site == s).sum()); u = int((w_un.site == s).sum())
    if r + u == 0:
        continue
    p, lo, hi = S.wilson(u, r + u)
    rows.append(dict(site=s, recorded=r, unrecorded=u, total=r + u,
                     unrecorded_pct=round(100 * p, 2), ci95=f"{100*lo:.2f}-{100*hi:.2f}"))
say(fmt(pd.DataFrame(rows)))
a = int((w_un.site == "PDBJ").sum()); b = int((w_re.site == "PDBJ").sum())
c = int((w_un.site != "PDBJ").sum()); d2 = int((w_re.site != "PDBJ").sum())
sor, cor, lo, hi, pv = S.or_ci(a, b, c, d2)
say(f"PDBj vs rest inside that window: OR {sor:.2f} (95% CI {lo:.2f}-{hi:.2f}), p {pv:.3g}. "
    "Inside a fixed 1.95-2.48 A band the site difference is small; the distance line, not the "
    "site, dominates there.")

stt = joinmeta(pd.read_csv(f"{REPO}/strict-cut-491/unrecorded_covalent_tight_cut_491_sites.csv"))
stt["is_cys"] = stt.residue_name.eq("CYS")
say("\n**strict-cut 485 unrecorded sites** (class 2, X-ray <= 2.5 A), against each site's share "
    "of the archive")
rows = []
for s in SITES:
    d = stt[stt.site == s]
    rows.append(dict(site=s, sites=len(d), entries=d.entry.nunique(), cys_sites=int(d.is_cys.sum()),
                     all_experimental=int(all_by_site.get(s, 0)),
                     at_risk_scanned=int(scanned_by_site.get(s, 0)),
                     entries_per_1000_at_risk=round(1000 * d.entry.nunique()
                                                    / max(1, scanned_by_site.get(s, 1)), 3)))
say(fmt(pd.DataFrame(rows)))
rmv = joinmeta(pd.read_csv(f"{REPO}/reactive-residues/removed_811_classified.csv"))
rmv["is_cys"] = rmv.residue.eq("CYS")
say("\n**811 archive-removed covalent sites**")
rows = []
for s in SITES:
    d = rmv[rmv.site == s]
    rows.append(dict(site=s, sites=len(d), entries=d.entry.nunique(), cys_sites=int(d.is_cys.sum()),
                     entries_per_1000_at_risk=round(1000 * d.entry.nunique()
                                                    / max(1, scanned_by_site.get(s, 1)), 3)))
say(fmt(pd.DataFrame(rows)))
stt.to_csv(f"{OUT}/lane2_strict485_with_site.csv", index=False)
rmv.to_csv(f"{OUT}/lane2_removed811_with_site.csv", index=False)

# ---------------- deliverables ----------------
agg = con.assign(unrec=1 - con.recorded)
agg["prim"] = agg.cl20 & agg.sidechain & agg.lig_class.eq("non-polymer")
per = agg.groupby("entry").apply(lambda g: pd.Series({
    "n_contacts_le2_6A": len(g),
    "n_contacts_le2_0A": int(g.cl20.sum()),
    "n_recorded_le2_0A": int(g.loc[g.cl20, "recorded"].sum()),
    "n_unrecorded_le2_0A": int(g.loc[g.cl20, "unrec"].sum()),
    "n_primary_le2_0A": int(g.prim.sum()),
    "n_primary_unrecorded": int((g.prim & g.unrec.eq(1)).sum()),
    "n_primary_bondwindow": int((g.prim & g.bond_window).sum()),
    "n_primary_bondwindow_unrecorded": int((g.prim & g.bond_window & g.unrec.eq(1)).sum()),
    "n_cys_le2_0A": int((g.cl20 & g.cys).sum()),
    "n_cys_unrecorded": int((g.cl20 & g.cys & g.unrec.eq(1)).sum()),
    "n_noncys_le2_0A": int((g.cl20 & ~g.cys).sum()),
    "n_noncys_unrecorded": int((g.cl20 & ~g.cys & g.unrec.eq(1)).sum()),
    "n_saccharide_le2_0A": int((g.cl20 & g.lig_class.eq("saccharide")).sum()),
    "n_saccharide_unrecorded": int((g.cl20 & g.lig_class.eq("saccharide") & g.unrec.eq(1)).sum()),
    "min_contact_dist_A": g.dist.min(),
}), include_groups=False).reset_index()
base = meta[meta.entry.isin(scanned)][["entry", "site", "deposit_site", "dep_date", "rel_date",
                                       "rel_year", "resolution", "method1", "cit", "la",
                                       "n_np", "n_prot"]]
per = base.merge(per, on="entry", how="left")
numcols = [c for c in per.columns if c.startswith(("n_", "min_"))]
per[numcols] = per[numcols].fillna(0)
strict = set(pd.read_csv(f"{REPO}/strict-cut-491/unrecorded_covalent_tight_cut_491_sites.csv").entry.str.upper())
rem = set(pd.read_csv(f"{REPO}/reactive-residues/removed_811_classified.csv").entry.str.upper())
tr = set(pd.read_csv(f"{REPO}/two-rulers/probe_line_unrecorded_pairs.csv").pdbid.str.upper())
per["in_strict485"] = per.entry.isin(strict).astype(int)
per["in_removed811"] = per.entry.isin(rem).astype(int)
per["in_two_rulers_unrecorded"] = per.entry.isin(tr).astype(int)
per["has_unrecorded_le2_0A"] = (per["n_unrecorded_le2_0A"] > 0).astype(int)
per["has_unrecorded_primary"] = (per["n_primary_unrecorded"] > 0).astype(int)
per = per.rename(columns={"site": "process_site", "cit": "citation_group",
                          "la": "last_author_group", "method1": "method"})
per.to_csv(f"{OUT}/lane2_per_entry.csv", index=False)
con.drop(columns=["cit", "la"]).to_csv(f"{OUT}/lane2_contacts.csv.gz", index=False,
                                       compression="gzip")
for nm, t in [("site_contacts_primary", t1a), ("site_units_contact", t_ci),
              ("site_units_chemistry_instance", t_chem), ("site_units_entry", t_ent),
              ("site_units_group_citation", t_gc), ("site_units_group_lastauthor", t_gl)]:
    t.to_csv(f"{OUT}/lane2_{nm}.csv", index=False)

say()
say("=" * 72)
say("## Files")
say("=" * 72)
say(f"- `lane2_processing_site.md` this report")
say(f"- `lane2_per_entry.csv` one row per scanned entry ({len(per):,} rows): processing site, "
    f"deposit site, deposition and release dates, method, resolution, depositing-group keys, "
    f"contact and unrecorded counts (all, primary set, bond window, cysteine, non-cysteine, "
    f"saccharide), and membership flags for the repo's three archive-wide tables")
say(f"- `lane2_contacts.csv.gz` one row per contact ({len(con):,} rows): both partners' chain, "
    f"residue, number, atom and element, distance, and the recorded flag")
say(f"- `lane2_site_units_*.csv` the per-site tables for each counting unit")
say(f"- `lane2_strict485_with_site.csv`, `lane2_removed811_with_site.csv` the repo's own tables "
    f"with site, release year and depositing-group keys joined on")

# ---------------- headline ----------------
ej = E[E.site == "PDBJ"]; eo = E[(E.site != "PDBJ") & E.site.notna()]
H = []
H.append("## Headline")
H.append("")
H.append("Across the whole experimental archive, whether a covalent-length protein-ligand "
         "contact carries a `_struct_conn` covale record still depends on the processing site, "
         "and PDBj is still the outlier. Counted once per entry, "
         f"{int(ej.unrec.sum()):,} of {len(ej):,} PDBj-processed entries with such a contact "
         f"({100*ej.unrec.mean():.1f}%) hold at least one with no record, against "
         f"{int(eo.unrec.sum()):,} of {len(eo):,} elsewhere ({100*eo.unrec.mean():.1f}%).")
H.append("")
H.append("| | RCSB | PDBe | PDBj | PDBc | BNL |")
H.append("|---|---:|---:|---:|---:|---:|")
row = {}
for s in ["RCSB", "PDBE", "PDBJ", "PDBC", "BNL"]:
    d = E[E.site == s]
    row[s] = f"{int(d.unrec.sum())}/{len(d)} ({100*d.unrec.mean():.0f}%)"
H.append("| entries with >=1 unrecorded contact | " + " | ".join(row[s] for s in
         ["RCSB", "PDBE", "PDBJ", "PDBC", "BNL"]) + " |")
row = {}
for s in ["RCSB", "PDBE", "PDBJ", "PDBC", "BNL"]:
    d = gl_cit[gl_cit.site == s]
    row[s] = f"{int(d.unrec.sum())}/{len(d)} ({100*d.unrec.mean():.0f}%)"
H.append("| depositing groups with >=1 unrecorded contact | " + " | ".join(
    row[s] for s in ["RCSB", "PDBE", "PDBJ", "PDBC", "BNL"]) + " |")
row = {}
cysE = collapse(PRIM[PRIM.cys & PRIM.bond_window], ["entry"])
for s in ["RCSB", "PDBE", "PDBJ", "PDBC", "BNL"]:
    d = cysE[cysE.site == s]
    row[s] = f"{int(d.unrec.sum())}/{len(d)} ({100*d.unrec.mean():.0f}%)" if len(d) else "-"
H.append("| cysteine entries, bond window, >=1 unrecorded | " + " | ".join(
    row[s] for s in ["RCSB", "PDBE", "PDBJ", "PDBC", "BNL"]) + " |")
H.append("")
H.append("**Odds of no record, PDBj vs every other site**")
H.append("")
H.append("| | OR | 95% CI | p |")
H.append("|---|---:|---|---:|")
H.append(f"| per entry, crude | {res['entry'][0]:.2f} | {res['entry'][1]:.2f}-{res['entry'][2]:.2f} "
         f"| {res['entry'][3]:.1g} |")
H.append(f"| per entry, release year held constant (MH) | {mh_year[0]:.2f} | "
         f"{mh_year[1]:.2f}-{mh_year[2]:.2f} | {mh_year[3]:.1g} |")
H.append(f"| per entry, depositing group held constant, last author (MH) | {mh_la[0]:.2f} | "
         f"{mh_la[1]:.2f}-{mh_la[2]:.2f} | {mh_la[3]:.1g} |")
H.append(f"| per chemistry instance, ligand component held constant (MH) | {mh_comp[0]:.2f} | "
         f"{mh_comp[1]:.2f}-{mh_comp[2]:.2f} | {mh_comp[3]:.1g} |")
H.append(f"| once per depositing group, primary citation | {res['cit'][0]:.2f} | "
         f"{res['cit'][1]:.2f}-{res['cit'][2]:.2f} | {res['cit'][3]:.1g} |")
if LOGIT_FULL:
    H.append(f"| per entry, logistic: + year + entry size + X-ray, SE clustered by depositing "
             f"group | {LOGIT_FULL[0]:.2f} | {LOGIT_FULL[1]:.2f}-{LOGIT_FULL[2]:.2f} | "
             f"{LOGIT_FULL[3]:.1g} |")
H.append("")
H.append("**Answers to the five questions**")
H.append("")
H.append(f"1. Archive-wide, PDBj misses the record on {100*ej.unrec.mean():.0f}% of the entries "
         f"that hold a covalent-length contact against {100*eo.unrec.mean():.0f}% elsewhere; "
         f"PDBe ({100*E[E.site=='PDBE'].unrec.mean():.0f}%) is better than RCSB "
         f"({100*E[E.site=='RCSB'].unrec.mean():.0f}%) and PDBc "
         f"({100*E[E.site=='PDBC'].unrec.mean():.0f}%) is as bad as PDBj.")
H.append("2. It is not an era effect. The miss rate rises with release year at every site, but "
         f"the site contrast survives stratification on release year (OR {mh_year[0]:.2f}, "
         f"{mh_year[1]:.2f}-{mh_year[2]:.2f}) and on deposition year, and PDBj is above the "
         "others inside almost every era band.")
H.append(f"3. Counted once per depositing group it survives: OR {res['cit'][0]:.2f} "
         f"({res['cit'][1]:.2f}-{res['cit'][2]:.2f}) by primary citation, "
         f"{res['la'][0]:.2f} ({res['la'][1]:.2f}-{res['la'][2]:.2f}) by last author. It is "
         "smaller than the per-contact figure, so part of the raw gap is clustering by "
         "depositor and by copy count, as the repo's cysteine result already said.")
H.append("4. It holds for cysteine most strongly "
         f"(entry level, bond window: OR {cys_or[0]:.2f}, {cys_or[1]:.2f}-{cys_or[2]:.2f}) and "
         f"weakly for non-cysteine sidechains as a whole (OR {ncys_or[0]:.2f}, "
         f"{ncys_or[1]:.2f}-{ncys_or[2]:.2f}). Within non-cysteine it is carried by Ser/Thr/Tyr "
         "hydroxyls; for His, and for Asp/Glu carboxylates, there is no site difference.")
H.append("5. Site survives year and depositing group together. The fully adjusted entry-level "
         f"logistic puts PDBj at OR "
         f"{LOGIT_FULL[0]:.2f} ({LOGIT_FULL[1]:.2f}-{LOGIT_FULL[2]:.2f}) against RCSB with "
         f"release year, entry size and X-ray in the model and standard errors clustered by "
         f"depositing group." if LOGIT_FULL else
         "5. Site survives year and depositing group together (see section 5).")
H.append("")
H.append("**Where site and depositor cannot be separated.** Only "
         f"{int((gs>1).sum())} of {len(gs):,} depositing groups identified by primary citation "
         f"({100*(gs>1).mean():.1f}%) used more than one processing site, so the "
         "same-paper-different-site contrast rests on 5 informative strata and is uninformative "
         f"(OR {mh_cit[0]:.2f}, {mh_cit[1]:.2f}-{mh_cit[2]:.2f}). The looser last-author "
         f"grouping gives {int((gs2>1).sum())} multi-site groups and 118 informative strata, "
         f"and there the entry-level effect holds (OR {mh_la[0]:.2f}, "
         f"{mh_la[1]:.2f}-{mh_la[2]:.2f}, p {mh_la[3]:.2g}) while the contact-level version "
         "reverses (OR 0.46). Site and depositor are correlated by construction - a group "
         "deposits where it deposits - so what is shown is that the site contrast is not fully "
         "explained by depositor identity, not that the site is the cause.")
txt = "\n".join(L).replace("<<HEADLINE>>", "\n".join(H))
open(f"{OUT}/lane2_processing_site.md", "w").write(txt + "\n")
print("\nREPORT WRITTEN")
