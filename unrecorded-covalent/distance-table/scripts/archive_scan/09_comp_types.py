"""Chemical-component types for every ligand component seen in the scan (RCSB Data API)."""
import glob, json, os, sys
import pandas as pd
sys.path.insert(0, "/home/claude/work")
import rcsb

D = "/home/claude/work/data"
comps = set()
for p in sorted(glob.glob(f"{D}/scan/contacts_*.tsv")):
    df = pd.read_csv(p, sep="\t", usecols=["lig_comp"])
    comps |= set(df.lig_comp.astype(str))
comps = sorted(c for c in comps if c and c != "nan")
print("distinct ligand components in scan:", len(comps))

Q = """{chem_comps(comp_ids:%s){chem_comp{id type formula name}
 rcsb_chem_comp_info{atom_count_heavy}}}"""
rows = []
for i in range(0, len(comps), 200):
    d = rcsb.gql(Q % json.dumps(comps[i:i + 200]))
    for c in (d.get("chem_comps") or []):
        if not c or not c.get("chem_comp"):
            continue
        cc = c["chem_comp"]
        rows.append({"comp": cc["id"], "comp_type": cc.get("type"),
                     "formula": cc.get("formula"), "comp_name": cc.get("name"),
                     "heavy": (c.get("rcsb_chem_comp_info") or {}).get("atom_count_heavy")})
    if (i // 200) % 5 == 0:
        print(" ", i, flush=True)
t = pd.DataFrame(rows)
print("types:", t.comp_type.value_counts().head(20).to_string())
t.to_pickle(f"{D}/comp_types.pkl")
t.to_csv(f"{D}/comp_types.csv", index=False)
print("missing:", len(set(comps) - set(t.comp)))
