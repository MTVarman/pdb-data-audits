"""Cross-check this scan against the audit repo's own published site tables."""
import glob, os, sys
import pandas as pd

D = "/home/claude/work/data"
REPO = "/home/claude/dl/pdb-data-audits/unrecorded-covalent"
con = pd.concat([pd.read_csv(p, sep="\t") for p in sorted(glob.glob(f"{D}/scan/contacts_*.tsv"))],
                ignore_index=True)
done = pd.concat([pd.read_csv(p, sep="\t", header=None, names=["entry", "n", "status"])
                  for p in sorted(glob.glob(f"{D}/scan/done_*.tsv"))], ignore_index=True)
scanned = set(done.loc[done.status == "ok", "entry"])
L = []
def say(s):
    print(s); L.append(str(s))

key = con.assign(k=con.entry + "|" + con.prot_chain + "|" + con.prot_res + con.prot_seq.astype(str)
                 + "|" + con.prot_atom + "|" + con.lig_comp + "|" + con.lig_atom)
rec_map = dict(zip(key.k, key.recorded))
dist_map = dict(zip(key.k, key.dist))

# ---- strict 485 (should all be UNRECORDED in this scan)
st = pd.read_csv(f"{REPO}/strict-cut-491/unrecorded_covalent_tight_cut_491_sites.csv")
st["k"] = (st.entry + "|" + st.chain + "|" + st.residue_name + st.residue_number.astype(str)
           + "|" + st.protein_atom + "|" + st.ligand_comp + "|" + st.ligand_atom)
st["in_scan_entry"] = st.entry.isin(scanned)
st["found"] = st.k.isin(rec_map)
st["recorded_now"] = st.k.map(rec_map)
say(f"strict-485: {len(st)} sites, {int(st.in_scan_entry.sum())} in entries this scan covered, "
    f"{int(st.found.sum())} matched at atom-pair level")
sub = st[st.found]
say(f"  of matched: {int((sub.recorded_now==0).sum())} still unrecorded, "
    f"{int((sub.recorded_now==1).sum())} now carry a covale record")
d = (sub.k.map(dist_map) - sub.distance_A).abs()
say(f"  distance agreement: max |delta| {d.max():.3f} A, median {d.median():.4f} A")

# ---- two-rulers recorded pairs (should all be RECORDED)
re_ = pd.read_csv(f"{REPO}/two-rulers/probe_line_recorded_pairs.csv")
re_["entry"] = re_.pdbid.str.upper()
pa = re_.protein_atom.str.split(":", expand=True)
lg = re_.ligand.str.rsplit("_", n=2, expand=True)
re_["k"] = (re_.entry + "|" + pa[0] + "|" + pa[1] + "|" + pa[2] + "|" + lg[0] + "|" + re_.ligand_atom)
re_["found"] = re_.k.isin(rec_map)
re_["recorded_now"] = re_.k.map(rec_map)
say(f"two-rulers recorded pairs: {len(re_)} rows, "
    f"{int(re_.entry.isin(scanned).sum())} in covered entries, {int(re_.found.sum())} matched")
s2 = re_[re_.found]
say(f"  of matched: {int((s2.recorded_now==1).sum())} recorded in this scan, "
    f"{int((s2.recorded_now==0).sum())} not")

# ---- two-rulers unrecorded pairs (should all be UNRECORDED)
un = pd.read_csv(f"{REPO}/two-rulers/probe_line_unrecorded_pairs.csv")
un["entry"] = un.pdbid.str.upper()
pa = un.protein_atom.str.split(":", expand=True)
lg = un.ligand.str.rsplit("_", n=2, expand=True)
un["k"] = (un.entry + "|" + pa[0] + "|" + pa[1] + "|" + pa[2] + "|" + lg[0] + "|" + un.ligand_atom)
un["found"] = un.k.isin(rec_map)
un["recorded_now"] = un.k.map(rec_map)
say(f"two-rulers unrecorded pairs: {len(un)} rows, "
    f"{int(un.entry.isin(scanned).sum())} in covered entries, {int(un.found.sum())} matched")
s3 = un[un.found]
flag = s3.residue_pair_bonded_elsewhere.astype(str).str.lower().eq("true")
say(f"  rows whose residue pair is NOT bonded elsewhere ({int((~flag).sum())} matched): "
    f"{int((s3.recorded_now[~flag]==0).sum())} unrecorded in this scan, "
    f"{int((s3.recorded_now[~flag]==1).sum())} recorded")
say(f"  rows flagged residue_pair_bonded_elsewhere ({int(flag.sum())} matched): "
    f"{int((s3.recorded_now[flag]==1).sum())} recorded in this scan, "
    f"{int((s3.recorded_now[flag]==0).sum())} not "
    f"(expected recorded: this scan's 'recorded' test is at residue-pair level, as in "
    f"cys-ligand-sites, so a pair bonded through another atom counts as recorded)")
open("/home/claude/work/validation.txt", "w").write("\n".join(L) + "\n")
