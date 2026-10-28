#!/usr/bin/env python3
"""Covalent KRAS temporal-holdout benchmark: Chai-1 and Boltz-2, two arms each.

Method follows chai_rnp.py / run_rnp_covalent.py (the six-entry RNP run) so the
numbers stay comparable: seed 42, 5 diffusion samples, MSAs from the ColabFold
server, arms

  smiles : ligand as SMILES, no bond declared   (what a benchmark normally feeds)
  bond   : covalent bond declared               (Boltz: CCD + bond constraint;
                                                 Chai: SMILES + bond.restraints)

Everything about each site -- chain, nucleophile atom, warhead atom, deposited
bond length, cofactors -- is read from the deposited mmCIF's own _struct_conn
record rather than trusted from the manifest. The manifest supplies the entry
list and is cross-checked against what the file says.

Two naming traps this handles, both of which silently corrupt a covalent run:

  1. Boltz-2 bond constraints reference ATOM NAMES. With a 'smiles:' ligand the
     atom names are generated, so a CCD name like C17 does not resolve and the
     constraint is quietly useless. The bond arm therefore uses 'ccd:'.
  2. Chai-1 renames SMILES ligand atoms element+counter in RDKit order. The
     restraint must use Chai's name for the warhead, not the CCD name. That
     mapping is recovered by a graph match between the CCD and SMILES molecules.

Usage:
  python kras_holdout.py prepare WORKDIR [--tier A] [--manifest run_manifest.csv]
  python kras_holdout.py predict-chai WORKDIR [--tier A]
  # Boltz runs from its own CLI, see the notebook cell
  python kras_holdout.py score WORKDIR
"""
import argparse, csv, json, sys
from collections import defaultdict
from pathlib import Path

RCSB_CIF = "https://files.rcsb.org/download/{}.cif"
RCSB_CCD = "https://files.rcsb.org/ligands/download/{}.cif"
NUCLEOTIDES = {"GDP", "GNP", "GTP", "GSP", "GCP"}
IONS = {"MG"}
AA3 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
       "MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","SEC","PYL"}
ONE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
       "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
       "THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M","SEC":"U","PYL":"O"}


def fetch(url, dest):
    import requests
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    r = requests.get(url, timeout=180)
    r.raise_for_status()
    dest.write_bytes(r.content)
    return dest


# ---------- RDKit / CCD helpers (same approach as chai_rnp.py) ----------
def ccd_mol(ccd_cif):
    """Heavy-atom RDKit mol from a CCD definition; atoms carry their CCD names."""
    import gemmi
    from rdkit import Chem
    blk = gemmi.cif.read(str(ccd_cif)).sole_block()
    atoms = blk.find("_chem_comp_atom.", ["atom_id", "type_symbol"])
    bonds = blk.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"])
    m, idx = Chem.RWMol(), {}
    for row in atoms:
        name, el = row[0].strip("'\""), row[1].capitalize()
        if el == "H":
            continue
        a = Chem.Atom(el)
        a.SetProp("name", name)
        idx[name] = m.AddAtom(a)
    order = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
             "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    for row in bonds:
        a1, a2 = row[0].strip("'\""), row[1].strip("'\"")
        if a1 in idx and a2 in idx:
            m.AddBond(idx[a1], idx[a2], order.get(row[2], Chem.BondType.SINGLE))
    m = m.GetMol()
    m.UpdatePropertyCache(strict=False)
    return m


def ccd_smiles(ccd_cif):
    import gemmi
    blk = gemmi.cif.read(str(ccd_cif)).sole_block()
    t = blk.find("_pdbx_chem_comp_descriptor.", ["type", "program", "descriptor"])
    cands = [(r[0].strip("'\""), r[1].strip("'\""), r[2].strip("'\"")) for r in t]
    for want_prog in ("OpenEye", None):
        for typ, prog, desc in cands:
            if "SMILES_CANONICAL" in typ.upper() and (want_prog is None or want_prog.lower() in prog.lower()):
                return desc
    for typ, _, desc in cands:
        if "SMILES" in typ.upper():
            return desc
    return None


def skeleton(mol):
    from rdkit import Chem
    rw = Chem.RWMol()
    for a in mol.GetAtoms():
        na = Chem.Atom(a.GetAtomicNum()); na.SetNoImplicit(True); rw.AddAtom(na)
    for b in mol.GetBonds():
        rw.AddBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), Chem.BondType.SINGLE)
    m = rw.GetMol(); m.UpdatePropertyCache(strict=False)
    return m


def smiles_mol(smi):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smi)
    return Chem.RemoveHs(m) if m else None


def chai_names(smi):
    """Chai-1 names SMILES ligand atoms element+counter in RDKit heavy-atom order."""
    m = smiles_mol(smi)
    cnt, names = defaultdict(int), []
    for a in m.GetAtoms():
        cnt[a.GetSymbol()] += 1
        names.append((a.GetSymbol() + str(cnt[a.GetSymbol()])).upper())
    return names


def warhead_in_smiles_order(ref, smi, ccd_atom):
    """Map a CCD atom name onto its index in SMILES order. None if ambiguous."""
    probe = smiles_mol(smi)
    if probe is None:
        return None
    matches = skeleton(ref).GetSubstructMatches(skeleton(probe), uniquify=False,
                                                useChirality=False, maxMatches=20000)
    if not matches:
        return None
    names = [a.GetProp("name") for a in ref.GetAtoms()]
    if ccd_atom not in names:
        return None
    ref_i = names.index(ccd_atom)
    hits = sorted({m.index(ref_i) for m in matches if ref_i in m})
    return hits[0] if len(hits) == 1 else None


# ---------- deposited-structure parsing ----------
def read_site(cif_path, want_res, want_resnum, want_lig):
    """Pull the covalent site straight from the entry's own _struct_conn."""
    import gemmi
    blk = gemmi.cif.read(str(cif_path)).sole_block()
    t = blk.find("_struct_conn.", ["conn_type_id",
                                   "ptnr1_label_comp_id", "ptnr1_label_seq_id", "ptnr1_label_atom_id",
                                   "ptnr1_auth_asym_id", "ptnr1_auth_seq_id",
                                   "ptnr2_label_comp_id", "ptnr2_label_seq_id", "ptnr2_label_atom_id",
                                   "ptnr2_auth_asym_id", "ptnr2_auth_seq_id", "pdbx_dist_value"])
    for r in t:
        if r[0].strip("'\"") != "covale":
            continue
        a = dict(comp=r[1].strip("'\""), seq=r[2], atom=r[3].strip("'\""),
                 chain=r[4].strip("'\""), auth=r[5], )
        b = dict(comp=r[6].strip("'\""), seq=r[7], atom=r[8].strip("'\""),
                 chain=r[9].strip("'\""), auth=r[10])
        for p, q in ((a, b), (b, a)):
            if p["comp"] in AA3 and q["comp"] not in AA3 and q["seq"] in (".", "?"):
                if p["comp"] == want_res and str(p["auth"]) == str(want_resnum) and q["comp"] == want_lig:
                    return dict(chain=p["chain"], res=p["comp"], resnum=int(p["auth"]),
                                patom=p["atom"], lig=q["comp"], latom=q["atom"],
                                lig_auth=q["auth"], lig_chain=q["chain"],
                                dist=float(r[11]) if r[11] not in (".", "?") else None)
    return None


def entity_seq(cif_path, auth_chain):
    """One-letter modelled sequence and auth->index map for the given auth chain."""
    import gemmi
    blk = gemmi.cif.read(str(cif_path)).sole_block()
    sch = blk.find("_pdbx_poly_seq_scheme.", ["asym_id", "seq_id", "auth_seq_num",
                                              "mon_id", "pdb_strand_id"])
    seq, amap = [], {}
    for row in sch:
        if row[4].strip("'\"") != auth_chain:
            continue
        mon = row[3].strip("'\"")
        seq.append(ONE.get(mon, "X"))
        a = row[2].strip("'\"")
        if a not in (".", "?"):
            amap[int(a)] = len(seq)
    return "".join(seq), amap


def cofactors_present(cif_path, auth_chain):
    """Nucleotide + Mg bound to this chain, as CCD codes."""
    import gemmi
    blk = gemmi.cif.read(str(cif_path)).sole_block()
    t = blk.find("_pdbx_nonpoly_scheme.", ["pdb_strand_id", "mon_id"])
    found = []
    for r in t:
        mon = r[1].strip("'\"")
        if mon in NUCLEOTIDES or mon in IONS:
            if mon not in found:
                found.append(mon)
    return found


# ---------- prepare ----------
def prepare(work, manifest, tier):
    work = Path(work)
    struct, inp = work / "structures", work / "inputs"
    struct.mkdir(parents=True, exist_ok=True)
    inp.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(manifest, newline="")))
    if tier:
        rows = [r for r in rows if r["tier"].startswith(tier)]
    out, problems = [], []

    for r in rows:
        e = r["pdb_id"]
        try:
            fetch(RCSB_CIF.format(e), struct / f"{e}.cif")
            site = read_site(struct / f"{e}.cif", r["residue_name"],
                             r["residue_num"], r["ligand_comp_id"])
            if site is None:
                problems.append(f"{e}: no matching covale row in deposited file"); continue
            if site["latom"] != r["ligand_atom"]:
                problems.append(f"{e}: ligand atom {site['latom']} != manifest {r['ligand_atom']}")
            fetch(RCSB_CCD.format(site["lig"]), struct / f"{site['lig']}.cif")
            ref = ccd_mol(struct / f"{site['lig']}.cif")
            smi = ccd_smiles(struct / f"{site['lig']}.cif")
            if not smi:
                problems.append(f"{e}: no SMILES for {site['lig']}"); continue
            seq, amap = entity_seq(struct / f"{e}.cif", site["chain"])
            k = amap.get(site["resnum"])
            if k is None or seq[k - 1] != ONE.get(site["res"], "X"):
                problems.append(f"{e}: residue {site['res']}{site['resnum']} not found in chain {site['chain']}")
                continue
            names = chai_names(smi)
            wi = warhead_in_smiles_order(ref, smi, site["latom"])
            chai_latom = names[wi] if wi is not None else None
            cof = cofactors_present(struct / f"{e}.cif", site["chain"])

            # ---- Chai-1 inputs ----
            fasta = [f">protein|name={e}-protein", seq, f">ligand|name={site['lig']}", smi]
            for c in cof:
                fetch(RCSB_CCD.format(c), struct / f"{c}.cif")
                csmi = ccd_smiles(struct / f"{c}.cif")
                if not csmi:
                    problems.append(f"{e}: no SMILES for cofactor {c}"); continue
                fasta += [f">ligand|name={c}", csmi]
            for arm in ("smiles", "bond"):
                d = inp / "chai" / f"{e}_{arm}"
                d.mkdir(parents=True, exist_ok=True)
                (d / "input.fasta").write_text("\n".join(fasta) + "\n")
                if arm == "bond":
                    if chai_latom is None:
                        problems.append(f"{e}: warhead {site['latom']} ambiguous in SMILES order; "
                                        f"chai bond arm skipped")
                        continue
                    (d / "bond.restraints").write_text(
                        "chainA,res_idxA,chainB,res_idxB,connection_type,confidence,"
                        "min_distance_angstrom,max_distance_angstrom,comment,restraint_id\n"
                        f"A,{ONE.get(site['res'],'X')}{k}@{site['patom']},B,@{chai_latom},"
                        f"covalent,1.0,0.0,0.0,protein-ligand,bond1\n")

            # ---- Boltz-2 inputs ----
            for arm in ("smiles", "bond"):
                d = inp / "boltz" / f"{e}_{arm}"
                d.mkdir(parents=True, exist_ok=True)
                L = ["version: 1", "sequences:",
                     "  - protein:", "      id: A", f"      sequence: {seq}"]
                if arm == "smiles":
                    L += ["  - ligand:", "      id: B", f"      smiles: '{smi}'"]
                else:
                    L += ["  - ligand:", "      id: B", f"      ccd: {site['lig']}"]
                for n, c in enumerate(cof):
                    L += ["  - ligand:", f"      id: {chr(ord('C') + n)}", f"      ccd: {c}"]
                if arm == "bond":
                    L += ["constraints:", "  - bond:",
                          f"      atom1: [A, {k}, {site['patom']}]",
                          f"      atom2: [B, 1, {site['latom']}]"]
                (d / f"{e}.yaml").write_text("\n".join(L) + "\n")

            out.append(dict(entry=e, tier=r["tier"], chain=site["chain"],
                            res=site["res"], resnum=site["resnum"], seq_index=k,
                            patom=site["patom"], lig=site["lig"],
                            latom_ccd=site["latom"], latom_chai=chai_latom,
                            deposited_bond_A=site["dist"],
                            manifest_bond_A=float(r["deposited_bond_A"]),
                            cofactors=cof, seq_len=len(seq), smiles=smi))
            print(f"{e} {site['res']}{site['resnum']}.{site['patom']} -> {site['lig']}.{site['latom']}"
                  f" (chai {chai_latom})  {site['dist']} A  cof={','.join(cof) or '-'}  seq={len(seq)}")
        except Exception as ex:
            problems.append(f"{e}: {type(ex).__name__}: {ex}")

    (work / "manifest.json").write_text(json.dumps(out, indent=2))
    print(f"\nprepared {len(out)} entries x 2 arms x 2 models")
    if problems:
        print("\nPROBLEMS:")
        for p in problems:
            print("  ", p)
    return out


# ---------- chai predict ----------
def predict_chai(work, tier, samples=5, seed=42, arm=None):
    import torch
    from chai_lab.chai1 import run_inference
    work = Path(work)
    man = json.loads((work / "manifest.json").read_text())
    keep = {m["entry"] for m in man if not tier or m["tier"].startswith(tier)}
    for d in sorted((work / "inputs" / "chai").iterdir()):
        if d.name.split("_")[0] not in keep:
            continue
        if arm and not d.name.endswith("_" + arm):
            continue
        out = work / "out" / "chai" / d.name
        if (out / "pred.model_idx_0.cif").exists():
            print("skip", d.name); continue
        out.mkdir(parents=True, exist_ok=True)
        kw = dict(fasta_file=d / "input.fasta", output_dir=out, num_trunk_recycles=3,
                  num_diffn_timesteps=200, seed=seed, device="cuda",
                  use_esm_embeddings=True, num_diffn_samples=samples)
        rp = d / "bond.restraints"
        if rp.exists():
            kw["constraint_path"] = rp
        print("running", d.name, "restraints" if rp.exists() else "")
        run_inference(**kw)
        torch.cuda.empty_cache()


# ---------- score ----------
def atoms_of(cif_path):
    import gemmi
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    return st


def score(work):
    import gemmi, math
    work = Path(work)
    man = {m["entry"]: m for m in json.loads((work / "manifest.json").read_text())}
    rows = []
    for model in ("chai", "boltz", "protenix"):
        root = work / "out" / model
        if not root.exists():
            continue
        for d in sorted(root.rglob("*.cif")):
            name = None
            for part in d.parts:
                cand = part.split("_")[0]
                if cand in man:
                    name, arm = cand, ("bond" if "_bond" in part else "smiles"); break
            if name is None:
                continue
            m = man[name]
            try:
                st = atoms_of(d)
            except Exception:
                continue
            nuc, lig_atoms, warhead = None, [], None
            for ch in st[0]:
                for res in ch:
                    if res.name in AA3:
                        if res.seqid.num in (m["seq_index"], m["resnum"]):
                            for a in res:
                                if a.name == m["patom"]:
                                    nuc = a.pos
                    elif res.name not in ("HOH",) and res.name not in NUCLEOTIDES and res.name not in IONS:
                        for a in res:
                            if a.element.name != "H":
                                lig_atoms.append((a.name, a.pos))
                                if a.name in (m["latom_ccd"], m["latom_chai"]):
                                    warhead = a.pos
            if nuc is None or not lig_atoms:
                continue
            def dist(p, q):
                return math.sqrt((p.x-q.x)**2 + (p.y-q.y)**2 + (p.z-q.z)**2)
            nearest = min(dist(nuc, p) for _, p in lig_atoms)
            rows.append(dict(model=model, entry=name, tier=m["tier"], arm=arm,
                             sample=d.name,
                             dist_warhead_A=round(dist(nuc, warhead), 3) if warhead else None,
                             dist_nearest_lig_atom_A=round(nearest, 3),
                             deposited_bond_A=m["deposited_bond_A"],
                             at_bond_distance=nearest <= 2.0))
    if not rows:
        print("no predictions found under out/ yet"); return
    with open(work / "results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    agg = defaultdict(list)
    for r in rows:
        agg[(r["model"], r["tier"], r["arm"])].append(r["dist_nearest_lig_atom_A"])
    print(f"{'model':<7}{'tier':<22}{'arm':<9}{'n':>4}{'median':>9}{'<=2A':>8}")
    for k in sorted(agg):
        v = sorted(agg[k]); n = len(v)
        med = v[n//2] if n % 2 else (v[n//2-1]+v[n//2])/2
        print(f"{k[0]:<7}{k[1]:<22}{k[2]:<9}{n:>4}{med:>9.2f}{str(sum(1 for x in v if x<=2.0))+'/'+str(n):>8}")
    print(f"\nwrote {work/'results.csv'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["prepare", "predict-chai", "score"])
    ap.add_argument("work")
    ap.add_argument("--manifest", default="run_manifest.csv")
    ap.add_argument("--tier", default=None)
    ap.add_argument("--arm", default=None, choices=["smiles", "bond"],
                    help="restrict predict-chai to one arm; default runs both")
    a = ap.parse_args()
    if a.cmd == "prepare":
        prepare(a.work, a.manifest, a.tier)
    elif a.cmd == "predict-chai":
        predict_chai(a.work, a.tier, arm=a.arm)
    else:
        score(a.work)


if __name__ == "__main__":
    main()
