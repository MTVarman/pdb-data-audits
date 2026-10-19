#"""
#Runs N' Poses covalent re-run.
#"""
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np
import requests

ROOT = Path(__file__).resolve().parent
INPUTS = ROOT / "inputs"
STRUCT = ROOT / "structures"
OUT = ROOT / "out"

SITES = [
    dict(entry="7EI4", ccd="J3U", protein_chain="A", residue=90,  res_name="CYS", protein_atom="SG", ligand_atom="CAT", cofactors=[]),
    dict(entry="7F1V", ccd="7XF", protein_chain="A", residue=211, res_name="LYS", protein_atom="NZ", ligand_atom="CAL", cofactors=[]),
    dict(entry="7YCC", ccd="IQC", protein_chain="A", residue=12,  res_name="CYS", protein_atom="SG", ligand_atom="C36", cofactors=["GDP", "MG"]),
    dict(entry="7YCE", ccd="IQN", protein_chain="A", residue=12,  res_name="CYS", protein_atom="SG", ligand_atom="C31", cofactors=["GDP", "MG"]),
    dict(entry="8DVH", ccd="BO2", protein_chain="A", residue=445, res_name="SER", protein_atom="OG", ligand_atom="B26", cofactors=[]),
    dict(entry="8X6R", ccd="Y9D", protein_chain="A", residue=12,  res_name="CYS", protein_atom="SG", ligand_atom="C46", cofactors=["GDP", "MG"]),
]

RCSB_CIF = "https://files.rcsb.org/download/{}.cif"
RCSB_CCD = "https://files.rcsb.org/ligands/download/{}.cif"

def fetch(url, dest):
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    dest.write_bytes(r.content)
    return dest

def ccd_smiles_offline(ccd):
    from rdkit import Chem
    try:
        mol = ccd_mol(STRUCT / f"{ccd}.cif")
        return Chem.MolToSmiles(mol)
    except Exception as e:
        print(f"WARNING: Failed to parse SMILES offline for {ccd}: {e}")
        return ""

def entity_info(cif_path, auth_chain):
    import gemmi
    doc = gemmi.cif.read(str(cif_path))
    blk = doc.sole_block()
    scheme = blk.find("_pdbx_poly_seq_scheme.", ["asym_id", "entity_id", "seq_id", "auth_seq_num", "pdb_strand_id"])
    entity_id, amap = None, {}
    for row in scheme:
        if row[4] == auth_chain:
            entity_id = row[1]
            if row[3] not in (".", "?"):
                amap[int(row[3])] = int(row[2])
    if entity_id is None:
        raise RuntimeError(f"chain {auth_chain} not in poly_seq_scheme of {cif_path}")
    ent = blk.find("_entity_poly.", ["entity_id", "pdbx_seq_one_letter_code_can"])
    seq = None
    for row in ent:
        if row[0] == entity_id:
            seq = gemmi.cif.as_string(row[1]).replace("\n", "").replace(" ", "")
    return seq, amap

def write_yaml(path, seq, ligand_block, cofactors, bond=None):
    lines = ["version: 1", "sequences:",
             "  - protein:", "      id: A", f"      sequence: {seq}",
             "  - ligand:", "      id: L"] + ["      " + ligand_block]
    for i, c in enumerate(cofactors):
        lines += ["  - ligand:", f"      id: {'XYZW'[i]}", f"      ccd: {c}"]
    if bond:
        res_idx, patom, latom = bond
        lines += ["constraints:", "  - bond:",
                  f"      atom1: [A, {res_idx}, {patom}]",
                  f"      atom2: [L, 1, {latom}]"]
    Path(path).write_text("\n".join(lines) + "\n")

def prepare():
    INPUTS.mkdir(exist_ok=True)
    STRUCT.mkdir(exist_ok=True)
    manifest = []
    for s in SITES:
        e = s["entry"]
        cif = fetch(RCSB_CIF.format(e), STRUCT / f"{e}.cif")
        fetch(RCSB_CCD.format(s["ccd"]), STRUCT / f"{s['ccd']}.cif")
        seq, amap = entity_info(cif, s["protein_chain"])
        if s["residue"] not in amap:
            raise RuntimeError(f"{e}: residue {s['residue']} not mapped for chain {s['protein_chain']}")
        res_idx = amap[s["residue"]]
        assert seq[res_idx - 1] == {"CYS": "C", "LYS": "K", "SER": "S"}[s["res_name"]], \
            f"{e}: seq index {res_idx} is {seq[res_idx-1]}, expected {s['res_name']}"
        smi = ccd_smiles_offline(s["ccd"])
        if not smi:
            raise RuntimeError(f"Could not obtain a valid SMILES structure offline for ligand {s['ccd']}")
        write_yaml(INPUTS / f"{e}_smiles.yaml", seq, f"smiles: '{smi}'", s["cofactors"])
        write_yaml(INPUTS / f" {e}_bond.yaml", seq, f"ccd: {s['ccd']}", s["cofactors"],
                   bond=(res_idx, s["protein_atom"], s["ligand_atom"]))
        manifest.append(dict(s, seq_index=res_idx, smiles=smi, seq_len=len(seq)))
    (ROOT / "manifest.json").write_text(json.dumps(manifest, indent=2))

def ccd_mol(ccd_cif):
    import gemmi
    from rdkit import Chem
    blk = gemmi.cif.read(str(ccd_cif)).sole_block()
    atoms = blk.find("_chem_comp_atom.", ["atom_id", "type_symbol"])
    bonds = blk.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"])
    m = Chem.RWMol()
    idx = {}
    for row in atoms:
        name, el = row[0].strip("'\""), row[1].capitalize()
        if el == "H":
            continue
        a = Chem.Atom(el)
        a.SetProp("name", name)
        idx[name] = m.AddAtom(a)
    order = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE, "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    for row in bonds:
        a1, a2 = row[0].strip("'\""), row[1].strip("'\"")
        if a1 in idx and a2 in idx:
            m.AddBond(idx[a1], idx[a2], order.get(row[2], Chem.BondType.SINGLE))
    m = m.GetMol()
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    return m

def crystal_coords(cif_path, chain, ccd):
    import gemmi
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    model = st[0]
    lig, ca = {}, {}
    for ch in model:
        if ch.name != chain:
            continue
        for res in ch:
            if res.name == ccd and not lig:
                for at in res:
                    if at.element.name != "H":
                        lig[at.name] = np.array(at.pos.tolist())
            elif res.find_atom("CA", "*") is not None:
                a = res.find_atom("CA", "*")
                if a:
                    ca[res.seqid.num] = np.array(a.pos.tolist())
    if not lig:
        raise RuntimeError(f"ligand {ccd} not found in chain {chain} of {cif_path}")
    return lig, ca

def predicted_coords(pred_cif):
    import gemmi
    st = gemmi.read_structure(str(pred_cif))
    model = st[0]
    lig, ca = [], {}
    for ch in model:
        if ch.name == "L":
            for res in ch:
                for at in res:
                    if at.element.name != "H":
                        lig.append((at.name, at.element.name.capitalize(), np.array(at.pos.tolist())))
        elif ch.name == "A":
            for i, res in enumerate(ch, start=1):
                a = res.find_atom("CA", "*")
                if a:
                    ca[i] = np.array(a.pos.tolist())
    return lig, ca

def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, qc - pc @ R.T

def symmetric_rmsd(ref_mol, ref_xyz_by_name, probe_mol, probe_xyz_by_idx):
    from rdkit import Chem
    import numpy as np
    if ref_mol is None or probe_mol is None:
        return float("nan"), 0
    
    # Generate a standard canonical SMILES molecule for both to ensure identical ring/bond representations
    ref_smiles = Chem.MolToSmiles(ref_mol, canonical=True)
    probe_smiles = Chem.MolToSmiles(probe_mol, canonical=True)
    
    mol_ref_std = Chem.MolFromSmiles(ref_smiles)
    mol_probe_std = Chem.MolFromSmiles(probe_smiles)
    
    if mol_ref_std is None or mol_probe_std is None:
        return float("nan"), 0
        
    # Map standard molecules back to original to get correct mapping keys
    matches_ref = ref_mol.GetSubstructMatches(mol_ref_std, uniquify=False)
    matches_probe = probe_mol.GetSubstructMatches(mol_probe_std, uniquify=False)
    
    if not matches_ref or not matches_probe:
        # Fallback to direct mapping as backup
        mol_ref_std = ref_mol
        mol_probe_std = probe_mol
        matches_ref = [tuple(range(ref_mol.GetNumAtoms()))]
        matches_probe = [tuple(range(probe_mol.GetNumAtoms()))]

    # Find substructure match between standard representations
    matches = mol_ref_std.GetSubstructMatches(mol_probe_std, uniquify=False, useChirality=False, maxMatches=1000)
    if not matches:
        # Extremely relaxed SMARTS match fallback
        q = Chem.AdjustQueryParameters()
        q.makeDummiesQueries = True
        mol_ref_std_relaxed = Chem.AdjustQueryProperties(mol_ref_std, q)
        matches = mol_ref_std_relaxed.GetSubstructMatches(mol_probe_std, uniquify=False, maxMatches=1000)
        if not matches:
            return float("nan"), 0
            
    best = float("inf")
    ref_match_map = matches_ref[0]
    probe_match_map = matches_probe[0]
    
    for match in matches:
        d = []
        for pi_std, ri_std in enumerate(match):
            # Map std indexes back to original indices to get names and coordinates
            ri = ref_match_map[ri_std]
            pi = probe_match_map[pi_std]
            name = ref_mol.GetAtomWithIdx(ri).GetProp("name")
            if name in ref_xyz_by_name and pi in probe_xyz_by_idx:
                d.append(np.sum((ref_xyz_by_name[name] - probe_xyz_by_idx[pi]) ** 2))
        if d:
            best = min(best, float(np.sqrt(np.mean(d))))
    return best, len(matches)

def score():
    from rdkit import Chem
    manifest = json.loads((ROOT / "manifest.json").read_text())
    rows = []
    for s in manifest:
        e, ccd = s["entry"], s["ccd"]
        ref_mol = ccd_mol(STRUCT / f"{ccd}.cif")
        lig_xtal, ca_xtal = crystal_coords(STRUCT / f"{e}.cif", s["protein_chain"], ccd)
        _, amap = entity_info(STRUCT / f"{e}.cif", s["protein_chain"])
        lig_center = np.mean(list(lig_xtal.values()), axis=0)
        pocket = [r for r, xyz in ca_xtal.items() if np.linalg.norm(xyz - lig_center) < 10.0 and r in amap]
        for cond in ("smiles", "bond"):
            name = f"{e}_{cond}"
            pdir = OUT / "boltz_results_inputs" / "predictions" / name
            if not pdir.exists():
                continue
            for model_cif in sorted(pdir.glob(f"{name}_model_*.cif")):
                k = int(model_cif.stem.split("_model_")[1])
                conf_file = pdir / f"confidence_{name}_model_{k}.json"
                conf = json.loads(conf_file.read_text()) if conf_file.exists() else {}
                lig_pred, ca_pred = predicted_coords(model_cif)
                P = np.array([ca_pred[amap[r]] for r in pocket if amap[r] in ca_pred])
                Q = np.array([ca_xtal[r] for r in pocket if amap[r] in ca_pred])
                if len(P) < 3:
                    continue
                R, t = kabsch(P, Q)
                pred_xyz = {i: xyz @ R.T + t for i, (_, _, xyz) in enumerate(lig_pred)}
                
                # Fallback to topology mapping via SMILES molecular graphs
                if cond == "bond":
                    probe = ref_mol
                    name_to_idx = {ref_mol.GetAtomWithIdx(i).GetProp("name"): i for i in range(ref_mol.GetNumAtoms())}
                    pred_by_idx = {name_to_idx[n]: xyz for n, xyz in ((nm, pred_xyz[i]) for i, (nm, _, _) in enumerate(lig_pred)) if n in name_to_idx}
                else:
                    probe = Chem.MolFromSmiles(s["smiles"])
                    # Strip hydrogens to match predicted heavy atom list
                    probe = Chem.RemoveHs(probe)
                    pred_by_idx = pred_xyz
                
                rmsd, nmatch = symmetric_rmsd(ref_mol, lig_xtal, probe, pred_by_idx)
                rows.append(dict(entry=e, ligand=ccd, condition=cond, model=k,
                                 confidence=conf.get("confidence_score"), iptm=conf.get("iptm"),
                                 ligand_rmsd_A=round(rmsd, 3), symmetry_matches=nmatch, pocket_ca=len(P)))
                print(f"{name} model {k}: RMSD {rmsd:.2f} Å  (conf {conf.get('confidence_score')})")
    if not rows:
        print("no predictions found under out/"); return
    with open(ROOT / "results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    lines = ["| Entry | Ligand | SMILES-only top / best | Bond declared top / best | Gap (top) |", "|---|---|---|---|---|"]
    for s in manifest:
        e = s["entry"]
        def pick(cond):
            r = [x for x in rows if x["entry"] == e and x["condition"] == cond]
            if not r: return None, None
            top = next((x["ligand_rmsd_A"], x) for x in r if x["model"] == 0)[0] if any(x["model"] == 0 for x in r) else None
            return top, min(x["ligand_rmsd_A"] for x in r)
        st, sb = pick("smiles"); bt, bb = pick("bond")
        gap = None if st is None or bt is None else round(st - bt, 2)
        lines.append(f"| {e} | {s['ccd']} | {st} / {sb} | {bt} / {bb} | {gap} |")
    (ROOT / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))

if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    {"prepare": prepare, "score": score}.get(cmd, lambda: print(__doc__))()
