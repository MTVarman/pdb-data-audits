"""Ligand pose RMSD vs deposited, symmetry-corrected, for co-folding predictions.

Superposes the predicted protein onto the deposited protein using CA atoms of the
residues present in both, applies that transform to the predicted ligand, then
reports the minimum heavy-atom RMSD over the ligand's automorphism group.

Atom correspondence:
  predicted ligand atoms are named element+per-element-counter in RDKit order of
  the SMILES (verified identical for chai_lab 0.6.1 and Protenix v1.0.0), so the
  predicted name maps directly onto a SMILES mol atom index.
  deposited ligand atoms carry CCD names, mapped onto SMILES indices by a graph
  match between the CCD mol and the SMILES mol.

Emits one row per (model, entry, arm, sample) with pose RMSD alongside the
nucleophile-to-nearest-ligand-heavy-atom distance, so the two metrics can be
compared directly.
"""
import csv, json, math, sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import gemmi
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

sys.path.insert(0, "/content")
from kras_holdout import ccd_mol, ccd_smiles, chai_names, skeleton, smiles_mol, entity_seq


def boltz_names(smi):
    """Boltz-2 names SMILES ligand atoms element + canonical rank over the H-added mol
    (boltz/data/parse/schema.py). Heavy atoms keep their heavy-atom indices because AddHs appends."""
    from rdkit.Chem import AllChem
    m = AllChem.AddHs(AllChem.MolFromSmiles(smi))
    rank = list(AllChem.CanonicalRankAtoms(m))
    nheavy = smiles_mol(smi).GetNumAtoms()
    return [(m.GetAtomWithIdx(i).GetSymbol().upper() + str(rank[i] + 1)) for i in range(nheavy)]

AA3 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
       "MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","SEC","PYL"}
COF = {"HOH","GDP","GNP","GTP","GSP","GCP","MG"}


def kabsch(P, Q):
    """Rotation+translation taking P onto Q (both (N,3)); returns R, t."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


def ccd_to_smiles_map(ref, smi):
    """Full map {CCD atom name -> SMILES mol atom index}, or None if ambiguous."""
    probe = smiles_mol(smi)
    if probe is None:
        return None
    matches = skeleton(ref).GetSubstructMatches(skeleton(probe), uniquify=False,
                                                useChirality=False, maxMatches=20000)
    if not matches:
        return None
    names = [a.GetProp("name") for a in ref.GetAtoms()]
    # matches[k][smiles_idx] = ref_idx ; invert, keep only unambiguous assignments
    per_ref = defaultdict(set)
    for m in matches:
        for smi_idx, ref_idx in enumerate(m):
            per_ref[ref_idx].add(smi_idx)
    out = {}
    for ref_idx, cands in per_ref.items():
        if len(cands) == 1:
            out[names[ref_idx]] = next(iter(cands))
    return out or None


def automorphisms(mol, cap=5000):
    ms = mol.GetSubstructMatches(mol, uniquify=False, useChirality=False, maxMatches=cap)
    return ms if ms else [tuple(range(mol.GetNumAtoms()))]


def protein_ca(st, want_chain=None):
    out = {}
    for ch in st[0]:
        for res in ch:
            if res.name in AA3:
                for a in res:
                    if a.name == "CA":
                        out[res.seqid.num] = np.array([a.pos.x, a.pos.y, a.pos.z])
    return out


def ligand_atoms(st, exclude=COF):
    out = {}
    for ch in st[0]:
        for res in ch:
            if res.name in AA3 or res.name in exclude:
                continue
            for a in res:
                if a.element.name != "H":
                    out[a.name] = np.array([a.pos.x, a.pos.y, a.pos.z])
    return out


def nucleophile(st, seq_index, resnum, patom):
    for ch in st[0]:
        for res in ch:
            if res.name in AA3 and res.seqid.num in (seq_index, resnum):
                for a in res:
                    if a.name == patom:
                        return np.array([a.pos.x, a.pos.y, a.pos.z])
    return None


def main(workdir):
    W = Path(workdir)
    man = {m["entry"]: m for m in json.loads((W / "manifest.json").read_text())}
    rows, problems = [], []

    for model in ("chai", "boltz", "protenix"):
        root = W / "out" / model
        if not root.exists():
            continue
        for jobdir in sorted(p for p in root.iterdir() if p.is_dir()):
            stem = jobdir.name
            entry = stem.split("_")[0]
            if entry not in man:
                continue
            arm = "bond" if stem.endswith("_bond") else "smiles"
            m = man[entry]

            dep_path = W / "structures" / f"{entry}.cif"
            lig_cif = W / "structures" / f"{m['lig']}.cif"
            if not dep_path.exists() or not lig_cif.exists():
                problems.append(f"{stem}: missing deposited files"); continue

            dep = gemmi.read_structure(str(dep_path))
            dep.setup_entities()
            dep.remove_alternative_conformations()
            smi = m["smiles"]
            mol = smiles_mol(smi)
            if mol is None:
                problems.append(f"{stem}: SMILES unparseable"); continue
            names = boltz_names(smi) if model == "boltz" else chai_names(smi)
            nmap = {n: i for i, n in enumerate(names)}
            cmap = ccd_to_smiles_map(ccd_mol(lig_cif), smi)
            if not cmap:
                problems.append(f"{stem}: no CCD->SMILES map"); continue

            dep_lig_named = ligand_atoms(dep)
            dep_coords = {}
            for nm, pos in dep_lig_named.items():
                if nm in cmap:
                    dep_coords[cmap[nm]] = pos
            if len(dep_coords) < 4:
                problems.append(f"{stem}: deposited ligand atoms unmapped"); continue

            _, amap = entity_seq(dep_path, m["chain"])
            dep_ca_auth = protein_ca(dep)
            dep_ca_idx = {amap[a]: xyz for a, xyz in dep_ca_auth.items() if a in amap}

            autos = automorphisms(mol)

            for cif in sorted(jobdir.glob("*.cif")):
                try:
                    pred = gemmi.read_structure(str(cif))
                    pred.setup_entities()
                except Exception as ex:
                    problems.append(f"{cif.name}: {ex}"); continue

                pred_ca = protein_ca(pred)
                common = sorted(set(pred_ca) & set(dep_ca_idx))
                if len(common) < 20:
                    problems.append(f"{stem}/{cif.name}: only {len(common)} common CA"); continue
                P = np.array([pred_ca[i] for i in common])
                Q = np.array([dep_ca_idx[i] for i in common])
                R, t = kabsch(P, Q)
                ca_rmsd = float(np.sqrt((((P @ R.T) + t - Q) ** 2).sum(1).mean()))

                pred_lig = ligand_atoms(pred)
                pred_coords = {}
                for nm, pos in pred_lig.items():
                    if nm in nmap:
                        pred_coords[nmap[nm]] = (R @ pos) + t
                shared = sorted(set(pred_coords) & set(dep_coords))
                if len(shared) < 4:
                    problems.append(f"{stem}/{cif.name}: {len(shared)} shared ligand atoms"); continue

                best = None
                for auto in autos:
                    # auto maps canonical index -> permuted index
                    pts, ref = [], []
                    ok = True
                    for i in shared:
                        j = auto[i]
                        if j not in pred_coords:
                            ok = False; break
                        pts.append(pred_coords[j]); ref.append(dep_coords[i])
                    if not ok or not pts:
                        continue
                    d = np.sqrt((((np.array(pts) - np.array(ref)) ** 2).sum(1)).mean())
                    if best is None or d < best:
                        best = float(d)
                if best is None:
                    problems.append(f"{stem}/{cif.name}: no valid automorphism"); continue

                nuc = nucleophile(pred, m["seq_index"], m["resnum"], m["patom"])
                near = None
                if nuc is not None and pred_lig:
                    near = min(float(np.linalg.norm(nuc - p)) for p in pred_lig.values())

                rows.append(dict(model=model, entry=entry, tier=m["tier"], arm=arm,
                                 sample=cif.name,
                                 pose_rmsd_A=round(best, 3),
                                 ca_rmsd_A=round(ca_rmsd, 3),
                                 dist_nearest_lig_atom_A=round(near, 3) if near is not None else None,
                                 n_lig_atoms_matched=len(shared),
                                 deposited_bond_A=m["deposited_bond_A"]))

    if not rows:
        print("no rows produced"); return
    outp = W / "pose_rmsd.csv"
    with open(outp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    agg = defaultdict(list)
    for r in rows:
        agg[(r["model"], r["tier"], r["arm"])].append(r)
    print(f"{'model':<9}{'tier':<22}{'arm':<8}{'n':>4}{'poseRMSD':>10}{'<=2A pose':>11}{'bond<=2A':>10}")
    for k in sorted(agg):
        v = agg[k]
        pr = sorted(x["pose_rmsd_A"] for x in v)
        med = pr[len(pr)//2] if len(pr) % 2 else (pr[len(pr)//2-1]+pr[len(pr)//2])/2
        npose = sum(1 for x in v if x["pose_rmsd_A"] <= 2.0)
        nbond = sum(1 for x in v if x["dist_nearest_lig_atom_A"] is not None
                    and x["dist_nearest_lig_atom_A"] <= 2.0)
        print(f"{k[0]:<9}{k[1]:<22}{k[2]:<8}{len(v):>4}{med:>10.2f}"
              f"{str(npose)+'/'+str(len(v)):>11}{str(nbond)+'/'+str(len(v)):>10}")

    disagree = [r for r in rows if r["dist_nearest_lig_atom_A"] is not None
                and r["pose_rmsd_A"] <= 2.0 and r["dist_nearest_lig_atom_A"] > 2.0]
    print(f"\nsamples a pose cut accepts but the bond metric rejects: {len(disagree)} / {len(rows)}")
    for r in sorted(disagree, key=lambda x: x["pose_rmsd_A"])[:15]:
        print(f"  {r['model']:<9}{r['entry']} {r['arm']:<7} pose {r['pose_rmsd_A']:>5.2f}  "
              f"nearest {r['dist_nearest_lig_atom_A']:>6.2f}")
    print(f"\nwrote {outp}  ({len(rows)} rows)")
    if problems:
        print(f"\n{len(problems)} problems:")
        for p in problems[:15]:
            print("  ", p)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "work_px")
