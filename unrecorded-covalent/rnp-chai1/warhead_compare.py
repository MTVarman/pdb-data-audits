"""Warhead-to-nucleophile distance for Boltz-2 and Chai-1 predictions, measured one way.

Neither model names its predicted ligand atoms the way the CCD does, and the two models do not
name them the same way as each other.  So instead of trusting any naming scheme, the predicted
ligand's bond graph is inferred from its own coordinates (covalent radii), graph-matched to the
CCD definition, and the warhead is whichever predicted atom maps to the CCD's warhead atom.
Identical treatment for both models, and for the crystal (where it must return the bond length).

Usage: python warhead_compare.py BOLTZ_PRED_ROOT CHAI_WORK STRUCT_DIR OUT_CSV
"""
import csv
import json
import sys
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem

TOL = 0.45  # A added to the sum of covalent radii when inferring a bond


def ccd_mol(ccd_cif):
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
    for row in bonds:
        a1, a2 = row[0].strip("'\""), row[1].strip("'\"")
        if a1 in idx and a2 in idx:
            m.AddBond(idx[a1], idx[a2], Chem.BondType.SINGLE)
    m = m.GetMol()
    m.UpdatePropertyCache(strict=False)
    return m


def geom_mol(elems, xyz):
    """Skeleton of a predicted ligand, bonds inferred from its own geometry."""
    rw = Chem.RWMol()
    for e in elems:
        a = Chem.Atom(e)
        a.SetNoImplicit(True)
        rw.AddAtom(a)
    r = [gemmi.Element(e).covalent_r for e in elems]
    n = len(elems)
    for i in range(n):
        for j in range(i + 1, n):
            if np.linalg.norm(xyz[i] - xyz[j]) < r[i] + r[j] + TOL:
                rw.AddBond(i, j, Chem.BondType.SINGLE)
    m = rw.GetMol()
    m.UpdatePropertyCache(strict=False)
    return m


def flat(mol):
    """Element + connectivity only, so CCD and geometry-derived graphs compare cleanly."""
    rw = Chem.RWMol()
    for a in mol.GetAtoms():
        na = Chem.Atom(a.GetAtomicNum())
        na.SetNoImplicit(True)
        rw.AddAtom(na)
    for b in mol.GetBonds():
        rw.AddBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), Chem.BondType.SINGLE)
    m = rw.GetMol()
    m.UpdatePropertyCache(strict=False)
    return m


def warhead_candidates(ref, elems, xyz, latom):
    """Indices in the predicted ligand that the CCD warhead atom can map to."""
    probe = geom_mol(elems, xyz)
    ref_names = [a.GetProp("name") for a in ref.GetAtoms()]
    wh_ref = ref_names.index(latom)
    matches = flat(ref).GetSubstructMatches(flat(probe), uniquify=False, useChirality=False,
                                            maxMatches=20000)
    if not matches:
        return [], 0, probe.GetNumBonds()
    cands = sorted({i for m in matches for i in range(len(m)) if m[i] == wh_ref})
    return cands, len(matches), probe.GetNumBonds()


def read_pred(cif, lig_chain, prot_chain="A", renumber=False):
    st = gemmi.read_structure(str(cif))
    names, elems, xyz, prot = [], [], [], {}
    for ch in st[0]:
        if ch.name == lig_chain:
            for res in ch:
                for a in res:
                    if a.element.name != "H":
                        names.append(a.name)
                        elems.append(a.element.name.capitalize())
                        xyz.append(np.array(a.pos.tolist()))
        elif ch.name == prot_chain:
            for i, res in enumerate(ch, start=1):
                k = i if renumber else res.seqid.num
                prot[k] = {a.name: np.array(a.pos.tolist()) for a in res}
    return names, elems, np.array(xyz), prot


def read_crystal(cif, chain, ccd):
    st = gemmi.read_structure(str(cif))
    names, elems, xyz, prot = [], [], [], {}
    for ch in st[0]:
        if ch.name != chain:
            continue
        for res in ch:
            if res.name == ccd and not names:
                for a in res:
                    if a.element.name != "H":
                        names.append(a.name)
                        elems.append(a.element.name.capitalize())
                        xyz.append(np.array(a.pos.tolist()))
            prot[res.seqid.num] = {a.name: np.array(a.pos.tolist()) for a in res}
    return names, elems, np.array(xyz), prot


def main():
    boltz_root, chai_work, struct, out_csv = (Path(p) for p in sys.argv[1:5])
    manifest = json.loads((chai_work / "manifest.json").read_text())
    rows = []
    for s in manifest:
        e, ccd, latom = s["entry"], s["ccd"], s["latom_ccd"]
        ref = ccd_mol(struct / f"{ccd}.cif")

        # crystal: the measurement must return the deposited bond length
        cn, ce, cx, cp = read_crystal(struct / f"{e}.cif", s["chain"], ccd)
        nuc_x = cp[s["residue"]][s["patom"]]
        cand, nm, nb = warhead_candidates(ref, ce, cx, latom)
        dx = [float(np.linalg.norm(cx[i] - nuc_x)) for i in cand]
        near_x = float(min(np.linalg.norm(v - nuc_x) for v in cx))
        print(f"{e} CRYSTAL warhead {min(dx):.2f} (cands {[cn[i] for i in cand]}, "
              f"matches {nm}) nearest {near_x:.2f}")
        rows.append(dict(entry=e, ligand=ccd, model="crystal", condition="-", sample="-",
                         warhead_A=round(min(dx), 2), nearest_A=round(near_x, 2),
                         warhead_is_nearest=abs(min(dx) - near_x) < 0.01,
                         n_warhead_candidates=len(cand)))

        for tag, root, lig_ch, pat, renum in (
                ("Boltz-2", boltz_root, "L", "{e}_{c}/{e}_{c}_model_{k}.cif", True),
                ("Chai-1", chai_work / "out", "B", "{e}_{c}/pred.model_idx_{k}.cif", False)):
            for cond in ("smiles", "bond"):
                for k in range(5):
                    f = root / pat.format(e=e, c=cond, k=k)
                    if not f.exists():
                        continue
                    pn, pe, px, pp = read_pred(f, lig_ch, renumber=renum)
                    nuc = pp[s["seq_index"]][s["patom"]]
                    cand, nm, nb = warhead_candidates(ref, pe, px, latom)
                    if not cand:
                        print(f"  !! {tag} {e} {cond} {k}: no graph match "
                              f"({len(pe)} atoms, {nb} inferred bonds)")
                        continue
                    d = min(float(np.linalg.norm(px[i] - nuc)) for i in cand)
                    near = float(min(np.linalg.norm(v - nuc) for v in px))
                    rows.append(dict(entry=e, ligand=ccd, model=tag, condition=cond, sample=k,
                                     warhead_A=round(d, 2), nearest_A=round(near, 2),
                                     warhead_is_nearest=abs(d - near) < 0.01,
                                     n_warhead_candidates=len(cand)))
                    if k == 0:
                        print(f"  {tag} {e} {cond} 0: warhead {d:.2f} (pred atom "
                              f"{[pn[i] for i in cand]}) nearest {near:.2f}")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote", out_csv, len(rows), "rows")


if __name__ == "__main__":
    main()
