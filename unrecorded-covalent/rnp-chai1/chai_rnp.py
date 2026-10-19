"""Chai-1 on the six Runs N' Poses entries with an unrecorded covalent bond.

Same design as the Boltz-2 run: one protein chain (+ GDP/Mg for KRAS), ligand given
  smiles : SMILES only, no bond (the benchmark's input)
  bond   : same SMILES plus a covalent restraint from the nucleophile to the warhead atom
Scored the same way: ligand heavy-atom RMSD after pocket-CA superposition (CA within 10 A of the
crystal ligand), minimised over graph automorphisms; plus the nucleophile-to-warhead and
nucleophile-to-nearest-ligand-atom distances in each prediction.

Usage:
  python chai_rnp.py prepare  WORKDIR      # download crystal + CCD files, write FASTA/restraints
  python chai_rnp.py predict  WORKDIR      # run Chai-1 (GPU)
  python chai_rnp.py score    WORKDIR      # results.csv + summary.md
"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

# Same six entries and the same SMILES the Boltz-2 run used (manifest.json of that run).
SITES = [
    dict(entry="7EI4", ccd="J3U", chain="A", residue=90, res="C", patom="SG", latom="CAT", cof=[],
         smiles="C=CC#CC#CC=CCC(O)C(O)CCCCC(=O)O"),
    dict(entry="7F1V", ccd="7XF", chain="A", residue=211, res="K", patom="NZ", latom="CAL", cof=[], smiles="CC1=C(O)C(CNC(CCS)C(=O)O)=C(COP(=O)(O)O)C=N1"),
    dict(entry="7YCC", ccd="IQC", chain="A", residue=12, res="C", patom="SG", latom="C36", cof=["GDP", "MG"], smiles="CCC(=O)N1CC2(CCN(C3=NC(NC4CCN(C)CC4)=NC4=C(F)C(C5=C(C)C=CC6=C5C=NN6)=C(Cl)C=C34)CC2)C1"),
    dict(entry="7YCE", ccd="IQN", chain="A", residue=12, res="C", patom="SG", latom="C31", cof=["GDP", "MG"], smiles="CCC(=O)N1CC2(CCN(C3=NC(OC4CCN(CC)CC4)=NC4=C3C=C(Cl)C(C3=C(C)C=CC5=C3C=NN5)=C4F)CC2)C1"),
    dict(entry="8DVH", ccd="BO2", chain="A", residue=445, res="S", patom="OG", latom="B26", cof=[], smiles="CC(C)CC(NC(=O)C(CC1=CC=CC=C1)NC(=O)C1=CN=CC=N1)B(O)O"),
    dict(entry="8X6R", ccd="Y9D", chain="A", residue=12, res="C", patom="SG", latom="C46", cof=["GDP", "MG"], smiles="C=CC1=CC2=C(N3CCC4(CC3)CN(C(=O)CC)C4)N=C(OC3CCN(C)CC3)N=C2C(OCC)=C1C1=C(C)C=CC2=C1C=NN2"),
]
COFACTOR_SMILES = {
    "GDP": "Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O",
    "MG": "[Mg+2]",
}
RCSB_CIF = "https://files.rcsb.org/download/{}.cif"
RCSB_CCD = "https://files.rcsb.org/ligands/download/{}.cif"


def fetch(url, dest):
    import requests
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    dest.write_bytes(r.content)
    return dest


def ccd_mol(ccd_cif):
    """Heavy-atom RDKit mol from a CCD definition, atoms carry their CCD names."""
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
    order = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE, "TRIP": Chem.BondType.TRIPLE,
             "AROM": Chem.BondType.AROMATIC}
    for row in bonds:
        a1, a2 = row[0].strip("'\""), row[1].strip("'\"")
        if a1 in idx and a2 in idx:
            m.AddBond(idx[a1], idx[a2], order.get(row[2], Chem.BondType.SINGLE))
    m = m.GetMol()
    m.UpdatePropertyCache(strict=False)
    return m


def skeleton(mol):
    """Element + connectivity only (all single, non-aromatic) so CCD and SMILES graphs compare cleanly."""
    from rdkit import Chem
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


def smiles_mol(smi):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smi)
    return Chem.RemoveHs(m)


def chai_names(smi):
    """Chai-1 names SMILES ligand atoms element + per-element counter, in RDKit atom order (heavy atoms first)."""
    m = smiles_mol(smi)
    cnt, names = defaultdict(int), []
    for a in m.GetAtoms():
        cnt[a.GetSymbol()] += 1
        names.append((a.GetSymbol() + str(cnt[a.GetSymbol()])).upper())
    return names


def ccd_to_smiles_index(ref, probe):
    """All graph matches: tuple t with t[probe_idx] = ref_idx."""
    return skeleton(ref).GetSubstructMatches(skeleton(probe), uniquify=False, useChirality=False, maxMatches=20000)


def entity_info(cif_path, auth_chain):
    import gemmi
    blk = gemmi.cif.read(str(cif_path)).sole_block()
    scheme = blk.find("_pdbx_poly_seq_scheme.", ["asym_id", "entity_id", "seq_id", "auth_seq_num", "pdb_strand_id"])
    entity_id, amap = None, {}
    for row in scheme:
        if row[4] == auth_chain:
            entity_id = row[1]
            if row[3] not in (".", "?"):
                amap[int(row[3])] = int(row[2])
    ent = blk.find("_entity_poly.", ["entity_id", "pdbx_seq_one_letter_code_can"])
    seq = next(gemmi.cif.as_string(r[1]).replace("\n", "").replace(" ", "") for r in ent if r[0] == entity_id)
    return seq, amap


def ensure_smiles(s, struct):
    if s["smiles"]:
        return s["smiles"]
    from rdkit import Chem
    return Chem.MolToSmiles(ccd_mol(struct / f"{s['ccd']}.cif"))


def prepare(work):
    work = Path(work)
    struct, inp = work / "structures", work / "inputs"
    struct.mkdir(parents=True, exist_ok=True)
    inp.mkdir(parents=True, exist_ok=True)
    manifest = []
    for s in SITES:
        e = s["entry"]
        fetch(RCSB_CIF.format(e), struct / f"{e}.cif")
        fetch(RCSB_CCD.format(s["ccd"]), struct / f"{s['ccd']}.cif")
        seq, amap = entity_info(struct / f"{e}.cif", s["chain"])
        k = amap[s["residue"]]
        assert seq[k - 1] == s["res"], f"{e}: index {k} is {seq[k-1]}, expected {s['res']}"
        smi = ensure_smiles(s, struct)
        ref = ccd_mol(struct / f"{s['ccd']}.cif")
        probe = smiles_mol(smi)
        matches = ccd_to_smiles_index(ref, probe)
        assert matches, f"{e}: SMILES does not match CCD graph"
        ref_names = [a.GetProp("name") for a in ref.GetAtoms()]
        warhead_ref = ref_names.index(s["latom"])
        # warhead position in SMILES order (take every match; the warhead atom must map to one index)
        hits = sorted({m.index(warhead_ref) for m in matches if warhead_ref in m})
        assert len(hits) == 1, f"{e}: warhead maps to {hits} in SMILES order"
        names = chai_names(smi)
        try:  # cross-check against Chai's own naming when chai_lab is installed
            from chai_lab.data.sources.rdkit import RefConformerGenerator
            got = [n.upper() for n in RefConformerGenerator().generate(smi).atom_names]
            heavy = [n for n in got if not n.startswith("H")]
            assert heavy[: len(names)] == names, f"{e}: chai naming differs: {heavy[:8]} vs {names[:8]}"
        except ImportError:
            pass
        lat = names[hits[0]]
        fasta = [f">protein|name={e}-protein", seq, f">ligand|name={s['ccd']}", smi]
        for c in s["cof"]:
            fasta += [f">ligand|name={c}", COFACTOR_SMILES[c]]
        for cond in ("smiles", "bond"):
            d = inp / f"{e}_{cond}"
            d.mkdir(exist_ok=True)
            (d / "input.fasta").write_text("\n".join(fasta) + "\n")
            if cond == "bond":
                (d / "bond.restraints").write_text(
                    "chainA,res_idxA,chainB,res_idxB,connection_type,confidence,min_distance_angstrom,"
                    "max_distance_angstrom,comment,restraint_id\n"
                    f"A,{s['res']}{k}@{s['patom']},B,@{lat},covalent,1.0,0.0,0.0,protein-ligand,bond1\n")
        manifest.append(dict(entry=e, ccd=s["ccd"], chain=s["chain"], residue=s["residue"], seq_index=k,
                             patom=s["patom"], latom_ccd=s["latom"], latom_chai=lat, smiles=smi,
                             cofactors=s["cof"], seq_len=len(seq)))
        print(f"{e}: seq index {k}, warhead {s['latom']} -> {lat}, {len(names)} heavy atoms")
    (work / "manifest.json").write_text(json.dumps(manifest, indent=2))


def predict(work, samples=5, seed=42):
    import shutil
    import torch
    from chai_lab.chai1 import run_inference
    work = Path(work)
    for s in json.loads((work / "manifest.json").read_text()):
        for cond in ("smiles", "bond"):
            d = work / "inputs" / f"{s['entry']}_{cond}"
            out = work / "out" / f"{s['entry']}_{cond}"
            if (out / "done").exists():
                print("skip", out.name)
                continue
            if out.exists():
                shutil.rmtree(out)
            out.mkdir(parents=True)
            kw = dict(fasta_file=d / "input.fasta", output_dir=out, num_trunk_recycles=3,
                      num_diffn_timesteps=200, num_diffn_samples=samples, seed=seed,
                      device="cuda:0", use_esm_embeddings=True, use_msa_server=True)
            if cond == "bond":
                kw["constraint_path"] = d / "bond.restraints"
            print(f"=== {out.name}", flush=True)
            run_inference(**kw)
            (out / "done").write_text("ok")
            torch.cuda.empty_cache()


def crystal(cif_path, chain, ccd):
    import gemmi
    st = gemmi.read_structure(str(cif_path))
    lig, ca, nuc = {}, {}, {}
    for ch in st[0]:
        if ch.name != chain:
            continue
        for res in ch:
            if res.name == ccd and not lig:
                lig = {a.name: np.array(a.pos.tolist()) for a in res if a.element.name != "H"}
            a = res.find_atom("CA", "*")
            tab = gemmi.find_tabulated_residue(res.name)
            if a and tab is not None and tab.is_amino_acid():
                ca[res.seqid.num] = np.array(a.pos.tolist())
                nuc[res.seqid.num] = {x.name: np.array(x.pos.tolist()) for x in res}
    return lig, ca, nuc


def predicted(cif_path):
    """Chai output: chain A protein (residues numbered from 1), chain B the ligand."""
    import gemmi
    st = gemmi.read_structure(str(cif_path))
    lig, ca, atoms = {}, {}, {}
    for ch in st[0]:
        if ch.name == "A":
            for res in ch:
                i = res.seqid.num  # Chai numbers residues 1..N along the input sequence
                a = res.find_atom("CA", "*")
                if a:
                    ca[i] = np.array(a.pos.tolist())
                atoms[i] = {x.name: np.array(x.pos.tolist()) for x in res}
        elif ch.name == "B":
            for res in ch:
                for a in res:
                    if a.element.name == "H":
                        continue
                    # Chai-1 writes ligand atom names as <element><per-element counter>_<fragment>,
                    # e.g. C15_1. Strip the fragment suffix to get the bare name chai_names() builds.
                    n = a.name.upper()
                    base = n.rsplit("_", 1)[0] if "_" in n else n
                    assert base not in lig, f"duplicate ligand atom name after suffix strip: {n}"
                    lig[base] = np.array(a.pos.tolist())
    return lig, ca, atoms


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - pc @ R.T


def score(work):
    work = Path(work)
    manifest = json.loads((work / "manifest.json").read_text())
    rows = []
    for s in manifest:
        e, ccd = s["entry"], s["ccd"]
        struct = work / "structures"
        ref = ccd_mol(struct / f"{ccd}.cif")
        ref_names = [a.GetProp("name") for a in ref.GetAtoms()]
        probe = smiles_mol(s["smiles"])
        pnames = chai_names(s["smiles"])
        matches = ccd_to_smiles_index(ref, probe)
        lig_x, ca_x, nuc_x = crystal(struct / f"{e}.cif", s["chain"], ccd)
        _, amap = entity_info(struct / f"{e}.cif", s["chain"])
        cen = np.mean(list(lig_x.values()), axis=0)
        pocket = [r for r, xyz in ca_x.items() if np.linalg.norm(xyz - cen) < 10.0 and r in amap]
        for cond in ("smiles", "bond"):
            out = work / "out" / f"{e}_{cond}"
            cifs = sorted(out.glob("pred.model_idx_*.cif"))
            for cif in cifs:
                k = int(cif.stem.split("_")[-1])
                sc = out / f"scores.model_idx_{k}.npz"
                agg = None
                if sc.exists():
                    z = np.load(sc)
                    if "aggregate_score" in z.files:
                        agg = float(np.ravel(z["aggregate_score"])[0])
                lig_p, ca_p, atoms_p = predicted(cif)
                pk = [r for r in pocket if amap[r] in ca_p]
                R, t = kabsch(np.array([ca_p[amap[r]] for r in pk]), np.array([ca_x[r] for r in pk]))
                X = {n: v @ R.T + t for n, v in lig_p.items()}
                missing = [n for n in pnames if n not in X]
                if missing:
                    print(f"{e} {cond} {k}: ligand atoms missing from prediction: {missing[:5]}")
                    continue
                # bond-geometry check of the name mapping: every SMILES bond 0.9-2.1 A
                bad = sum(1 for b in probe.GetBonds()
                          if not 0.9 < np.linalg.norm(X[pnames[b.GetBeginAtomIdx()]] - X[pnames[b.GetEndAtomIdx()]]) < 2.1)
                best, best_m = 1e9, None
                for m in matches:
                    d = [np.sum((lig_x[ref_names[m[i]]] - X[pnames[i]]) ** 2) for i in range(len(pnames))
                         if ref_names[m[i]] in lig_x]
                    v = float(np.sqrt(np.mean(d)))
                    if v < best:
                        best, best_m = v, m
                # distances inside the prediction itself
                nucl = atoms_p[s["seq_index"]].get(s["patom"])
                wh_ref = ref_names.index(s["latom_ccd"])
                wh_idx = best_m.index(wh_ref)
                d_wh = float(np.linalg.norm(lig_p[pnames[wh_idx]] - nucl))
                d_near = float(min(np.linalg.norm(v - nucl) for v in lig_p.values()))
                cen_off = float(np.linalg.norm(np.mean([X[n] for n in pnames], axis=0) - cen))
                rows.append(dict(entry=e, ligand=ccd, condition=cond, model=k, aggregate_score=agg,
                                 ligand_rmsd_A=round(best, 3), warhead_to_nucleophile_A=round(d_wh, 2),
                                 nearest_ligand_atom_to_nucleophile_A=round(d_near, 2),
                                 centroid_offset_A=round(cen_off, 2), bad_bonds=bad,
                                 symmetry_matches=len(matches), pocket_ca=len(pk)))
                print(f"{e} {cond} model {k}: RMSD {best:.2f}  warhead {d_wh:.2f}  near {d_near:.2f}  "
                      f"agg {agg}  badbonds {bad}")
    if not rows:
        print("no predictions found")
        return
    with open(work / "results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    lines = ["| Entry | Ligand | SMILES-only top / best | Bond declared top / best | Warhead to nucleophile, SMILES-only (all samples) | Bond declared (all samples) |",
             "|---|---|---|---|---|---|"]
    for s in manifest:
        cells = []
        rng = []
        for cond in ("smiles", "bond"):
            r = [x for x in rows if x["entry"] == s["entry"] and x["condition"] == cond]
            if not r:
                cells.append("-")
                rng.append("-")
                continue
            top = (max(r, key=lambda x: x["aggregate_score"]) if all(x["aggregate_score"] is not None for x in r)
                   else min(r, key=lambda x: x["model"]))
            cells.append(f"{top['ligand_rmsd_A']:.2f} / {min(x['ligand_rmsd_A'] for x in r):.2f}")
            ws = [x["warhead_to_nucleophile_A"] for x in r]
            rng.append(f"{min(ws):.2f}-{max(ws):.2f}")
        lines.append(f"| {s['entry']} | {s['ccd']} | {cells[0]} | {cells[1]} | {rng[0]} | {rng[1]} |")
    (work / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    cmd, work = sys.argv[1], sys.argv[2]
    {"prepare": prepare, "predict": predict, "score": score}[cmd](work)
