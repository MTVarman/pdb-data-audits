#!/usr/bin/env python3
"""Protenix third arm for the covalent KRAS holdout.

Reads work/manifest.json written by `kras_holdout.py prepare`, so the sites,
sequences, SMILES and warhead atoms are the same objects the Chai-1 and Boltz-2
arms used. Nothing is re-derived and nothing can drift between arms.

Two facts established by reading protenix/data/inference/json_parser.py rather
than inferring from behaviour:

  1. Protenix names SMILES ligand atoms element symbol + per-element counter in
     RDKit atom order (rdkit_mol_to_atom_info). Chem.AddHs runs first but appends
     hydrogens after the heavy atoms, so heavy-atom numbering is unchanged. This
     is the SAME scheme chai_lab 0.6.1 uses -- verified identical on all 21
     ligands in this set. So manifest.json's `latom_chai` is also the Protenix
     name, and no new mapping is needed.

  2. Declaring the bond needs no CCD lookup when the ligand is given as SMILES:
     the bond references entity/position/atom directly. Protenix therefore does
     not inherit Boltz-2's frozen-dictionary failure, where 11 of 13 post-cutoff
     ligands could not be bond-declared at all.

Entity numbering in `covalent_bonds` is the 1-based index into `sequences`.

Models (configs/configs_model_type.py):
  protenix_base_default_v1.0.0    training cutoff 2021-09-30  (AF3-matched)
  protenix_base_20250630_v1.0.0   training cutoff 2025-06-30

Usage:
  python protenix_arm.py prepare WORKDIR [--tier A]
  python protenix_arm.py cmds    WORKDIR [--model protenix_base_default_v1.0.0]
  python protenix_arm.py collect WORKDIR
"""
import argparse, json, shutil
from pathlib import Path

NUCLEOTIDES = {"GDP", "GNP", "GTP", "GSP", "GCP"}
IONS = {"MG"}


def build_job(m, arm):
    """One Protenix inference job dict for an entry and arm."""
    seqs = [{"proteinChain": {"sequence": m["seq"], "count": 1}},
            {"ligand": {"ligand": m["smiles"], "count": 1}}]
    for c in m["cofactors"]:
        if c in IONS:
            seqs.append({"ion": {"ion": c, "count": 1}})
        else:
            seqs.append({"ligand": {"ligand": f"CCD_{c}", "count": 1}})

    job = {"name": f"{m['entry']}_{arm}", "sequences": seqs}

    if arm == "bond":
        job["covalent_bonds"] = [{
            "entity1": "1", "position1": str(m["seq_index"]), "atom1": m["patom"],
            "entity2": "2", "position2": "1",                 "atom2": m["latom_chai"],
        }]
    return job


def prepare(work, tier):
    work = Path(work)
    man = json.loads((work / "manifest.json").read_text())

    # kras_holdout.prepare does not store the sequence; recover it the same way.
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from kras_holdout import entity_seq

    out_dir = work / "inputs" / "protenix"
    out_dir.mkdir(parents=True, exist_ok=True)
    made, skipped = [], []

    for m in man:
        if tier and not m["tier"].startswith(tier):
            continue
        cif = work / "structures" / f"{m['entry']}.cif"
        seq, _ = entity_seq(cif, m["chain"])
        m = dict(m, seq=seq)
        for arm in ("smiles", "bond"):
            if arm == "bond" and not m.get("latom_chai"):
                skipped.append(f"{m['entry']}: warhead ambiguous in SMILES order, bond arm skipped")
                continue
            p = out_dir / f"{m['entry']}_{arm}.json"
            p.write_text(json.dumps([build_job(m, arm)], indent=2))
            made.append(p.name)
        print(f"{m['entry']}  {m['res']}{m['seq_index']}.{m['patom']} -> "
              f"{m['lig']}.{m['latom_chai']}  cof={','.join(m['cofactors']) or '-'}  seq={len(seq)}")

    print(f"\nwrote {len(made)} job files to {out_dir}")
    for s in skipped:
        print("  SKIP", s)
    return made


def cmds(work, model):
    work = Path(work)
    d = work / "inputs" / "protenix"
    lines = []
    for p in sorted(d.glob("*.json")):
        out = work / "out_raw" / "protenix" / p.stem
        lines.append(
            f'protenix pred -i "{p}" -o "{out}" -s 42 -n {model} '
            f'-e 5 -p 200 -c 10 --use_default_params true'
        )
    print("\n".join(lines))
    return lines


def collect(work):
    """Flatten Protenix output into work/out/protenix/<entry>_<arm>/ for scoring."""
    work = Path(work)
    raw = work / "out_raw" / "protenix"
    if not raw.exists():
        print("no output at", raw); return
    n = 0
    for job in sorted(raw.iterdir()):
        if not job.is_dir():
            continue
        dest = work / "out" / "protenix" / job.name
        dest.mkdir(parents=True, exist_ok=True)
        cifs = sorted(job.rglob("*.cif"))
        for i, c in enumerate(cifs):
            shutil.copy(c, dest / f"sample_{i}.cif")
            n += 1
        print(f"{job.name}: {len(cifs)} models")
    print(f"\ncollected {n} structures into {work/'out'/'protenix'}")
    print("now run:  python kras_holdout.py score", work)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["prepare", "cmds", "collect"])
    ap.add_argument("work")
    ap.add_argument("--tier", default=None)
    ap.add_argument("--model", default="protenix_base_default_v1.0.0")
    a = ap.parse_args()
    if a.cmd == "prepare":
        prepare(a.work, a.tier)
    elif a.cmd == "cmds":
        cmds(a.work, a.model)
    else:
        collect(a.work)


if __name__ == "__main__":
    main()
