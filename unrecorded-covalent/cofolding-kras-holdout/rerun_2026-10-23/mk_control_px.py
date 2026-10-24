import json, pathlib, sys
W = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "work_control")
IONS = {"MG"}
man = json.load(open(W/"manifest.json"))
out = W/"inputs"/"protenix"; out.mkdir(parents=True, exist_ok=True)
ONE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
       "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
       "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
sys.path.insert(0, "/content")
from kras_holdout import entity_seq
made, bad = 0, []
for m in man:
    e = m["entry"]
    seq, _ = entity_seq(W/"structures"/f"{e}.cif", m["chain"])
    pos = m["seq_index"]
    if seq[pos-1] != ONE.get(m["res"], "X"):
        bad.append(f"{e}: seq[{pos}]={seq[pos-1]} != {m['res']}"); continue
    seqs = [{"proteinChain": {"sequence": seq, "count": 1}},
            {"ligand": {"ligand": m["smiles"], "count": 1}}]
    for c in m["cofactors"]:
        seqs.append({"ion": {"ion": c, "count": 1}} if c in IONS
                    else {"ligand": {"ligand": f"CCD_{c}", "count": 1}})
    job = {"name": f"{e}_smiles", "sequences": seqs}
    (out/f"{e}_smiles.json").write_text(json.dumps([job], indent=2)); made += 1
print(f"wrote {made} control protenix job files (smiles arm only)")
for b in bad: print("  SKIP", b)
assert not bad, "manifest/sequence mismatch"
