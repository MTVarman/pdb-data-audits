import json, sys
ONE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
       "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
       "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
W = sys.argv[1]
man = {m["entry"]: m for m in json.load(open(f"{W}/manifest.json"))}
allok = True
for e, m in man.items():
    j = json.load(open(f"{W}/inputs/protenix/{e}_bond.json"))[0]
    seq = j["sequences"][0]["proteinChain"]["sequence"]
    bd = j["covalent_bonds"][0]
    pos, exp = int(bd["position1"]), ONE[m["res"]]
    ok = seq[pos-1] == exp and bd["atom2"] == m["latom_chai"] and bd["entity2"] == "2"
    allok &= ok
    print(f"{e:<7} auth {m['resnum']:>4} posn {pos:>4} expect {exp} got {seq[pos-1]} {'OK' if ok else 'FAIL'}  {bd['atom2']}")
assert allok, "position or atom mismatch - stop"
print("all bond declarations check out")
