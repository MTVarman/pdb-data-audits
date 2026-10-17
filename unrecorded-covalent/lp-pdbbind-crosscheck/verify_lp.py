import csv, hashlib, re, sys
from collections import Counter

lp_path, our_path, issue_path = sys.argv[1], sys.argv[2], sys.argv[3]

# A. hash
h = hashlib.sha256(open(lp_path, 'rb').read()).hexdigest()
print("SHA256 repo file :", h)
print("SHA256 cited     : 7bb4e54e66c0c58ab74263edb8c2a2d677c2474cfc4ae4dec1eb31a231b10f8a")
print("MATCH" if h == "7bb4e54e66c0c58ab74263edb8c2a2d677c2474cfc4ae4dec1eb31a231b10f8a" else "NO MATCH")
print()

# B. columns + counts
rows = list(csv.DictReader(open(lp_path, newline='')))
cols = list(rows[0].keys())
print("rows:", len(rows))
print("columns:", cols)
def pick(pat):
    for c in cols:
        if re.search(pat, c, re.I): return c
    return None
c_id   = pick(r'^(pdb|pdbid|id|Unnamed: 0|)$') or cols[0]
c_cov  = pick(r'covalent')
c_split= pick(r'split')
c_smi  = pick(r'smiles')
c_hdr  = pick(r'header|category|type')
print("using id:", repr(c_id), "| covalent:", repr(c_cov), "| split:", repr(c_split), "| smiles:", repr(c_smi), "| header:", repr(c_hdr))
lp = {}
for r in rows:
    lp[r[c_id].strip().upper()] = r
print("covalent counts:", Counter(r[c_cov] for r in rows))
print("split counts   :", Counter(r[c_split] for r in rows) if c_split else "n/a")
print()

# C. strict 289 overlap
our = list(csv.DictReader(open(our_path, newline='')))
strict = sorted(set(r['entry'].upper() for r in our))
print("our strict entries:", len(strict))
hit = [e for e in strict if e in lp]
print("strict ∩ LP:", len(hit), "(doc claims 34)")
cov_split = Counter(lp[e][c_cov] for e in hit)
print("  covalent labels among hits:", dict(cov_split), "(doc claims 27 False / 7 True)")
print("  hits:", hit)
print()

# D. the 27 named entries
claims = {
 '3O0U':('test','confirmed'),'5Z1D':('test','confirmed'),'2R9M':('val','confirmed'),'3N3G':('val','confirmed'),
 '1NXY':('train','confirmed'),'1NYM':('train','confirmed'),'2A4R':('train','confirmed'),'2VJ1':('train','confirmed'),
 '3FV7':('train','confirmed'),'3FYZ':('train','confirmed'),'3FZC':('train','confirmed'),'4QLQ':('train','confirmed'),
 '4YV2':('train','confirmed'),'5YKP':('train','confirmed'),'5YR4':('train','confirmed'),'5YR5':('train','confirmed'),
 '5YR6':('train','confirmed'),'6NZG':('train','confirmed'),
 '5YJO':('train','plausible'),'6IJL':('train','plausible'),
 '3AAV':('test','doubtful'),'3KIG':('test','doubtful'),'4KB7':('train','doubtful'),'4R1E':('train','doubtful'),'4U6W':('train','doubtful'),
 '3BMO':('test','excluded'),'3JQF':('test','excluded'),
}
print(f"{'pdb':6}{'status':11}{'claimed_split':14}{'LP_split':10}{'LP_covalent':12}{'in_strict':10}verdict")
bad = 0
for e,(sp,st) in claims.items():
    r = lp.get(e)
    if r is None:
        print(f"{e:6}{st:11}{sp:14}{'—':10}{'NOT IN LP':12}{str(e in strict):10}MISSING"); bad += 1; continue
    lsp = r[c_split] if c_split else '?'
    lcv = r[c_cov]
    ok = (lsp == sp) and (str(lcv).strip().lower() in ('false','0'))
    if not ok: bad += 1
    print(f"{e:6}{st:11}{sp:14}{lsp:10}{str(lcv):12}{str(e in strict):10}{'ok' if ok else '** MISMATCH **'}")
print("mismatches:", bad)
print()

# names in LP not claimed, among strict hits
unclaimed = [e for e in hit if e not in claims]
print("strict hits NOT in the doc's 27 (should be the 7 covalent=True):")
for e in unclaimed:
    print("  ", e, "covalent =", lp[e][c_cov], "| split =", lp[e][c_split] if c_split else '?')
print()

# E. SMILES
print("SMILES checks:")
for e in ['2R9M','3N3G','3O0U','5Z1D','2A4R','1NXY','4QLQ']:
    r = lp.get(e)
    if not r: print(e, "not in LP"); continue
    s = r[c_smi]
    flags = []
    if 'C#N' in s or 'N#C' in s: flags.append('NITRILE C#N present')
    if re.search(r'C\(=N\)|=N\)|C=N', s): flags.append('C=N present')
    if re.search(r'C=C', s): flags.append('C=C present')
    if 'B(O)O' in s or 'B(O)' in s or 'OB(O)' in s: flags.append('boronic B(O)O present')
    print(f"  {e} [{r[c_split] if c_split else '?'}] {'; '.join(flags) or 'no flagged motifs'}")
    print(f"     {s}")
print()

# F. issue #8 overlap
txt = open(issue_path).read()
ids8 = set(re.findall(r'\b[1-9][A-Za-z0-9]{3}\b', txt))
ids8 = {i.upper() for i in ids8 if re.match(r'^[1-9][A-Z0-9]{3}$', i.upper())}
overlap = sorted(set(claims) & ids8)
print("PDB-like tokens in our filed issue #8 text:", len(ids8))
print("overlap with the doc's 27:", overlap)
