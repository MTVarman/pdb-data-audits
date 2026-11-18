"""Unrecorded contacts on residues other than cysteine, from the validator's bond line to 0.40 A beyond it,
sorted by the fixed rule of sample_scan/analyze2.py (verdict) on the columns written by sample_scan/evidence.py.
Usage: python3 -I past_line_triage.py data/sample_window_evidence.csv.gz data/past_line_triage.csv"""
import csv, gzip, sys, collections
COV = {'C': 0.77, 'N': 0.70, 'O': 0.66, 'S': 1.04, 'P': 1.10, 'B': 0.86, 'Se': 1.17}
def op(path): return gzip.open(path, 'rt', newline='') if path.endswith('.gz') else open(path, newline='')
def verdict(r):
    if r['hbond_capable_pair'] == 'True': return 'O-O, N-O or N-N pair'
    if r['family'] == 'additive_cryoprotectant': return 'additive'
    if r['family'] == 'class1_cofactor_ptm': return 'cofactor with known covalent attachment'
    if r['family'] == 'class3_nucleotide_phosphate': return 'phosphoryl or nucleotide'
    if r['residue_already_bonded'] == 'True': return 'residue already bonded through another atom'
    if r['ligand_atom_tier'] == 'C': return 'unreactive ligand atom'
    if r['ligand_atom_tier'] in ('A', 'B') and r['lig_el'] in ('C', 'P', 'S', 'B'):
        if r['ligand_atom_kind'].startswith('C=') or r['ligand_atom_kind'] in ('C#N / C#C',): return 'bond or near-attack complex'
        return 'candidate bond'
    return 'other'
out = []
for r in csv.DictReader(op(sys.argv[1])):
    if not (r['occupancy_class'] == 'full' and float(r['prot_occ'] or 0) > 0 and float(r['lig_occ'] or 0) > 0): continue
    if r['residue'] == 'CYS' or r['recorded_covale'] == 'True': continue
    pl = round(COV[r['lig_el']] + COV[r['prot_el']] + 0.2, 2); d = float(r['dist_A'])
    if not (pl < d <= pl + 0.40): continue
    out.append({'pdbid': r['pdbid'], 'protein_atom': r['protein_atom'], 'ligand': r['ligand'], 'ligand_atom': r['ligand_atom'], 'pair': r['pair'], 'dist_A': r['dist_A'],
                'validator_line_A': pl, 'triage': verdict(r), 'density_class': r['density_class'], 'clash_listed_at_pair': r['clash_listed'],
                'ligand_atom_tier': r['ligand_atom_tier'], 'ligand_atom_kind': r['ligand_atom_kind'], 'family': r['family'],
                'resolution': r['resolution'], 'deposited': r['deposited'], 'title': r['title']})
with open(sys.argv[2], 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
print(len(out), 'contacts in', len({o['pdbid'] for o in out}), 'entries;', dict(collections.Counter(o['triage'] for o in out)))
