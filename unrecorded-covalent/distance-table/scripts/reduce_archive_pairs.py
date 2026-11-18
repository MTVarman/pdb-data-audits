"""Reduce the archive-wide atom-pair table to one row per protein residue / ligand residue pair:
the closest side-chain-nucleophile-to-ligand atom pair.  Non-polymer organic ligands only.
Usage: python3 -I reduce_archive_pairs.py data/archive_atom_pairs.csv.gz data/archive_residue_pairs.csv.gz"""
import csv, gzip, sys, collections
NUC = {('CYS','SG'),('SER','OG'),('THR','OG1'),('TYR','OH'),('LYS','NZ'),('HIS','ND1'),('HIS','NE2'),
       ('ASP','OD1'),('ASP','OD2'),('GLU','OE1'),('GLU','OE2'),('SEC','SE')}
best = {}
cls = collections.Counter()
for r in csv.DictReader(gzip.open(sys.argv[1], 'rt', newline='')):
    cls[r['lig_class']] += 1
    if (r['prot_res'], r['prot_atom']) not in NUC: continue
    if r['lig_class'] != 'non-polymer': continue
    k = (r['entry'], r['prot_chain'], r['prot_res'], r['prot_seq'], r['lig_chain'], r['lig_comp'], r['lig_seq'])
    d = float(r['dist'])
    if k not in best or d < float(best[k]['dist']): best[k] = r
print('ligand classes in source table:', dict(cls))
cols = ['entry','prot_chain','prot_res','prot_seq','prot_atom','prot_el','lig_chain','lig_comp','lig_seq','lig_atom','lig_el','dist','recorded','site','dep_year','rel_year','resolution','method1']
with (gzip.open(sys.argv[2], 'wt', newline='') if sys.argv[2].endswith('.gz') else open(sys.argv[2], 'w', newline='')) as f:
    w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore'); w.writeheader()
    for k in sorted(best): w.writerow(best[k])
print('residue pairs written:', len(best), 'entries:', len({k[0] for k in best}))
