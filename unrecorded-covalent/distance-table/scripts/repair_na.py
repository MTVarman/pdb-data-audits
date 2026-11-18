"""Restore the identifier value "NA" (a chain named NA, a ligand atom named NA) where a table export left the cell empty.
The scan's own output has no empty identifier cells; every empty cell in these three columns was the string NA.
Usage: python3 -I repair_na.py in.csv.gz out.csv.gz"""
import csv, gzip, sys
n = {'prot_chain': 0, 'lig_chain': 0, 'lig_atom': 0}
with gzip.open(sys.argv[1], 'rt', newline='') as f, gzip.open(sys.argv[2], 'wt', newline='') as g:
    r = csv.DictReader(f); w = csv.DictWriter(g, fieldnames=r.fieldnames); w.writeheader()
    for row in r:
        for k in n:
            if row[k] == '': row[k] = 'NA'; n[k] += 1
        for k in ('prot_res', 'prot_seq', 'prot_atom', 'prot_el', 'lig_comp', 'lig_seq', 'lig_el'):
            assert row[k] != '', (k, row)
        w.writerow(row)
print('restored', n)
