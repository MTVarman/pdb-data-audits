#!/usr/bin/env python3
"""Recount the strict cut of unrecorded protein-ligand covalent contacts.

Repository: github.com/MTVarman/pdb-data-audits  ·  Dataset: doi:10.5281/zenodo.23441847
Usage:  python3 verify_strict.py
Sits at the repository root, next to README.md and unrecorded-covalent/. Standard library only.

Predicates on unrecorded-covalent/strict-cut-491/unrecorded_covalent_tight_cut_491_sites.csv:
  site              : distinct site_id (the folder keeps the name 491 so links resolve; six duplicated
                      rows were collapsed, see that folder's README)
  strict cut        : triage class 2, X-ray, resolution <= 2.5 A
  density-supported : both partners RSCC >= 0.8, RSR <= 0.3, occupancy > 0 in the wwPDB validation
                      report. Recomputed here from the six raw columns and compared with the
                      density_class column.
  clash outlier at the flagged pair: clash_class outside {'no_clash_recorded_for_pair',
                      'clash_between_same_residue_and_ligand_other_atoms', 'no_clash_data_in_report'}
"""
import csv, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(HERE, 'unrecorded-covalent', 'strict-cut-491', 'unrecorded_covalent_tight_cut_491_sites.csv')

results = []
def check(label, got, expected):
    ok = got == expected
    results.append(ok)
    print(f'{"PASS" if ok else "FAIL"}  {label}: got {got}, expected {expected}')

def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

def supported(r):
    v = [num(r[k]) for k in ('protein_rscc', 'ligand_rscc', 'protein_rsr', 'ligand_rsr',
                             'protein_occupancy', 'ligand_occupancy')]
    if any(x is None for x in v):
        return False
    prs, lrs, prr, lrr, po, lo = v
    return prs >= 0.8 and lrs >= 0.8 and prr <= 0.3 and lrr <= 0.3 and po > 0 and lo > 0

if not os.path.exists(P):
    print(f'FAIL  missing {P}: verify_strict.py belongs at the root of the repository / unzipped dataset')
    sys.exit(1)

rows = list(csv.DictReader(open(P, newline='', encoding='utf-8')))
sites = {r['site_id']: r for r in rows}
S = list(sites.values())

check('rows', len(rows), 485)
check('unique site_ids', len(sites), 485)
check('entries', len({r['entry'] for r in S}), 289)
check('all class 2, X-ray, resolution <= 2.5 A',
      all(r['triage_class'] == '2' and r['method'] == 'X-ray' and num(r['resolution_A']) is not None
          and num(r['resolution_A']) <= 2.5 for r in S), True)

ds_col = [r for r in S if r['density_class'] == 'density_supported']
ds_raw = [r for r in S if supported(r)]
check('density-supported sites (density_class column)', len(ds_col), 397)
check('density-supported entries', len({r['entry'] for r in ds_col}), 233)
check('density-supported sites recomputed from RSCC / RSR / occupancy', len(ds_raw), 397)
check('recomputed set equals the density_class column',
      {r['site_id'] for r in ds_raw} == {r['site_id'] for r in ds_col}, True)

known = {'no_clash_recorded_for_pair', 'clash_between_same_residue_and_ligand_other_atoms', 'no_clash_data_in_report'}
check('sites reported as a clash outlier at the flagged pair', sum(r['clash_class'] not in known for r in S), 0)

first = {}
for r in S:
    first.setdefault(r['entry'], r['release_date'])
since = sum(d >= '2020-01-01' for d in first.values())
check('entries released since 2020 (README: more than half of 289)', since, 160)

print(f'\n{sum(results)}/{len(results)} checks passed')
sys.exit(0 if all(results) else 1)
