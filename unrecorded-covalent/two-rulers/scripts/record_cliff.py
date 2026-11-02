"""Same bond, same entry: is it recorded in every copy?

Usage: python3 -I record_cliff.py index_all.jsonl out.csv
Groups: same ligand component + same ligand atom to Cys SG within one entry; both
partners full occupancy, no crossing altlocs, 1.6-2.5 A. A group is mixed when at least
one copy has a _struct_conn covale record and at least one copy has none (and no record
on that residue pair either).
"""
import json, collections, csv, sys
COF = set('PEB PUB PCB CYC PVB BLA BLV LBV BLR RBV HEC HEM HEB HAS HEO FAD FMN FDA FAE PLP PMP BTN LPA RET PLM MYR FAR GER PNS GSH ACE FOR'.split())
rows = []
for line in open(sys.argv[1]):
    if '"Covale"' not in line:
        continue
    r = json.loads(line)
    grp = collections.defaultdict(list)
    for c in r.get('contacts', []):
        if ':CYS' in c['prot'] and c['prot'].endswith(':SG') and c['lig_el'] == 'C' and c['prot_occ'] >= 0.99 and c['lig_occ'] >= 0.99 \
           and not (c['prot_altloc'] and c['lig_altloc'] and c['prot_altloc'] != c['lig_altloc']):
            grp[(c['lig'].split('_')[0], c['lig_atom'])].append(c)
    for (comp, atom), cs in grp.items():
        recd = [c for c in cs if 'Covale' in c['conn_pair'] and 1.6 <= c['dist'] < 2.5]
        un = [c for c in cs if not c['conn_pair'] and not any('Covale' in x for x in c['conn_residues']) and 1.6 <= c['dist'] < 2.5]
        if recd and un:
            rows.append({'pdbid': r['pdbid'], 'ligand': comp, 'ligand_atom': atom, 'class': 'cofactor' if comp in COF else 'ligand',
                         'method': r.get('method'), 'resolution': r.get('resolution'),
                         'recorded_copies': ' '.join(f"{c['prot'].split(':')[0]}:{c['dist']:.3f}" for c in sorted(recd, key=lambda c: c['dist'])),
                         'unrecorded_copies': ' '.join(f"{c['prot'].split(':')[0]}:{c['dist']:.3f}" for c in sorted(un, key=lambda c: c['dist'])),
                         'max_recorded_A': max(c['dist'] for c in recd), 'min_unrecorded_A': min(c['dist'] for c in un), 'title': r.get('title', '')})
rows.sort(key=lambda x: (x['class'] != 'ligand', x['pdbid']))
with open(sys.argv[2], 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(len(rows), 'groups,', len({x['pdbid'] for x in rows}), 'entries')
