"""Within-entry record split for non-cysteine nucleophiles (second ruler of the published two-rulers audit: other residues).

Usage: python3 -I cliff.py contacts_all.csv out.csv
A group is: one entry, one ligand component, one ligand atom name, one nucleophile
residue type + atom name, full occupancy, no altloc crossing. A group is mixed when
at least one copy carries a _struct_conn covale record for that atom pair and at
least one carries none, with no record on that residue pair either.
"""
import collections, csv, sys


def main():
    src, dst = sys.argv[1], sys.argv[2]
    rows = [r for r in csv.DictReader(open(src)) if r['occupancy_class'] == 'full']
    g = collections.defaultdict(list)
    for r in rows:
        if r['residue_already_bonded'] == 'True':
            continue
        k = (r['pdbid'], r['ligand_comp'], r['ligand_atom'], r['residue'],
             r['protein_atom'].split(':')[2])
        g[k].append(r)
    out = []
    for k, v in g.items():
        if len(v) < 2:
            continue
        rec = [x for x in v if x['recorded_covale'] == 'True']
        un = [x for x in v if x['recorded_covale'] == 'False']
        if not rec or not un:
            continue
        rd = sorted(float(x['dist_A']) for x in rec)
        ud = sorted(float(x['dist_A']) for x in un)
        out.append({'pdbid': k[0], 'residue': k[3], 'nucleophile_atom': k[4],
                    'ligand_comp': k[1], 'ligand_atom': k[2], 'pair': v[0]['pair'],
                    'probe_line_A': v[0]['probe_line_A'],
                    'n_copies': len(v), 'n_recorded': len(rec), 'n_unrecorded': len(un),
                    'recorded_min': rd[0], 'recorded_max': rd[-1],
                    'unrecorded_min': ud[0], 'unrecorded_max': ud[-1],
                    'split_clean': rd[-1] < ud[0],
                    'resolution': v[0]['resolution'], 'title': v[0]['title'][:90],
                    'recorded_chains': ' '.join(x['protein_atom'].split(':')[0] for x in rec),
                    'unrecorded_chains': ' '.join(x['protein_atom'].split(':')[0] for x in un)})
    out.sort(key=lambda x: (x['residue'] == 'CYS', x['residue'], x['pdbid']))
    with open(dst, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    nc = [x for x in out if x['residue'] != 'CYS']
    print('mixed groups (some copies recorded, some not):', len(out),
          'in', len({x['pdbid'] for x in out}), 'entries')
    print('  non-Cys:', len(nc), 'in', len({x['pdbid'] for x in nc}), 'entries')
    print('  non-Cys with a clean split (every recorded copy shorter than every unrecorded copy):',
          sum(1 for x in nc if x['split_clean']))
    print('  by residue:', dict(collections.Counter(x['residue'] for x in nc)))
    print('  by pair:', dict(collections.Counter(x['pair'] for x in nc)))
    print('\n  non-Cys groups with a clean split:')
    for x in nc:
        if x['split_clean']:
            print(f"    {x['pdbid']} {x['residue']}{'':2s} {x['nucleophile_atom']:4s} {x['ligand_comp']:7s} "
                  f"{x['ligand_atom']:5s} {x['pair']:4s} line={x['probe_line_A']} "
                  f"recorded {x['recorded_min']:.3f}-{x['recorded_max']:.3f} ({x['n_recorded']})  "
                  f"NOT recorded {x['unrecorded_min']:.3f}-{x['unrecorded_max']:.3f} ({x['n_unrecorded']})  "
                  f"| {x['title'][:50]}")


if __name__ == '__main__':
    main()
