"""Does the wwPDB clash list ever contain a boron atom?

Usage: python3 -I boron.py val_dir contacts_all.csv out.csv
For every entry that holds an in-band unrecorded Ser/Thr-boronate contact, read the
cached validation report and count (a) clashes whose atom is a boron atom of the
ligand component, (b) clashes on any other atom of the same ligand residue.
Also gets clash status for every B-pair contact at any distance.
"""
import collections, csv, glob, gzip, json, os, sys
import xml.etree.ElementTree as ET

BORON_NAMES = set()


def main():
    vd, src, dst = sys.argv[1], sys.argv[2], sys.argv[3]
    rows = [r for r in csv.DictReader(open(src)) if r['lig_el'] == 'B' and r['occupancy_class'] == 'full']
    print('B-pair contacts, full occupancy, any distance:', len(rows),
          'entries', len({r['pdbid'] for r in rows}))
    by = collections.defaultdict(list)
    for r in rows:
        by[r['pdbid']].append(r)
    out, lig_clash_total, b_clash_total, n_rep = [], 0, 0, 0
    for pdbid, rs in sorted(by.items()):
        e = pdbid.lower()
        p = os.path.join(vd, f'{e}_validation.xml.gz')
        if not os.path.exists(p):
            continue
        try:
            subs = ET.fromstring(gzip.decompress(open(p, 'rb').read())).findall('ModelledSubgroup')
        except Exception:
            continue
        n_rep += 1
        # every clash in the whole report, by atom name, for the ligand residues of interest
        ligcomps = {r['ligand'].rsplit('_', 2)[0] for r in rs}
        batoms = {r['ligand_atom'] for r in rs}
        for s in subs:
            if s.get('resname') in ligcomps:
                for cl in s.findall('clash'):
                    lig_clash_total += 1
                    if cl.get('atom') in batoms:
                        b_clash_total += 1
        for r in rs:
            lname, lch, lnum = r['ligand'].rsplit('_', 2)
            lnum = ''.join(c for c in lnum if c.isdigit() or c == '-')
            lc = [s for s in subs if s.get('model', '1') == '1' and s.get('chain') == lch
                  and s.get('resnum') == lnum and s.get('resname') == lname]
            nb = sum(1 for s in lc for cl in s.findall('clash') if cl.get('atom') == r['ligand_atom'])
            nall = sum(1 for s in lc for cl in s.findall('clash'))
            out.append({'pdbid': pdbid, 'residue': r['residue'], 'pair': r['pair'],
                        'protein_atom': r['protein_atom'], 'ligand': r['ligand'],
                        'ligand_atom': r['ligand_atom'], 'dist_A': r['dist_A'],
                        'probe_line_A': r['probe_line_A'], 'recorded_covale': r['recorded_covale'],
                        'clashes_on_this_boron_atom': nb,
                        'clashes_anywhere_on_this_ligand_residue': nall})
    with open(dst, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print('reports read', n_rep)
    print('clash records on any atom of these boronate ligand residues:', lig_clash_total)
    print('of those, on the boron atom itself:', b_clash_total)
    print('rows with a clash on the boron atom:', sum(1 for x in out if x['clashes_on_this_boron_atom']))
    print('rows with a clash anywhere on the ligand residue:',
          sum(1 for x in out if x['clashes_anywhere_on_this_ligand_residue']))
    d = [x for x in out if x['recorded_covale'] == 'False']
    print('\nunrecorded B-pair contacts by distance, clash on boron:')
    lo = 1.6
    while lo < 2.65:
        g = [x for x in d if lo <= float(x['dist_A']) < lo + 0.1]
        if g:
            print(f'  {lo:.2f}-{lo+0.1:.2f}  n={len(g):4d}  clash-on-B={sum(1 for x in g if x["clashes_on_this_boron_atom"]):3d}')
        lo = round(lo + 0.1, 2)


if __name__ == '__main__':
    main()
