"""Is each in-band unrecorded contact in the entry's own _pdbx_validate_close_contact table?

Usage: python3 -I closecheck.py scan_dir window_evidence.csv out.csv
"""
import csv, glob, json, re, sys


def main():
    scan_dir, win, dst = sys.argv[1], sys.argv[2], sys.argv[3]
    rows = [r for r in csv.DictReader(open(win))
            if r['in_band'] == 'True' and r['recorded_covale'] == 'False'
            and r['occupancy_class'] == 'full']
    need = {r['pdbid'] for r in rows}
    close = {}
    for p in sorted(glob.glob(f'{scan_dir}/shard_*.jsonl')):
        for line in open(p):
            if '"pdbid":"' not in line:
                continue
            r = json.loads(line)
            if r['pdbid'] in need:
                close[r['pdbid']] = r.get('close') or []
    out = []
    for r in rows:
        pch, pres, patom = r['protein_atom'].split(':')
        m = re.match(r'([A-Z]+)(-?\d+)', pres)
        pname, pnum = m.group(1), m.group(2)
        lname, lch, lnum = r['ligand'].rsplit('_', 2)
        lnum = re.match(r'-?\d+', lnum).group(0)
        hit, hitd = False, ''
        for c in close.get(r['pdbid'], []):
            # [model, ch1, comp1, seq1, atom1, ch2, comp2, seq2, atom2, dist]
            a = (c[1], c[2], c[3], c[4])
            b = (c[5], c[6], c[7], c[8])
            want_p = (pch, pname, pnum, patom)
            want_l = (lch, lname, lnum, r['ligand_atom'])
            if (a == want_p and b == want_l) or (a == want_l and b == want_p):
                hit, hitd = True, c[9]
        out.append({'pdbid': r['pdbid'], 'residue': r['residue'], 'pair': r['pair'],
                    'protein_atom': r['protein_atom'], 'ligand': r['ligand'],
                    'ligand_atom': r['ligand_atom'], 'dist_A': r['dist_A'],
                    'probe_line_A': r['probe_line_A'], 'verdict_clash_listed': r['clash_listed'],
                    'in_close_contact_table': hit, 'close_contact_dist': hitd,
                    'n_close_rows_in_entry': len(close.get(r['pdbid'], []))})
    with open(dst, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    import collections
    print('in-band unrecorded rows', len(out))
    print('in the PDB close-contact table:', sum(1 for x in out if x['in_close_contact_table']))
    for pair in sorted({x['pair'] for x in out}):
        b = [x for x in out if x['pair'] == pair]
        print(f'  {pair:6s} {len(b):4d} rows, in close-contact table '
              f'{sum(1 for x in b if x["in_close_contact_table"]):4d}, clash-listed '
              f'{sum(1 for x in b if x["verdict_clash_listed"] == "True"):4d}')


if __name__ == '__main__':
    main()
