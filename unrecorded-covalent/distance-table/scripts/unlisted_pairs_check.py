"""Unrecorded pairs under 2.195 A that are in neither contact table of their file: why not?
For each, read the two atoms from the deposited coordinates: occupancies, alternate-conformer labels, and any
_struct_conn row of another type that joins them.
Usage: python3 -I unlisted_pairs_check.py archive_residue_pairs.csv.gz boundary_exact_distances.csv cif_dir out.csv
(downloads missing files into cif_dir)"""
import csv, gzip, os, sys, urllib.request, collections
import gemmi
src, bnd, cifdir, dst = sys.argv[1:5]
BX = {}
for b in csv.DictReader(open(bnd, newline='')):
    if b['source'] == 'archive':
        BX[(b['entry'], b['protein_chain'], b['protein_residue'], b['protein_seq'], b['protein_atom'], b['ligand_chain'], b['ligand_comp'], b['ligand_seq'], b['ligand_atom'])] = float(b['dist_exact'])
BOND = ('C-S', 'C-O', 'C-N', 'S-S', 'O-P', 'O-S', 'N-S', 'N-P')
def canon(a, b): return '-'.join(sorted([a.capitalize(), b.capitalize()]))
rows = []
for r in csv.DictReader(gzip.open(src, 'rt', newline='')):
    if r['header_read'] != '1': continue
    rec = r['recorded'] == '1' or r['disulf'] == '1'
    d = BX.get((r['entry'], r['prot_chain'], r['prot_res'], r['prot_seq'], r['prot_atom'], r['lig_chain'], r['lig_comp'], r['lig_seq'], r['lig_atom']), float(r['dist']))
    if canon(r['prot_el'], r['lig_el']) in BOND and not rec and d < 2.195 and r['close_listed'] != '1' and r['symm_listed'] != '1':
        rows.append(r)
out = []
for e in sorted({r['entry'] for r in rows}):
    p = os.path.join(cifdir, f'{e}.cif.gz')
    if not os.path.exists(p):
        req = urllib.request.Request(f'https://files.rcsb.org/download/{e}.cif.gz', headers={'User-Agent': 'pdb-struct-audit/1.0'})
        with urllib.request.urlopen(req, timeout=180) as u, open(p, 'wb') as f: f.write(u.read())
    st = gemmi.read_structure(p); st.setup_entities(); m = st[0]
    ns = gemmi.NeighborSearch(m, st.cell, 4).populate()
    for r in [x for x in rows if x['entry'] == e]:
        best = None
        for ch in m:
            if ch.name != r['prot_chain']: continue
            for res in ch:
                if res.name != r['prot_res'] or str(res.seqid.num) != r['prot_seq']: continue
                for at in res:
                    if at.name != r['prot_atom']: continue
                    for mk in ns.find_atoms(at.pos, '\0', radius=2.7):
                        c = mk.to_cra(m)
                        if c.chain.name != r['lig_chain'] or c.residue.name != r['lig_comp'] or str(c.residue.seqid.num) != r['lig_seq'] or c.atom.name != r['lig_atom']: continue
                        if at.altloc != '\0' and c.atom.altloc != '\0' and at.altloc != c.atom.altloc: continue
                        d = st.cell.find_nearest_pbc_image(at.pos, c.atom.pos, mk.image_idx).dist()
                        if best is None or d < best[0]: best = (d, at, c.atom, mk.image_idx)
        if best is None:
            out.append({'entry': e, 'prot': f"{r['prot_chain']}:{r['prot_res']}{r['prot_seq']}:{r['prot_atom']}", 'lig': f"{r['lig_comp']}_{r['lig_chain']}_{r['lig_seq']}:{r['lig_atom']}",
                        'scan_dist': r['dist'], 'dist': '', 'prot_occ': '', 'lig_occ': '', 'prot_alt': '', 'lig_alt': '', 'symmetry_mate': '', 'conn_types': r['conn_types'], 'reason': 'pair not found'})
            continue
        d, a, b, img = best
        both_partial = a.occ < 1.0 and b.occ < 1.0
        reason = ('both atoms below full occupancy' if both_partial else 'symmetry mate' if img != 0 else
                  'joined by a record of another type' if r['conn_types'] else 'other')
        out.append({'entry': e, 'prot': f"{r['prot_chain']}:{r['prot_res']}{r['prot_seq']}:{r['prot_atom']}", 'lig': f"{r['lig_comp']}_{r['lig_chain']}_{r['lig_seq']}:{r['lig_atom']}",
                    'scan_dist': r['dist'], 'dist': f'{d:.3f}', 'prot_occ': round(a.occ, 2), 'lig_occ': round(b.occ, 2), 'prot_alt': a.altloc.strip('\0'), 'lig_alt': b.altloc.strip('\0'),
                    'symmetry_mate': img != 0, 'conn_types': r['conn_types'], 'method': r['method1'], 'resolution': r['resolution'], 'dep_year': r['dep_year'], 'reason': reason})
cols = ['entry', 'prot', 'lig', 'scan_dist', 'dist', 'prot_occ', 'lig_occ', 'prot_alt', 'lig_alt', 'symmetry_mate', 'conn_types', 'method', 'resolution', 'dep_year', 'reason']
with open(dst, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(out)
print(len(out), 'pairs in', len({r['entry'] for r in out}), 'entries;', dict(collections.Counter(r['reason'] for r in out)))
