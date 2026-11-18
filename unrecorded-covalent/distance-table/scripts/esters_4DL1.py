"""The sixteen heme ester links of entry 4DL1: exact distance from the deposited coordinates, whether a covale
record joins the two atoms, and whether the file's own close-contact table lists them.
Usage: python3 -I esters_4DL1.py 4DL1.cif data/esters_4DL1.csv   (4DL1.cif from https://files.rcsb.org/download/4DL1.cif)"""
import csv, sys, gemmi
st = gemmi.read_structure(sys.argv[1]); st.setup_entities(); m = st[0]
blk = gemmi.cif.read(sys.argv[1]).sole_block()
rec = set(); close = {}
for r in blk.find('_struct_conn.', ['conn_type_id', 'ptnr1_auth_asym_id', 'ptnr1_auth_seq_id', 'ptnr1_label_atom_id', 'ptnr2_auth_asym_id', 'ptnr2_auth_seq_id', 'ptnr2_label_atom_id']):
    if r[0] == 'covale': rec.add(frozenset([(r[1], r[2], r[3]), (r[4], r[5], r[6])]))
for i, r in enumerate(blk.find('_pdbx_validate_close_contact.', ['auth_asym_id_1', 'auth_seq_id_1', 'auth_atom_id_1', 'auth_asym_id_2', 'auth_seq_id_2', 'auth_atom_id_2', 'dist'])):
    close[frozenset([(r[0], r[1], r[2]), (r[3], r[4], r[5])])] = (i + 1, r[6])
ns = gemmi.NeighborSearch(m, st.cell, 4).populate()
pts = []
for ch in m:
    for res in ch:
        if res.name in ('ASP', 'GLU'):
            for at in res:
                if at.name in ('OD1', 'OD2', 'OE1', 'OE2'):
                    for mk in ns.find_atoms(at.pos, '\0', radius=2.1):
                        c = mk.to_cra(m)
                        if c.residue.name == 'HEM' and c.atom.element.name == 'C':
                            d = c.atom.pos.dist(at.pos)
                            key = frozenset([(ch.name, str(res.seqid.num), at.name), (c.chain.name, str(c.residue.seqid.num), c.atom.name)])
                            row, cd = close.get(key, ('', ''))
                            pts.append((d, key in rec, f'{res.name.title()}{res.seqid.num} {ch.name}', key in close, row, cd))
pts.sort()
with open(sys.argv[2], 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['distance_A', 'recorded', 'residue_chain', 'listed_as_close_contact', 'close_contact_row', 'close_contact_dist'])
    for d, r, lab, cl, row, cd in pts: w.writerow([f'{d:.5f}', r, lab, cl, row, cd])
print(len(pts), 'esters,', sum(p[1] for p in pts), 'recorded,', sum(p[3] for p in pts), 'listed as close contacts; close-contact rows in file:', len(close))
