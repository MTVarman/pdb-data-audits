"""Every serine-to-boron contact of vaborbactam (component 4D6) in the given entries, alternate conformers
counted separately: distance from the deposited coordinates, whether a covale record joins the two atoms in
that conformer, and whether the file's own close-contact table lists them.
Usage: python3 -I vaborbactam_entries.py out.csv 4XUX.cif 4XUZ.cif ...   (files from https://files.rcsb.org/download/)"""
import csv, sys, gemmi
out = []
def val(x): return '' if x in ('?', '.') else gemmi.cif.as_string(x)
for path in sys.argv[2:]:
    blk = gemmi.cif.read(path).sole_block()
    st = gemmi.make_structure_from_block(blk); st.setup_entities(); m = st[0]
    rec = [[val(x) for x in r] for r in blk.find('_struct_conn.', ['conn_type_id', 'ptnr1_auth_asym_id', 'ptnr1_auth_seq_id', 'ptnr1_label_atom_id', '?pdbx_ptnr1_label_alt_id',
                                                                    'ptnr2_auth_asym_id', 'ptnr2_auth_seq_id', 'ptnr2_label_atom_id', '?pdbx_ptnr2_label_alt_id'])]
    rec = [r[1:] for r in rec if r[0] == 'covale']
    close = [[val(x) for x in r] for r in blk.find('_pdbx_validate_close_contact.', ['auth_asym_id_1', 'auth_seq_id_1', 'auth_atom_id_1', '?label_alt_id_1',
                                                                                    'auth_asym_id_2', 'auth_seq_id_2', 'auth_atom_id_2', '?label_alt_id_2'])]
    def has(rows_, a, b):
        """a, b: (chain, residue number, atom name, altloc).  A row with no alt id stands for every conformer."""
        for r in rows_:
            for x, y in ((r[:4], r[4:]), (r[4:], r[:4])):
                if x[:3] == list(a[:3]) and y[:3] == list(b[:3]) and x[3] in ('', a[3]) and y[3] in ('', b[3]): return True
        return False
    meta = [gemmi.cif.as_string(blk.find_value(t) or '?') for t in ('_exptl.method', '_refine.ls_d_res_high', '_pdbx_database_status.recvd_initial_deposition_date', '_pdbx_database_status.process_site')]
    ns = gemmi.NeighborSearch(m, st.cell, 4).populate()
    for ch in m:
        for res in ch:
            if res.name != '4D6': continue
            for at in res:
                if at.element.name != 'B': continue
                for mk in ns.find_atoms(at.pos, '\0', radius=2.6):
                    if mk.image_idx != 0: continue
                    c = mk.to_cra(m)
                    if c.residue.name != 'SER' or c.atom.name != 'OG': continue
                    la, pa = at.altloc.strip('\0'), c.atom.altloc.strip('\0')
                    if la and pa and la != pa: continue
                    a = (c.chain.name, str(c.residue.seqid.num), 'OG', pa); b = (ch.name, str(res.seqid.num), at.name, la)
                    out.append({'distance_A': f'{c.atom.pos.dist(at.pos):.4f}', 'pdbid': st.name.upper(), 'protein_atom': f'SER{c.residue.seqid.num} {c.chain.name} OG',
                                'ligand_atom': f'4D6 {ch.name}{res.seqid.num} {at.name}', 'ligand_conformer': la, 'recorded_covale': has(rec, a, b),
                                'listed_as_close_contact': has(close, a, b), 'ligand_atom_occupancy': round(at.occ, 2), 'protein_atom_occupancy': round(c.atom.occ, 2),
                                'method': meta[0], 'resolution_A': meta[1], 'deposited': meta[2], 'process_site': meta[3]})
out.sort(key=lambda r: float(r['distance_A']))
with open(sys.argv[1], 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
print(len(out), 'contacts in', len({r['pdbid'] for r in out}), 'entries')
