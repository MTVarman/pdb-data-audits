#!/usr/bin/env python3
"""Anchor-bond scan. For each entry, pair every protein nucleophile atom with the
nearest heavy atom of a target component within 2.3 A and report whether the
entry carries a covale record for that chain, residue and component.

usage: scan.py <family> <entry-list.txt> <cif-dir> <out.csv>
families: hemec phycobilin biliverdin pyp retinal
"""
import sys, os, csv, gzip, urllib.request
import gemmi
FAM = {'hemec': (['HEC'], 'CYS', 'SG'), 'phycobilin': (['CYC', 'PEB', 'PUB'], 'CYS', 'SG'),
       'biliverdin': (['3Q8', 'EL5', 'BLA', 'LBV', 'LBW'], 'CYS', 'SG'),
       'pyp': (['HC4'], 'CYS', 'SG'), 'retinal': (['RET'], 'LYS', 'NZ')}
def fetch(e, d):
    p = os.path.join(d, e + '.cif.gz')
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        urllib.request.urlretrieve(f'https://files.rcsb.org/download/{e}.cif.gz', p)
    return p
def scan(e, path, comps, resname, atom):
    b = gemmi.cif.read(path).sole_block()
    st = gemmi.make_structure_from_block(b); st.setup_entities(); st.remove_hydrogens()
    rec = set()
    for r in b.find('_struct_conn.', ['conn_type_id', 'ptnr1_auth_asym_id', 'ptnr1_auth_comp_id', 'ptnr1_auth_seq_id',
                                      'ptnr2_auth_asym_id', 'ptnr2_auth_comp_id', 'ptnr2_auth_seq_id']):
        if not r[0].lower().startswith('covale'): continue
        for (c1, n1, s1), (c2, n2, s2) in (((r[1], r[2], r[3]), (r[4], r[5], r[6])), ((r[4], r[5], r[6]), (r[1], r[2], r[3]))):
            if n1 == resname and n2 in comps: rec.add((c1, s1, n2))
    m = st[0]; ns = gemmi.NeighborSearch(m, st.cell, 3).populate(); out = []
    for ch in m:
        for res in ch:
            if res.name != resname: continue
            best = None                       # one contact per residue: nearest pair over alternate conformations
            for a in res:
                if a.name != atom or a.occ <= 0: continue
                for mk in ns.find_atoms(a.pos, '\0', radius=2.3):
                    c = mk.to_cra(m)
                    if c.residue.name not in comps or c.atom.occ <= 0: continue
                    if a.has_altloc() and c.atom.has_altloc() and a.altloc != c.atom.altloc: continue
                    d = a.pos.dist(c.atom.pos) if (mk.image_idx == 0 or not st.cell.is_crystal()) else \
                        a.pos.dist(st.cell.find_nearest_pbc_position(a.pos, c.atom.pos, mk.image_idx))
                    if best is None or d < best[0]: best = (d, c.residue.name)
            if best:
                    out.append([best[1], e, ch.name, res.seqid.num,
                                'recorded' if (ch.name, str(res.seqid.num), best[1]) in rec else 'NOT recorded', round(best[0], 3)])
    return out
if __name__ == '__main__':
    fam, lst, cifdir, outp = sys.argv[1:5]
    comps, resname, atom = FAM[fam]; rows = []
    for e in [x.strip().upper() for x in open(lst) if x.strip()]:
        rows += scan(e, fetch(e, cifdir), comps, resname, atom)
    with open(outp, 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['comp', 'entry', 'chain', 'res', 'record', 'dist_A']); w.writerows(rows)
    print(fam, 'entries', len({r[1] for r in rows}), 'contacts', len(rows), 'unrecorded', sum(r[4] != 'recorded' for r in rows))
