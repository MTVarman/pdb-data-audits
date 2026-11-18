"""Third, independent re-scan of a random sample of entries, compared row by row with the scan table.
Usage: python3 -I rescan_sample.py sample_contacts.csv cif_dir n seed   (first argument uncompressed; downloads the entries it needs)"""
import csv, sys, os, random, collections, urllib.request, gemmi
src, cdir, n, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
NUC = {('CYS','SG'),('SER','OG'),('THR','OG1'),('TYR','OH'),('LYS','NZ'),('HIS','ND1'),('HIS','NE2'),
       ('ASP','OD1'),('ASP','OD2'),('GLU','OE1'),('GLU','OE2'),('SEC','SE')}
LIGEL = {'C','N','O','S','P','B','Se'}
rows = list(csv.DictReader(open(src)))
by = collections.defaultdict(list)
for r in rows: by[r['pdbid']].append(r)
ids = sorted(by); random.Random(seed).shuffle(ids); ids = ids[:n]
os.makedirs(cdir, exist_ok=True)
tot = collections.Counter(); bad = []
for pid in ids:
    p = os.path.join(cdir, pid + '.cif')
    if not os.path.exists(p):
        with urllib.request.urlopen(f'https://files.rcsb.org/download/{pid}.cif', timeout=120) as r: open(p, 'wb').write(r.read())
    doc = gemmi.cif.read(p); b = doc.sole_block()
    cov = set()
    t = b.find('_struct_conn.', ['conn_type_id','ptnr1_auth_asym_id','ptnr1_auth_comp_id','ptnr1_auth_seq_id','ptnr1_label_atom_id',
                                 'ptnr2_auth_asym_id','ptnr2_auth_comp_id','ptnr2_auth_seq_id','ptnr2_label_atom_id'])
    for r in t:
        v = [gemmi.cif.as_string(r[i]) for i in range(9)]
        if v[0].lower().startswith('covale'):
            cov.add(frozenset([(v[1], v[2], v[3], v[4]), (v[5], v[6], v[7], v[8])]))
    st = gemmi.read_structure(p); st.setup_entities()
    m = st[0]
    ns = gemmi.NeighborSearch(m, st.cell, 5).populate()
    mine = {}
    for ch in m:
        for res in ch:
            if res.entity_type != gemmi.EntityType.Polymer: continue
            for at in res:
                if (res.name, at.name) not in NUC: continue
                best = {}
                for mk in ns.find_atoms(at.pos, '\0', radius=2.6):
                    if mk.image_idx != 0: continue
                    c = mk.to_cra(m)
                    lr = c.residue
                    if lr.entity_type not in (gemmi.EntityType.NonPolymer, gemmi.EntityType.Branched): continue
                    if lr.name in ('HOH','DOD','WAT','H2O'): continue
                    if c.atom.element.name not in LIGEL: continue
                    if not any(a.element.name == 'C' for a in lr): continue
                    d = c.atom.pos.dist(at.pos)
                    if d >= 2.6: continue
                    k = (c.chain.name, lr.name, str(lr.seqid.num) + lr.seqid.icode.strip(), at.altloc, c.atom.altloc)
                    kk = (c.chain.name, lr.name, str(lr.seqid.num) + lr.seqid.icode.strip())
                    rec = frozenset([(ch.name, res.name, str(res.seqid.num), at.name), (c.chain.name, lr.name, str(lr.seqid.num), c.atom.name)]) in cov
                    if kk not in best or d < best[kk][0]: best[kk] = (d, c.atom.name, rec)
                for kk, (d, an, rec) in best.items():
                    mine[(f'{ch.name}:{res.name}{res.seqid.num}{res.seqid.icode.strip()}:{at.name}', f'{kk[1]}_{kk[0]}_{kk[2]}')] = (d, an, rec)
    theirs = {(r['protein_atom'], r['ligand']): (float(r['dist_A']), r['ligand_atom'], r['recorded_covale'] == 'True') for r in by[pid]}
    for k, (d, an, rec) in theirs.items():
        tot['table rows'] += 1
        if k not in mine: tot['row not reproduced'] += 1; bad.append((pid, k, 'absent')); continue
        d2, an2, rec2 = mine[k]
        if abs(d - d2) > 0.0015 or an != an2: tot['distance or atom differs'] += 1; bad.append((pid, k, d, an, round(d2, 3), an2))
        else: tot['distance and atom agree'] += 1
        if rec == rec2: tot['recorded flag agrees'] += 1
        else: tot['recorded flag differs'] += 1; bad.append((pid, k, 'rec', rec, rec2))
    for k in mine:
        if k not in theirs: tot['rescan row absent from table'] += 1; bad.append((pid, k, 'extra', round(mine[k][0], 3)))
print(len(ids), 'entries re-scanned;', dict(tot))
for x in bad[:25]: print('  ', x)
