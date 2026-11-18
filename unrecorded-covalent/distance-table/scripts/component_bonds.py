"""Every PDB entry containing a given ligand component: the distance from its named atom to the nearest protein
nucleophile, and whether _struct_conn records a covale row for that atom pair.
Usage: python3 -I component_bonds.py COMP ELEMENT cif_dir [out.csv]"""
import sys, os, csv, json, urllib.request, gemmi
comp, el, cdir = sys.argv[1], sys.argv[2], sys.argv[3]
q = {"query": {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id", "operator": "exact_match", "value": comp}},
     "return_type": "entry", "request_options": {"return_all_hits": True}}
req = urllib.request.Request('https://search.rcsb.org/rcsbsearch/v2/query', data=json.dumps(q).encode(), headers={'Content-Type': 'application/json'})
ids = sorted(h['identifier'] for h in json.load(urllib.request.urlopen(req, timeout=60))['result_set'])
print(comp, 'entries:', len(ids), ids)
rows = []
for pid in ids:
    p = f'{cdir}/{pid}.cif'
    if not os.path.exists(p) or os.path.getsize(p) < 1000:
        open(p, 'wb').write(urllib.request.urlopen(f'https://files.rcsb.org/download/{pid}.cif', timeout=120).read())
    b = gemmi.cif.read(p).sole_block()
    cov = set()
    for r in b.find('_struct_conn.', ['conn_type_id', 'ptnr1_auth_asym_id', 'ptnr1_auth_comp_id', 'ptnr1_auth_seq_id', 'ptnr1_label_atom_id', 'ptnr2_auth_asym_id', 'ptnr2_auth_comp_id', 'ptnr2_auth_seq_id', 'ptnr2_label_atom_id']):
        v = [gemmi.cif.as_string(r[i]) for i in range(9)]
        if v[0].lower().startswith('covale'): cov.add(frozenset([(v[1], v[2], v[3], v[4]), (v[5], v[6], v[7], v[8])]))
    st = gemmi.read_structure(p); st.setup_entities(); m = st[0]
    ns = gemmi.NeighborSearch(m, st.cell, 4).populate()
    meth = gemmi.cif.as_string(b.find_value('_exptl.method')); res = b.find_value('_refine.ls_d_res_high') or b.find_value('_em_3d_reconstruction.resolution')
    for ch in m:
        for r_ in ch:
            if r_.name != comp: continue
            for a in r_:
                if a.element.name != el: continue
                best = None
                for mk in ns.find_atoms(a.pos, '\0', radius=2.6):
                    c = mk.to_cra(m)
                    if c.residue.entity_type != gemmi.EntityType.Polymer or c.atom.element.name not in ('O', 'N', 'S'): continue
                    if c.atom.name in ('O', 'N', 'OXT'): continue
                    d = c.atom.pos.dist(a.pos)
                    if best is None or d < best[0]: best = (d, c)
                if best:
                    d, c = best
                    rec = frozenset([(c.chain.name, c.residue.name, str(c.residue.seqid.num), c.atom.name), (ch.name, comp, str(r_.seqid.num), a.name)]) in cov
                    rows.append((d, pid, f'{c.residue.name}{c.residue.seqid.num} {c.chain.name} {c.atom.name}', f'{comp} {ch.name}{r_.seqid.num} {a.name}', rec, a.occ, c.atom.occ, meth, res,
                                 gemmi.cif.as_string(b.find_value('_pdbx_database_status.recvd_initial_deposition_date')), gemmi.cif.as_string(b.find_value('_pdbx_database_status.process_site'))))
for r in sorted(rows): print(f"  {r[0]:.4f}  {r[1]}  {r[2]:16s} {r[3]:16s} {'recorded' if r[4] else 'NOT recorded':12s} occ {r[5]:.2f}/{r[6]:.2f}  {r[7][:5]} {r[8]}  dep {r[9]} {r[10]}")

if len(sys.argv) > 4:
    with open(sys.argv[4], 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['distance_A', 'pdbid', 'protein_atom', 'ligand_atom', 'recorded_covale', 'ligand_atom_occupancy', 'protein_atom_occupancy', 'method', 'resolution_A', 'deposited', 'process_site'])
        for r in sorted(rows): w.writerow([f'{r[0]:.4f}', r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10]])
