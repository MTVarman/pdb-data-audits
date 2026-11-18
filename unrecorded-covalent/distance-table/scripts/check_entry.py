import sys, gemmi
path, resname, resnum = sys.argv[1], sys.argv[2], int(sys.argv[3])
cutoff = float(sys.argv[4]) if len(sys.argv) > 4 else 2.6
doc = gemmi.cif.read(path)
b = doc.sole_block()
print("title:", b.find_value('_struct.title'))
print("resolution:", b.find_value('_refine.ls_d_res_high'), "method:", b.find_value('_exptl.method'))
print("process_site:", b.find_value('_pdbx_database_status.process_site'),
      "deposit:", b.find_value('_pdbx_database_status.recvd_initial_deposition_date'))
print("nonpoly:", [(r[0], r[1]) for r in b.find('_pdbx_entity_nonpoly.', ['comp_id', 'name'])])
sc = b.find('_struct_conn.', ['id','conn_type_id','ptnr1_auth_asym_id','ptnr1_auth_comp_id','ptnr1_auth_seq_id','ptnr1_label_atom_id',
                              'ptnr2_auth_asym_id','ptnr2_auth_comp_id','ptnr2_auth_seq_id','ptnr2_label_atom_id','pdbx_dist_value'])
print("struct_conn rows:", len(sc))
for r in sc:
    row = [r[i] for i in range(11)]
    if row[1] == 'covale' and (resname in (row[3], row[7])) and (str(resnum) in (row[4], row[8])):
        print("  CONN", row)
print("all covale rows not involving sugars:")
for r in sc:
    row = [r[i] for i in range(11)]
    if row[1] == 'covale' and not {'NAG','MAN','BMA','FUC'} & {row[3], row[7]}:
        print("  ", row)
st = gemmi.read_structure(path)
st.setup_entities()
ns = gemmi.NeighborSearch(st[0], st.cell, 4).populate()
for ch in st[0]:
    for res in ch:
        if res.name == resname and res.seqid.num == resnum:
            for at in res:
                if at.element.name in ('O','N','S') and at.name not in ('O','N'):
                    for m in ns.find_atoms(at.pos, '\0', radius=cutoff):
                        cra = m.to_cra(st[0])
                        if cra.residue is res: continue
                        if cra.residue.het_flag != 'H': continue
                        d = cra.atom.pos.dist(at.pos)
                        print(f"  chain {ch.name} {res.name}{res.seqid.num} {at.name} alt[{at.altloc!r}] occ {at.occ:.2f}  --  "
                              f"{cra.residue.name} {cra.chain.name}{cra.residue.seqid.num} {cra.atom.name} ({cra.atom.element.name}) "
                              f"alt[{cra.atom.altloc!r}] occ {cra.atom.occ:.2f}  d={d:.3f}")
