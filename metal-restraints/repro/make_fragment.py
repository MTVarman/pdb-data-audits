import gemmi,sys
src,out=sys.argv[1],sys.argv[2]
st=gemmi.read_structure(src); st.setup_entities()
m=st[0]
# find a CLA whose MG is within 2.0 A of a HIS NE2
ns=gemmi.NeighborSearch(m,st.cell,4).populate()
best=None
for ch in m:
    for res in ch:
        if res.name!='CLA': continue
        mg=res.find_atom('MG','*')
        if not mg: continue
        for mk in ns.find_atoms(mg.pos,'\0',radius=2.1):
            cra=mk.to_cra(m)
            if cra.residue.name=='HIS' and cra.atom.name=='NE2':
                best=(ch.name,res,cra.chain.name,cra.residue, mg.pos.dist(cra.atom.pos)); break
        if best: break
    if best: break
print('pair: CLA',best[0],best[1].seqid.num,'<-> HIS',best[2],best[3].seqid.num,'dist %.3f'%best[4])
hch=m[best[2]]; hnum=best[3].seqid.num
ns2=gemmi.Structure(); ns2.cell=gemmi.UnitCell(200,200,200,90,90,90); ns2.spacegroup_hm='P 1'
nm=gemmi.Model('1'); pc=gemmi.Chain('A'); lc=gemmi.Chain('B')
for res in hch:
    if res.name in ('CLA',): continue
    if hnum-3<=res.seqid.num<=hnum+3 and res.het_flag=='A': pc.add_residue(res)
lc.add_residue(best[1])
nm.add_chain(pc); nm.add_chain(lc); ns2.add_model(nm); ns2.setup_entities()
ns2.write_pdb(out)
print('wrote',out,'protein residues',[r.name+str(r.seqid.num) for r in pc],'ligand',best[1].name)
