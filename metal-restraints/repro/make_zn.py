import gemmi,sys
src,out=sys.argv[1],sys.argv[2]
st=gemmi.read_structure(src); st.setup_entities(); m=st[0]
ns=gemmi.NeighborSearch(m,st.cell,4).populate()
best=None
for ch in m:
    for res in ch:
        if res.name!='ZN': continue
        zn=res[0]; lig=[]
        for mk in ns.find_atoms(zn.pos,'\0',radius=2.7):
            cra=mk.to_cra(m)
            if cra.residue is res: continue
            lig.append((cra.chain.name,cra.residue.name,cra.residue.seqid.num,cra.atom.name,round(zn.pos.dist(cra.atom.pos),3),cra.residue))
        his=[l for l in lig if l[1]=='HIS' and l[3] in('NE2','ND1')]
        rec=[l for l in lig if l[3] in ('SG','NE2','ND1')]
        if len(his)>=2 and len(rec)<4 and all(abs(l[4]-2.30)<0.02 for l in his):
            best=(ch.name,res,lig); break
    if best: break
print('ZN',best[0],best[1].seqid.num,'ligand atoms:',[(l[1],l[2],l[3],l[4]) for l in best[2]])
ns2=gemmi.Structure(); ns2.cell=gemmi.UnitCell(200,200,200,90,90,90); ns2.spacegroup_hm='P 1'
nm=gemmi.Model('1'); pc=gemmi.Chain('A'); lc=gemmi.Chain('B'); seen=set()
for l in sorted(best[2],key=lambda l:l[2]):
    if (l[0],l[2]) in seen or l[1]=='HOH': continue
    seen.add((l[0],l[2])); pc.add_residue(l[5])
lc.add_residue(best[1]); nm.add_chain(pc); nm.add_chain(lc); ns2.add_model(nm); ns2.setup_entities(); ns2.write_pdb(out)
print('wrote',out,[r.name+str(r.seqid.num) for r in pc])
