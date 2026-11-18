"""ASAP-Polaris-OpenADMET ligand-pose challenge structures (Zenodo 15582067, structures.tar.gz; OESpruce-prepared PDBs).
For each structure: minimum distance from any Cys145 SG to a heavy atom of the ligand (resname UNK), and whether the
file carries any LINK/CONECT record. Protein label (SARS-CoV-2 vs MERS) is inferred by sequence at 145+-3."""
import gemmi, glob, os, pandas as pd
D='/workspace/mpro_covalent_census'
rows=[]
for p in sorted(glob.glob(f'{D}/data/asap/structures/*.pdb')):
    txt=open(p).read(); st=gemmi.read_structure(p); st.remove_hydrogens(); m=st[0]
    sg=[a for ch in m for r in ch if r.name=='CYS' and r.seqid.num==145 for a in r if a.name=='SG']
    lig=[a for ch in m for r in ch if r.name=='UNK' for a in r]
    d=min((x.pos.dist(y.pos),y.name,y.element.name) for x in sg for y in lig) if sg and lig else (None,None,None)
    seq=gemmi.one_letter_code([r.name for r in m[0]])
    rows.append({'file':os.path.basename(p),'set':os.path.basename(p).split('_')[0],'n_cys145_sg':len(sg),'n_lig_heavy':len(lig),'min_sg_lig':None if d[0] is None else round(d[0],3),'lig_atom':d[1],'lig_el':d[2],
                 'has_LINK':'\nLINK' in txt,'has_CONECT':'\nCONECT' in txt,'chainA_len':len(seq),'is_sars2_like':'SGFRKMAF' in seq})
X=pd.DataFrame(rows); X.to_csv(f'{D}/out/asap_polaris_cys145.tsv',sep='\t',index=False)
print(len(X),X.set.value_counts().to_dict(),'SARS2-like',X.is_sars2_like.sum())
print('min SG-ligand distance bins',pd.cut(X.min_sg_lig,[0,1.6,2.0,2.4,3.0,99]).value_counts().sort_index().to_dict())
print(X[X.min_sg_lig<2.4].to_string())
