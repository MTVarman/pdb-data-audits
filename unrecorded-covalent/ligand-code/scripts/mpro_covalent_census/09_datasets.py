"""Join every main-set Cys145 contact (recorded + unrecorded) to datasets already on the box.
PLINDER 2026-09 annotation table (entry + ligand label_asym + CCD; same key as plinder_covalent_audit/02_join.py),
Synthyra/PLINDER legacy split mirror (third-party), Runs N' Poses v5, BioLiP2 weekly (entry + CCD + ligand chain + seq),
LP-PDBBind / PDBbind v2020 (entry), DiffDock time split, MISATO, PoseBusters v1/v2 (entry + CCD)."""
import pandas as pd, numpy as np
D='/workspace/mpro_covalent_census'; DL='/workspace/dataset_covalent_labels/data'
M=pd.read_csv(f'{D}/out/contacts_main_cys145_C_1.6-2.4.tsv',sep='\t',keep_default_na=False,low_memory=False)
M['recorded']=M.recorded.astype(str)=='True'; M['pdb']=M.entry.str.lower()
ents=sorted(M.pdb.unique())
cols=['entry_pdb_id','system_id','ligand_ccd_code','ligand_unique_ccd_code','ligand_asym_id','ligand_is_covalent','ligand_covalent_linkages','ligand_member_asym_ids','system_pass_validation_criteria']
P=pd.read_parquet('/workspace/dsets/plinder/annotation_table.parquet',columns=cols,filters=[('entry_pdb_id','in',ents)])
P['members']=P.ligand_member_asym_ids.apply(lambda x:set(x) if x is not None and len(x) else set())
P['links']=P.ligand_covalent_linkages.apply(lambda x:list(x) if x is not None else [])
Pg={k:g for k,g in P.groupby('entry_pdb_id')}
S=pd.concat([pd.read_parquet(f'/workspace/plinder_covalent_audit/data/synthyra_{s}.parquet',columns=['system_id','split','entry_pdb_id','ligand_ccd_code','ligand_instance_chain','ligand_is_covalent']) for s in ['train','valid','test']]).drop_duplicates()
S['asym']=S.ligand_instance_chain.str.split('.',n=1).str[1]; Sg={k:g for k,g in S.groupby('entry_pdb_id')}
R=pd.read_csv('/workspace/dsets/rnp/annotations.csv',low_memory=False,usecols=['system_id','entry_pdb_id','ligand_instance_chain','ligand_ccd_code'])
R['asym']=R.ligand_instance_chain.str.split('.',n=1).str[1]; Rg={k:g for k,g in R.groupby('entry_pdb_id')}
B=pd.read_csv(f'{DL}/biolip/biolip_ligand_rows.tsv.gz',sep='\t',dtype=str); B=B[B.entry.isin(set(M.entry))]
Bk=set(zip(B.entry,B.lig_comp,B.lig_chain,B.lig_seq)); Be=set(B.entry)
LP=pd.read_csv(f'{DL}/lists/lp_pdbbind.tsv',sep='\t',dtype=str); LPd=LP.set_index('entry')
DD=pd.read_csv(f'{DL}/lists/diffdock_timesplit.tsv',sep='\t',dtype=str).set_index('entry')
MI=pd.read_csv(f'{DL}/lists/misato_md.tsv',sep='\t',dtype=str).set_index('entry')
PB1=pd.read_csv(f'{DL}/lists/posebusters_v1_428.tsv',sep='\t',dtype=str); PB2=pd.read_csv(f'{DL}/lists/posebusters_v2_308.tsv',sep='\t',dtype=str)
j=lambda xs:';'.join(sorted(set(map(str,xs))))
out=[]
for _,r in M.iterrows():
    e,a,comp=r.pdb,r.lig_label_asym,r.lig_comp; o={}
    g=Pg.get(e); o['plinder_entry']=g is not None
    if g is not None:
        m=g[g.apply(lambda x:(x.ligand_asym_id==a or a in x.members),axis=1)]
        mc=m[(m.ligand_ccd_code==comp)|m.ligand_unique_ccd_code.fillna('').str.contains(comp,regex=False)]; m=mc if len(mc) else m
    else: m=pd.DataFrame()
    o['plinder_ligand']=len(m)>0; o['plinder_is_covalent']=j(m.ligand_is_covalent) if len(m) else ''
    o['plinder_n_systems']=m.system_id.nunique() if len(m) else 0
    o['plinder_pass_validation']=j(m.system_pass_validation_criteria) if len(m) else ''
    sg=Sg.get(e); sm=sg[(sg.asym==a)&(sg.ligand_ccd_code==comp)] if sg is not None else pd.DataFrame()
    o['legacy_split']=j(sm.split) if len(sm) else ''; o['legacy_is_covalent']=j(sm.ligand_is_covalent) if len(sm) else ''
    rg=Rg.get(e); o['rnp_entry']=rg is not None; o['rnp_ligand']=bool(rg is not None and ((rg.asym==a)&(rg.ligand_ccd_code==comp)).any())
    o['biolip_entry']=r.entry in Be; o['biolip_ligand_site']=(r.entry,comp,r.lig_chain,str(r.lig_seq)) in Bk
    o['lp_pdbbind']=r.entry in LPd.index
    o['lp_split']=LPd.new_split.get(r.entry,'') if r.entry in LPd.index else ''; o['lp_covalent']=LPd.covalent.get(r.entry,'') if r.entry in LPd.index else ''
    o['pdbbind_v2020_category']=LPd.category.get(r.entry,'') if r.entry in LPd.index else ''
    o['diffdock_split']=DD.diffdock_split.get(r.entry,'') if r.entry in DD.index else ''
    o['misato_split']=MI.misato_split.get(r.entry,'') if r.entry in MI.index else ''
    o['posebusters_v1']=((PB1.entry==r.entry)&(PB1.lig_comp==comp)).any(); o['posebusters_v2']=((PB2.entry==r.entry)&(PB2.lig_comp==comp)).any()
    out.append(o)
X=pd.concat([M.reset_index(drop=True),pd.DataFrame(out)],axis=1)
X.to_csv(f'{D}/out/contacts_main_datasets.tsv',sep='\t',index=False)
for b in ['<2.0','2.0-2.4']:
    for rec in [False,True]:
        x=X[(X.band==b)&(X.recorded==rec)]
        print(f'\n== {b} recorded={rec}: contacts {len(x)} entries {x.entry.nunique()}')
        print(' PLINDER entry',x.plinder_entry.sum(),'ligand',x.plinder_ligand.sum(),'is_covalent',x[x.plinder_ligand].plinder_is_covalent.value_counts().to_dict())
        print(' legacy split mirror',x.legacy_split.value_counts().to_dict(),'| legacy is_cov',x[x.legacy_split!=''].legacy_is_covalent.value_counts().to_dict())
        print(' RnP entry',x.rnp_entry.sum(),'ligand',x.rnp_ligand.sum(),x[x.rnp_entry].entry.unique().tolist())
        print(' BioLiP entry',x.biolip_entry.sum(),'ligand site',x.biolip_ligand_site.sum(),'entries',x[x.biolip_ligand_site].entry.nunique())
        print(' LP/PDBbind v2020',x.lp_pdbbind.sum(),'DiffDock',(x.diffdock_split!='').sum(),'MISATO',(x.misato_split!='').sum(),'PB v1',x.posebusters_v1.sum(),'v2',x.posebusters_v2.sum())
