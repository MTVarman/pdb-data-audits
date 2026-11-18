"""For every covalent-complex entry, pull annotation fields that could explain missing covale rows:
_pdbx_validate_close_contact rows naming the Cys145 SG and the ligand atom (did validation flag the 'bond' as a clash?),
all _pdbx_database_status fields, _pdbx_audit_revision_category entries for struct_conn, chem_comp of the ligand
(_chem_comp.pdbx_type / mon_nstd), and _pdbx_entry_details fields."""
import gemmi, pandas as pd, json
D='/workspace/mpro_covalent_census'
M=pd.read_csv(f'{D}/out/contacts_main_cys145_C_1.6-2.4.tsv',sep='\t',keep_default_na=False,low_memory=False)
M['recorded']=M.recorded.astype(str)=='True'
rows=[]
for e,g in M.groupby('entry'):
    b=gemmi.cif.read(f'{D}/data/cif/{e}.cif.gz').sole_block()
    cc=b.get_mmcif_category('_pdbx_validate_close_contact.')
    close=set()
    if cc:
        for i in range(len(cc['id'])):
            a=(cc['auth_asym_id_1'][i],cc['auth_seq_id_1'][i],cc['auth_comp_id_1'][i],cc['auth_atom_id_1'][i]); bb=(cc['auth_asym_id_2'][i],cc['auth_seq_id_2'][i],cc['auth_comp_id_2'][i],cc['auth_atom_id_2'][i])
            close.add((a,bb)); close.add((bb,a))
    st=b.get_mmcif_category('_pdbx_database_status.')
    stat={k:(v[0] if v else None) for k,v in (st or {}).items()}
    rc=b.get_mmcif_category('_pdbx_audit_revision_category.')
    sc_rev=sorted({str(x) for x in (rc or {}).get('category',[]) if x and 'struct_conn' in str(x)})
    ri=b.get_mmcif_category('_pdbx_audit_revision_details.')
    for _,r in g.iterrows():
        k=((r.nuc_chain,str(r.nuc_seq),'CYS','SG'),(r.lig_chain,str(r.lig_seq),r.lig_comp,r.lig_atom))
        rows.append({'entry':e,'nuc_chain':r.nuc_chain,'lig_comp':r.lig_comp,'dist':r.dist,'band':r.band,'recorded':r.recorded,'process_site':r.process_site,
                     'flagged_close_contact':k in close,'n_close_contacts_entry':len(cc['id']) if cc else 0,
                     'status_fields':json.dumps({k2:v for k2,v in stat.items() if k2 not in ('entry_id',)}),
                     'struct_conn_revised':';'.join(sc_rev)})
X=pd.DataFrame(rows); X.to_csv(f'{D}/out/annotation_fields.tsv',sep='\t',index=False)
print(pd.crosstab([X.band,X.recorded],X.flagged_close_contact))
print(pd.crosstab([X.process_site,X.recorded],X.flagged_close_contact))
print('struct_conn revised among unrecorded',(X[~X.recorded].struct_conn_revised!='').sum(),'recorded',(X[X.recorded].struct_conn_revised!='').sum())
keys=set()
for s in X.status_fields: keys|=set(json.loads(s))
print(sorted(keys))
