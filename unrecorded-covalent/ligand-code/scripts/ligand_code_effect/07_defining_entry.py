"""For each CCD code, read its defining entry (_chem_comp.pdbx_model_coordinates_db_code) and record whether that
entry has a _struct_conn covale row joining the code to a polymer residue ('code born linked').
Local copies are used where present (cys_covalent_scan/cif, mpro_covalent_census/data/cif); others are downloaded."""
import pandas as pd, os, gemmi, urllib.request, concurrent.futures as cf, json
O='/workspace/ligand_code_effect'
f=pd.read_csv(f'{O}/out/contacts_features.tsv',sep='\t',low_memory=False)
pairs=f[['lig_comp','ccd_defining_entry']].drop_duplicates().dropna()
def path(e):
    for d in ('/workspace/cys_covalent_scan/cif','/workspace/mpro_covalent_census/data/cif',f'{O}/data/cif',f'{O}/data/defcif'):
        p=f'{d}/{e}.cif.gz'
        if os.path.exists(p): return p
    p=f'{O}/data/defcif/{e}.cif.gz'
    try: open(p,'wb').write(urllib.request.urlopen(f'https://files.rcsb.org/download/{e}.cif.gz',timeout=90).read()); return p
    except Exception as ex: return None
def v(x):
    s=gemmi.cif.as_string(x); return None if s in ('?','.','') else s
def job(row):
    c,e=row; e=e.upper(); p=path(e)
    if not p: return c,e,None,None,None,None
    try: b=gemmi.cif.read(p).sole_block()
    except Exception: return c,e,None,None,None,None
    t=b.find('_struct_conn.',['conn_type_id','ptnr1_label_comp_id','ptnr2_label_comp_id'])
    n=sum(1 for r in t if (v(r[0]) or '').startswith('covale') and c in (v(r[1]),v(r[2])))
    present=any(v(x)==c for x in b.find_values('_pdbx_nonpoly_scheme.mon_id')) or any(v(x)==c for x in b.find_values('_chem_comp.id'))
    site=v(b.find_value('_pdbx_database_status.process_site')) if b.find_value('_pdbx_database_status.process_site') else None
    dep=v(b.find_value('_pdbx_database_status.recvd_initial_deposition_date')) if b.find_value('_pdbx_database_status.recvd_initial_deposition_date') else None
    return c,e,n,present,site,dep
with cf.ThreadPoolExecutor(12) as ex: res=list(ex.map(job,pairs.itertuples(index=False,name=None)))
d=pd.DataFrame(res,columns=['lig_comp','defining_entry','def_n_covale_to_comp','def_comp_present','def_process_site','def_deposit_date'])
d.to_csv(f'{O}/out/defining_entry_links.tsv',sep='\t',index=False)
print(len(d), d.def_n_covale_to_comp.isna().sum(), (d.def_n_covale_to_comp>0).sum(), (d.def_n_covale_to_comp==0).sum())
