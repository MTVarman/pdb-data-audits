"""For every unrecorded main-set Cys145 contact, check every archived major version (files-versioned.wwpdb.org)
for a covale row joining that Cys to that ligand residue. Reuses files already fetched by the PLINDER audit."""
import os, re, shutil, urllib.request, pandas as pd, gemmi
from concurrent.futures import ThreadPoolExecutor
D='/workspace/mpro_covalent_census'; V=f'{D}/data/versioned'; OLD='/workspace/plinder_covalent_audit/data/versioned'
BASE='https://files-versioned.wwpdb.org/pdb_versioned/data/entries'
M=pd.read_csv(f'{D}/out/contacts_main_cys145_C_1.6-2.4.tsv',sep='\t',keep_default_na=False,low_memory=False)
M['recorded']=M.recorded.astype(str)=='True'; u=M[~M.recorded].copy()
def listing(e):
    e=e.lower(); url=f'{BASE}/{e[1:3]}/pdb_0000{e}/'
    try: html=urllib.request.urlopen(url,timeout=60).read().decode()
    except Exception: return e,None
    files=sorted(set(re.findall(rf'(pdb_0000{e}_xyz_v(\d+)-(\d+)\.cif\.gz)',html)),key=lambda x:(int(x[1]),int(x[2])))
    out=[]
    for f,ma,mi in files:
        p=f'{V}/{f}'
        if not os.path.exists(p):
            if os.path.exists(f'{OLD}/{f}'): shutil.copy(f'{OLD}/{f}',p)
            else: urllib.request.urlretrieve(url+f,p)
        out.append((f'{ma}.{mi}',p))
    return e,out
with ThreadPoolExecutor(6) as ex: L=dict(ex.map(listing,sorted(u.entry.unique())))
def cov(p):
    b=gemmi.cif.read(p).sole_block()
    t=b.find('_struct_conn.',['conn_type_id','ptnr1_auth_asym_id','ptnr1_auth_comp_id','ptnr1_auth_seq_id','ptnr2_auth_asym_id','ptnr2_auth_comp_id','ptnr2_auth_seq_id'])
    s=set()
    for r in t:
        if not r.str(0).startswith('covale'): continue
        a=(r.str(1),r.str(2),r.str(3)); bb=(r.str(4),r.str(5),r.str(6)); s.add((a,bb)); s.add((bb,a))
    return s
cache={}; rows=[]
for _,r in u.iterrows():
    vers=L.get(r.entry.lower()); cys=(str(r.nuc_chain),'CYS',str(r.nuc_seq)); lig=(str(r.lig_chain),r.lig_comp,str(r.lig_seq)); hist=[]
    for v,p in vers or []:
        if p not in cache: cache[p]=cov(p)
        hist.append(f"{v}:{'covale' if (cys,lig) in cache[p] else 'none'}")
    rows.append({'archived_versions':';'.join(v for v,_ in vers) if vers else 'LISTING_FAILED','covale_by_version':';'.join(hist),'covale_in_any_earlier_version':any(h.endswith(':covale') for h in hist)})
res=pd.concat([u.reset_index(drop=True),pd.DataFrame(rows)],axis=1)
res[['entry','process_site','deposit_year','lig_comp','drug','lig_chain','lig_seq','nuc_chain','nuc_seq','dist','band','archived_versions','covale_by_version','covale_in_any_earlier_version']].to_csv(f'{D}/out/unrecorded_version_history.tsv',sep='\t',index=False)
print(res.groupby('band').covale_in_any_earlier_version.agg(['sum','size'])); print('listing failed',sum(v is None for v in L.values()))
print(res[res.covale_in_any_earlier_version][['entry','process_site','lig_comp','covale_by_version']].to_string())
