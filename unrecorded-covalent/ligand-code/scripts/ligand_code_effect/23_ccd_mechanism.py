"""Step 4 (3): CCD-record features of carried-forward codes vs controls. Uses current CCD files (downloaded 16 Nov 2026 AKST)."""
import json, pandas as pd, numpy as np, sys
from scipy.stats import fisher_exact
sys.argv=['x']; O='/workspace/ligand_code_effect'; OUT=f'{O}/out/step4'
C=json.load(open(f'{O}/out/ccd_parsed_archive.json'))
src=open(f'{O}/scripts/04_features.py').read(); ns={}
exec(src[src.index('HAL='):src.index('rows=[]')].replace("C=json.load","#"),{'C':C},ns) if False else None
HAL={'Cl','Br','I'}
def feat(c,a):
    r=C.get(c) or {}; nb=r.get('nb',{}); el=r.get('el',{}); lv=set(r.get('leaving_atoms',[]))
    if a not in nb: return dict(atom_in_ccd=False)
    heavy=[(n,o) for n,o in nb[a] if el.get(n)!='H']; orders=[o for n,o in heavy]
    lv_nb=any(n in lv for n,o in nb[a])
    if 'TRIP' in orders: k='nitrile/alkyne (pre-reaction)'
    elif any(n in lv and el.get(n)!='H' for n,o in heavy) or any(el.get(n) in HAL for n,o in heavy): k='leaving group attached (pre-reaction)'
    elif a in r.get('arom',{}): k='aromatic'
    elif any(o=='DOUB' for n,o in heavy): k='sp2 C=X (pre-reaction or acyl)'
    else: k='saturated (post-reaction-like)'
    pcm=bool(r.get('pcm_rows')) or r.get('pdbx_pcm')=='Y'
    pcm_atom=any((p[1]==a) for p in (r.get('pcm_rows') or []))
    return dict(atom_in_ccd=True,link_ready=lv_nb,any_leaving=bool(lv),attach_class=k,pcm=pcm,pcm_on_attacked_atom=pcm_atom)
S=pd.read_csv(f'{OUT}/strict_units.tsv',sep='\t',low_memory=False); S['dep']=pd.to_datetime(S.dep)
rows=[]
for c,g in S.sort_values(['dep','rel','entry']).groupby('lig_comp'):
    f0=g.iloc[0]; later=g[g.dep>f0.dep]
    if not len(later): continue
    cls=('carried (first missed, later missed)' if later.miss.sum()>0 else 'escaped (first missed, later all recorded)') if f0.miss==1 else ('first recorded, later missed' if later.miss.sum()>0 else 'first recorded, later all recorded')
    r=C.get(c) or {}; fe=feat(c,f0.lig_atom)
    au=r.get('audit') or []
    mods_after_first=[f'{a[0]}@{a[1]}' for a in au if a[1] and a[1]>str(f0.dep.date()) and a[0] not in ('Initial release','Create component')]
    lm=later[later.miss==1]
    rows.append(dict(lig_comp=c,cls=cls,first_entry=f0.entry,first_dep=str(f0.dep.date()),first_site=f0.site,first_lab=f0.lab,attacked_atom=f0.lig_atom,nuc=f0.nuc_atom,
        n_later=len(later),later_miss=int(later.miss.sum()),later_labs=later.lab.nunique(),created=r.get('pdbx_initial_date'),modified=r.get('pdbx_modified_date'),create_site=r.get('create_site'),
        def_entry=r.get('pdbx_model_coordinates_db_code'),mods_after_first=';'.join(mods_after_first),
        modified_after_last_miss=bool(len(lm)) and bool(r.get('pdbx_modified_date')) and r.get('pdbx_modified_date')>str(lm.dep.max().date()),
        pcm_dates=';'.join(r.get('pcm_modify_dates') or []),name=(r.get('name') or '')[:40],**fe))
D=pd.DataFrame(rows); D.to_csv(f'{OUT}/ccd_features_by_code.tsv',sep='\t',index=False)
L=[];P=lambda *a:L.append(' '.join(map(str,a)))
P('codes with >=1 later strict use:',len(D)); P(D.cls.value_counts().to_string())
for k in ['link_ready','any_leaving','pcm','pcm_on_attacked_atom']:
    P(f'\n{k}:'); t=pd.crosstab(D.cls,D[k].fillna(False)); P(t.to_string())
P('\nattach class:'); P(pd.crosstab(D.cls,D.attach_class).to_string())
fm=D[D.cls.str.startswith('carried')|D.cls.str.startswith('escaped')]
P('\nfirst-missed codes (carried + escaped):'); pd.set_option('display.width',260)
P(fm[['lig_comp','cls','first_entry','first_dep','first_site','first_lab','n_later','later_miss','later_labs','attacked_atom','attach_class','link_ready','pcm','created','modified','mods_after_first','name']].sort_values('cls').to_string(index=False))
# carried vs first-recorded-all-recorded: feature contrasts
a=D[D.cls.str.startswith('carried')]; b=D[D.cls=='first recorded, later all recorded']
for k in ['link_ready','pcm']:
    x=[[int(a[k].fillna(False).sum()),int((~a[k].fillna(False).astype(bool)).sum())],[int(b[k].fillna(False).sum()),int((~b[k].fillna(False).astype(bool)).sum())]]
    P(f'{k}: carried {x[0][0]}/{len(a)} vs always-recorded {x[1][0]}/{len(b)}  Fisher p={fisher_exact(x)[1]:.2g}')
# unit-level: later units, CCD link_ready / pcm vs miss (circularity caveat)
S2=S[S.dep>S.groupby('lig_comp').dep.transform('min')].copy()
F=pd.DataFrame([dict(lig_comp=c,**feat(c,a)) for c,a in S2[['lig_comp','lig_atom']].drop_duplicates().values])
S2=S2.merge(F,on=['lig_comp'],how='left').drop_duplicates(['lig_comp','entry'])
for k in ['link_ready','pcm']:
    t=S2.groupby(S2[k].fillna(False)).miss.agg(['sum','size']); P(f'unit level later uses by {k}:',t.to_dict())
open(f'{OUT}/23_stdout.txt','w').write('\n'.join(L)); print('\n'.join(L))
