"""Step 4 (3b): is PCM annotation a cause or a consequence? Compare PCM date to deposition date of each later unit;
and test pre-reaction attachment chemistry for carried codes."""
import json, pandas as pd, numpy as np
from scipy.stats import fisher_exact
O='/workspace/ligand_code_effect'; OUT=f'{O}/out/step4'
C=json.load(open(f'{O}/out/ccd_parsed_archive.json'))
S=pd.read_csv(f'{OUT}/strict_units.tsv',sep='\t',low_memory=False); S['dep']=pd.to_datetime(S.dep)
D=pd.read_csv(f'{OUT}/ccd_features_by_code.tsv',sep='\t')
L=[];P=lambda *a:L.append(' '.join(map(str,a)))
allpcm=[d for c in C.values() for d in (c.get('pcm_modify_dates') or [])]
P('PCM modify dates across all parsed CCDs, top:',pd.Series(allpcm).value_counts().head(5).to_dict())
def pcm_first(c):
    r=C.get(c) or {}; ds=sorted(r.get('pcm_modify_dates') or [])
    if not (r.get('pcm_rows') or r.get('pdbx_pcm')=='Y'): return None
    return ds[0] if ds else (r.get('pdbx_initial_date') or None)
S['pcm_date']=S.lig_comp.map(pcm_first)
S['pcm_before_dep']=S.pcm_date.notna()&(S.pcm_date<S.dep.dt.strftime('%Y-%m-%d'))
first=S.groupby('lig_comp').dep.transform('min'); Lt=S[S.dep>first]
P('later strict units:',len(Lt))
for lab,x in [('PCM existed before deposition',Lt[Lt.pcm_before_dep]),('no PCM at deposition (added later or never)',Lt[~Lt.pcm_before_dep])]:
    P(f'  {lab}: {int(x.miss.sum())}/{len(x)} missed')
R=Lt[Lt.dep>='2024-10-01']
P('later units deposited from 2024-10-01 (after the 2024-09-27 PCM bulk update):',len(R))
for lab,x in [('code had PCM',R[R.pcm_before_dep]),('code had no PCM',R[~R.pcm_before_dep])]:
    P(f'  {lab}: {int(x.miss.sum())}/{len(x)} missed; entries missed: {",".join(x[x.miss==1].entry)[:200]}')
a=R[R.pcm_before_dep]; b=R[~R.pcm_before_dep]
if len(a) and len(b):
    t=[[int(b.miss.sum()),len(b)-int(b.miss.sum())],[int(a.miss.sum()),len(a)-int(a.miss.sum())]]; o,p=fisher_exact(t); P('  Fisher (no PCM vs PCM): OR',round(o,1) if np.isfinite(o) else o,'p',f'{p:.2g}')
    P('  no-PCM units, by code:',b.groupby('lig_comp').miss.agg(['sum','size']).to_dict('index'))
    P('  no-PCM units at PDBj:',int((b.site=="PDBJ").sum()),' misses outside PDBj',int(b[b.site!="PDBJ"].miss.sum()),'/',int((b.site!="PDBJ").sum()))
# chemistry of attack atom, carried vs always-recorded
pre=lambda s: s.str.contains('pre-reaction')
a=D[D.cls.str.startswith('carried')]; b=D[D.cls=='first recorded, later all recorded']
t=[[int(pre(a.attach_class).sum()),len(a)-int(pre(a.attach_class).sum())],[int(pre(b.attach_class).sum()),len(b)-int(pre(b.attach_class).sum())]]
P('pre-reaction attachment atom in CCD: carried',t[0][0],'/',len(a),'vs always-recorded',t[1][0],'/',len(b),'Fisher p',f'{fisher_exact(t)[1]:.2g}')
fm=D[D.cls.str.startswith(('carried','escaped'))].copy(); fm['pcm']=fm.pcm.astype(str).eq('True')
P('first-missed codes lacking PCM today:',fm[~fm.pcm].lig_comp.tolist(),'| all are carried:',bool((fm[~fm.pcm].cls.str.startswith('carried')).all()))
P('first-missed codes with PCM today: PCM dates',fm[fm.pcm][['lig_comp','pcm_dates']].values.tolist())
open(f'{OUT}/24_stdout.txt','w').write('\n'.join(L)); print('\n'.join(L))
