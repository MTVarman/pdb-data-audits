"""Step 3 (b) same-drug cross-code splits (same lab+day, same lab, same site+year) and (c) codes that get 'fixed'.
Same drug = identical skeleton key (scripts/03) OR Morgan-2 Tanimoto >= 0.8 on CCD SMILES (leaving atoms kept)."""
import json, pandas as pd, numpy as np, itertools, collections
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog('rdApp.*')
O='/workspace/ligand_code_effect'; OUT=f'{O}/out/step3'
C=json.load(open(f'{O}/out/ccd_parsed_archive.json'))
E=pd.read_csv(f'{OUT}/entry_level_bonded.tsv',sep='\t',low_memory=False)
E['dep']=pd.to_datetime(E.dep)
E=E[E.cat=='drug-like'].copy()
B=pd.read_csv(f'{O}/out/archive_bonds_lt2.tsv',sep='\t',low_memory=False)
FP={}
for c in E.lig_comp.unique():
    s=(C.get(c) or {}).get('smiles'); m=Chem.MolFromSmiles(s) if s else None
    if m is None and s: m=Chem.MolFromSmiles(s,sanitize=False)
    if m is not None:
        try: m.UpdatePropertyCache(strict=False); FP[c]=AllChem.GetMorganFingerprintAsBitVect(m,2,2048)
        except Exception: pass
def sim(a,b):
    if C.get(a,{}).get('skeleton') and C.get(a,{}).get('skeleton')==C.get(b,{}).get('skeleton'): return 1.0,'skeleton'
    if a in FP and b in FP: return DataStructs.TanimotoSimilarity(FP[a],FP[b]),'tanimoto'
    return 0,'none'
# per code summary within a grouping
def codestat(g):
    return dict(n=len(g),miss=int(g.miss.sum()),entries=','.join(g.entry),dist=','.join(f'{d:.2f}' for d in g.dmin),sites=','.join(sorted(set(g.site))))
rows=[]
for key,cols in [('lab+day',['au_last','dep']),('lab',['au_last']),('site+year',['site','year'])]:
    for k,g in E[E.au_last.notna()].groupby(cols):
        per={c:x for c,x in g.groupby('lig_comp')}
        rec=[c for c,x in per.items() if x.miss.sum()==0]; mis=[c for c,x in per.items() if x.miss.mean()==1]
        for a in rec:
            for b in mis:
                s,how=sim(a,b)
                if s<0.8: continue
                xa,xb=per[a],per[b]
                same_entry=bool(set(xa.entry)&set(xb.entry))
                rows.append(dict(level=key,group=' | '.join(map(str,k if isinstance(k,tuple) else (k,))),code_rec=a,code_miss=b,similarity=round(s,2),how=how,same_entry=same_entry,
                    rec_n=len(xa),rec_entries=','.join(xa.entry),rec_dist=','.join(f'{d:.2f}' for d in xa.dmin),rec_sites=','.join(sorted(set(xa.site))),
                    miss_n=len(xb),miss_entries=','.join(xb.entry),miss_dist=','.join(f'{d:.2f}' for d in xb.dmin),miss_sites=','.join(sorted(set(xb.site))),
                    rec_birth=xa.birth_batch.iloc[0],miss_birth=xb.birth_batch.iloc[0],
                    rec_ccd_created=C.get(a,{}).get('pdbx_initial_date'),miss_ccd_created=C.get(b,{}).get('pdbx_initial_date'),
                    rec_def=C.get(a,{}).get('pdbx_model_coordinates_db_code'),miss_def=C.get(b,{}).get('pdbx_model_coordinates_db_code'),
                    name_rec=(C.get(a,{}).get('name') or '')[:45],name_miss=(C.get(b,{}).get('name') or '')[:45]))
S=pd.DataFrame(rows).sort_values(['level','similarity'],ascending=[True,False]); S.to_csv(f'{OUT}/same_drug_code_splits.tsv',sep='\t',index=False)
L=[]; P=lambda *a: L.append(' '.join(map(str,a)))
for lev in ['lab+day','lab','site+year']:
    s=S[S.level==lev]; P(f'(b) {lev}: discordant same-drug code pairs {len(s)} in {s.group.nunique()} groups; same entry {int(s.same_entry.sum())}; recorded side {int(s.rec_n.sum())} entries, missed side {int(s.miss_n.sum())} entries; missed-side site PDBJ {int(s.miss_sites.str.contains("PDBJ").sum())}; missed-side code born: {s.miss_birth.value_counts().to_dict()}')
# baseline: how many lab+day groups have >=2 same-drug codes at all (concordant or not)
cnt=collections.Counter()
for k,g in E[E.au_last.notna()].groupby(['au_last','dep']):
    cs=list(g.lig_comp.unique())
    for a,b in itertools.combinations(cs,2):
        s,_=sim(a,b)
        if s>=0.8:
            ma,mb=g[g.lig_comp==a].miss.mean(),g[g.lig_comp==b].miss.mean()
            cnt['pairs']+=1; cnt['concordant_both_rec']+= (ma==0 and mb==0); cnt['concordant_both_miss']+=(ma==1 and mb==1); cnt['discordant_full']+=((ma==0 and mb==1) or (ma==1 and mb==0))
P('(b) all same-drug code pairs within a lab+day group:',dict(cnt))
# ---- (c) fixes ----
H=[]
for c,g in E.sort_values(['dep','rel','entry']).groupby('lig_comp'):
    if len(g)<3: continue
    m=list(g.miss); 
    if m[0]==0: continue
    k=m.index(0) if 0 in m else len(m)
    if k==len(m): H.append(dict(lig_comp=c,type='never recorded',n=len(g),pre_misses=k,sites=','.join(sorted(set(g.site))))); continue
    post=g.iloc[k:]; rr=1-post.miss.mean()
    typ='fixed' if (rr>=0.8 and (post.miss==0).sum()>=2) else 'flip-flop'
    pre=g.iloc[:k]; fx=g.iloc[k]
    au=C.get(c,{}).get('audit') or []
    lo=str(pre.dep.iloc[-1].date()); hi=str(fx.dep.date())
    ev=[f'{a[0]}@{a[1]}' for a in au if a[1] and lo<=a[1]<=hi]
    H.append(dict(lig_comp=c,type=typ,n=len(g),pre_misses=k,post_n=len(post),post_recorded=int((post.miss==0).sum()),first_entry=g.entry.iloc[0],first_dep=str(g.dep.iloc[0].date()),
        last_miss_before=pre.entry.iloc[-1],fix_entry=fx.entry,fix_dep=hi,lag_days=(fx.dep-pre.dep.iloc[0]).days,pre_sites=','.join(sorted(set(pre.site))),fix_site=fx.site,
        same_lab=fx.au_last in set(pre.au_last),pre_dist=','.join(f'{d:.2f}' for d in pre.dmin),fix_dist=f'{fx.dmin:.2f}',birth=g.birth_batch.iloc[0],ccd_events_between=';'.join(ev),
        ccd_created=C.get(c,{}).get('pdbx_initial_date'),ccd_modified=C.get(c,{}).get('pdbx_modified_date'),name=(C.get(c,{}).get('name') or '')[:45]))
H=pd.DataFrame(H); H.to_csv(f'{OUT}/code_fixes.tsv',sep='\t',index=False)
allc=E.groupby('lig_comp').size(); P(f'(c) drug-like codes with >=3 bonded entries: {int((allc>=3).sum())}; starting with a miss: {len(H)}'); P(H.type.value_counts().to_dict())
F=H[H.type=='fixed']
if len(F):
    P('fixed: pre-fix misses',int(F.pre_misses.sum()),'post entries',int(F.post_n.sum()),'post recorded',int(F.post_recorded.sum()))
    P('fix year',F.fix_dep.str[:4].value_counts().sort_index().to_dict())
    P('pre-fix sites include PDBJ',int(F.pre_sites.str.contains('PDBJ').sum()),'| fix site == PDBJ',int((F.fix_site=='PDBJ').sum()),'| site changed at fix',int(sum(f not in p.split(',') for f,p in zip(F.fix_site,F.pre_sites))),'/',len(F))
    P('same lab at fix',int(F.same_lab.sum()),'/',len(F),'| CCD audit event between last miss and fix',int((F.ccd_events_between!='').sum()))
    P('median lag days',F.lag_days.median())
    pd.set_option('display.width',250); P(F[['lig_comp','name','n','pre_misses','post_n','post_recorded','first_entry','first_dep','fix_entry','fix_dep','pre_sites','fix_site','same_lab','pre_dist','fix_dist','ccd_events_between']].to_string())
open(f'{OUT}/18_stdout.txt','w').write('\n'.join(L)); print('\n'.join(L))
