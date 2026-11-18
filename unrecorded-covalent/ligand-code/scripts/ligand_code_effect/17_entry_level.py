"""Step 3 (entry level). Unit = (code, entry) with >=1 nucleophile bond < 2.0 A (recorded / missed / mixed).
Miss = entry has >=1 unrecorded bond <2.0 A for that code (missed or mixed).
Birth = status of the code's birth batch (all entries sharing its first deposition date in the archive).
(a) logistic miss ~ birth + site + year (+ nucleophile), clustered by code; attribution of misses.
Writes out/step3/*.tsv"""
import json, pandas as pd, numpy as np, statsmodels.formula.api as smf, statsmodels.api as sm
O='/workspace/ligand_code_effect'; A=f'{O}/archive'; OUT=f'{O}/out/step3'
TL=pd.read_csv(f'{A}/code_timelines.tsv',sep='\t',low_memory=False)
CD=pd.read_csv(f'{O}/out/archive_codes.tsv',sep='\t',index_col=0)
B=pd.read_csv(f'{O}/out/archive_bonds_lt2.tsv',sep='\t',low_memory=False)
TL=TL.join(CD[['cat','birth_batch','birth_dep']],on='lig_comp')
TL['dep']=pd.to_datetime(TL.dep,errors='coerce'); TL['birth_dep']=pd.to_datetime(TL.birth_dep,errors='coerce')
TL['bonded']=TL.status.isin(['recorded','missed','mixed'])
TL['miss']=TL.status.isin(['missed','mixed']).astype(int)
TL['phase']=np.where(TL.dep>TL.birth_dep,'later','birth_batch')
TL['year']=TL.dep.dt.year; TL['yearc']=TL.year-2015
TL['site']=TL.ps.fillna('NA').replace({'BNL':'RCSB','NDB':'RCSB'}); TL['site']=TL.site.where(TL.site.isin(['RCSB','PDBE','PDBJ','PDBC']),'other')
# nucleophile of the entry's closest bond
nb=B.sort_values('dist').groupby(['lig_comp','entry']).nuc.first()
TL=TL.join(nb.rename('nuc'),on=['lig_comp','entry'])
E=TL[TL.bonded].copy(); E.to_csv(f'{OUT}/entry_level_bonded.tsv',sep='\t',index=False)
lines=[]; P=lambda *a: lines.append(' '.join(map(str,a)))
P('bonded (code,entry) units:',len(E),' misses:',int(E.miss.sum()),' codes:',E.lig_comp.nunique(),' entries:',E.entry.nunique())
for sub,x0 in [('all',E),('drug-like',E[E.cat=='drug-like']),('non-bilin/heme',E[~E.cat.isin(['bilin','heme'])])]:
    P(f'\n### {sub}: units {len(x0)}, misses {int(x0.miss.sum())}')
    t=x0.groupby(['phase','birth_batch']).agg(units=('miss','size'),misses=('miss','sum'),codes=('lig_comp','nunique')); t['rate']=(t.misses/t.units).round(3); P(t.to_string())
    x=x0[(x0.phase=='later')&x0.birth_batch.isin(['recorded','missed','unbonded','linked_other','mixed'])].copy()
    t2=x.groupby(['birth_batch',np.where(x.site=='PDBJ','PDBj','other')]).agg(units=('miss','size'),misses=('miss','sum')); t2['rate']=(t2.misses/t2.units).round(3); P(t2.to_string())
    x['code']=x.lig_comp.astype('category').cat.codes
    try:
        m=smf.glm('miss ~ C(birth_batch, Treatment("recorded")) + C(site, Treatment("RCSB")) + yearc',x,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':x.code})
        ci=m.conf_int(); R=[]
        for k in m.params.index: R.append(dict(subset=sub,term=k,OR=round(np.exp(m.params[k]),2),lo=round(np.exp(ci.loc[k,0]),2),hi=round(np.exp(ci.loc[k,1]),2),p=f'{m.pvalues[k]:.2g}',n=int(m.nobs),codes=x.code.nunique()))
        R=pd.DataFrame(R); P(R.to_string()); R.to_csv(f'{OUT}/logit_{sub.replace("/","_")}.tsv',sep='\t',index=False)
        # attribution
        tot=int(x0.miss.sum()); AT=[]
        AT.append(dict(bucket='birth batch itself missed (first use)',misses=int(x0[x0.phase=='birth_batch'].miss.sum())))
        for g in ['unbonded','linked_other','missed','recorded','mixed']:
            s=x[x.birth_batch==g]
            if not len(s): continue
            s2=s.copy(); s2['birth_batch']='recorded'; ex=float(m.predict(s2).sum())
            AT.append(dict(bucket=f'later use, code born {g}',units=len(s),misses=int(s.miss.sum()),expected_if_born_recorded=round(ex,1),excess=round(s.miss.sum()-ex,1)))
        AT=pd.DataFrame(AT); AT['share_of_all_misses']=(AT.misses/tot).round(3); AT['excess_share']=(AT.excess/tot).round(3)
        P('total misses',tot); P(AT.to_string()); AT.to_csv(f'{OUT}/attribution_{sub.replace("/","_")}.tsv',sep='\t',index=False)
        # PDBj-free sensitivity
        y=x[x.site!='PDBJ'].copy(); y['code']=y.lig_comp.astype('category').cat.codes
        m2=smf.glm('miss ~ C(birth_batch, Treatment("recorded")) + C(site, Treatment("RCSB")) + yearc',y,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':y.code})
        k='C(birth_batch, Treatment("recorded"))[T.unbonded]'; ci=m2.conf_int()
        P('excluding PDBj: born-unbonded OR',round(np.exp(m2.params[k]),2),round(np.exp(ci.loc[k,0]),2),round(np.exp(ci.loc[k,1]),2),'p',f'{m2.pvalues[k]:.2g}','n',int(m2.nobs))
    except Exception as ex: P('model failed',ex)
# strict "inheritance": codes born unbonded; first bonded entry status -> later bonded entries
D=E[E.cat=='drug-like'].sort_values(['lig_comp','dep','rel','entry'])
rows=[]
for c,g in D.groupby('lig_comp'):
    if CD.loc[c,'birth_batch']!='unbonded' or len(g)<2: continue
    f=g.iloc[0]; rest=g.iloc[1:]
    rows.append(dict(lig_comp=c,first_bonded=f.entry,first_status=f.status,first_site=f.site,n_later=len(rest),later_miss=int(rest.miss.sum()),later_pdbj=int((rest.site=='PDBJ').sum())))
R=pd.DataFrame(rows); R.to_csv(f'{OUT}/born_unbonded_first_bonded.tsv',sep='\t',index=False)
if len(R):
    R['fmiss']=R.first_status!='recorded'
    P('\n### drug-like codes born unbonded with >=2 bonded entries:',len(R))
    P(R.groupby('fmiss').agg(codes=('lig_comp','size'),later=('n_later','sum'),later_miss=('later_miss','sum')).to_string())
open(f'{OUT}/17_stdout.txt','w').write('\n'.join(lines)); print('\n'.join(lines))
