"""Step 3: does the status of a code's FIRST bonded use predict later bonded uses? Drug-like, specific (<=10 unbonded uses),
bond-like distances. Adjusted for site and year, clustered by code; also split by same lab as the first bonded use."""
import pandas as pd, numpy as np, statsmodels.formula.api as smf, statsmodels.api as sm
O='/workspace/ligand_code_effect'; OUT=f'{O}/out/step3'
E=pd.read_csv(f'{OUT}/entry_level_bonded.tsv',sep='\t',low_memory=False)
TL=pd.read_csv(f'{O}/archive/code_timelines.tsv',sep='\t',low_memory=False)
B=pd.read_csv(f'{O}/out/archive_bonds_lt2.tsv',sep='\t',low_memory=False)
nun=TL[TL.status=='unbonded'].groupby('lig_comp').size(); E['nun']=E.lig_comp.map(nun).fillna(0)
cb=B.sort_values('dist').groupby(['lig_comp','entry']).first()[['nuc_atom','dist']]; E=E.join(cb,on=['lig_comp','entry'])
E['bondlike']=np.where(E.nuc_atom=='SG',E.dist.between(1.65,1.95),E.dist.between(1.20,1.65))
L=[];P=lambda *a:L.append(' '.join(map(str,a)))
for name,D in [('drug-like',E[E.cat=='drug-like']),('drug-like specific+bond-like',E[(E.cat=='drug-like')&(E.nun<=10)&E.bondlike])]:
    D=D.sort_values(['lig_comp','dep','rel','entry']).copy()
    f=D.groupby('lig_comp').head(1).set_index('lig_comp')
    D['first_miss']=D.lig_comp.map(f.miss); D['first_dep']=D.lig_comp.map(f.dep); D['first_lab']=D.lig_comp.map(f.au_last)
    x=D[D.dep>D.first_dep].copy(); x['same_lab']=x.au_last==x.first_lab; x['code']=x.lig_comp.astype('category').cat.codes
    t=x.groupby(['first_miss','same_lab']).agg(units=('miss','size'),misses=('miss','sum'),codes=('lig_comp','nunique')); t['rate']=(t.misses/t.units).round(3)
    P(f'\n## {name}: later bonded uses {len(x)} (deposited after the first bonded use)'); P(t.to_string())
    m=smf.glm('miss ~ first_miss + same_lab + C(site, Treatment("RCSB")) + yearc',x,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':x.code})
    ci=m.conf_int()
    for k in ['first_miss','same_lab[T.True]','C(site, Treatment("RCSB"))[T.PDBJ]']:
        P(k,'OR',round(np.exp(m.params[k]),2),f'({np.exp(ci.loc[k,0]):.2f}-{np.exp(ci.loc[k,1]):.2f})','p',f'{m.pvalues[k]:.2g}')
    y=x[~x.same_lab].copy(); y['code']=y.lig_comp.astype('category').cat.codes
    m2=smf.glm('miss ~ first_miss + C(site, Treatment("RCSB")) + yearc',y,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':y.code}); ci=m2.conf_int()
    P('different-lab later uses only: first_miss OR',round(np.exp(m2.params['first_miss']),2),f'({np.exp(ci.loc["first_miss",0]):.2f}-{np.exp(ci.loc["first_miss",1]):.2f})','p',f'{m2.pvalues["first_miss"]:.2g}','n',len(y),'misses',int(y.miss.sum()))
    z=x[x.first_miss==1]; P('later uses of codes whose first bonded use missed: PDBj',int((z.site=='PDBJ').sum()),'of',len(z),'; misses outside PDBj',int(z[z.site!='PDBJ'].miss.sum()),'of',int((z.site!='PDBJ').sum()))
    z.to_csv(f'{OUT}/later_uses_after_missed_first_{name.split()[0]}{"_strict" if "specific" in name else ""}.tsv',sep='\t',index=False)
open(f'{OUT}/20_stdout.txt','w').write('\n'.join(L)); print('\n'.join(L))
