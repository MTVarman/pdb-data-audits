"""Test 'code history': for contacts whose CCD code was defined by an EARLIER, different entry, does the defining
entry lacking a covale row for that code predict a missing record here?  Not circular: the defining entry is a
different, earlier deposition."""
import pandas as pd, numpy as np, statsmodels.formula.api as smf, statsmodels.api as sm
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import StratifiedTable
O='/workspace/ligand_code_effect/out'
s=pd.read_csv(f'{O}/contacts_lt2_model_input.tsv',sep='\t',low_memory=False)
d=pd.read_csv(f'{O}/defining_entry_links.tsv',sep='\t')
s=s.merge(d,on='lig_comp',how='left')
s['reused']=(s.ccd_defining_entry.str.upper()!=s.entry)&s.def_n_covale_to_comp.notna()
s['earlier']=pd.to_datetime(s.def_deposit_date,errors='coerce')<pd.to_datetime(s.deposit_date,errors='coerce')
r=s[s.reused&s.earlier].copy(); r['born_unlinked']=(r.def_n_covale_to_comp==0).astype(int)
L=[]
def fish(x,tag):
    a=((x.born_unlinked==1)&(x.miss==1)).sum(); b=((x.born_unlinked==1)&(x.miss==0)).sum()
    c=((x.born_unlinked==0)&(x.miss==1)).sum(); dd=((x.born_unlinked==0)&(x.miss==0)).sum()
    o,p=fisher_exact([[a,b],[c,dd]]); L.append(dict(subset=tag,miss_born_unlinked=f'{a}/{a+b}',miss_born_linked=f'{c}/{c+dd}',OR=round(o,1),p=f'{p:.2g}',
        codes_born_unlinked=x[x.born_unlinked==1].lig_comp.nunique(),entries_born_unlinked=x[x.born_unlinked==1].entry.nunique()))
fish(r,'all reused codes (defining entry earlier)')
for st,x in r.groupby('site'): fish(x,st)
fish(r[r.site!='PDBJ'],'non-PDBj')
tabs=[np.array([[((x.born_unlinked==1)&(x.miss==1)).sum(),((x.born_unlinked==1)&(x.miss==0)).sum()],[((x.born_unlinked==0)&(x.miss==1)).sum(),((x.born_unlinked==0)&(x.miss==0)).sum()]])+0.0 for st,x in r.groupby('site')]
tabs=[t for t in tabs if t.sum(1).min()>0]; st_=StratifiedTable(tabs)
L.append(dict(subset='Mantel-Haenszel by site',OR=round(st_.oddsratio_pooled,1),p=f'{st_.test_null_odds().pvalue:.2g}'))
T=pd.DataFrame(L); T.to_csv(f'{O}/code_history_tests.tsv',sep='\t',index=False); print(T.to_string())
r['code']=r.lig_comp.astype('category').cat.codes
R=[]
for name,form in [('site+year','miss ~ C(site, Treatment("RCSB")) + yearc'),('site+year+born_unlinked','miss ~ C(site, Treatment("RCSB")) + yearc + born_unlinked'),
                  ('site+year+born_unlinked+unreacted+leaving','miss ~ C(site, Treatment("RCSB")) + yearc + born_unlinked + unreacted + leaving')]:
    m=smf.glm(form,r,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':r.code})
    for k in m.params.index: R.append(dict(model=name,term=k,OR=round(np.exp(m.params[k]),2),ci_lo=round(np.exp(m.conf_int().loc[k,0]),2),ci_hi=round(np.exp(m.conf_int().loc[k,1]),2),p=f'{m.pvalues[k]:.2g}',n=int(m.nobs),aic=round(m.aic,1)))
R=pd.DataFrame(R); R.to_csv(f'{O}/code_history_logistic.tsv',sep='\t',index=False); print(R.to_string())
# explain the misses: share of all misses < 2.0 A by category
s['bucket']=np.select([s.newcode==1,s.reused&s.earlier&(s.def_n_covale_to_comp==0),s.reused&s.earlier&(s.def_n_covale_to_comp>0)],['code new with this entry','reused code, born unlinked','reused code, born linked'],'other/unknown')
b=s.groupby(['bucket','pdbj']).agg(n=('miss','size'),miss=('miss','sum')).reset_index(); print(b.to_string())
b.to_csv(f'{O}/miss_buckets.tsv',sep='\t',index=False)
r.to_csv(f'{O}/reused_contacts.tsv',sep='\t',index=False)
print(r[(r.born_unlinked==1)].groupby(['lig_comp','site','ccd_defining_entry']).agg(n=('miss','size'),miss=('miss','sum'),entries=('entry',lambda x:','.join(sorted(set(x)))[:80])).sort_values('n',ascending=False).to_string())
