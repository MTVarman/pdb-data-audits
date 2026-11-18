"""Do CCD features predict a missing covale record for ligand-C / Cys-SG contacts < 2.0 A, beyond site and year?"""
import pandas as pd, numpy as np, statsmodels.formula.api as smf, statsmodels.api as sm
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import StratifiedTable
O='/workspace/ligand_code_effect/out'
f=pd.read_csv(f'{O}/contacts_features.tsv',sep='\t',low_memory=False)
s=f[f.lt2].copy(); s['miss']=(~s.recorded).astype(int)
s['site']=s.process_site.replace({'BNL':'RCSB'}).fillna('NA')
s['pdbj']=(s.site=='PDBJ').astype(int)
s['year']=pd.to_numeric(s.deposit_year,errors='coerce')
s['yearc']=s.year-2020
s['unreacted']=(s.attach_state=='unreacted').astype(int)
s['reacted']=(s.attach_state=='reacted').astype(int)
s['leaving']=s.ccd_has_leaving.astype(bool).astype(int)
s['leaving_at_attach']=s.attach_has_leaving_nb.astype(bool).astype(int)
s['newcode']=s.code_defined_by_this_entry.astype(bool).astype(int)
s['new_at_dep']=s.code_new_at_dep.fillna(False).astype(bool).astype(int)
s['ccd_year']=pd.to_datetime(s.ccd_initial_date).dt.year
s['pcm_blank']=(s.ccd_pcm_flag!='Y').astype(int)
s=s[s.year.notna()]
L=[]
def fish(col,lab,sub=s,tag='all'):
    t=pd.crosstab(sub[col],sub.miss).reindex(index=[1,0],columns=[1,0]).fillna(0).astype(int)
    a,b=t.loc[1,1],t.loc[1,0]; c,d=t.loc[0,1],t.loc[0,0]
    orr,p=fisher_exact([[a,b],[c,d]])
    L.append(dict(test='fisher',subset=tag,feature=lab,miss_with=f'{a}/{a+b}',miss_without=f'{c}/{c+d}',OR=round(orr,2),p=f'{p:.2g}'))
feats=[('unreacted','CCD attachment atom in pre-reaction form'),('reacted','CCD attachment atom in post-reaction form'),
       ('leaving','CCD defines any leaving atom'),('leaving_at_attach','leaving atom on attachment atom'),
       ('newcode','code first defined by this entry'),('new_at_dep','code created on/after deposition date'),('pcm_blank','CCD has no PCM annotation (circular)')]
for c,l in feats:
    fish(c,l)
    fish(c,l,s[s.site!='PDBJ'],'non-PDBj'); fish(c,l,s[s.site=='PDBJ'],'PDBj')
# Mantel-Haenszel across site strata
for c,l in feats:
    tabs=[]
    for st,x in s.groupby('site'):
        t=np.array([[((x[c]==1)&(x.miss==1)).sum(),((x[c]==1)&(x.miss==0)).sum()],[((x[c]==0)&(x.miss==1)).sum(),((x[c]==0)&(x.miss==0)).sum()]])+0.0
        if t.sum(1).min()>0: tabs.append(t)
    st_=StratifiedTable(tabs); L.append(dict(test='Mantel-Haenszel by site',subset='all',feature=l,OR=round(st_.oddsratio_pooled,2),p=f'{st_.test_null_odds().pvalue:.2g}'))
pd.DataFrame(L).to_csv(f'{O}/feature_tests.tsv',sep='\t',index=False)
print(pd.DataFrame(L).to_string())
# logistic regressions, SE clustered by CCD code
R=[]
s['code']=s.lig_comp.astype('category').cat.codes
for name,form in [('site+year','miss ~ C(site, Treatment("RCSB")) + yearc'),
    ('site+year+unreacted','miss ~ C(site, Treatment("RCSB")) + yearc + unreacted'),
    ('site+year+leaving','miss ~ C(site, Treatment("RCSB")) + yearc + leaving'),
    ('site+year+newcode','miss ~ C(site, Treatment("RCSB")) + yearc + newcode'),
    ('site+year+all3','miss ~ C(site, Treatment("RCSB")) + yearc + unreacted + leaving + newcode'),
    ('site+year+all3+pcm','miss ~ C(site, Treatment("RCSB")) + yearc + unreacted + leaving + newcode + pcm_blank')]:
    m=smf.glm(form,s,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':s.code})
    for k in m.params.index:
        R.append(dict(model=name,term=k,OR=round(np.exp(m.params[k]),2),ci_lo=round(np.exp(m.conf_int().loc[k,0]),2),ci_hi=round(np.exp(m.conf_int().loc[k,1]),2),p=f'{m.pvalues[k]:.2g}',n=int(m.nobs),aic=round(m.aic,1)))
pd.DataFrame(R).to_csv(f'{O}/logistic_models.tsv',sep='\t',index=False)
print(pd.DataFrame(R).to_string())
s.to_csv(f'{O}/contacts_lt2_model_input.tsv',sep='\t',index=False)
print(pd.crosstab([s.attach_class],[s.site,s.miss]))
# CCD age at deposition and CCD creation year (continuous), adjusted for site and year
A=[]
x=s.dropna(subset=['ccd_age_at_dep_y','ccd_year']).copy(); x['age']=x.ccd_age_at_dep_y.clip(lower=0); x['ccdy']=x.ccd_year-2010
for name,form,term in [('site+year+ccd_age_at_deposition_years','miss ~ C(site, Treatment("RCSB")) + yearc + age','age'),('site+year+ccd_creation_year','miss ~ C(site, Treatment("RCSB")) + yearc + ccdy','ccdy')]:
    m=smf.glm(form,x,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':x.code})
    A.append(dict(model=name,term=term,OR_per_year=round(np.exp(m.params[term]),3),ci_lo=round(np.exp(m.conf_int().loc[term,0]),3),ci_hi=round(np.exp(m.conf_int().loc[term,1]),3),p=f'{m.pvalues[term]:.2g}',n=int(m.nobs)))
pd.DataFrame(A).to_csv(f'{O}/ccd_age_models.tsv',sep='\t',index=False); print(pd.DataFrame(A))
# leaving atom bonded to the attachment atom (CCD written with the linkage in mind)
R2=[]
for name,form in [('site+year+leaving_at_attach','miss ~ C(site, Treatment("RCSB")) + yearc + leaving_at_attach'),
                  ('site+year+leaving_at_attach+unreacted+newcode','miss ~ C(site, Treatment("RCSB")) + yearc + leaving_at_attach + unreacted + newcode')]:
    m=smf.glm(form,s,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':s.code})
    for k in m.params.index: R2.append(dict(model=name,term=k,OR=round(np.exp(m.params[k]),2),ci_lo=round(np.exp(m.conf_int().loc[k,0]),2),ci_hi=round(np.exp(m.conf_int().loc[k,1]),2),p=f'{m.pvalues[k]:.2g}',n=int(m.nobs),aic=round(m.aic,1)))
pd.concat([pd.read_csv(f'{O}/logistic_models.tsv',sep='\t'),pd.DataFrame(R2)]).to_csv(f'{O}/logistic_models.tsv',sep='\t',index=False)
print(pd.DataFrame(R2).to_string())
print(pd.crosstab(s[s.leaving_at_attach==1].attach_class,s[s.leaving_at_attach==1].miss))
print(pd.crosstab(s[s.leaving_at_attach==1].site,s[s.leaving_at_attach==1].miss))
