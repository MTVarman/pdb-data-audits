"""Step 5: test copy-forward of a code's link template.
(1) template fields of recorded units vs the code's first recorded unit, against an era/site-matched baseline of other codes;
(2) misses vs refinement software, release era, revision activity on struct_conn, independent of code;
(3) template oddities and whether they propagate."""
import json, pandas as pd, numpy as np, collections, random, statsmodels.formula.api as smf, statsmodels.api as sm
from scipy.stats import fisher_exact, wilcoxon
O='/workspace/ligand_code_effect'; OUT=f'{O}/out/step5'
CC=json.load(open(f'{O}/out/ccd_parsed_archive.json'))
S=pd.read_csv(f'{O}/out/step4/strict_units.tsv',sep='\t',low_memory=False); S['dep']=pd.to_datetime(S.dep)
J={}
for l in open(f'{OUT}/templates.jsonl'):
    r=json.loads(l); J[r['e']]=r
L=[];P=lambda *a:L.append(' '.join(map(str,a)))
def mode(x):
    x=[y for y in x]; return collections.Counter(x).most_common(1)[0][0] if x else None
rows=[]
for _,u in S.iterrows():
    r=J[u.entry]; cv=[x for x in r['cov'] if x['code']==u.lig_comp and (x['nuc_res'],x['nuc_atom']) in {('CYS','SG'),('SER','OG'),('LYS','NZ'),('HIS','NE2'),('HIS','ND1')}]
    cc=[x for x in r['cc'] if x['code']==u.lig_comp]
    ref=r['refine'][0].split('|') if r['refine'] else ['none','']
    revs=sorted(r['rev'],key=lambda x:int(x[0])) if r['rev'] else []
    rel=revs[0][3] if revs else None
    sc_rev=[x for x in r['rev_cat'] if x[1]=='struct_conn' and x[0]!='1']
    sc_rev_dates=[d for o,_,_,d in revs if o in {x[0] for x in sc_rev}]
    pcm_atoms={p[1] for p in (CC.get(u.lig_comp,{}).get('pcm_rows') or []) if p[1]}
    ccd_atoms=set((CC.get(u.lig_comp,{}).get('el') or {}).keys())
    lv=set(CC.get(u.lig_comp,{}).get('leaving_atoms') or [])
    d=dict(entry=u.entry,code=u.lig_comp,miss=u.miss,status=u.status,dep=u.dep,site=u.site,lab=u.lab,year=u.year,
           refine=ref[0].upper().replace('PHENIX.REAL_SPACE_REFINE','PHENIX'),refine_ver=ref[1] if len(ref)>1 else '',release=rel,
           n_cov=len(cv),struct_conn_revised=len(sc_rev)>0,struct_conn_rev_dates=';'.join(sc_rev_dates),hpm=r.get('pdbx_entry_hpm'))
    if cv:
        d.update(pair=mode([f"{x['nuc_res']}.{x['nuc_atom']}-{x['lig_atom']}" for x in cv]),lig_atom=mode([x['lig_atom'] for x in cv]),
                 leaving=mode([x['leaving'] for x in cv]),prot_first=mode([x['prot_first'] for x in cv]),value_order=mode([x['value_order'] for x in cv]),
                 details=mode([x['details'] for x in cv]),covale_type=mode([x['conn_type'] for x in cv]),dist_max=max(float(x['dist_str']) for x in cv if x['dist_str']),
                 n_pairs=len({(x['nuc_atom'],x['lig_atom']) for x in cv}),
                 atom_not_in_ccd=any(x['lig_atom'] not in ccd_atoms for x in cv) if ccd_atoms else None,
                 atom_is_leaving=any(x['lig_atom'] in lv for x in cv),
                 atom_not_pcm=(mode([x['lig_atom'] for x in cv]) not in pcm_atoms) if pcm_atoms else None,
                 # oddity: same nucleophile listed in close contacts with a DIFFERENT atom of the code, closer than the covale distance
                 shorter_other_atom=any(c['lig_atom']!=mode([x['lig_atom'] for x in cv]) and c['dist']<min(float(x['dist_str']) for x in cv if x['dist_str']) for c in cc),
                 long_link=max(float(x['dist_str']) for x in cv if x['dist_str'])>2.3)
    rows.append(d)
U=pd.DataFrame(rows); U.to_csv(f'{OUT}/unit_templates.tsv',sep='\t',index=False)
P('units',len(U),'recorded with covale rows',int(U.n_cov.gt(0).sum()),'misses',int(U.miss.sum()))
# sanity: recorded units with no covale rows (mixed definition / other nucleophile)
P('recorded/mixed units lacking covale rows for code:',int(((U.miss==0)&(U.n_cov==0)).sum()))
# ---------- (1) template copy test ----------
R=U[U.n_cov>0].sort_values(['dep','release','entry'])
first=R.groupby('code').head(1).set_index('code')
FIELDS=['pair','leaving','prot_first','value_order','details','covale_type','n_pairs']
later=R[R.dep>R.code.map(first.dep)].copy()
P('\n(1) codes with >=2 recorded units (later recorded after first recorded):',later.code.nunique(),'later recorded units',len(later))
first['fyear']=first.dep.dt.year
rng=random.Random(0); res=[]
for f in FIELDS:
    own=[];base=[]
    for _,u in later.iterrows():
        fc=first.loc[u.code]; own.append(str(u[f])==str(fc[f]))
        if f=='pair': continue
        pool=first[(first.fyear==fc.fyear)&(first.site==fc.site)&(first.index!=u.code)]
        if len(pool)<3: pool=first[(abs(first.fyear-fc.fyear)<=1)&(first.index!=u.code)]
        base.append(np.mean([str(u[f])==str(x) for x in pool[f]]) if len(pool) else np.nan)
    own=np.array(own,float); 
    if f!='pair':
        base=np.array(base,float); ok=~np.isnan(base)
        diff=own[ok]-base[ok]
        try: p=wilcoxon(diff).pvalue
        except Exception: p=np.nan
        res.append(dict(field=f,later_units=int(ok.sum()),match_own_first=round(own[ok].mean(),3),match_era_site_baseline=round(base[ok].mean(),3),excess=round(diff.mean(),3),wilcoxon_p=f'{p:.2g}'))
    else: res.append(dict(field=f,later_units=len(own),match_own_first=round(own.mean(),3)))
RES=pd.DataFrame(res); RES.to_csv(f'{OUT}/template_match.tsv',sep='\t',index=False); P(RES.to_string(index=False))
# informative subset: codes whose first template has a non-modal value (an "oddity") in a field
P('\n(1b) non-default first-template values, and whether later units repeat them vs the era baseline')
modal={f:R[f].astype(str).mode()[0] for f in FIELDS if f!='pair'}
P('archive modal values:',modal)
ODD=[]
for f in ['leaving','prot_first','value_order','details','covale_type','n_pairs']:
    fc=first[first[f].astype(str)!=modal[f]]
    lu=later[later.code.isin(fc.index)]
    rep=int(sum(str(u[f])==str(fc.loc[u.code,f]) for _,u in lu.iterrows()))
    lu2=later[~later.code.isin(fc.index)]
    rep2=int(sum(str(u[f])!=modal[f] for _,u in lu2.iterrows()))
    if len(lu): 
        o,p=fisher_exact([[rep,len(lu)-rep],[rep2,len(lu2)-rep2]])
        ODD.append(dict(field=f,codes_with_odd_first=len(fc),their_later_units=len(lu),later_repeat_odd=rep,other_codes_later_units=len(lu2),other_later_nondefault=rep2,OR=round(o,2),p=f'{p:.2g}'))
ODD=pd.DataFrame(ODD); ODD.to_csv(f'{OUT}/oddity_first_propagation.tsv',sep='\t',index=False); P(ODD.to_string(index=False))
# ---------- (3) link-atom oddities ----------
P('\n(3) link-atom oddities among recorded units')
for k in ['shorter_other_atom','long_link','atom_not_in_ccd','atom_is_leaving','atom_not_pcm']:
    x=R[k].dropna().astype(bool); P(f'  {k}: {int(x.sum())}/{len(x)} units, codes {R[R[k].fillna(False).astype(bool)].code.nunique()}')
    fo=first[first[k].fillna(False).astype(bool)].index; lu=later[later.code.isin(fo)]; lo=later[~later.code.isin(fo)]
    a=int(lu[k].fillna(False).astype(bool).sum()); b=int(lo[k].fillna(False).astype(bool).sum())
    if len(fo):
        o,p=fisher_exact([[a,len(lu)-a],[b,len(lo)-b]])
        P(f'     first recorded unit odd: {len(fo)} codes; their later units odd {a}/{len(lu)} vs other codes later odd {b}/{len(lo)}  OR {o:.2f} p={p:.2g}')
R.to_csv(f'{OUT}/recorded_templates.tsv',sep='\t',index=False)
# ---------- (2) era / software / revisions ----------
U['era']=pd.cut(U.dep.dt.year,[0,2007,2013,2019,2030],labels=['<=2007','2008-13','2014-19','2020+'])
U['ref']=U.refine.where(U.refine.isin(['PHENIX','REFMAC','BUSTER','CNS']),'other')
U['rel_lag']=(pd.to_datetime(U.release)-U.dep).dt.days
P('\n(2) misses by refinement program:'); P(U.groupby('ref').miss.agg(['sum','size']).to_string())
P('misses by deposition era:'); P(U.groupby('era',observed=True).miss.agg(['sum','size']).to_string())
P('struct_conn touched in a later revision: recorded',int(U[(U.miss==0)].struct_conn_revised.sum()),'/',int((U.miss==0).sum()),' missed',int(U[U.miss==1].struct_conn_revised.sum()),'/',int(U.miss.sum()))
# code history from step 4
U=U.merge(S[['entry','lig_comp','c_first_miss','c_first_dep','later','yearc']].rename(columns={'lig_comp':'code'}),on=['entry','code'])
U['pdbj']=(U.site=='PDBJ').astype(int)
fe=U.groupby('code').head(0)
cf=S.sort_values('dep').groupby('lig_comp').head(1).set_index('lig_comp')
U['c_first_era']=U.code.map(pd.cut(cf.dep.dt.year,[0,2007,2013,2019,2030],labels=['<=2007','2008-13','2014-19','2020+']).astype(str))
U['c_first_site']=U.code.map(cf.site); U['c_first_ref']=U.code.map(U.sort_values('dep').groupby('code').head(1).set_index('code').ref)
U.to_csv(f'{OUT}/unit_templates.tsv',sep='\t',index=False)
M=U[U.later].copy(); M['codeid']=M.code.astype('category').cat.codes
P('\nlater units',len(M),'misses',int(M.miss.sum()))
for form in ['miss ~ c_first_miss + C(ref, Treatment("PHENIX")) + C(era) + pdbj',
             'miss ~ c_first_miss + C(c_first_era) + C(ref, Treatment("PHENIX")) + C(era) + pdbj',
             'miss ~ c_first_miss + C(c_first_site, Treatment("RCSB")) + C(ref, Treatment("PHENIX")) + C(era) + pdbj']:
    m=smf.glm(form,M,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':M.codeid}); ci=m.conf_int()
    P('\n'+form)
    for k in m.params.index:
        if k=='Intercept': continue
        P(f'   {k}: OR {np.exp(m.params[k]):.2f} ({np.exp(ci.loc[k,0]):.2f}-{np.exp(ci.loc[k,1]):.2f}) p={m.pvalues[k]:.2g}')
# release-order test for carried misses: was the later entry deposited after the first entry's release?
cm=M[M.c_first_miss==1].copy(); cm['first_release']=cm.code.map(U.sort_values('dep').groupby('code').head(1).set_index('code').release)
cm['after_first_release']=cm.dep>pd.to_datetime(cm.first_release)
P('\nlater uses of first-missed codes: deposited after the first entry was released:',int(cm.after_first_release.sum()),'/',len(cm),
  '| missed among after:',int(cm[cm.after_first_release].miss.sum()),'/',int(cm.after_first_release.sum()),' among before:',int(cm[~cm.after_first_release].miss.sum()),'/',int((~cm.after_first_release).sum()))
cm.to_csv(f'{OUT}/first_missed_later_units.tsv',sep='\t',index=False)
open(f'{OUT}/27_stdout.txt','w').write('\n'.join(L)); print('\n'.join(L))
