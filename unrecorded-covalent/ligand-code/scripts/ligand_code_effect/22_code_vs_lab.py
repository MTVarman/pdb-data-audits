"""Step 4: does the carried-forward miss follow the CCD code, the lab (last author) or the site?
Strict set = drug-like codes with <=10 unbonded uses, closest bond at a bond-like distance (SG 1.65-1.95; O/N 1.20-1.65 A).
Unit = (code, entry). Miss = >=1 unrecorded bond <2.0 A for that code in that entry."""
import json, pandas as pd, numpy as np, statsmodels.formula.api as smf, statsmodels.api as sm, itertools, sys
from scipy.stats import fisher_exact
from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM
from statsmodels.discrete.conditional_models import ConditionalLogit
from rdkit import Chem, DataStructs, RDLogger; from rdkit.Chem import AllChem; RDLogger.DisableLog('rdApp.*')
O='/workspace/ligand_code_effect'; OUT=f'{O}/out/step4'
C=json.load(open(f'{O}/out/ccd_parsed_archive.json'))
E=pd.read_csv(f'{O}/out/step3/entry_level_bonded.tsv',sep='\t',low_memory=False)
TL=pd.read_csv(f'{O}/archive/code_timelines.tsv',sep='\t',low_memory=False)
B=pd.read_csv(f'{O}/out/archive_bonds_lt2.tsv',sep='\t',low_memory=False)
nun=TL[TL.status=='unbonded'].groupby('lig_comp').size(); E['nun']=E.lig_comp.map(nun).fillna(0)
cb=B.sort_values('dist').groupby(['lig_comp','entry']).first()[['nuc_atom','dist','lig_atom']]; E=E.join(cb,on=['lig_comp','entry'])
E['bondlike']=np.where(E.nuc_atom=='SG',E.dist.between(1.65,1.95),E.dist.between(1.20,1.65)); E['dep']=pd.to_datetime(E.dep)
S=E[(E.cat=='drug-like')&(E.nun<=10)&E.bondlike&E.au_last.notna()].sort_values(['dep','rel','entry']).copy().reset_index(drop=True)
S['lab']=S.au_last; S['pdbj']=(S.site=='PDBJ').astype(int)
S['ik14']=[((C.get(c) or {}).get('inchikey') or c)[:14] for c in S.lig_comp]
S['skel']=[(C.get(c) or {}).get('skeleton') or c for c in S.lig_comp]
L=[]; P=lambda *a: L.append(' '.join(map(str,a)))
P('strict units',len(S),'misses',int(S.miss.sum()),'codes',S.lig_comp.nunique(),'labs',S.lab.nunique())
# code history
f=S.groupby('lig_comp').head(1).set_index('lig_comp')
S['c_first_miss']=S.lig_comp.map(f.miss); S['c_first_dep']=S.lig_comp.map(f.dep); S['c_first_lab']=S.lig_comp.map(f.lab); S['c_first_site']=S.lig_comp.map(f.site); S['c_first_entry']=S.lig_comp.map(f.entry)
S['later']=S.dep>S.c_first_dep
S['same_lab']=S.lab==S.c_first_lab; S['same_site']=S.site==S.c_first_site
# lab history (any code, strict set): earlier units strictly before this deposition date
S['lab_prior_n']=0; S['lab_prior_miss']=0
for lab,g in S.groupby('lab'):
    for i,r in g.iterrows():
        prev=g[g.dep<r.dep]; S.at[i,'lab_prior_n']=len(prev); S.at[i,'lab_prior_miss']=int(prev.miss.sum()>0)
S.to_csv(f'{OUT}/strict_units.tsv',sep='\t',index=False)
def rate(x): return f'{int(x.miss.sum())}/{len(x)} ({100*x.miss.mean():.0f}%)' if len(x) else '0/0'
def orci(a,b,c,d):  # a miss/b rec in group1, c miss/d rec in group2 ; Haldane-corrected Wald CI
    o,p=fisher_exact([[a,b],[c,d]]); A,Bq,Cc,D=[v+.5 for v in (a,b,c,d)]; lo=np.log(A*D/(Bq*Cc)); se=np.sqrt(1/A+1/Bq+1/Cc+1/D)
    return f'OR {np.exp(lo):.1f} ({np.exp(lo-1.96*se):.1f}-{np.exp(lo+1.96*se):.1f}), Fisher p={p:.2g}'
def compare(x1,x2): return orci(int(x1.miss.sum()),int((1-x1.miss).sum()),int(x2.miss.sum()),int((1-x2.miss).sum()))
SUBS={'strict':S,'drop HUR/80I':S[~S.lig_comp.isin(['HUR','80I'])],'drop Li J. lab':S[S.lab!='Li, J.'],'drop PDBj':S[S.site!='PDBJ']}
RES=[]
for name,X in SUBS.items():
    P(f'\n==================== {name}: units {len(X)}, misses {int(X.miss.sum())}')
    # recompute code-first within subset (so dropped entries cannot define history)
    X=X.copy(); f=X.groupby('lig_comp').head(1).set_index('lig_comp')
    for k,col in [('c_first_miss','miss'),('c_first_dep','dep'),('c_first_lab','lab'),('c_first_site','site')]: X[k]=X.lig_comp.map(f[col])
    X['later']=X.dep>X.c_first_dep; X['same_lab']=X.lab==X.c_first_lab; X['same_site']=X.site==X.c_first_site
    Lt=X[X.later]
    # ---- 1a same code, different labs ----
    P('[1a] later uses of a code, by first-use status x lab x site')
    for fm in [1,0]:
        for sl,ss in [(True,True),(True,False),(False,True),(False,False)]:
            x=Lt[(Lt.c_first_miss==fm)&(Lt.same_lab==sl)&(Lt.same_site==ss)]
            P(f'   first {"MISSED" if fm else "recorded"} | {"same" if sl else "other"} lab | {"same" if ss else "other"} site : {rate(x)}  codes={x.lig_comp.nunique()}')
    a=Lt[(Lt.c_first_miss==1)&~Lt.same_lab]; b=Lt[(Lt.c_first_miss==0)&~Lt.same_lab]
    P('   other-lab uses: first-missed vs first-recorded codes:',rate(a),'vs',rate(b),compare(a,b))
    a2=Lt[(Lt.c_first_miss==1)&~Lt.same_lab&~Lt.same_site]; b2=Lt[(Lt.c_first_miss==0)&~Lt.same_lab&~Lt.same_site]
    P('   other-lab AND other-site uses:',rate(a2),'vs',rate(b2),compare(a2,b2))
    a3=Lt[(Lt.c_first_miss==1)&Lt.same_lab]; P('   same-lab reuse of first-missed codes:',rate(a3))
    RES.append(dict(subset=name,design='same code, other lab',fm=rate(a),fr=rate(b),stat=compare(a,b)))
    RES.append(dict(subset=name,design='same code, other lab + other site',fm=rate(a2),fr=rate(b2),stat=compare(a2,b2)))
    # ---- 1b same lab, different codes ----
    P('[1b] labs with an early miss: their LATER units, same code vs other codes')
    rows=[]
    for lab,g in X.groupby('lab'):
        m=g[g.miss==1]
        if not len(m): continue
        t0=m.dep.min(); mc=set(m[m.dep==t0].lig_comp); mik=set(X[X.lig_comp.isin(mc)].ik14); msk=set(X[X.lig_comp.isin(mc)].skel)
        for _,r in g[g.dep>t0].iterrows():
            rel='same code' if r.lig_comp in mc else ('other code, same drug (InChIKey-14 or skeleton)' if (r.ik14 in mik or r.skel in msk) else 'other code, different drug')
            rows.append(dict(lab=lab,entry=r.entry,code=r.lig_comp,missed_codes=','.join(sorted(mc)),relation=rel,miss=r.miss,site=r.site,dep=r.dep,first_miss_dep=t0))
    Rb=pd.DataFrame(rows)
    if name=='strict': Rb.to_csv(f'{OUT}/lab_after_miss_units.tsv',sep='\t',index=False)
    # baseline: labs whose first unit was recorded and had no miss before -> later units
    base=X[(X.lab_prior_n>0)&(X.lab_prior_miss==0)]
    for rel in ['same code','other code, same drug (InChIKey-14 or skeleton)','other code, different drug']:
        x=Rb[Rb.relation==rel] if len(Rb) else Rb; P(f'   {rel}: {rate(x) if len(x) else "0/0"}  labs={x.lab.nunique() if len(x) else 0}')
    P('   baseline (later units of labs with earlier units, none missed):',rate(base))
    if len(Rb):
        oc=Rb[Rb.relation!='same code']; sc=Rb[Rb.relation=='same code']
        P('   same code vs other codes (same lab after its first miss):',compare(sc,oc))
        P('   other codes after a lab miss vs baseline labs:',compare(oc,base))
        RES.append(dict(subset=name,design='same lab after miss: same code vs other codes',fm=rate(sc),fr=rate(oc),stat=compare(sc,oc)))
        RES.append(dict(subset=name,design='same lab after miss, other codes vs labs with no prior miss',fm=rate(oc),fr=rate(base),stat=compare(oc,base)))
    # ---- 1c joint models ----
    M=Lt.copy(); M['lab_prior_miss_other']=0
    # lab prior miss on a DIFFERENT code (so it is not the code history)
    for lab,g in X.groupby('lab'):
        for i,r in M[M.lab==lab].iterrows():
            prev=g[(g.dep<r.dep)&(g.lig_comp!=r.lig_comp)]; M.at[i,'lab_prior_miss_other']=int(prev.miss.sum()>0)
    M['code']=M.lig_comp.astype('category').cat.codes; M['labid']=M.lab.astype('category').cat.codes
    site_term='+ pdbj' if name!='drop PDBj' else ''
    try:
        m=smf.glm(f'miss ~ c_first_miss + lab_prior_miss_other + same_lab {site_term} + yearc',M,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':M.code})
        ci=m.conf_int()
        for k in m.params.index:
            if k=='Intercept': continue
            P(f'   [1c GLM, cluster=code] {k}: OR {np.exp(m.params[k]):.2f} ({np.exp(ci.loc[k,0]):.2f}-{np.exp(ci.loc[k,1]):.2f}) p={m.pvalues[k]:.2g}')
            RES.append(dict(subset=name,design=f'GLM {k}',stat=f'OR {np.exp(m.params[k]):.2f} ({np.exp(ci.loc[k,0]):.2f}-{np.exp(ci.loc[k,1]):.2f}) p={m.pvalues[k]:.2g}',n=len(M)))
        m2=smf.glm(f'miss ~ c_first_miss + lab_prior_miss_other + same_lab {site_term} + yearc',M,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':M.labid})
        k='c_first_miss'; ci=m2.conf_int(); P(f'   [1c GLM, cluster=lab] c_first_miss OR {np.exp(m2.params[k]):.2f} ({np.exp(ci.loc[k,0]):.2f}-{np.exp(ci.loc[k,1]):.2f}) p={m2.pvalues[k]:.2g}')
        k='lab_prior_miss_other'; P(f'   [1c GLM, cluster=lab] lab_prior_miss_other OR {np.exp(m2.params[k]):.2f} ({np.exp(ci.loc[k,0]):.2f}-{np.exp(ci.loc[k,1]):.2f}) p={m2.pvalues[k]:.2g}')
    except Exception as ex: P('   GLM failed',ex)
    # fixed-effect swap: conditional logit within CODE (does lab history matter holding code fixed?) and within LAB (does code history matter holding lab fixed?)
    for gk,cov in [('lig_comp',['lab_prior_miss_other','same_lab']),('lab',['c_first_miss'])]:
        G=M[M.groupby(gk).miss.transform('nunique')==2]
        G=G[[c for c in cov if G[c].astype(float).nunique()>1]+[gk,'miss']]
        cov2=[c for c in cov if c in G.columns]
        if len(G)>5 and cov2:
            try:
                cl=ConditionalLogit(G.miss,G[cov2].astype(float),groups=G[gk]).fit(disp=0); ci=cl.conf_int()
                for k in cov2:
                    P(f'   [1c conditional logit within {gk}] {k}: OR {np.exp(cl.params[k]):.2f} ({np.exp(ci.loc[k,0]):.2f}-{np.exp(ci.loc[k,1]):.2f}) p={cl.pvalues[k]:.2g} (informative groups {G[gk].nunique()}, units {len(G)}, misses {int(G.miss.sum())})')
                    RES.append(dict(subset=name,design=f'clogit within {gk}: {k}',stat=f'OR {np.exp(cl.params[k]):.2f} ({np.exp(ci.loc[k,0]):.2f}-{np.exp(ci.loc[k,1]):.2f}) p={cl.pvalues[k]:.2g}',n=len(G)))
            except Exception as ex: P('   clogit failed',gk,ex)
        else: P(f'   [1c conditional logit within {gk}] not estimable: informative units {len(G)}')
    print(name,'done',flush=True)
    # conditional logit within drug (InChIKey-14 group): only groups with >=2 codes and outcome variation inform
    for gk in ['ik14','skel']:
        G=M[M.groupby(gk).lig_comp.transform('nunique')>=2]
        G=G[G.groupby(gk).miss.transform('nunique')==2]
        if len(G)>5 and G.c_first_miss.nunique()>1:
            try:
                cv=[k for k in ['c_first_miss','lab_prior_miss_other'] if G[k].nunique()>1]
                cl=ConditionalLogit(G.miss,G[cv].astype(float),groups=G[gk]).fit(disp=0)
                ci=cl.conf_int()
                for k in cv:
                    P(f'   [1c conditional logit within drug={gk}] {k}: OR {np.exp(cl.params[k]):.2f} ({np.exp(ci.loc[k,0]):.2f}-{np.exp(ci.loc[k,1]):.2f}) p={cl.pvalues[k]:.2g}  (informative groups {G[gk].nunique()}, units {len(G)})')
            except Exception as ex: P('   clogit failed',gk,ex)
        else: P(f'   [1c conditional logit within drug={gk}] not estimable: informative units {len(G)}')
    # 2x2 within drug groups having >=2 codes
    D=M[M.groupby('ik14').lig_comp.transform('nunique')>=2]
    P('   [1c 2x2, drugs with >=2 codes (InChIKey-14)] units',len(D))
    for cf in [1,0]:
        for lp in [1,0]:
            x=D[(D.c_first_miss==cf)&(D.lab_prior_miss_other==lp)]; P(f'      code first {"missed" if cf else "recorded"}, lab prior miss on other code={lp}: {rate(x)}')
pd.DataFrame(RES).to_csv(f'{OUT}/design_results.tsv',sep='\t',index=False)
open(f'{OUT}/22_stdout.txt','w').write('\n'.join(L)); print('\n'.join(L))
