"""Step 6b. Rerun the ligand-code carry-forward tests on contacts inside the MAXIT windows only
(all ten windows of Varman distance-table note v1.1, Table 1, that an S/O/N nucleophile can reach; see script 29).
Unit = (code, entry) with >=1 in-window contact < 2.0 A; miss = >=1 of those in-window contacts unrecorded.
Out-of-window and no-window contacts are ignored (a unit whose only miss is out of window counts as recorded).
In-window strict set = drug-like codes with <=10 unbonded uses, last author known, >=1 in-window contact
(the in-window test replaces script 22's bond-like distance test). Code history is recomputed on the in-window units.
Design code follows scripts 20 (step-3 model) and 22 (step-4 designs) unchanged; a validation pass re-runs the same
functions on the original strict set and must reproduce the published step-3/4 numbers.
Writes out/step6/30_stdout.txt, window_units.tsv, window_strict_units.tsv, window_results.tsv, window_attribution.tsv."""
import json, pandas as pd, numpy as np, statsmodels.formula.api as smf, statsmodels.api as sm, warnings
from scipy.stats import fisher_exact
warnings.filterwarnings('ignore')
O='/workspace/ligand_code_effect'; OUT=f'{O}/out/step6'
CCD=json.load(open(f'{O}/out/ccd_parsed_archive.json'))
E=pd.read_csv(f'{O}/out/step3/entry_level_bonded.tsv',sep='\t',dtype={'entry':str},low_memory=False)
TL=pd.read_csv(f'{O}/archive/code_timelines.tsv',sep='\t',dtype={'entry':str},low_memory=False)
B=pd.read_csv(f'{O}/out/archive_bonds_lt2.tsv',sep='\t',dtype={'entry':str},low_memory=False)
W=pd.read_csv(f'{OUT}/contacts_windowed.tsv',sep='\t',dtype={'entry':str},low_memory=False)
nun=TL[TL.status=='unbonded'].groupby('lig_comp').size(); E['nun']=E.lig_comp.map(nun).fillna(0)
cb=B.sort_values('dist').groupby(['lig_comp','entry']).first()[['nuc_atom','dist','lig_atom']]; E=E.join(cb,on=['lig_comp','entry'])
E['bondlike']=np.where(E.nuc_atom=='SG',E.dist.between(1.65,1.95),E.dist.between(1.20,1.65)); E['dep']=pd.to_datetime(E.dep)
I=W[W.inside&(W.dist<2.0)]
iu=I.groupby(['lig_comp','entry']).agg(n_in=('recorded','size'),n_in_unrec=('recorded',lambda x:int((~x).sum())))
ku=set(zip(E.lig_comp,E.entry)); notin=[k not in ku for k in zip(I.lig_comp,I.entry)]
E=E.join(iu,on=['lig_comp','entry']); E['n_in']=E.n_in.fillna(0).astype(int); E['n_in_unrec']=E.n_in_unrec.fillna(0).astype(int)
E['miss_orig']=E.miss
L=[]; P=lambda *a: L.append(' '.join(map(str,a))); RES=[]
def rate(x): return f'{int(x.miss.sum())}/{len(x)} ({100*x.miss.mean():.0f}%)' if len(x) else '0/0'
def orci(a,b,c,d):
    o,p=fisher_exact([[a,b],[c,d]]); A,Bq,Cc,D=[v+.5 for v in (a,b,c,d)]; lo=np.log(A*D/(Bq*Cc)); se=np.sqrt(1/A+1/Bq+1/Cc+1/D)
    return f'OR {np.exp(lo):.1f} ({np.exp(lo-1.96*se):.1f}-{np.exp(lo+1.96*se):.1f}), Fisher p={p:.2g}'
def compare(x1,x2): return orci(int(x1.miss.sum()),int((1-x1.miss).sum()),int(x2.miss.sum()),int((1-x2.miss).sum()))
def fmt(m,k):
    ci=m.conf_int(); return f'OR {np.exp(m.params[k]):.2f} ({np.exp(ci.loc[k,0]):.2f}-{np.exp(ci.loc[k,1]):.2f}) p={m.pvalues[k]:.2g}'
def prep(S):
    S=S.sort_values(['dep','rel','entry']).copy().reset_index(drop=True)
    S['lab']=S.au_last; S['pdbj']=(S.site=='PDBJ').astype(int)
    S['ik14']=[((CCD.get(c) or {}).get('inchikey') or c)[:14] for c in S.lig_comp]
    S['skel']=[(CCD.get(c) or {}).get('skeleton') or c for c in S.lig_comp]
    S['lab_prior_n']=0; S['lab_prior_miss']=0
    for lab,g in S.groupby('lab'):
        for i,r in g.iterrows():
            prev=g[g.dep<r.dep]; S.at[i,'lab_prior_n']=len(prev); S.at[i,'lab_prior_miss']=int(prev.miss.sum()>0)
    return S
def step3(S,tag):
    D=S.sort_values(['lig_comp','dep','rel','entry']).copy(); f=D.groupby('lig_comp').head(1).set_index('lig_comp')
    D['first_miss']=D.lig_comp.map(f.miss); D['first_dep']=D.lig_comp.map(f.dep); D['first_lab']=D.lig_comp.map(f.au_last)
    x=D[D.dep>D.first_dep].copy(); x['same_lab']=x.au_last==x.first_lab; x['code']=x.lig_comp.astype('category').cat.codes
    a=x[x.first_miss==1]; b=x[x.first_miss==0]
    P(f'[{tag}] step-3 design: later uses {len(x)}; first missed {rate(a)} codes={a.lig_comp.nunique()}; first recorded {rate(b)}')
    m=smf.glm('miss ~ first_miss + same_lab + C(site, Treatment("RCSB")) + yearc',x,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':x.code})
    for k in ['first_miss','same_lab[T.True]','C(site, Treatment("RCSB"))[T.PDBJ]']:
        P(f'   [{tag}] step-3 GLM {k}: {fmt(m,k)}'); RES.append(dict(set=tag,subset='all',design=f'step-3 GLM {k}',stat=fmt(m,k),n=len(x)))
    y=x[~x.same_lab].copy(); y['code']=y.lig_comp.astype('category').cat.codes
    m2=smf.glm('miss ~ first_miss + C(site, Treatment("RCSB")) + yearc',y,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':y.code})
    P(f'   [{tag}] step-3 GLM other-lab later uses only: first_miss {fmt(m2,"first_miss")} n={len(y)} misses={int(y.miss.sum())}')
    RES.append(dict(set=tag,subset='all',design='step-3 GLM other-lab only first_miss',stat=fmt(m2,'first_miss'),n=len(y)))
    return D,x,m
def step4(S,tag):
    SUBS={'strict':S,'drop HUR/80I':S[~S.lig_comp.isin(['HUR','80I'])],'drop Li J. lab':S[S.lab!='Li, J.'],'drop PDBj':S[S.site!='PDBJ']}
    for name,X in SUBS.items():
        X=X.copy(); f=X.groupby('lig_comp').head(1).set_index('lig_comp')
        for k,col in [('c_first_miss','miss'),('c_first_dep','dep'),('c_first_lab','lab'),('c_first_site','site')]: X[k]=X.lig_comp.map(f[col])
        X['later']=X.dep>X.c_first_dep; X['same_lab']=X.lab==X.c_first_lab; X['same_site']=X.site==X.c_first_site
        Lt=X[X.later]
        P(f'[{tag}] == {name}: units {len(X)}, misses {int(X.miss.sum())}, first-missed codes with a later use {Lt[Lt.c_first_miss==1].lig_comp.nunique()}')
        for fm in [1,0]:
            for sl,ss in [(True,True),(True,False),(False,True),(False,False)]:
                x=Lt[(Lt.c_first_miss==fm)&(Lt.same_lab==sl)&(Lt.same_site==ss)]
                P(f'   first {"MISSED" if fm else "recorded"} | {"same" if sl else "other"} lab | {"same" if ss else "other"} site : {rate(x)}  codes={x.lig_comp.nunique()}')
        if True:
            fmc=Lt[Lt.c_first_miss==1].groupby('lig_comp').agg(later=('miss','size'),later_missed=('miss','sum'),first_entry=('c_first_entry','first') if 'c_first_entry' in Lt else ('entry','first'))
            P('   first-missed codes with a later use (code: later missed/later):',', '.join(f'{c} {int(r.later_missed)}/{int(r.later)}' for c,r in fmc.iterrows()))
        a=Lt[(Lt.c_first_miss==1)&~Lt.same_lab]; b=Lt[(Lt.c_first_miss==0)&~Lt.same_lab]
        P(f'   other-lab uses: {rate(a)} vs {rate(b)} {compare(a,b)}'); RES.append(dict(set=tag,subset=name,design='same code, other lab',fm=rate(a),fr=rate(b),stat=compare(a,b)))
        a2=Lt[(Lt.c_first_miss==1)&~Lt.same_lab&~Lt.same_site]; b2=Lt[(Lt.c_first_miss==0)&~Lt.same_lab&~Lt.same_site]
        P(f'   other-lab AND other-site uses: {rate(a2)} vs {rate(b2)} {compare(a2,b2)}'); RES.append(dict(set=tag,subset=name,design='same code, other lab + other site',fm=rate(a2),fr=rate(b2),stat=compare(a2,b2)))
        rows=[]
        for lab,g in X.groupby('lab'):
            m=g[g.miss==1]
            if not len(m): continue
            t0=m.dep.min(); mc=set(m[m.dep==t0].lig_comp); mik=set(X[X.lig_comp.isin(mc)].ik14); msk=set(X[X.lig_comp.isin(mc)].skel)
            for _,r in g[g.dep>t0].iterrows():
                rel='same code' if r.lig_comp in mc else ('other code, same drug' if (r.ik14 in mik or r.skel in msk) else 'other code, different drug')
                rows.append(dict(lab=lab,entry=r.entry,code=r.lig_comp,missed_codes=','.join(sorted(mc)),relation=rel,miss=r.miss,site=r.site))
        Rb=pd.DataFrame(rows); base=X[(X.lab_prior_n>0)&(X.lab_prior_miss==0)]
        for rel in ['same code','other code, same drug','other code, different drug']:
            x=Rb[Rb.relation==rel] if len(Rb) else Rb
            P(f'   lab after first miss, {rel}: {rate(x) if len(x) else "0/0"} labs={x.lab.nunique() if len(x) else 0}')
            RES.append(dict(set=tag,subset=name,design=f'lab after first miss: {rel}',fm=rate(x) if len(x) else '0/0',stat=(';'.join(sorted(f'{r.missed_codes}->{r.code} {r.entry} miss={r.miss}' for r in x.itertuples())) if rel=='other code, same drug' and len(x) else '')))
        P(f'   baseline labs (earlier units, none missed): {rate(base)}'); RES.append(dict(set=tag,subset=name,design='baseline labs',fm=rate(base)))
        if name=='strict' and len(Rb): Rb.to_csv(f'{OUT}/{tag}_lab_after_miss_units.tsv',sep='\t',index=False)
        M=Lt.copy(); M['lab_prior_miss_other']=0
        for lab,g in X.groupby('lab'):
            for i,r in M[M.lab==lab].iterrows():
                prev=g[(g.dep<r.dep)&(g.lig_comp!=r.lig_comp)]; M.at[i,'lab_prior_miss_other']=int(prev.miss.sum()>0)
        M['code']=M.lig_comp.astype('category').cat.codes
        st='+ pdbj' if name!='drop PDBj' else ''
        try:
            m=smf.glm(f'miss ~ c_first_miss + lab_prior_miss_other + same_lab {st} + yearc',M,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':M.code})
            for k in m.params.index:
                if k=='Intercept': continue
                P(f'   joint GLM (cluster=code) {k}: {fmt(m,k)}'); RES.append(dict(set=tag,subset=name,design=f'joint GLM {k}',stat=fmt(m,k),n=len(M)))
        except Exception as ex: P('   joint GLM failed',ex)
# ---------------- validation on the original strict set ----------------
S0=prep(E[(E.cat=='drug-like')&(E.nun<=10)&E.bondlike&E.au_last.notna()])
P('########## VALIDATION: original strict set (bond-like distance test, misses at any element pair < 2.0 A)')
P('strict units',len(S0),'misses',int(S0.miss.sum()))
step3(S0,'orig'); step4(S0,'orig')
# ---------------- in-window ----------------
U=E[E.n_in>0].copy(); U['miss']=(U.n_in_unrec>0).astype(int)
U[['lig_comp','entry','cat','site','ps','dep','au_last','n_in','n_in_unrec','miss','miss_orig','nun','bondlike']].to_csv(f'{OUT}/window_units.tsv',sep='\t',index=False)
P('\n########## IN-WINDOW')
P('in-window contacts < 2.0 A',len(I),'unrecorded',int((~I.recorded).sum()),'; not in any bonded unit (left out):',int(sum(notin)),'unrecorded',int((~I.recorded[notin]).sum()))
P('bonded units',len(E),'-> units with >=1 in-window contact < 2.0 A',len(U),'misses',int(U.miss.sum()),'codes',U.lig_comp.nunique())
P('units dropped (no in-window contact)',len(E)-len(U),'of which missed',int(E[E.n_in==0].miss.sum()))
P('in-window units recorded in window but with an out-of-window miss',int(((U.miss==0)&(U.miss_orig==1)).sum()))
dl=U[U.cat=='drug-like']; P('drug-like in-window units',len(dl),'misses',int(dl.miss.sum()))
S1=prep(U[(U.cat=='drug-like')&(U.nun<=10)&U.au_last.notna()])
P('in-window strict units',len(S1),'misses',int(S1.miss.sum()),'codes',S1.lig_comp.nunique(),'labs',S1.lab.nunique())
x=S1.merge(S0[['lig_comp','entry']],on=['lig_comp','entry'],how='left',indicator=True)
P('   overlap with the original strict set:',int((x._merge=='both').sum()),'units; new to the in-window set:',int((x._merge=='left_only').sum()))
S1.to_csv(f'{OUT}/window_strict_units.tsv',sep='\t',index=False)
D1,x1,m1=step3(S1,'win'); step4(S1,'win')
S2=prep(U[(U.cat=='drug-like')&(U.nun<=10)&U.au_last.notna()&U.bondlike])
P('\n[sensitivity win+bondlike] units',len(S2),'misses',int(S2.miss.sum()))
step3(S2,'win+bondlike'); step4(S2,'win+bondlike')
# ---------------- how much of the in-window misses does code history account for? ----------------
P('\n########## ATTRIBUTION of in-window unit misses (code history recomputed on each set)')
AT=[]
for name,X0 in [('all in-window units',U),('drug-like in-window units',U[U.cat=='drug-like']),('in-window strict set',S1)]:
    D=X0.sort_values(['lig_comp','dep','rel','entry']).copy(); f=D.groupby('lig_comp').head(1).set_index('lig_comp')
    D['first_miss']=D.lig_comp.map(f.miss); D['first_dep']=D.lig_comp.map(f.dep)
    D['phase']=np.where(D.dep>D.first_dep,'later','first day')
    tot=int(D.miss.sum()); fd=D[D.phase=='first day']; lt=D[D.phase=='later']
    a=lt[lt.first_miss==1]; b=lt[lt.first_miss==0]
    x=lt.copy(); x['code']=x.lig_comp.astype('category').cat.codes
    m=smf.glm('miss ~ first_miss + C(site, Treatment("RCSB")) + yearc',x,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':x.code})
    a2=a.copy(); a2['first_miss']=0; ex=float(m.predict(a2).sum())
    bm=int(X0.n_in_unrec[X0.index.isin(a.index)].sum()) if False else None
    row=dict(set=name,units=len(D),misses=tot,first_day_units=len(fd),first_day_misses=int(fd.miss.sum()),
             later_units_first_missed=len(a),later_misses_first_missed=int(a.miss.sum()),codes_first_missed_with_later=a.lig_comp.nunique(),
             later_units_first_recorded=len(b),later_misses_first_recorded=int(b.miss.sum()),
             expected_if_first_recorded=round(ex,1),excess=round(a.miss.sum()-ex,1),
             share_later_first_missed=round(a.miss.sum()/tot,3),share_excess=round((a.miss.sum()-ex)/tot,3),
             first_miss_OR=fmt(m,'first_miss'))
    for ph,sub in [('first day',fd),('later, first recorded',b),('later, first missed',a)]:
        mm=sub[sub.miss==1]; row[f'{ph}: misses at PDBj or PDBc']=int(mm.site.isin(['PDBJ','PDBC']).sum())
    AT.append(row)
    P(f'## {name}: units {len(D)}, misses {tot}'); [P(f'   {k}: {v}') for k,v in row.items() if k!='set']
AT=pd.DataFrame(AT); AT.to_csv(f'{OUT}/window_attribution.tsv',sep='\t',index=False)
# bond-level share: in-window unrecorded contacts <2.0 A that sit in later uses of first-missed codes (all units)
D=U.sort_values(['lig_comp','dep','rel','entry']).copy(); f=D.groupby('lig_comp').head(1).set_index('lig_comp')
D['first_miss']=D.lig_comp.map(f.miss); D['first_dep']=D.lig_comp.map(f.dep)
lt=D[(D.dep>D.first_dep)&(D.first_miss==1)]
P(f'\nbond level (all in-window units): in-window unrecorded contacts < 2.0 A {int(U.n_in_unrec.sum())}; in later uses of first-missed codes {int(lt.n_in_unrec.sum())}; of which PEB {int(lt[lt.lig_comp=="PEB"].n_in_unrec.sum())}')
P('first-day (first deposition date of the code) in-window unrecorded contacts:',int(D[D.dep<=D.first_dep].n_in_unrec.sum()))
pd.DataFrame(RES).to_csv(f'{OUT}/window_results.tsv',sep='\t',index=False)
open(f'{OUT}/30_stdout.txt','w').write('\n'.join(L)); print('\n'.join(L))
