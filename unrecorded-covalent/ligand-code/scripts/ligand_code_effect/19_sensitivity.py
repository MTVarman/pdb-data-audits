"""Step 3 sensitivity for (a): many 'born unbonded' codes are commodity ligands (lipids, nucleotides, cholesterol) whose
<2.0 A contacts are likely clashes, not bonds. Restrict to (i) codes with <=10 unbonded uses in the archive ('specific'),
(ii) bond-like distances at the closest bond (SG 1.65-1.95 A; O/N 1.20-1.65 A), and (iii) both."""
import pandas as pd, numpy as np, statsmodels.formula.api as smf, statsmodels.api as sm
O='/workspace/ligand_code_effect'; OUT=f'{O}/out/step3'
E=pd.read_csv(f'{OUT}/entry_level_bonded.tsv',sep='\t',low_memory=False)
TL=pd.read_csv(f'{O}/archive/code_timelines.tsv',sep='\t',low_memory=False)
B=pd.read_csv(f'{O}/out/archive_bonds_lt2.tsv',sep='\t',low_memory=False)
nun=TL[TL.status=='unbonded'].groupby('lig_comp').size(); E['n_unbonded_uses']=E.lig_comp.map(nun).fillna(0).astype(int)
cb=B.sort_values('dist').groupby(['lig_comp','entry']).first()[['nuc_atom','dist']]
E=E.join(cb,on=['lig_comp','entry'])
E['bondlike']=np.where(E.nuc_atom=='SG',E.dist.between(1.65,1.95),E.dist.between(1.20,1.65))
L=[]; P=lambda *a: L.append(' '.join(map(str,a))); R=[]
D=E[E.cat=='drug-like']
for name,x0 in [('drug-like',D),('drug-like, specific codes (<=10 unbonded uses)',D[D.n_unbonded_uses<=10]),
                ('drug-like, bond-like distance',D[D.bondlike]),('drug-like, specific + bond-like',D[(D.n_unbonded_uses<=10)&D.bondlike])]:
    tot=int(x0.miss.sum()); first=int(x0[x0.phase=='birth_batch'].miss.sum())
    x=x0[(x0.phase=='later')&x0.birth_batch.isin(['recorded','missed','unbonded','linked_other'])].copy(); x['code']=x.lig_comp.astype('category').cat.codes
    m=smf.glm('miss ~ C(birth_batch, Treatment("recorded")) + C(site, Treatment("RCSB")) + yearc',x,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':x.code})
    ci=m.conf_int(); row=dict(subset=name,units=len(x0),misses=tot,first_use_misses=first)
    for g in ['unbonded','missed']:
        k=f'C(birth_batch, Treatment("recorded"))[T.{g}]'; s=x[x.birth_batch==g]; s2=s.copy(); s2['birth_batch']='recorded'; ex=float(m.predict(s2).sum())
        row.update({f'{g}_later_units':len(s),f'{g}_later_misses':int(s.miss.sum()),f'{g}_excess':round(s.miss.sum()-ex,1),f'{g}_OR':round(np.exp(m.params[k]),2),f'{g}_CI':f'{np.exp(ci.loc[k,0]):.2f}-{np.exp(ci.loc[k,1]):.2f}',f'{g}_p':f'{m.pvalues[k]:.2g}'})
    r=x[x.birth_batch=='recorded']; row.update(recorded_later_units=len(r),recorded_later_misses=int(r.miss.sum()))
    k='C(site, Treatment("RCSB"))[T.PDBJ]'; row.update(PDBJ_OR=round(np.exp(m.params[k]),2),PDBJ_p=f'{m.pvalues[k]:.2g}',codes=x.code.nunique())
    R.append(row)
R=pd.DataFrame(R); R.to_csv(f'{OUT}/attribution_sensitivity.tsv',sep='\t',index=False); P(R.T.to_string())
open(f'{OUT}/19_stdout.txt','w').write('\n'.join(L)); print('\n'.join(L))
