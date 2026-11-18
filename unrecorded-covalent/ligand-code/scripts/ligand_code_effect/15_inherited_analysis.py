"""(a) misses explained by inheriting a code first used for an unbonded model; (b) same lab + same day splits across
codes; (c) codes whose recording starts partway through their history. Site and year controls throughout."""
import json, pandas as pd, numpy as np, statsmodels.formula.api as smf, statsmodels.api as sm, sys
from scipy.stats import fisher_exact
O='/workspace/ligand_code_effect'; A=f'{O}/archive'; OUT=f'{O}/out'
CCD=json.load(open(f'{OUT}/ccd_parsed_archive.json'))
B=pd.read_csv(f'{A}/bonds_all.tsv',sep='\t',dtype={'nuc_seq':str,'lig_seq':str},low_memory=False)
TL=pd.read_csv(f'{A}/code_timelines.tsv',sep='\t',low_memory=False)
HEME={'HEM','HEC','HEA','HEB','HEO','HDD','HDM','HAS','HEV','HE5','DHE','SRM','COH','CCH','HCO','HNI','ZNH','1CP','FEC','HIF','VEA','VER','HMG','HE6','MHM','PNI','HKL','6CO','6CQ','CLN','ISW','MD9','HEG','HES','HP5','UFE','T2Z','HCB','HNN','HME','HRM','HFM','HEW'}
BILIN={'BLA','LBV','LBW','CYC','PEB','PUB','3Q8','EL5','BLR','RBF','PCB','URO','DBV','15P','BV','PVN','BPB'}
FLAVIN={'FAD','FMN','FDA','FAE','FMA','6FA','8FM','FNR','FAS','FSH','FAO','C1L','LFN'}
ADD={'BME','DTT','DTU','DTV','GSH','CYS','2ME','MRY','CME','TCE','DTD','MPO','CSD','DCY','GTT','GDS','GSF','GSN','MTL','ACM','MPD','EDO','GOL','PEG','PGE','1PE','PG4','ACT','FMT','DMS','EPE','MES','TRS','BU3','CIT','SO4','PO4','IMD','NHE','MLI','IPA','EOH','MOH','BCN','B3P','ACY','AKG','TLA','TAR','PGO','P6G','SIN','12P','CAC','ACE','SCN','CO3','NO3','BCT','OXY','UNX','UNL','HED','DIO','ETF','FLC','BEZ','PYR','ETX'}
LIPID={'PLM','MYR','STE','OLA','PAM','HXA','DAO','DKA','UND','LAU','PLC','COA','ACO','SCA','MLC','CAA','HXC','DCC','MCA','OCA','DCR','LPP','PEE','PTY','CDL'}
OTHERCOF={'PLP','PMP','RET','NAD','NAP','NDP','NAI','NAJ','ATP','ADP','AMP','ANP','GTP','GDP','GMP','UMP','UDP','CTP','TPP','TDP','SAM','SAH','PQQ','TPQ','LPA','BTN','FFO','THF','MTX','PPV','FOL','H4B','MOO','MGD','PCD','SF4','FES','F3S'}
METALS={'Fe','Zn','Cu','Co','Ni','Mn','Mg','Ca','Hg','Pt','Ru','Rh','Pd','Au','Ag','Cd','Mo','W','V','Cr','Os','Ir','Re','Ga','Al','Sn','Pb','Na','K','Bi','Sb','Te','Tl','Gd','Eu','Tb','Yb','Lu','U'}
def cat(c):
    if c in HEME: return 'heme'
    if c in BILIN: return 'bilin'
    if c in FLAVIN: return 'flavin'
    if c in ADD: return 'additive'
    if c in LIPID: return 'lipid/acyl'
    if c in OTHERCOF: return 'cofactor'
    r=CCD.get(c)
    if r is None: return 'unknown'
    els=set(r['el'].values())
    if els & METALS: return 'metal-containing'
    if r['n_heavy']<6: return 'small(<6 heavy)'
    return 'drug-like'
codes=sorted(set(TL.lig_comp)); CAT={c:cat(c) for c in codes}
TL['cat']=TL.lig_comp.map(CAT); B['cat']=B.lig_comp.map(CAT)
TL['dep']=pd.to_datetime(TL.dep,errors='coerce')
first=TL.sort_values(['lig_comp','dep','rel','entry']).groupby('lig_comp').head(1).set_index('lig_comp')
birth=first.status.rename('birth'); bdate=first.dep.rename('birth_dep')
# birth batch: all entries sharing the first deposition date; status of the batch = combined
bb=TL.merge(bdate,left_on='lig_comp',right_index=True); bb=bb[bb.dep==bb.birth_dep]
def comb(s):
    s=set(s)
    if s&{'recorded','mixed'} and s&{'missed','mixed'}: return 'mixed'
    if 'recorded' in s or 'mixed' in s: return 'recorded'
    if 'missed' in s: return 'missed'
    if 'linked_other' in s: return 'linked_other'
    return 'unbonded'
bbs=bb.groupby('lig_comp').status.agg(comb).rename('birth_batch')
codesdf=pd.DataFrame({'cat':pd.Series(CAT)}).join(birth).join(bbs).join(bdate).join(first[['entry','ps']].rename(columns={'entry':'first_entry','ps':'first_site'}))
codesdf['ccd_defining_entry']=[ (CCD.get(c) or {}).get('pdbx_model_coordinates_db_code') for c in codesdf.index]
codesdf['ccd_initial_date']=[ (CCD.get(c) or {}).get('pdbx_initial_date') for c in codesdf.index]
codesdf['n_entries']=TL.groupby('lig_comp').size()
codesdf['skeleton']=[(CCD.get(c) or {}).get('skeleton') for c in codesdf.index]
codesdf.to_csv(f'{OUT}/archive_codes.tsv',sep='\t')
# bonds < 2.0 A in later entries (deposited after the birth batch)
b=B[B.lt2].copy(); b['dep']=pd.to_datetime(b.dep,errors='coerce')
b=b.join(codesdf[['birth_batch','birth_dep']],on='lig_comp')
b['phase']=np.where(b.dep>b.birth_dep,'later','birth_batch')
b['miss']=(~b.recorded).astype(int); b['year']=b.dep.dt.year; b['yearc']=b.year-2015
b['site']=b.ps.fillna('NA').replace({'BNL':'RCSB','NDB':'RCSB'}); b['site']=b.site.where(b.site.isin(['RCSB','PDBE','PDBJ','PDBC']),'other')
b['nuc']=b.nuc_res
b.to_csv(f'{OUT}/archive_bonds_lt2.tsv',sep='\t',index=False)
res={}
def summ(x,tag):
    t=x.groupby(['birth_batch']).agg(n=('miss','size'),miss=('miss','sum'),codes=('lig_comp','nunique'),entries=('entry','nunique')); t['rate']=(t.miss/t.n).round(4)
    t2=x.groupby(['birth_batch',x.site=='PDBJ']).agg(n=('miss','size'),miss=('miss','sum')).unstack()
    print('\n==',tag,'\n',t.to_string(),'\n',t2.to_string()); return t
L=[]; R=[]
for tag,x in [('all codes',b[b.phase=='later']),('drug-like codes',b[(b.phase=='later')&(b.cat=='drug-like')])]:
    t=summ(x,tag); t['subset']=tag; L.append(t.reset_index())
    x=x[x.birth_batch.isin(['recorded','missed','mixed','unbonded','linked_other'])].copy()
    x['code']=x.lig_comp.astype('category').cat.codes
    m=smf.glm('miss ~ C(birth_batch, Treatment("recorded")) + C(site, Treatment("RCSB")) + yearc + C(nuc, Treatment("CYS"))',x,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':x.code})
    ci=m.conf_int()
    for k in m.params.index: R.append(dict(subset=tag,term=k,OR=round(np.exp(m.params[k]),2),ci_lo=round(np.exp(ci.loc[k,0]),2),ci_hi=round(np.exp(ci.loc[k,1]),2),p=f'{m.pvalues[k]:.2g}',n=int(m.nobs)))
    # attribution: expected misses for born-unbonded / born-missed bonds had their code been born recorded
    for grp in ['unbonded','missed','linked_other']:
        s=x[x.birth_batch==grp].copy()
        if len(s)==0: continue
        s2=s.copy(); s2['birth_batch']='recorded'
        exp=float(m.predict(s2).sum()); obs=int(s.miss.sum())
        res[(tag,grp)]=dict(subset=tag,birth=grp,bonds=len(s),observed_misses=obs,expected_if_born_recorded=round(exp,1),excess=round(obs-exp,1),
                            all_later_misses=int(x.miss.sum()),share_of_all_later_misses=round(obs/x.miss.sum(),3),excess_share=round((obs-exp)/x.miss.sum(),3))
pd.concat(L).to_csv(f'{OUT}/archive_birth_summary.tsv',sep='\t',index=False)
pd.DataFrame(R).to_csv(f'{OUT}/archive_birth_logistic.tsv',sep='\t',index=False); print(pd.DataFrame(R).to_string())
pd.DataFrame(res.values()).to_csv(f'{OUT}/archive_attribution.tsv',sep='\t',index=False); print(pd.DataFrame(res.values()).to_string())
