"""Main census numbers. Unit definitions:
contact = Mpro Cys145 SG -> nearest drug-like ligand carbon (non-polymer, not additive/thiol, not metal-containing) at 1.6-<2.4 A
          (one per Cys145 per entry chain; same convention as the Cys scan)
entry   = covalent complex if it has >=1 such contact; 'has_unrecorded_lt2' if >=1 unrecorded contact under 2.0 A;
          'all_unrecorded' if every contact is unrecorded."""
import pandas as pd, numpy as np, json, sys
sys.path.insert(0,'/workspace/cys_covalent_scan/scripts'); from fisher import fisher
D='/workspace/mpro_covalent_census'; O=f'{D}/out'
C=pd.read_csv(f'{O}/contacts_all.tsv',sep='\t',low_memory=False,dtype={'deposit_year':str,'release_year':str})
E=pd.read_csv(f'{O}/entries.tsv',sep='\t',dtype={'deposit_year':str,'release_year':str})
L=[]; P=lambda *a: (L.append(' '.join(str(x) for x in a)), print(*a))
c145=C[(C.nsp5_pos==145)&(C.nuc_atom=='SG')]
M=c145[(c145.lig_el=='C')&c145.category.isin(['other_ligand','named_drug'])&(c145.dist>=1.6)&(c145.dist<2.4)&c145.nearest_for_nuc].copy()
M['band']=np.where(M.dist<2.0,'<2.0','2.0-2.4')
M['fragment']=M.lig_heavy_atoms<=20
M['xchem']=M.xchem_like_id|M.deposit_group.notna()
M.to_csv(f'{O}/contacts_main_cys145_C_1.6-2.4.tsv',sep='\t',index=False)
P('Mpro entries scanned',len(E),'| sites',E.process_site.value_counts().to_dict(),'| methods',E.method.value_counts().to_dict())
P('entries with any Cys145 SG -> nonpolymer/peptide partner <=2.6 A',c145.entry.nunique())
P('MAIN contacts',len(M),'entries',M.entry.nunique())
def frac(g): return f'{int(g.sum())}/{len(g)} ({100*g.mean():.1f}%)'
for b in ['<2.0','2.0-2.4']:
    P(f'recorded {b}:',frac(M[M.band==b].recorded))
# entry level
g=M.groupby('entry')
EN=pd.DataFrame({'n_contacts':g.size(),'n_rec':g.recorded.sum(),'n_unrec_lt2':g.apply(lambda x:((~x.recorded)&(x.band=='<2.0')).sum()),
                 'min_dist':g.dist.min(),'ligs':g.lig_comp.apply(lambda x:';'.join(sorted(set(x)))),'drug':g.drug.apply(lambda x:';'.join(sorted(set(x.dropna()))))}).reset_index()
EN=EN.merge(E,on='entry',how='left')
EN['has_unrec_lt2']=EN.n_unrec_lt2>0; EN['all_unrec']=EN.n_rec==0; EN['has_unrec_any']=EN.n_rec<EN.n_contacts
EN['fragment_entry']=EN.entry.isin(M[M.fragment].entry); EN['xchem']=EN.xchem_like_id|EN.deposit_group.notna()
EN.to_csv(f'{O}/entries_covalent.tsv',sep='\t',index=False)
P('covalent complex entries',len(EN),'| with >=1 unrecorded <2.0 A contact',frac(EN.has_unrecorded_lt2 if 'has_unrecorded_lt2' in EN else EN.has_unrec_lt2),
  '| with >=1 unrecorded contact (<2.4)',frac(EN.has_unrec_any),'| no contact recorded at all',frac(EN.all_unrec))
# by site
rows=[]
for s,x in M.groupby('process_site'):
    for b in ['<2.0','2.0-2.4']:
        y=x[x.band==b]; rows.append({'site':s,'band':b,'contacts':len(y),'recorded':int(y.recorded.sum()),'unrecorded':int((~y.recorded).sum()),'pct_recorded':round(100*y.recorded.mean(),1) if len(y) else None})
    en=EN[EN.process_site==s]
    rows.append({'site':s,'band':'ENTRIES','contacts':len(en),'recorded':int((~en.has_unrec_lt2).sum()),'unrecorded':int(en.has_unrec_lt2.sum()),'pct_recorded':round(100*(~en.has_unrec_lt2).mean(),1)})
S=pd.DataFrame(rows); S.to_csv(f'{O}/by_site.tsv',sep='\t',index=False); P(S.to_string())
lt=M[M.band=='<2.0']; pj=lt[lt.process_site=='PDBJ']; ot=lt[lt.process_site!='PDBJ']
a,b_,c_,d_=int((~pj.recorded).sum()),int(pj.recorded.sum()),int((~ot.recorded).sum()),int(ot.recorded.sum())
P('Fisher <2.0 contacts PDBj unrec/rec vs others',a,b_,c_,d_,'OR=%.1f p=%.2g'%fisher(a,b_,c_,d_))
ep=EN[EN.process_site=='PDBJ']; eo=EN[EN.process_site!='PDBJ']
a,b_,c_,d_=int(ep.has_unrec_lt2.sum()),int((~ep.has_unrec_lt2).sum()),int(eo.has_unrec_lt2.sum()),int((~eo.has_unrec_lt2).sum())
P('Fisher entries PDBj',a,b_,'others',c_,d_,'OR=%.1f p=%.2g'%fisher(a,b_,c_,d_))
# by year (deposition) x site
Y=lt.assign(site2=np.where(lt.process_site=='PDBJ','PDBJ','other')).groupby(['deposit_year','site2']).recorded.agg(['sum','count']).unstack()
Y.to_csv(f'{O}/by_deposit_year_site_lt2.tsv',sep='\t'); P('by deposit year (<2.0 contacts; sum=recorded,count=all)\n'+Y.to_string())
Yr=lt.groupby('release_year').recorded.agg(['sum','count']); Yr.to_csv(f'{O}/by_release_year_lt2.tsv',sep='\t'); P('by release year\n'+Yr.to_string())
Ye=EN.assign(site2=np.where(EN.process_site=='PDBJ','PDBJ','other')).groupby(['deposit_year','site2']).has_unrec_lt2.agg(['sum','count']).unstack()
Ye.to_csv(f'{O}/by_deposit_year_site_entries.tsv',sep='\t'); P('entries with unrecorded <2.0 by deposit year\n'+Ye.to_string())
# method
P('by method <2.0', lt.groupby('method').recorded.agg(['sum','count']).to_dict())
# fragments / xchem
for nm,msk in [('fragment(<=20 heavy)',lt.fragment),('drug-like(>20)',~lt.fragment),('XChem-like/group dep',lt.xchem),('not group dep',~lt.xchem)]:
    P('<2.0',nm,frac(lt[msk].recorded),'entries',lt[msk].entry.nunique())
P('group deposition IDs among covalent entries:',EN[EN.deposit_group.notna()].deposit_group.value_counts().to_dict())
P('XChem-like ID prefixes among all Mpro entries',E[E.xchem_like_id].entry.str[:2].value_counts().to_dict(),'with group dep',E.deposit_group.notna().sum())
# per drug
DR=['4WI','ZGW','OW1','7YY','K36','UED','U5G','SV6','V2M','80I','YDL','7ON','A1AA0','NNA','FHR','ZQB','VR4','A1EN0','PQL','HUR']
import importlib.util; _s=importlib.util.spec_from_file_location('t',f'{D}/scripts/06_tables.py')
DRN={'4WI':'nirmatrelvir (bound, thioimidate)','ZGW':'nirmatrelvir (nitrile)','OW1':'nirmatrelvir (reacted form)','7YY':'ensitrelvir (non-covalent control)','K36':'GC376 (bisulfite adduct)','UED':'GC376/GC373 (bound form)','U5G':'boceprevir (bound)','SV6':'telaprevir (bound)','V2M':'PF-00835231 (bound)','80I':'lufotrelvir (PF-07304814)','YDL':'ibuzatrelvir (bound)','7ON':'leritrelvir (bound)','A1AA0':'zevotrelvir/EDP-235 (bound)','NNA':'narlaprevir (bound)','FHR':'bofutrelvir/11a','ZQB':'pomotrelvir (bound)','VR4':'GC376 alpha-ketoamide analog','A1EN0':'leritrelvir (unreacted-form CCD)','PQL':'simnotrelvir (nitrile CCD)','HUR':'YH-53'}
rows=[]
allpairs=C[(C.nsp5_pos==145)&(C.nuc_atom=='SG')]
for comp in sorted(set(M.lig_comp)|set(DR)):
    ents=E[E.nonpoly_comps.fillna('').str.split(';').apply(lambda s:comp in s)]
    if comp not in DR and len(ents)<5: continue
    x=M[M.lig_comp==comp]
    for s,es in ents.groupby('process_site'):
        xs=x[x.process_site==s]; xl=xs[xs.band=='<2.0']
        rows.append({'comp':comp,'drug':DRN.get(comp,''),'site':s,'entries_with_comp':len(es),
            'entries_with_cys145_contact':xs.entry.nunique(),'contacts_lt2':len(xl),'recorded_lt2':int(xl.recorded.sum()),
            'contacts_2.0-2.4':int((xs.band=='2.0-2.4').sum()),'recorded_2.0-2.4':int(xs[xs.band=='2.0-2.4'].recorded.sum()),
            'entries_with_unrec_lt2':int(EN[(EN.process_site==s)&EN.entry.isin(xl[~xl.recorded].entry)].shape[0]),
            'entries_covale_row_to_comp':int(es.covale_pos145_partners.fillna('').str.split(';').apply(lambda t:comp in t).sum()),
            'dist_range':f'{xs.dist.min():.3f}-{xs.dist.max():.3f}' if len(xs) else ''})
DT=pd.DataFrame(rows); DT.to_csv(f'{O}/per_drug_by_site.tsv',sep='\t',index=False)
P(DT.to_string())
# ensitrelvir control: min SG distance
P('ensitrelvir 7YY: entries',E.nonpoly_comps.fillna('').str.contains(r'(?:^|;)7YY(?:;|$)').sum(),'any SG contact <=2.6:',len(allpairs[allpairs.lig_comp=='7YY']),'covale rows to 7YY:',E.covale_pos145_partners.fillna('').str.contains('7YY').sum())
# comparison to cys-scan baseline
P('Cys-scan baseline (drug-like C, <2.0, all Cys sites, 24,144-entry sample): 197/4094 unrecorded (4.8%); PDBj 97/526 (18.4%); others 100/3568 (2.8%); entries: PDBj 47/324 (14.5%), others 58/2118 (2.7%)')
P('Mpro <2.0 unrecorded',int((~lt.recorded).sum()),'/',len(lt),f'({100*(~lt.recorded).mean():.1f}%)','PDBj',int((~pj.recorded).sum()),'/',len(pj),'others',int((~ot.recorded).sum()),'/',len(ot))
# overlap with cys scan
U=pd.read_csv('/workspace/cys_covalent_scan/out/unrecorded_under_2.0A_druglike.tsv',sep='\t')
cs=pd.read_csv('/workspace/cys_covalent_scan/out/entries_scanned.tsv',sep='\t')
un_m=set(lt[~lt.recorded].entry)
P('Mpro entries with unrecorded <2.0:',len(un_m),'| also in Cys-scan sample:',len(un_m&set(cs.entry)),'| already among Cys-scan 197 misses:',len(un_m&set(U.entry)),'| new:',len(un_m-set(U.entry)))
P('Covalent Mpro entries in Cys-scan sample',len(set(EN.entry)&set(cs.entry)),'of',len(EN))
open(f'{O}/summary_numbers.txt','w').write('\n'.join(L)+'\n')
