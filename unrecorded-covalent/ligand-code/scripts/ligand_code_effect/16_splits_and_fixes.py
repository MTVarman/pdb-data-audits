"""(b) same lab + same day splits across codes; (c) codes whose recording starts partway through their history."""
import json, pandas as pd, numpy as np, collections
O='/workspace/ligand_code_effect'; A=f'{O}/archive'; OUT=f'{O}/out'
CCD=json.load(open(f'{OUT}/ccd_parsed_archive.json'))
TL=pd.read_csv(f'{A}/code_timelines.tsv',sep='\t',low_memory=False)
CD=pd.read_csv(f'{OUT}/archive_codes.tsv',sep='\t',index_col=0)
TL=TL.join(CD[['cat','skeleton','birth_batch']],on='lig_comp')
bt=TL[TL.status.isin(['recorded','missed','mixed'])].copy()
# ---------- (b) ----------
bt['lab_day']=bt.au_last.fillna('?')+' | '+bt.dep.astype(str)
rows=[]; same_code=[]
for k,g in bt[bt.au_last.notna()].groupby('lab_day'):
    per=g.groupby('lig_comp').status.agg(lambda s:'recorded' if set(s)=={'recorded'} else 'missed' if set(s)=={'missed'} else 'split')
    for c,s in per.items():
        if s=='split' or (g[g.lig_comp==c].status=='mixed').any():
            x=g[g.lig_comp==c]; same_code.append(dict(lab_day=k,lig_comp=c,cat=x.cat.iloc[0],entries_recorded=','.join(x[x.status=='recorded'].entry),entries_missed=','.join(x[x.status=='missed'].entry),entries_mixed=','.join(x[x.status=='mixed'].entry),site=','.join(sorted(set(x.ps.astype(str))))))
    rec=[c for c,s in per.items() if s=='recorded']; mis=[c for c,s in per.items() if s=='missed']
    for c1 in rec:
        for c2 in mis:
            s1=CCD.get(c1,{}).get('skeleton'); s2=CCD.get(c2,{}).get('skeleton')
            x1=g[g.lig_comp==c1]; x2=g[g.lig_comp==c2]
            rows.append(dict(lab_day=k,code_recorded=c1,code_missed=c2,same_skeleton=bool(s1 and s1==s2),cat1=x1.cat.iloc[0],cat2=x2.cat.iloc[0],
                entries_recorded=','.join(x1.entry),entries_missed=','.join(x2.entry),n_rec=len(x1),n_mis=len(x2),site=','.join(sorted(set(g.ps.astype(str)))),
                name1=(CCD.get(c1,{}).get('name') or '')[:50],name2=(CCD.get(c2,{}).get('name') or '')[:50]))
S=pd.DataFrame(rows); S.to_csv(f'{OUT}/archive_labday_code_splits.tsv',sep='\t',index=False)
SC=pd.DataFrame(same_code); SC.to_csv(f'{OUT}/archive_labday_same_code_splits.tsv',sep='\t',index=False)
print('(b) lab-day groups with bonded entries:',bt[bt.au_last.notna()].lab_day.nunique())
print('    discordant code pairs (one code all recorded, another all missed, same lab+day):',len(S),'in',S.lab_day.nunique() if len(S) else 0,'groups')
if len(S): print('    same skeleton (same drug):',int(S.same_skeleton.sum()),'pairs in',S[S.same_skeleton].lab_day.nunique(),'groups'); print(S[S.same_skeleton].to_string())
print('    same code split within a lab-day group:',len(SC))
# ---------- (c) ----------
H=[]
for c,g in bt.sort_values(['dep','rel','entry']).groupby('lig_comp'):
    if len(g)<2: continue
    st=list(g.status); d=list(g.dep); s=list(g.ps)
    isrec=[x=='recorded' for x in st]; ismis=[x in ('missed','mixed') for x in st]
    if all(isrec): typ='always recorded'
    elif all(x=='missed' for x in st): typ='never recorded'
    else:
        fm=ismis.index(True) if any(ismis) else None; fr=isrec.index(True) if any(isrec) else None
        if fr is None: typ='never fully recorded'
        elif fm is not None and fm<fr:
            typ='fixed (recording starts later)'
        else: typ='recorded first, misses later'
    rec=dict(lig_comp=c,cat=g.cat.iloc[0],birth_batch=g.birth_batch.iloc[0],n_bonded_entries=len(g),n_recorded=sum(isrec),n_missed=sum(ismis),type=typ,
             first_bonded=d[0],first_bonded_status=st[0],first_bonded_site=s[0])
    if typ.startswith('fixed'):
        fr=isrec.index(True); pre=g.iloc[:fr]; post=g.iloc[fr+1:]
        rec.update(fix_entry=g.entry.iloc[fr],fix_dep=d[fr],fix_site=s[fr],misses_before_fix=int(sum(ismis[:fr])),pre_fix_sites=','.join(sorted(set(map(str,pre.ps)))),
                   pre_fix_entries=','.join(pre.entry)[:120],misses_after_fix=int(sum(ismis[fr+1:])),recorded_after_fix=int(sum(isrec[fr+1:])),
                   post_fix_miss_sites=','.join(sorted(set(map(str,post[post.status!='recorded'].ps)))))
        # CCD audit events between the last pre-fix entry and the fixing entry
        au=(CCD.get(c) or {}).get('audit') or []
        lo=str(pre.dep.iloc[-1]); hi=str(d[fr])
        rec['ccd_events_between']=';'.join(f'{a[0]}@{a[1]}' for a in au if a[1] and lo<=a[1]<=hi)
        rec['ccd_events_all']=';'.join(f'{a[0]}@{a[1]}' for a in au)
    H.append(rec)
H=pd.DataFrame(H); H.to_csv(f'{OUT}/archive_code_histories.tsv',sep='\t',index=False)
print('\n(c) codes with >=2 bonded entries:',len(H)); print(pd.crosstab(H.type,H.cat,margins=True).to_string())
F=H[H.type.str.startswith('fixed')]
if len(F):
    F=F.copy(); F['fix_year']=F.fix_dep.str[:4]
    print('fix year:',F.fix_year.value_counts().sort_index().to_dict())
    print('fix site:',F.fix_site.value_counts().to_dict(),' pre-fix sites:',F.pre_fix_sites.value_counts().head(10).to_dict())
    print('site changed at fix (fix site not among pre-fix sites):',int(sum(fs not in ps.split(',') for fs,ps in zip(F.fix_site.astype(str),F.pre_fix_sites))),'/',len(F))
    print('misses after fix: codes',int((F.misses_after_fix>0).sum()),' total entries',int(F.misses_after_fix.sum()),' recorded after fix',int(F.recorded_after_fix.sum()))
    print('CCD event between last miss and fix:',int((F.ccd_events_between.fillna('')!='').sum()),'/',len(F), collections.Counter(e.split('@')[0] for x in F.ccd_events_between.fillna('') for e in x.split(';') if e).most_common(8))
