"""Archive-wide bonds and per-code timelines.
Bond = (entry, nucleophile residue, ligand residue). Recorded = a covale row joins them. Unrecorded = no covale row
but the file's own _pdbx_validate_close_contact lists a nucleophile atom within 2.0 A of a ligand atom of element
C, B, P, S, Si, Se or As (element from the CCD). Nucleophiles: CYS SG, SER OG, LYS NZ, HIS NE2/ND1.
Per (code, entry) status: recorded / missed / mixed / unbonded (no nucleophile bond <2.0 A and no covale row) /
linked_other (unbonded to a nucleophile but covale-linked to another polymer atom)."""
import json, pandas as pd, numpy as np, collections, sys
O='/workspace/ligand_code_effect'; A=f'{O}/archive'
PARSED=sys.argv[1] if len(sys.argv)>1 else f'{A}/parsed.jsonl'
C=json.load(open(f'{O}/out/ccd_parsed_archive.json'))
MAIN={('CYS','SG'),('SER','OG'),('LYS','NZ'),('HIS','NE2'),('HIS','ND1')}
EL={'C','B','P','S','Si','Se','As'}
def elem(c,a):
    r=C.get(c); 
    if r and a in r['el']: return r['el'][a]
    return None
meta=[]; bonds=[]; usage=collections.defaultdict(list)
for l in open(PARSED):
    r=json.loads(l)
    if 'err' in r: continue
    e=r['e']; meta.append({k:r.get(k) for k in ('e','ps','ds','dep','rel','au_last','au_first','meth','grp','hpm','has_cc_table')})
    for c in r['np']: usage[c].append(e)
    rec={}
    for x in r['cov']:
        if (x[0],x[1]) not in MAIN: continue
        k=(x[0],x[2],x[3],x[4] or '',x[5],x[7],x[8],x[9] or '')
        try: d=float(x[10])
        except: d=None
        el=elem(x[5],x[6])
        rec.setdefault(k,[]).append((d,x[1],x[6],el))
    unrec={}
    for x in r['cc']:
        if (x[0],x[1]) not in MAIN or x[10]>=2.0: continue
        el=elem(x[5],x[6])
        if el not in EL: continue
        k=(x[0],x[2],x[3],x[4] or '',x[5],x[7],x[8],x[9] or '')
        if k in rec: continue
        if k not in unrec or x[10]<unrec[k][0]: unrec[k]=(x[10],x[1],x[6],el)
    nuc_rec={(k[0],k[1],k[2],k[3]) for k in rec}
    for k,v in rec.items():
        ds=[z for z in v if z[0] is not None]; best=min(ds) if ds else v[0]
        bonds.append(dict(entry=e,nuc_res=k[0],nuc_chain=k[1],nuc_seq=k[2],nuc_ins=k[3],lig_comp=k[4],lig_chain=k[5],lig_seq=k[6],lig_ins=k[7],
                          dist=best[0],nuc_atom=best[1],lig_atom=best[2],lig_el=best[3],recorded=True,nuc_cov_other=False))
    for k,v in unrec.items():
        bonds.append(dict(entry=e,nuc_res=k[0],nuc_chain=k[1],nuc_seq=k[2],nuc_ins=k[3],lig_comp=k[4],lig_chain=k[5],lig_seq=k[6],lig_ins=k[7],
                          dist=v[0],nuc_atom=v[1],lig_atom=v[2],lig_el=v[3],recorded=False,nuc_cov_other=(k[0],k[1],k[2],k[3]) in nuc_rec))
M=pd.DataFrame(meta).set_index('e'); B=pd.DataFrame(bonds)
B=B.join(M[['ps','dep','rel','au_last','grp','meth']],on='entry')
B['lt2']=B.dist.notna()&(B.dist<2.0)
B.to_csv(f'{A}/bonds_all.tsv',sep='\t',index=False); M.to_csv(f'{A}/entries_meta.tsv',sep='\t')
link_other=collections.defaultdict(set)
for l in open(PARSED):
    r=json.loads(l)
    for c in r.get('link_any',[]): link_other[c].add(r['e'])
# target codes: >=1 bond < 2.0 A (recorded or not) to a main nucleophile
T=sorted(set(B[B.lt2].lig_comp))
L=[]
g=B.groupby(['lig_comp','entry'])
agg=g.agg(n_lt2_rec=('recorded',lambda x:int((x&B.loc[x.index,'lt2']).sum())),n_lt2_unrec=('recorded',lambda x:int(((~x)&B.loc[x.index,'lt2']).sum())),
          n_rec_any=('recorded','sum'),dmin=('dist','min')).reset_index()
aggd={(a.lig_comp,a.entry):a for a in agg.itertuples()}
for c in T:
    for e in usage.get(c,[]):
        a=aggd.get((c,e))
        if a is None: st='linked_other' if e in link_other[c] else 'unbonded'; nr=nu=0; dm=None; ra=0
        else:
            nr,nu,ra,dm=a.n_lt2_rec,a.n_lt2_unrec,a.n_rec_any,a.dmin
            st='mixed' if (nu>0 and ra>0) else 'missed' if nu>0 else 'recorded' if ra>0 else ('linked_other' if e in link_other[c] else 'unbonded')
        m=M.loc[e]
        L.append(dict(lig_comp=c,entry=e,status=st,n_lt2_rec=nr,n_lt2_unrec=nu,n_rec_any=ra,dmin=dm,ps=m.ps,dep=m.dep,rel=m.rel,au_last=m.au_last,grp=m.grp,meth=m.meth,has_cc_table=m.has_cc_table))
TL=pd.DataFrame(L).sort_values(['lig_comp','dep','rel','entry'])
TL['order']=TL.groupby('lig_comp').cumcount()
TL.to_csv(f'{A}/code_timelines.tsv',sep='\t',index=False)
print('entries',len(M),'bonds',len(B),'lt2',int(B.lt2.sum()),'unrec lt2',int((B.lt2&~B.recorded).sum()),'target codes',len(T),'timeline rows',len(TL))
print(TL.status.value_counts().to_dict())
