"""Step 5 summary numbers + figure out/copy_forward_test.png"""
import json, re, pandas as pd, numpy as np, statsmodels.formula.api as smf, statsmodels.api as sm, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O='/workspace/ligand_code_effect'; OUT=f'{O}/out'; S5=f'{OUT}/step5'
CC=json.load(open(f'{OUT}/ccd_parsed_archive.json'))
R=pd.read_csv(f'{S5}/recorded_templates.tsv',sep='\t',low_memory=False); R['dep']=pd.to_datetime(R.dep)
U=pd.read_csv(f'{S5}/unit_templates.tsv',sep='\t',low_memory=False)
TM=pd.read_csv(f'{S5}/template_match.tsv',sep='\t')
L=[];P=lambda *a:L.append(' '.join(map(str,a)))
def lr(c,a):
    r=CC.get(c) or {}; nb=r.get('nb',{}); lv=set(r.get('leaving_atoms',[]))
    return (any(n in lv for n,o in nb.get(a,[])) if a in nb else None)
R['ccd_link_ready']=[lr(c,a) for c,a in zip(R.code,R.lig_atom)]
ct=pd.crosstab(R.leaving,R.ccd_link_ready.astype(str)); P('leaving flag vs CCD link-ready atom:\n'+ct.to_string())
P('codes with >1 leaving-flag value:',int((R.groupby('code').leaving.nunique()>1).sum()),'of',int((R.groupby('code').size()>1).sum()),'codes with >=2 recorded units')
f=R.sort_values(['dep','release','entry']).groupby('code').head(1).set_index('code')
Lt=R[R.dep>R.code.map(f.dep)].copy()
P('later recorded units',len(Lt),'codes',Lt.code.nunique(),'; atom pair differs from first recorded:',int((Lt.pair!=Lt.code.map(f.pair)).sum()),'units in',Lt[Lt.pair!=Lt.code.map(f.pair)].code.nunique(),'codes')
odd={'prot_first':lambda x:str(x)!='True','value_order':lambda x:isinstance(x,str) and x not in ('?',''),'details':lambda x:isinstance(x,str) and x not in ('?','')}
OD=[]
for fld,isodd in odd.items():
    oc=[c for c in f.index if isodd(f.loc[c,fld])]
    x=Lt[Lt.code.isin(oc)].copy(); x['rep']=[str(u[fld])==str(f.loc[u.code,fld]) for _,u in x.iterrows()]; x['same_lab']=[u.lab==f.loc[u.code,'lab'] for _,u in x.iterrows()]
    base=Lt[~Lt.code.isin(oc)]; bodd=int(sum(isodd(v) for v in base[fld]))
    d=dict(field=fld,codes_odd_first=len(oc),later_units=len(x),same_lab_repeat=f"{int(x[x.same_lab].rep.sum())}/{int(x.same_lab.sum())}",other_lab_repeat=f"{int(x[~x.same_lab].rep.sum())}/{int((~x.same_lab).sum())}",
           baseline_odd_other_codes=f"{bodd}/{len(base)}")
    OD.append(d)
OD=pd.DataFrame(OD); OD.to_csv(f'{S5}/oddity_same_vs_other_lab.tsv',sep='\t',index=False); P(OD.to_string(index=False))
# oddities: link atom not the shortest listed, atom name not in current CCD
P('link atom differs from a closer listed atom of the same ligand:',int(R.shorter_other_atom.sum()),'units,',R[R.shorter_other_atom].code.nunique(),'codes; first-unit odd codes',int(f.shorter_other_atom.sum()),
  '; their later units odd',int(Lt[Lt.code.isin(f[f.shorter_other_atom].index)].shorter_other_atom.sum()),'/',int(Lt.code.isin(f[f.shorter_other_atom].index).sum()))
P('link atom name absent from current CCD:',int(R.atom_not_in_ccd.fillna(False).astype(bool).sum()),'units:',R[R.atom_not_in_ccd.fillna(False).astype(bool)].groupby('code').entry.apply(lambda s:','.join(s)).to_dict())
# struct_conn revised
for m in [0,1]:
    x=U[U.miss==m]; P(f'miss={m}: struct_conn changed in a later revision {int(x.struct_conn_revised.sum())}/{len(x)}')
open(f'{S5}/28_stdout.txt','w').write('\n'.join(L)); print('\n'.join(L))
# ---------- figure ----------
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5}); GREEN,RED,GREY,BLUE='#1b7f5b','#c0392b','#7f8c8d','#2c6fbb'
fig=plt.figure(figsize=(13.5,9)); gs=fig.add_gridspec(2,3,hspace=0.6,wspace=0.45)
# A leaving flag vs CCD
ax=fig.add_subplot(gs[0,0])
a=int(((R.leaving=='one')&(R.ccd_link_ready==True)).sum()); b=int(((R.leaving=='one')&(R.ccd_link_ready==False)).sum())
c=int(((R.leaving=='none')&(R.ccd_link_ready==True)).sum()); d=int(((R.leaving=='none')&(R.ccd_link_ready==False)).sum())
M=np.array([[a,b],[c,d]]); ax.imshow(np.log10(M+1),cmap='Greens')
for i in range(2):
    for j in range(2): ax.text(j,i,f'{M[i,j]:,}',ha='center',va='center',fontsize=11,color='k')
ax.set_xticks([0,1]); ax.set_xticklabels(['CCD: leaving atom\non link atom','CCD: none']); ax.set_yticks([0,1]); ax.set_yticklabels(["file flag\n'one'","file flag\n'none'"])
ax.set_title("A. The covale leaving-atom flag is set by\nthe dictionary entry, not by earlier files",fontsize=9,loc='left')
# B template match vs baseline
ax=fig.add_subplot(gs[0,1]); T=TM[TM.field.isin(['leaving','n_pairs','prot_first','value_order','details'])]
xs=np.arange(len(T)); ax.bar(xs-.18,T.match_own_first*100,.36,color=BLUE,label="same as the code's first recorded entry")
ax.bar(xs+.18,T.match_era_site_baseline*100,.36,color=GREY,label='same as other codes first seen in same year + site')
ax.set_xticks(xs); ax.set_xticklabels(['leaving\nflag','no. of\natom pairs','atom\norder','value\norder','details'],fontsize=7.5); ax.set_ylim(50,103); ax.set_ylabel('% of 1,591 later recorded entries')
ax.legend(fontsize=6.8,frameon=False,loc='upper center',bbox_to_anchor=(.5,-.28)); ax.set_title('B. Later records match the first one, but only\nwhere chemistry or the dictionary fixes the value',fontsize=9,loc='left')
for s in ['top','right']: ax.spines[s].set_visible(False)
# C oddities same lab vs other lab
ax=fig.add_subplot(gs[0,2]); items=[]
for _,r in OD.iterrows():
    for lab,col in [('same_lab_repeat','same lab'),('other_lab_repeat','other lab')]:
        m,n=map(int,r[lab].split('/')); items.append((f"{ {'prot_first':'order','value_order':'v_order','details':'details'}[r.field]}\n{col.replace(' lab','')}",m,n,GREEN if col=='same lab' else RED))
for i,(lab,m,n,cl) in enumerate(items):
    ax.bar(i,100*m/n if n else 0,color=cl); ax.text(i,(100*m/n if n else 0)+2,f'{m}/{n}',ha='center',fontsize=8)
ax.set_xticks(range(len(items))); ax.set_xticklabels([x[0] for x in items],fontsize=6.8); ax.set_ylabel('% of later entries repeating the oddity')
ax.set_title("C. Quirks with no chemical reason repeat within\na lab (green), rarely across labs (red)",fontsize=9,loc='left')
for s in ['top','right']: ax.spines[s].set_visible(False)
# D forest: misses
ax=fig.add_subplot(gs[1,0:2])
txt=open(f'{S5}/27_stdout.txt').read()
blk=txt.split('miss ~ c_first_miss + C(c_first_era)')[1].split('miss ~')[0]
rows=re.findall(r'   (.+?): OR ([\d.]+) \(([\d.]+)-([\d.]+)\) p=(\S+)',blk)
lab={'c_first_miss':"code's first bonded use missed",'pdbj':'processed at PDBj'}
names=[];k=0
for t,o,lo,hi,p in rows:
    nm=lab.get(t,t.replace('C(c_first_era)[T.',"code's first use era ").replace('C(ref, Treatment("PHENIX"))[T.','refinement: ').replace('C(era)[T.','entry era ').replace(']',''))
    c=RED if t=='c_first_miss' else ('#2c3e50' if t=='pdbj' else GREY)
    ax.plot([float(lo),float(hi)],[k,k],color=c,lw=1.6); ax.scatter(float(o),k,color=c,s=24,zorder=3); names.append(f'{nm}  (p={p})'); k+=1
ax.axvline(1,color='k',lw=.8,ls=':'); ax.set_xscale('log'); ax.set_yticks(range(k)); ax.set_yticklabels(names,fontsize=7.5); ax.invert_yaxis()
ax.set_xlabel('OR for a later miss (1,703 later uses, 115 missed; CI clustered by code). Reference levels: PHENIX; entry era ≤2007; first-use era 2008-13')
ax.set_title('D. Misses do not track refinement software, deposition era or the era of the\ncode\'s first entry; the code\'s first-use status and PDBj remain',fontsize=9,loc='left')
for s in ['top','right']: ax.spines[s].set_visible(False)
# E struct_conn revisions
ax=fig.add_subplot(gs[1,2])
vals=[(int(U[U.miss==m].struct_conn_revised.sum()),int((U.miss==m).sum()),lab_,c) for m,lab_,c in [(0,'recorded',GREEN),(1,'missed',RED)]]
for i,(m,n,lab_,c) in enumerate(vals): ax.bar(i,100*m/n,color=c,width=.6); ax.text(i,100*m/n+1.5,f'{m}/{n}',ha='center')
ax.set_xticks([0,1]); ax.set_xticklabels(['bond\nrecorded','bond\nmissed']); ax.set_ylabel('% of entries whose _struct_conn was\nchanged in a later revision')
ax.set_title('E. Missed entries were re-touched just as\noften, yet the bond is still absent',fontsize=9,loc='left')
for s in ['top','right']: ax.spines[s].set_visible(False)
fig.text(.01,.005,'Strict set from step 4 (5,355 code–entry units; 5,054 with covale rows). Header mmCIFs from files.rcsb.org (16 Nov 2026 AKST). Annotator identity and annotation-software version are not in public files. Scripts 26–28.',fontsize=7,color='#555')
for ext in ['png','pdf','svg']: fig.savefig(f'{OUT}/copy_forward_test.{ext}',dpi=200,bbox_inches='tight')
print('fig ok')
