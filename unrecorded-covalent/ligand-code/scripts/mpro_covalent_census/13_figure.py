"""One figure: (A) Cys145 bond records by processing site (contacts <2.0 A), (B) entries with a missing record by
deposition year, PDBj vs other sites, (C) recorded/total by drug and site."""
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
D='/workspace/mpro_covalent_census'; O=f'{D}/out'
M=pd.read_csv(f'{O}/contacts_main_cys145_C_1.6-2.4.tsv',sep='\t',dtype={'deposit_year':str},low_memory=False)
EN=pd.read_csv(f'{O}/entries_covalent.tsv',sep='\t',dtype={'deposit_year':str})
lt=M[M.dist<2.0]
fig,ax=plt.subplots(1,3,figsize=(17,5.6),gridspec_kw={'width_ratios':[1,1.25,1.6]})
# A
sites=['RCSB','PDBE','PDBC','PDBJ']; lab={'RCSB':'RCSB PDB','PDBE':'PDBe','PDBC':'PDBc','PDBJ':'PDBj'}
rec=[int(lt[lt.process_site==s].recorded.sum()) for s in sites]; tot=[int((lt.process_site==s).sum()) for s in sites]
un=[t-r for r,t in zip(rec,tot)]
a=ax[0]; x=np.arange(len(sites))
a.bar(x,[100*u/t for u,t in zip(un,tot)],color=['#9aa5b1','#9aa5b1','#e09f3e','#c0392b'])
for i,(u,t) in enumerate(zip(un,tot)): a.text(i,100*u/t+1,f'{u}/{t}',ha='center',fontsize=10)
a.set_xticks(x,[lab[s] for s in sites]); a.set_ylabel('% of Cys145–ligand bonds (<2.0 Å)\nwith no _struct_conn record'); a.set_ylim(0,40)
a.set_title('A  Missing bond records by processing site',loc='left',fontsize=11,fontweight='bold')
# B
EN['grp']=np.where(EN.process_site=='PDBJ','PDBj','RCSB/PDBe/PDBc')
yrs=[str(y) for y in range(2020,2027)]
a=ax[1]; w=0.38
for k,(g,c) in enumerate([('PDBj','#c0392b'),('RCSB/PDBe/PDBc','#9aa5b1')]):
    v=[];lb=[]
    for y in yrs:
        z=EN[(EN.grp==g)&(EN.deposit_year==y)]; n=len(z); m=int(z.has_unrec_lt2.sum())
        v.append(100*m/n if n else np.nan); lb.append(f'{m}/{n}' if n else '–')
    xs=np.arange(len(yrs))+(k-0.5)*w
    a.bar(xs,np.nan_to_num(v),w,color=c,label=g)
    for xx,vv,l in zip(xs,v,lb): a.text(xx,(0 if np.isnan(vv) else vv)+1.2,l,ha='center',fontsize=7.5,rotation=90)
a.set_xticks(np.arange(len(yrs)),yrs); a.set_xlabel('deposition year'); a.set_ylabel('% of covalent Mpro entries with ≥1\nunrecorded bond (<2.0 Å)'); a.set_ylim(0,80); a.legend(frameon=False,fontsize=9,loc='upper right')
a.set_title('B  PDBj misses cluster in 2021–2023 depositions',loc='left',fontsize=11,fontweight='bold')
# C
drugs=[('4WI','nirmatrelvir'),('80I','lufotrelvir (PF-07304814)'),('HUR','YH-53'),('7ON','leritrelvir, bound-form CCD 7ON'),('A1EN0','leritrelvir, CCD A1EN0'),('YDL','ibuzatrelvir'),('V2M','PF-00835231'),('UED','GC373/GC376 (UED)'),('K36','GC376 (K36)'),('U5G','boceprevir'),('SV6','telaprevir'),('FHR','bofutrelvir (11a)')]
a=ax[2]
for i,(c,n) in enumerate(drugs):
    for j,s in enumerate(sites):
        z=lt[(lt.lig_comp==c)&(lt.process_site==s)]
        if not len(z): continue
        r=int(z.recorded.sum()); t=len(z); f=r/t
        a.scatter(j,i,s=60+14*t,c=[plt.cm.RdYlGn(f)],edgecolor='k',lw=0.5,vmin=0,vmax=1)
        a.text(j+0.22,i,f'{r}/{t}',va='center',fontsize=8.5)
a.set_yticks(range(len(drugs)),[n for _,n in drugs],fontsize=9); a.invert_yaxis(); a.set_xticks(range(len(sites)),[lab[s] for s in sites]); a.set_xlim(-0.5,3.8)
a.set_title('C  Recorded / modelled Cys145 bonds, by drug and site',loc='left',fontsize=11,fontweight='bold')
a.text(-0.5,len(drugs)+0.35,'Ensitrelvir (non-covalent control): 28 entries, 0 SG contacts ≤2.6 Å, 0 covale rows.',fontsize=8.5,va='top')
fig.suptitle(f'SARS-CoV-2 Mpro: {len(EN)} covalent complexes among 1,957 entries. All {int((~lt.recorded).sum())} unrecorded Cys145 bonds (<2.0 Å) are listed as clashes in the file\'s own close-contact table, and PLINDER labels every one non-covalent.',fontsize=10.5,y=1.0)
plt.tight_layout()
for ext in ['png','pdf','svg']: plt.savefig(f'{O}/mpro_covalent_census.{ext}',dpi=200,bbox_inches='tight')
