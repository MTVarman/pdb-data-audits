"""Figure: out/inherited_code_timeline.png (+ .pdf/.svg)."""
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.dates as mdates
O='/workspace/ligand_code_effect'; OUT=f'{O}/out'
TL=pd.read_csv(f'{O}/archive/code_timelines.tsv',sep='\t',low_memory=False); TL['dep']=pd.to_datetime(TL.dep)
SEN=pd.read_csv(f'{OUT}/step3/attribution_sensitivity.tsv',sep='\t')
PAIRS=[('Leritrelvir (prior step)','7ON','A1EN0','Xiong lab, PDBc'),
       ('N4-hydroxy-dCMP, Cys adduct','7K8','NOH','Rode lab'),
       ('Open cephalosporin acyl-Ser','CEF','FZS','CEF: many labs; FZS: Schmidt lab, one day'),
       ('Diol fatty-acid adduct','3KI','J3U','Yang lab, PDBj'),
       ('Bilin, phytochrome Cys','M1V','BVX','Lagarias lab; 6XHH vs 6XHG same day'),
       ('Benzoxaborole, Ser','J00','H9O','Yan lab, PDBj')]
COL={'recorded':'#1b7f5b','missed':'#c0392b','mixed':'#e67e22','unbonded':'#9aa5b1','linked_other':'#9aa5b1'}
MK={'recorded':'o','missed':'X','mixed':'D','unbonded':'o','linked_other':'s'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5})
fig=plt.figure(figsize=(13,9.6))
gs=fig.add_gridspec(3,3,height_ratios=[1,1,1.3],hspace=0.7,wspace=0.42)
for i,(title,a,b,lab) in enumerate(PAIRS):
    ax=fig.add_subplot(gs[i//3,i%3])
    sub=TL[TL.lig_comp.isin([a,b])]
    for y,c in [(1,a),(0,b)]:
        t=TL[TL.lig_comp==c]
        for st,g in t.groupby('status'):
            ax.scatter(g.dep,np.full(len(g),y)+np.random.default_rng(1).uniform(-.12,.12,len(g)),marker=MK[st],s=46 if st!='unbonded' else 30,
                       c=COL[st] if st!='unbonded' else 'white',edgecolors=COL[st],linewidths=1.1,zorder=3)
        nb=t.status.isin(['recorded','missed','mixed']); nr=(t.status=='recorded').sum()
        nu=(t.status=='unbonded').sum()
        ax.text(1.01,y,f'{nr}/{nb.sum()} rec.'+(f'\n+{nu} unbonded' if nu else ''),transform=ax.get_yaxis_transform(),va='center',fontsize=7.5)
    lo,hi=sub.dep.min(),sub.dep.max(); pad=max((hi-lo)*0.08,pd.Timedelta(days=20))
    ax.set_xlim(lo-pad,hi+pad); ax.set_ylim(-.6,1.6); ax.set_yticks([1,0]); ax.set_yticklabels([a,b],fontweight='bold')
    loc=mdates.AutoDateLocator(maxticks=4); ax.xaxis.set_major_locator(loc); ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
    for s in ['top','right']: ax.spines[s].set_visible(False)
    ax.set_title(f'{"ABCDEF"[i]}. {title}\n{lab}',fontsize=9,loc='left')
    ax.grid(axis='x',alpha=.25)
# legend
from matplotlib.lines import Line2D
h=[Line2D([],[],marker='o',ls='',mfc=COL['recorded'],mec=COL['recorded'],label='bond <2.0 Å, covale recorded'),
   Line2D([],[],marker='X',ls='',mfc=COL['missed'],mec=COL['missed'],label='bond <2.0 Å, not recorded (file lists it as a close contact)'),
   Line2D([],[],marker='D',ls='',mfc=COL['mixed'],mec=COL['mixed'],label='mixed (some copies recorded)'),
   Line2D([],[],marker='o',ls='',mfc='white',mec=COL['unbonded'],label='code used with no nucleophile bond')]
fig.legend(handles=h,loc='upper center',ncol=4,frameon=False,bbox_to_anchor=(.5,.985),fontsize=8.5)
# panel G: first bonded use status -> later uses
ax=fig.add_subplot(gs[2,0:2])
E=pd.read_csv(f'{OUT}/step3/entry_level_bonded.tsv',sep='\t',low_memory=False)
import re
txt=open(f'{OUT}/step3/20_stdout.txt').read()
# recompute strict table for plotting
B=pd.read_csv(f'{OUT}/archive_bonds_lt2.tsv',sep='\t',low_memory=False)
nun=TL[TL.status=='unbonded'].groupby('lig_comp').size(); E['nun']=E.lig_comp.map(nun).fillna(0)
cb=B.sort_values('dist').groupby(['lig_comp','entry']).first()[['nuc_atom','dist']]; E=E.join(cb,on=['lig_comp','entry'])
E['bondlike']=np.where(E.nuc_atom=='SG',E.dist.between(1.65,1.95),E.dist.between(1.20,1.65)); E['dep']=pd.to_datetime(E.dep)
D=E[(E.cat=='drug-like')&(E.nun<=10)&E.bondlike].sort_values(['lig_comp','dep','rel','entry']).copy()
f=D.groupby('lig_comp').head(1).set_index('lig_comp'); D['fm']=D.lig_comp.map(f.miss); D['fd']=D.lig_comp.map(f.dep)
x=D[D.dep>D.fd]; x=x.assign(pj=np.where(x.site=='PDBJ','PDBj','other sites'))
t=x.groupby(['fm','pj']).agg(n=('miss','size'),m=('miss','sum')).reset_index()
xs=[];labels=[];k=0
for fm in [0,1]:
    for pj in ['other sites','PDBj']:
        r=t[(t.fm==fm)&(t.pj==pj)].iloc[0]; rate=r.m/r.n
        ax.bar(k,rate*100,color='#c0392b' if fm else '#1b7f5b',alpha=.55 if pj=='other sites' else .9,width=.7)
        ax.text(k,rate*100+1.5,f'{int(r.m)}/{int(r.n)}',ha='center',fontsize=8)
        labels.append(f'{pj}'); xs.append(k); k+=1
    k+=.6
ax.set_xticks(xs); ax.set_xticklabels(labels); ax.set_ylabel('later bonded uses not recorded (%)')
ax.text(0.5,-0.2,"code's first bonded use recorded",ha='center',transform=ax.get_xaxis_transform(),fontsize=8.5)
ax.text(3.1,-0.2,"code's first bonded use missed",ha='center',transform=ax.get_xaxis_transform(),fontsize=8.5)
ax.set_title('G. A miss on a code\'s first bonded use predicts later misses (drug-like, specific codes, bond-like distances)\n'
             'adjusted OR 14.2 (4.3–47), site/year/same-lab controlled; different-lab reuses only: OR 8.7 (2.3–33)',fontsize=9,loc='left')
for s in ['top','right']: ax.spines[s].set_visible(False)
# panel H: born-unbonded OR attenuation
ax=fig.add_subplot(gs[2,2])
lab=['all drug-like','specific codes','bond-like dist.','specific + bond-like']
for j,(g,c,off) in enumerate([('unbonded','#7f8c8d',-.13),('missed','#c0392b',.13)]):
    for i,r in SEN.iterrows():
        lo,hi=map(float,r[f'{g}_CI'].split('-')); o=r[f'{g}_OR']
        ax.plot([lo,hi],[i+off]*2,color=c,lw=1.6); ax.scatter(o,i+off,color=c,zorder=3,s=26,label=('code born on an unbonded model' if g=='unbonded' else 'code born on a missed bond') if i==0 else None)
ax.axvline(1,color='k',lw=.8,ls=':'); ax.set_xscale('log'); ax.set_yticks(range(4)); ax.set_yticklabels(lab); ax.invert_yaxis()
ax.set_xlabel('OR for a later miss vs code born recorded\n(site + year adjusted, CI clustered by code)')
ax.set_title('H. "Born unbonded" (grey) shrinks under stricter\nfilters; "born on a missed bond" (red) does not',fontsize=9,loc='left'); 
for s in ['top','right']: ax.spines[s].set_visible(False)
fig.text(.01,.005,'Source: header mmCIFs for all 199,326 PDB entries with a non-polymer entity (files.rcsb.org, fetched 16 Nov 2026 AKST). Bond = Cys SG / Ser OG / Lys NZ / His NE2,ND1 within 2.0 Å of a ligand C/B/P/S/Si/Se/As atom, '
         'from _struct_conn covale\nor the same file\'s _pdbx_validate_close_contact table. Panels A–F show every entry using each code (deposition date). Scripts 11–21, /workspace/ligand_code_effect.',fontsize=7,color='#555')
for ext in ['png','pdf','svg']: fig.savefig(f'{OUT}/inherited_code_timeline.{ext}',dpi=200,bbox_inches='tight')
print('ok')
