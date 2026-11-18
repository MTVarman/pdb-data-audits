import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
O='/workspace/ligand_code_effect/out'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
NAVY='#1f3a5f'; ORANGE='#d9822b'; GREY='#9aa5b1'; RED='#c0392b'; GREEN='#2e8b57'
fig=plt.figure(figsize=(15,10.5)); gs=fig.add_gridspec(2,2,hspace=0.45,wspace=0.62)
# A leritrelvir geometry
ax=fig.add_subplot(gs[0,0]); d=pd.read_csv(f'{O}/leritrelvir_copies.tsv',sep='\t')
d['batch']=d.deposit_date.map({'2026-02-11':'X','2025-03-17':'L','2023-02-21':'I'})
for (b,c),x in d.groupby(['batch','lig_comp']):
    col=NAVY if c=='7ON' else ORANGE; mk={'X':'o','L':'s','I':'D'}[b]
    ax.scatter(x.sg145_c7_dist,x.c7_angle_sum,s=60,marker=mk,facecolor=col if c=='7ON' else 'white',edgecolor=col,linewidth=1.8,zorder=3)
ax.axvline(2.0,color=GREY,ls='--',lw=1); ax.text(2.01,321,'2.0 Å',color=GREY,fontsize=9)
ax.axvspan(1.81,1.94,ymin=0.62,ymax=0.98,color='#f4e1c8',alpha=0.5,zorder=0)
ax.text(2.22,352,'Xiong lab, deposited 11 Feb 2026 (PDBc)\nsame geometry, different code:\n7ON 22/22 copies recorded\nA1EN0 0/10 recorded',fontsize=9,va='top')
ax.set_xlabel('Cys145 SG – ligand C7 distance (Å)'); ax.set_ylabel('Sum of bond angles at C7 (°)\n360 = flat ketone, ~328 = tetrahedral')
ax.set_title('A. Leritrelvir: the bond record follows the code, not the model',loc='left',fontweight='bold')
h=[Line2D([],[],marker='o',ls='',mfc=NAVY,mec=NAVY,label='7ON "bound form" (record present)'),Line2D([],[],marker='o',ls='',mfc='white',mec=ORANGE,mew=1.8,label='A1EN0 ketone form (no record)'),
   Line2D([],[],marker='o',ls='',color='k',mfc='none',label='Xiong 2026 batch'),Line2D([],[],marker='s',ls='',color='k',mfc='none',label='Li 2025 batch'),Line2D([],[],marker='D',ls='',color='k',mfc='none',label='8IGN (2023)')]
ax.legend(handles=h,fontsize=8,loc='lower right',frameon=False); ax.set_xlim(1.55,2.95); ax.set_ylim(315,362)
# B same drug, two codes
ax=fig.add_subplot(gs[0,1])
rows=[('Leritrelvir','7ON','PDBC',22,22,'post'),('Leritrelvir','A1EN0','PDBC',0,10,'pre'),
      ('YH-53','V7G','RCSB',5,5,'post'),('YH-53','HUR','PDBJ',0,12,'pre'),
      ('Oridonin','ODN','PDBC',2,2,'post'),('Oridonin','A1EP4','PDBC',0,1,'post'),
      ('Nirmatrelvir','4WI','RCSB+PDBe+PDBc',74,74,'post'),('Nirmatrelvir','4WI','PDBJ',7,30,'post')]
y=np.arange(len(rows))[::-1]
for yi,(drug,c,site,r,n,st) in zip(y,rows):
    col=NAVY if r==n else (RED if r==0 else ORANGE)
    ax.barh(yi,r/n*100,color=col,height=0.62); ax.text(max(r/n*100,0)+2,yi,f'{r}/{n}',va='center',fontsize=9)
ax.set_yticks(y); ax.set_yticklabels([f'{a} {c}\n{s}, {st}-reaction CCD' for a,c,s,_,_,st in rows],fontsize=8.5)
for k in (1.5,3.5,5.5): ax.axhline(y[0]-k,color='#dde3ea',lw=0.8)
ax.set_xlim(0,118); ax.set_xlabel('Bonds < 2.0 Å with a _struct_conn covale record (%)')
ax.set_title('B. Same drug, different CCD code',loc='left',fontweight='bold')
ax.set_xlabel('Bonds < 2.0 Å with a _struct_conn covale record (%)\nNirmatrelvir is the control: one code, recorded everywhere except PDBj')
# C adjusted odds ratios
ax=fig.add_subplot(gs[1,0])
L=pd.read_csv(f'{O}/logistic_models.tsv',sep='\t'); H=pd.read_csv(f'{O}/code_history_logistic.tsv',sep='\t')
def g(df,model,term): r=df[(df.model==model)&(df.term==term)].iloc[0]; return r.OR,r.ci_lo,r.ci_hi,r.p,int(r.n)
items=[('CCD is the pre-reaction form\n(ketone, nitrile, alkene…)',*g(L,'site+year+unreacted','unreacted')),
       ('CCD defines leaving atoms (anywhere)',*g(L,'site+year+leaving','leaving')),
       ('Leaving atom flagged on the\nattachment atom (link-ready CCD)',*g(L,'site+year+leaving_at_attach','leaving_at_attach')),
       ('Code created for this entry',*g(L,'site+year+newcode','newcode')),
       ('Reused code whose defining entry\nhad no bond record ("born unlinked")',*g(H,'site+year+born_unlinked','born_unlinked')),
       ('Processing site PDBj (vs RCSB)',*g(L,'site+year','C(site, Treatment("RCSB"))[T.PDBJ]'))]
for i,(lab,o,lo,hi,p,n) in enumerate(items[::-1]):
    sig=float(p)<0.01; col=(RED if o>1 else GREEN) if sig else GREY
    ax.plot([lo,hi],[i,i],color=col,lw=2.2); ax.plot(o,i,'o',color=col,ms=8)
    ax.text(60,i,f'OR {o:.2f}  p={p}',va='center',fontsize=8.5)
ax.axvline(1,color='k',lw=0.8); ax.set_xscale('log'); ax.set_xlim(0.02,400)
ax.set_yticks(range(len(items))); ax.set_yticklabels([x[0] for x in items[::-1]],fontsize=8.5)
ax.set_xlabel('Odds ratio for a MISSING record (log scale), adjusted for site and year;\n95% CI clustered by CCD code. n = 4,135 bonds < 2.0 Å (1,180 for "born unlinked")',fontsize=8.5)
ax.set_title('C. Pre/post-reaction form does not predict misses',loc='left',fontweight='bold')
# D buckets
ax=fig.add_subplot(gs[1,1]); B=pd.read_csv(f'{O}/miss_buckets.tsv',sep='\t')
order=['code new with this entry','reused code, born linked','reused code, born unlinked']
xs=np.arange(3); w=0.36
for k,(pj,col,lab) in enumerate([(0,NAVY,'RCSB, PDBe, PDBc, BNL'),(1,ORANGE,'PDBj')]):
    vals=[]; 
    for o in order:
        r=B[(B.bucket==o)&(B.pdbj==pj)].iloc[0]; vals.append((r.miss,r.n))
    ax.bar(xs+(k-0.5)*w,[m/n*100 for m,n in vals],w,color=col,label=lab)
    for xi,(m,n) in zip(xs,vals): ax.text(xi+(k-0.5)*w,m/n*100+1,f'{m}/{n}',ha='center',fontsize=8.5)
ax.set_xticks(xs); ax.set_xticklabels(['Code created\nfor this entry','Reused code,\ndefining entry had\nthe bond recorded','Reused code,\ndefining entry had\nno bond record'],fontsize=9)
ax.set_ylabel('Bonds < 2.0 Å with no record (%)'); ax.legend(frameon=False,fontsize=9,title='Processing site',title_fontsize=9)
ax.set_title('D. Misses by where the ligand code came from',loc='left',fontweight='bold')
fig.suptitle('Does the ligand code decide whether a covalent Cys bond gets recorded?\nPartly: the code\'s history and a link-ready definition matter; pre- vs post-reaction form does not',fontsize=13,fontweight='bold',x=0.06,ha='left',y=0.985)
fig.text(0.06,0.01,'Data: 4,135 ligand-carbon–Cys SG contacts < 2.0 Å from the Mpro census (Cys145) and the Cys covalent scan (24,144-entry text-selected sample), PDB files downloaded 14–16 Nov 2026; CCD files 16 Nov 2026. '
 '"Recorded" = _struct_conn covale row. M. T. Varman with Grok Bot.',fontsize=7.5,color='#555')
for e in ('png','pdf','svg'): fig.savefig(f'{O}/ligand_code_effect.{e}',dpi=200,bbox_inches='tight')
