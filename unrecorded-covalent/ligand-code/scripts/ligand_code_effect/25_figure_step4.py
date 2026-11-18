"""Figure: out/carryover_code_vs_lab.png"""
import pandas as pd, numpy as np, re, json, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O='/workspace/ligand_code_effect'; OUT=f'{O}/out'
S=pd.read_csv(f'{OUT}/step4/strict_units.tsv',sep='\t',low_memory=False); S['dep']=pd.to_datetime(S.dep)
R=pd.read_csv(f'{OUT}/step4/design_results.tsv',sep='\t')
D=pd.read_csv(f'{OUT}/step4/ccd_features_by_code.tsv',sep='\t')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5})
fig=plt.figure(figsize=(13.5,9.6)); gs=fig.add_gridspec(2,3,hspace=0.55,wspace=0.42)
GREEN,RED,GREY='#1b7f5b','#c0392b','#7f8c8d'
def bars(ax,items,title,ylab='later bonded uses not recorded (%)'):
    for i,(lab,m,n,c) in enumerate(items):
        ax.bar(i,100*m/n if n else 0,color=c,width=.7); ax.text(i,(100*m/n if n else 0)+1.5,f'{m}/{n}',ha='center',fontsize=8)
    ax.set_xticks(range(len(items))); ax.set_xticklabels([x[0] for x in items],fontsize=7); ax.set_ylabel(ylab); ax.set_title(title,fontsize=9,loc='left')
    for s in ['top','right']: ax.spines[s].set_visible(False)
# A: same code, by lab/site
Lt=S[S.later]
def cnt(q): x=Lt.query(q); return int(x.miss.sum()),len(x)
items=[]
for fm,c in [(0,GREEN),(1,RED)]:
    for lab,q in [('same\nlab',f'c_first_miss=={fm} and same_lab'),('new lab\nsame site',f'c_first_miss=={fm} and not same_lab and same_site'),('new lab\nnew site',f'c_first_miss=={fm} and not same_lab and not same_site')]:
        m,n=cnt(q); items.append((lab,m,n,c))
ax=fig.add_subplot(gs[0,0]); bars(ax,items,'A. Same code, new labs: the miss follows the code\n(green: code\'s first bonded use recorded; red: missed)')
ax.axvline(2.5,color='k',lw=.5,ls=':')
# B: same lab after first miss
X=pd.read_csv(f'{OUT}/step4/lab_after_miss_units.tsv',sep='\t')
base=S[(S.lab_prior_n>0)&(S.lab_prior_miss==0)]
it=[]
for lab,rel,c in [('same code','same code',RED),('other code,\nsame drug','other code, same drug (InChIKey-14 or skeleton)',GREEN),('other code,\nother drug','other code, different drug','#e67e22')]:
    x=X[X.relation==rel]; it.append((lab,int(x.miss.sum()),len(x),c))
it.append(('labs with no\nearlier miss',int(base.miss.sum()),len(base),GREY))
ax=fig.add_subplot(gs[0,1]); bars(ax,it,'B. Same lab after its first miss: switching code\nends the miss; other codes carry a smaller excess')
# C: forest
ax=fig.add_subplot(gs[0,2])
subs=['strict','drop HUR/80I','drop Li J. lab','drop PDBj']; terms=[('GLM c_first_miss','code: first bonded use missed',RED,-.2),('GLM lab_prior_miss_other','lab: earlier miss on another code',GREY,0),('GLM pdbj','site: PDBj',"#2c3e50",.2)]
for i,s in enumerate(subs):
    for t,lab,c,off in terms:
        r=R[(R.subset==s)&(R.design==t)]
        if not len(r): continue
        o,lo,hi=map(float,re.findall(r'OR ([\d.]+) \(([\d.]+)-([\d.]+)\)',r.stat.iloc[0])[0])
        ax.plot([lo,hi],[i+off]*2,color=c,lw=1.6); ax.scatter(o,i+off,color=c,s=24,zorder=3,label=lab if i==0 else None)
ax.axvline(1,color='k',lw=.8,ls=':'); ax.set_xscale('log'); ax.set_yticks(range(4)); ax.set_yticklabels(subs); ax.invert_yaxis()
ax.set_xlabel('adjusted OR for a later miss (95% CI, clustered by code)'); ax.legend(fontsize=7,frameon=False,loc='upper center',bbox_to_anchor=(.45,-.2),ncol=1)
ax.set_title('C. Joint model: code history survives every\nsensitivity; lab history does not survive dropping PDBj',fontsize=9,loc='left')
for s in ['top','right']: ax.spines[s].set_visible(False)
# D: case timelines
ax=fig.add_subplot(gs[1,0:2])
cases=[('MKZ','Phillips (PDBe) → Wang (RCSB)'),('XVE','Mittl (PDBe) → Gondry (RCSB)'),('59H','Owens (PDBe) → Davies (RCSB)'),('HUR','Senda → Zhang, Li (all PDBj)'),
       ('ALD','escaped: 1BP4 (RCSB 1998) → 4 labs'),('FFQ','escaped: 1UAE (RCSB 1996) → 3 labs')]
labs_=[]
for y,(c,lab) in enumerate(cases):
    g=S[S.lig_comp==c].sort_values('dep')
    for k,(_,r) in enumerate(g.iterrows()):
        ax.scatter(r.dep,y,marker='X' if r.miss else 'o',s=60 if k==0 else 40,c=RED if r.miss else GREEN,edgecolors='k' if k==0 else 'none',linewidths=1,zorder=3)
    labs_.append(f'{c}: {lab}\n{int(g.miss.sum())}/{len(g)} uses missed')
    ax.plot([g.dep.min(),g.dep.max()],[y,y],color='#ccc',lw=1,zorder=1)
ax.set_yticks(range(len(cases))); ax.set_yticklabels(labs_,fontsize=7.6); ax.invert_yaxis(); ax.set_xlim(pd.Timestamp('1994-01-01'),pd.Timestamp('2026-06-01'))
ax.set_title('D. Codes whose first bonded use was missed (black-edged marker = first use). Red X = later miss, green = recorded',fontsize=9,loc='left')
for s in ['top','right']: ax.spines[s].set_visible(False)
ax.grid(axis='x',alpha=.25)
# E: PCM
ax=fig.add_subplot(gs[1,2])
def k(cl): x=D[D.cls.str.startswith(cl)]; return int((x.pcm.astype(str)!='True').sum()),len(x)
it=[('carried\n(first missed,\nlater missed)',*k('carried'),RED),('escaped\n(first missed,\nlater recorded)',*k('escaped'),GREEN),('always\nrecorded',*k('first recorded, later all'),GREY)]
bars(ax,it,'E. CCD today: codes with NO PCM annotation\n(PCM added in bulk 2024-09-27: a consequence, not a cause)','% of codes without PCM today')
fig.text(.01,.005,'Strict set: drug-like CCD codes with ≤10 unbonded uses; closest nucleophile bond at bond-like distance (SG 1.65–1.95 Å; O/N 1.20–1.65 Å). 5,355 (code, entry) units, 334 missed. Lab = last author. '
         'After the 2024-09-27 PCM update, later uses of codes\nwithout PCM were missed 1/28 vs 6/271 with PCM (Fisher p = 0.5): no evidence that a missing PCM entry causes misses. Scripts 22–25, /workspace/ligand_code_effect.',fontsize=7,color='#555')
for ext in ['png','pdf','svg']: fig.savefig(f'{OUT}/carryover_code_vs_lab.{ext}',dpi=200,bbox_inches='tight')
print('ok')
