"""Same drug / chemotype under more than one CCD code. Grouping key: heavy-atom skeleton (all bonds single, no
stereo, leaving atoms removed) so pre- and post-reaction forms of one compound fall together. Contacts < 2.0 A only."""
import pandas as pd
O='/workspace/ligand_code_effect/out'
f=pd.read_csv(f'{O}/contacts_features.tsv',sep='\t',low_memory=False); s=f[f.lt2]
g=s.groupby('skeleton').lig_comp.nunique(); multi=g[g>1].index; s=s[s.skeleton.isin(multi)]
gid={k:i+1 for i,k in enumerate(sorted(multi))}
t=s.groupby(['skeleton','lig_comp','process_site']).agg(n=('recorded','size'),recorded=('recorded','sum'),entries=('entry',lambda x:','.join(sorted(set(x)))),
   dist_min=('dist','min'),dist_max=('dist','max'),deposit_first=('deposit_date','min'),deposit_last=('deposit_date','max'),
   attach_class=('attach_class','first'),attach_state=('attach_state','first'),ccd_pcm=('ccd_pcm_flag','first'),ccd_initial_date=('ccd_initial_date','first'),
   ccd_defining_entry=('ccd_defining_entry','first'),ccd_n_leaving=('ccd_n_leaving','first'),inchikey=('inchikey','first'),name=('ccd_name','first')).reset_index()
t.insert(0,'group',t.skeleton.map(gid)); t['rate']=(t.recorded/t.n).round(3)
t.to_csv(f'{O}/multicode_groups.tsv',sep='\t',index=False)
c=s.groupby(['skeleton','lig_comp']).agg(n=('recorded','size'),rec=('recorded','sum'),pcm=('ccd_pcm_flag','first'),state=('attach_state','first')).reset_index()
c['rate']=c.rec/c.n
sp=c.groupby('skeleton').rate.agg(lambda x:x.max()-x.min()); split=sp[sp>0].index
print('groups with >1 code (<2.0 A contacts):',len(multi),'contacts',len(s),'entries',s.entry.nunique())
print('groups whose codes differ in recorded rate:',len(split))
cs=c[c.skeleton.isin(split)]
low=cs.loc[cs.groupby('skeleton').rate.idxmin()]; high=cs.loc[cs.groupby('skeleton').rate.idxmax()]
print('lower-rate code PCM blank:',(low.pcm!='Y').sum(),'/',len(low),' higher-rate code PCM Y:',(high.pcm=='Y').sum(),'/',len(high))
print('lower-rate code pre-reaction form:',(low.state=='unreacted').sum(),'/',len(low),' higher-rate code pre-reaction form:',(high.state=='unreacted').sum())
# same-site splits
ss=s[s.skeleton.isin(split)].groupby(['skeleton','process_site']).lig_comp.nunique(); print('split groups with 2 codes at the same site:',ss[ss>1].reset_index()[['process_site']].value_counts().to_dict())
# groups where pre- and post-reaction codes coexist: compare
both=c.groupby('skeleton').state.agg(lambda x:{'unreacted','reacted'}<=set(x)); B=c[c.skeleton.isin(both[both].index)]
pre=B[B.state=='unreacted']; post=B[B.state=='reacted']
print('groups with both a pre- and a post-reaction code:',both.sum(),' pre-form contacts recorded',int(pre.rec.sum()),'/',int(pre.n.sum()),' post-form',int(post.rec.sum()),'/',int(post.n.sum()))
low.assign(role='lower').to_csv(f'{O}/multicode_split_groups.tsv',sep='\t',index=False)
pd.concat([low.assign(role='lower'),high.assign(role='higher')]).sort_values(['skeleton','role']).to_csv(f'{O}/multicode_split_groups.tsv',sep='\t',index=False)
