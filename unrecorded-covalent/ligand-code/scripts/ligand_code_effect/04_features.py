"""Attach CCD features to each contact; classify the ligand attachment atom's chemistry in the CCD definition."""
import json, pandas as pd, datetime as dt
O='/workspace/ligand_code_effect'
C=json.load(open(f'{O}/out/ccd_parsed.json'))
d=pd.read_csv(f'{O}/out/contacts_merged.tsv',sep='\t',dtype={'cys_seq':str,'lig_seq':str},keep_default_na=False)
HAL={'Cl','Br','I'}
def warhead(c,a):
    r=C.get(c)
    if not r or a not in r['nb']: return 'atom_not_in_ccd','unknown',False,False
    nb=r['nb'][a]; el=r['el']; lv=set(r['leaving_atoms'])
    heavy=[(n,o) for n,o in nb if el.get(n)!='H']; nH=sum(el.get(n)=='H' for n,o in nb)
    lv_nb=any(n in lv for n,o in nb); lv_heavy_nb=any(n in lv and el.get(n)!='H' for n,o in nb)
    hn=[n for n,o in heavy]; ring3=any(y in [z for z,_ in r['nb'].get(x,[])] for x in hn for y in hn if x<y)
    orders=[o for n,o in heavy]
    if 'TRIP' in orders: k='nitrile' if any(el.get(n)=='N' and o=='TRIP' for n,o in heavy) else 'alkyne'; s='unreacted'
    elif lv_heavy_nb or any(el.get(n) in HAL for n,o in heavy): k='leaving_group_attached'; s='unreacted'
    elif ring3 and any(el.get(n) in ('O','N','S') for n in hn): k='epoxide_aziridine'; s='unreacted'
    elif a in r['arom']: k='aromatic'; s='ambiguous'
    elif any(o=='DOUB' and el.get(n)=='O' for n,o in heavy):
        het=[n for n,o in heavy if o=='SING' and el.get(n) in ('N','O','S')]
        if het: k='acyl_C(=O)X'; s='ambiguous'
        else: k='ketone_aldehyde'; s='unreacted'
    elif any(o=='DOUB' and el.get(n)=='N' for n,o in heavy): k='imine_C=N'; s='reacted' if nH>0 else 'ambiguous'
    elif any(o=='DOUB' and el.get(n)=='C' for n,o in heavy): k='alkene_michael'; s='unreacted'
    elif any(o=='DOUB' for n,o in heavy): k='other_double'; s='ambiguous'
    elif all(o=='SING' for o in orders): k='sp3_saturated'; s='reacted'
    else: k='other'; s='ambiguous'
    return k,s,lv_nb,lv_heavy_nb
rows=[]
for i,x in d.iterrows():
    r=C.get(x.lig_comp,{})
    k,s,lvn,lvh=warhead(x.lig_comp,x.lig_atom)
    pcm=r.get('pcm_rows') or []
    dep=x.deposit_date or None; ini=r.get('pdbx_initial_date')
    age=None
    if dep and ini:
        age=(dt.date.fromisoformat(dep)-dt.date.fromisoformat(ini)).days/365.25
    rows.append(dict(ccd_type=r.get('type'),ccd_pdbx_type=r.get('pdbx_type'),ccd_formula=r.get('formula'),ccd_n_atoms=r.get('n_atoms'),
        ccd_n_heavy=r.get('n_heavy'),ccd_initial_date=ini,ccd_release_date=r.get('ccd_release_date'),ccd_create_site=r.get('create_site'),
        ccd_defining_entry=r.get('pdbx_model_coordinates_db_code'),ccd_pcm_flag=r.get('pdbx_pcm') or 'blank',
        ccd_pcm_link_atoms=';'.join(sorted({p[1] for p in pcm if p[1]})), ccd_pcm_partner_res=';'.join(sorted({p[0] for p in pcm if p[0]})),
        ccd_pcm_modify_dates=';'.join(r.get('pcm_modify_dates') or []),
        ccd_n_leaving=len(r.get('leaving_atoms',[])),ccd_n_leaving_heavy=len(r.get('leaving_heavy',[])),
        ccd_has_leaving=bool(r.get('leaving_atoms')),attach_has_leaving_nb=lvn,attach_class=k,attach_state=s,
        ccd_age_at_dep_y=None if age is None else round(age,3),code_new_at_dep=None if age is None else age<=0,
        code_defined_by_this_entry=r.get('pdbx_model_coordinates_db_code')==x.entry,
        skeleton=r.get('skeleton'),inchikey=r.get('inchikey'),ccd_name=(r.get('name') or '')[:120]))
f=pd.concat([d.reset_index(drop=True),pd.DataFrame(rows)],axis=1)
f['lt2']=f.dist<2.0
f.to_csv(f'{O}/out/contacts_features.tsv',sep='\t',index=False)
s=f[f.lt2]
print(len(f),len(s),(~s.recorded).sum())
print(pd.crosstab(s.attach_class,s.recorded,margins=True))
print(pd.crosstab(s.ccd_pcm_flag,s.recorded))
