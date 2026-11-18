"""Step 5: re-read header mmCIFs for every entry in the strict set and extract, for that entry's code:
full covale rows (atom pair, order, leaving flag, details, dist string, value_order, role, symmetry),
close contacts (nucleophile to that code), refinement software, revision history and categories touched per revision."""
import gemmi, gzip, re, json, os, multiprocessing as mp, pandas as pd
O='/workspace/ligand_code_effect'; A=f'{O}/archive'
RX=re.compile(r'^loop_\s*\n(?=(#|loop_|data_|\Z))',re.M)
NUC={('CYS','SG'),('SER','OG'),('LYS','NZ'),('HIS','NE2'),('HIS','ND1')}
def v(x):
    if x is None: return None
    s=gemmi.cif.as_string(x); return None if s in ('?','.','') else s
def raw(x): return None if x is None else str(x)
def un(x): return None if x is None else gemmi.cif.as_string(x)
def tab(b,cat,tags):
    c=b.find_mmcif_category(cat); have=set(c.tags)
    if not have: return []
    t=b.find(cat,[x if (cat+x) in have else '?'+x for x in tags])
    return [[((raw(r[i]) if tags[i] in RAWTAGS else un(r[i])) if r.has(i) else None) for i in range(len(tags))] for r in t]
RAWTAGS={'pdbx_dist_value'}
def work(args):
    e,codes=args; p=f'{A}/headers/{e}.cif.gz'
    try: b=gemmi.cif.read_string(RX.sub('',gzip.open(p,'rt').read())).sole_block()
    except Exception as ex: return dict(e=e,err=str(ex)[:80])
    T=['id','conn_type_id','pdbx_leaving_atom_flag','ptnr1_label_comp_id','ptnr1_label_atom_id','ptnr1_auth_asym_id','ptnr1_auth_seq_id',
       'ptnr2_label_comp_id','ptnr2_label_atom_id','ptnr2_auth_asym_id','ptnr2_auth_seq_id','details','pdbx_dist_value','pdbx_value_order','pdbx_role','ptnr1_symmetry','ptnr2_symmetry','pdbx_ptnr1_label_alt_id','pdbx_ptnr2_label_alt_id']
    cov=[]; cov_all_n=0
    for r in tab(b,'_struct_conn.',T):
        if not (r[1] or '').startswith('covale'): continue
        cov_all_n+=1
        for code in codes:
            if r[3]==code or r[7]==code:
                prot_first= (r[3],r[4]) in NUC
                nuc=(r[3],r[4]) if prot_first else (r[7],r[8]); lig=r[8] if prot_first else r[4]
                other=r[7] if prot_first else r[3]
                cov.append(dict(code=code,id=r[0],conn_type=r[1],leaving=r[2],prot_first=prot_first,nuc_res=nuc[0],nuc_atom=nuc[1],lig_atom=lig,partner_comp=other,
                                details=r[11],dist_str=r[12],value_order=r[13],role=r[14],sym=f'{r[15]}|{r[16]}',alt=f'{r[17]}|{r[18]}'))
    cc=[]
    for r in tab(b,'_pdbx_validate_close_contact.',['auth_atom_id_1','auth_comp_id_1','auth_asym_id_1','auth_seq_id_1','auth_atom_id_2','auth_comp_id_2','auth_asym_id_2','auth_seq_id_2','dist','PDB_model_num']):
        if r[9] not in (None,'1'): continue
        for a,bb in ((0,4),(4,0)):
            if (gemmi.cif.as_string(r[a+1]),gemmi.cif.as_string(r[a])) in NUC and gemmi.cif.as_string(r[bb+1]) in codes:
                cc.append(dict(code=gemmi.cif.as_string(r[bb+1]),nuc_res=gemmi.cif.as_string(r[a+1]),nuc_atom=gemmi.cif.as_string(r[a]),nuc_chain=r[a+2],nuc_seq=r[a+3],lig_atom=gemmi.cif.as_string(r[bb]),lig_chain=r[bb+2],dist=float(r[8])))
    sw=[(v(x[0]),v(x[1]),v(x[2])) for x in tab(b,'_software.',['classification','name','version'])]
    ref=[f'{n}|{ver}' for c,n,ver in sw if c and 'refine' in c.lower()]
    rh=[(v(x[0]),v(x[1]),v(x[2]),v(x[3])) for x in tab(b,'_pdbx_audit_revision_history.',['ordinal','major_revision','minor_revision','revision_date'])]
    rc=[(v(x[0]),v(x[1])) for x in tab(b,'_pdbx_audit_revision_category.',['revision_ordinal','category'])]
    rg=[(v(x[0]),v(x[1])) for x in tab(b,'_pdbx_audit_revision_group.',['revision_ordinal','group'])]
    rd=[(v(x[0]),v(x[1]),v(x[2])) for x in tab(b,'_pdbx_audit_revision_details.',['revision_ordinal','type','description'])]
    return dict(e=e,cov=cov,cov_all_n=cov_all_n,cc=cc,refine=ref,rev=rh,rev_cat=rc,rev_grp=rg,rev_det=rd,
                pdbx_entry_hpm=v(b.find_value('_pdbx_entry_details.has_protein_modification')) if b.find_value('_pdbx_entry_details.has_protein_modification') is not None else None)
if __name__=='__main__':
    S=pd.read_csv(f'{O}/out/step4/strict_units.tsv',sep='\t',low_memory=False)
    job=S.groupby('entry').lig_comp.agg(lambda s:sorted(set(s))).reset_index().values.tolist()
    with mp.Pool(7) as pool, open(f'{O}/out/step5/templates.jsonl','w') as fo:
        for r in pool.imap_unordered(work,job,chunksize=20): fo.write(json.dumps(r)+'\n')
    print('entries',len(job))
