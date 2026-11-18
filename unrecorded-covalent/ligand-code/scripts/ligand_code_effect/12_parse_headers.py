"""Parse every header mmCIF: entry metadata, non-polymer comps, covale rows joining a protein nucleophile atom to a
non-polymer residue, and _pdbx_validate_close_contact rows (<=2.2 A) between a protein nucleophile atom and a
non-polymer atom. Nucleophile atoms: CYS SG, SER OG, LYS NZ, HIS NE2/ND1 (THR OG1, TYR OH kept, flagged).
Output: archive/parsed.jsonl (one line per entry)."""
import gemmi, gzip, re, json, sys, os, glob, multiprocessing as mp
A='/workspace/ligand_code_effect/archive'
NUC={('CYS','SG'),('SER','OG'),('LYS','NZ'),('HIS','NE2'),('HIS','ND1'),('THR','OG1'),('TYR','OH')}
RX=re.compile(r'^loop_\s*\n(?=(#|loop_|data_|\Z))',re.M)
def v(x):
    if x is None: return None
    s=gemmi.cif.as_string(x); return None if s in ('?','.','') else s
def rows(b,cat,tags):
    t=b.find(cat,tags) if all(b.find_mmcif_category(cat).has_tag(cat+x) if False else True for x in tags) else []
    return [[v(r[i]) for i in range(len(tags))] for r in t]
def opt(b,cat,tags):
    c=b.find_mmcif_category(cat); have=set(c.tags)
    if not have: return []
    idx=[cat+x for x in tags]; tab=b.find(cat,[x if (cat+x) in have else '?'+x for x in tags])
    return [[v(r[i]) if r.has(i) else None for i in range(len(tags))] for r in tab]
def one(b,tag):
    x=b.find_value(tag); return v(x) if x is not None else None
def parse(p):
    e=os.path.basename(p).split('.')[0]
    try:
        t=gzip.open(p,'rt').read(); b=gemmi.cif.read_string(RX.sub('',t)).sole_block()
    except Exception as ex: return {'e':e,'err':str(ex)[:100]}
    cc=opt(b,'_chem_comp.',['id','type'])
    nonpoly={c for c,ty in cc if ty and ty.lower()=='non-polymer' and c not in ('HOH','DOD')}
    revs=[x[0] for x in opt(b,'_pdbx_audit_revision_history.',['revision_date']) if x[0]]
    au=[x[0] for x in opt(b,'_audit_author.',['name']) if x[0]]
    meth=';'.join(x[0] for x in opt(b,'_exptl.',['method']) if x[0])
    grp=';'.join(sorted({x[0] for x in opt(b,'_pdbx_deposit_group.',['group_id']) if x[0]}))
    r=dict(e=e,ps=one(b,'_pdbx_database_status.process_site'),ds=one(b,'_pdbx_database_status.deposit_site'),
           dep=one(b,'_pdbx_database_status.recvd_initial_deposition_date'),rel=min(revs) if revs else None,
           au_last=au[-1] if au else None,au_first=au[0] if au else None,meth=meth,grp=grp or None,
           hpm=one(b,'_pdbx_entry_details.has_protein_modification'),np=sorted(nonpoly))
    cov=[]
    for x in opt(b,'_struct_conn.',['conn_type_id','ptnr1_label_comp_id','ptnr1_label_atom_id','ptnr1_auth_asym_id','ptnr1_auth_seq_id','pdbx_ptnr1_PDB_ins_code',
                                     'ptnr2_label_comp_id','ptnr2_label_atom_id','ptnr2_auth_asym_id','ptnr2_auth_seq_id','pdbx_ptnr2_PDB_ins_code','pdbx_dist_value','ptnr1_auth_comp_id','ptnr2_auth_comp_id']):
        if not (x[0] or '').startswith('covale'): continue
        for a,bb in ((1,6),(6,1)):
            if (x[a],x[a+1]) in NUC and x[bb] in nonpoly:
                cov.append([x[a],x[a+1],x[a+2],x[a+3],x[a+4],x[bb],x[bb+1],x[bb+2],x[bb+3],x[bb+4],x[11]])
    # covale rows joining a non-polymer comp to any polymer residue (any atom), for 'linked elsewhere'
    anylink=set()
    for x in opt(b,'_struct_conn.',['conn_type_id','ptnr1_label_comp_id','ptnr2_label_comp_id']):
        if (x[0] or '').startswith('covale'):
            if x[1] in nonpoly and x[2] not in nonpoly: anylink.add(x[1])
            if x[2] in nonpoly and x[1] not in nonpoly: anylink.add(x[2])
    ccl=[]
    for x in opt(b,'_pdbx_validate_close_contact.',['auth_atom_id_1','auth_asym_id_1','auth_comp_id_1','auth_seq_id_1','PDB_ins_code_1','label_alt_id_1',
                                                    'auth_atom_id_2','auth_asym_id_2','auth_comp_id_2','auth_seq_id_2','PDB_ins_code_2','label_alt_id_2','dist','PDB_model_num']):
        try: d=float(x[12])
        except: continue
        if d>2.2 or (x[13] not in (None,'1')): continue
        for a,bb in ((0,6),(6,0)):
            if (x[a+2],x[a]) in NUC and x[bb+2] in nonpoly:
                ccl.append([x[a+2],x[a],x[a+1],x[a+3],x[a+4],x[bb+2],x[bb],x[bb+1],x[bb+3],x[bb+4],d,x[a+5],x[bb+5]])
    r['cov']=cov; r['cc']=ccl; r['link_any']=sorted(anylink)
    r['has_cc_table']=bool(b.find_mmcif_category('_pdbx_validate_close_contact.').tags)
    return r
if __name__=='__main__':
    fs=sorted(glob.glob(f'{A}/headers/*.cif.gz'))
    if len(sys.argv)>1: fs=fs[:int(sys.argv[1])]
    out=sys.argv[2] if len(sys.argv)>2 else f'{A}/parsed.jsonl'
    with mp.Pool(6) as pool, open(out,'w') as fo:
        for i,r in enumerate(pool.imap_unordered(parse,fs,chunksize=50)):
            fo.write(json.dumps(r,separators=(',',':'))+'\n')
    print('parsed',len(fs))
