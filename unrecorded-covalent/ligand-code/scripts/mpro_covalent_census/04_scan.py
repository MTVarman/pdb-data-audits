"""Scan every SARS-CoV-2 Mpro entry for ligand atoms near the catalytic nucleophiles, with _struct_conn match.
Conventions follow /workspace/cys_covalent_scan/scripts/03_scan.py (model 1, H removed, occ>0, altlocs must match
or be blank, nearest periodic image via gemmi NeighborSearch).
Nucleophile atoms scanned in Mpro chains (entity from out/mpro_entities.tsv; nsp5 position mapped via label_seq):
  every Cys SG (pos145 flagged), Ser OG at nsp5 pos 145 (C145S mutants) and all Ser OG, His41 NE2/ND1.
Partners: any heavy, occ>0 atom within 2.6 A belonging to (a) a non-polymer residue (not water, >1 heavy atom, not
all-metal) or (b) a residue of a polymer chain whose entity is NOT an Mpro entity (peptide/peptidomimetic inhibitors).
One row per (nucleophile residue, partner residue, nucleophile atom), nearest partner atom kept.
Recorded = a _struct_conn row with conn_type_id covale* joining that nucleophile residue to that partner residue
(auth chain/seq/ins/comp), either order, any symmetry. Atom-exact flag also kept.
Entry metadata: process site, dates, method, resolution, deposit group, citation, authors, software, title,
nonpolymer components with heavy-atom counts, and every covale row touching an Mpro Cys145."""
import sys, os, json, time, gzip
import multiprocessing as mp
import gemmi
D='/workspace/mpro_covalent_census'; CIF=f'{D}/data/cif'; R=2.6
WATER={'HOH','DOD','WAT','H2O'}
AA={'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','SEC','MSE'}
ENT={}
for l in open(f'{D}/out/mpro_entities.tsv'):
    if l.startswith('#') or l.startswith('entry\t'): continue
    f=l.rstrip('\n').split('\t'); ENT.setdefault(f[0],{})[f[1].split('_')[1]]=int(f[7])
def v(x):
    if x is None: return None
    x=gemmi.cif.as_string(x); return None if x in ('?','.','') else x
def first(b,tag):
    x=b.find_value(tag); return v(x) if x is not None else None
def col(b,tag): return [v(x) for x in b.find_values(tag)]
def meta(b):
    m={'process_site':first(b,'_pdbx_database_status.process_site'),
       'deposit_site':first(b,'_pdbx_database_status.deposit_site'),
       'deposit_date':first(b,'_pdbx_database_status.recvd_initial_deposition_date'),
       'status_code':first(b,'_pdbx_database_status.status_code'),
       'sg_entry':first(b,'_pdbx_database_status.SG_entry'),
       'title':first(b,'_struct.title'),
       'method':';'.join(x for x in col(b,'_exptl.method') if x) or None,'resolution':None}
    for tag in ('_refine.ls_d_res_high','_em_3d_reconstruction.resolution','_reflns.d_resolution_high'):
        x=first(b,tag)
        if x:
            try: m['resolution']=float(x); break
            except ValueError: pass
    revs=[x for x in col(b,'_pdbx_audit_revision_history.revision_date') if x]
    m['initial_release']=min(revs) if revs else None; m['last_revision']=max(revs) if revs else None
    m['n_revisions']=len(revs)
    m['deposit_group']=';'.join(sorted({x for x in col(b,'_pdbx_deposit_group.group_id') if x})) or None
    m['deposit_group_title']=(first(b,'_pdbx_deposit_group.group_title') or '')[:200] or None
    m['deposit_group_type']=first(b,'_pdbx_deposit_group.group_type')
    cid=col(b,'_citation.id'); ct=col(b,'_citation.title'); cj=col(b,'_citation.journal_abbrev'); cy=col(b,'_citation.year'); cd=col(b,'_citation.pdbx_database_id_DOI'); cp=col(b,'_citation.pdbx_database_id_PubMed')
    if cid:
        k=cid.index('primary') if 'primary' in cid else 0
        g=lambda L:(L[k] if k<len(L) else None)
        m.update(citation_title=(g(ct) or '')[:250],citation_journal=g(cj),citation_year=g(cy),citation_doi=g(cd),citation_pubmed=g(cp))
    m['authors']=';'.join(x for x in col(b,'_audit_author.name') if x)[:400]
    sn=col(b,'_software.name'); sv=col(b,'_software.version'); sc=col(b,'_software.classification')
    m['software']=';'.join(f'{a}|{v_ or ""}|{c or ""}' for a,v_,c in zip(sn,sv+[None]*len(sn),sc+[None]*len(sn)) if a)
    return m
def conns(b):
    cat=b.get_mmcif_category('_struct_conn.')
    if not cat: return []
    n=len(next(iter(cat.values()))); out=[]
    cols=['id','conn_type_id','ptnr1_auth_asym_id','ptnr1_auth_seq_id','pdbx_ptnr1_PDB_ins_code','ptnr1_label_comp_id','ptnr1_auth_comp_id','ptnr1_label_atom_id',
          'ptnr2_auth_asym_id','ptnr2_auth_seq_id','pdbx_ptnr2_PDB_ins_code','ptnr2_label_comp_id','ptnr2_auth_comp_id','ptnr2_label_atom_id','pdbx_dist_value','pdbx_leaving_atom_flag','pdbx_role','details']
    for i in range(n):
        d={}
        for c in cols:
            x=cat.get(c); x=x[i] if x is not None else None
            d[c]=None if x in (None,False,True,'?','.','') else str(x)
        out.append(d)
    return out
def rk(ch,seq,ins,comp): return (ch,str(seq),(ins or '').strip(),comp)
def scan(e):
    res={'entry':e}
    try:
        doc=gemmi.cif.read(f'{CIF}/{e}.cif.gz'); b=doc.sole_block()
        res.update(meta(b))
        cl=conns(b); rec={}; anyc={}
        for c in cl:
            k1=rk(c['ptnr1_auth_asym_id'],c['ptnr1_auth_seq_id'],c['pdbx_ptnr1_PDB_ins_code'],c['ptnr1_auth_comp_id'] or c['ptnr1_label_comp_id'])
            k2=rk(c['ptnr2_auth_asym_id'],c['ptnr2_auth_seq_id'],c['pdbx_ptnr2_PDB_ins_code'],c['ptnr2_auth_comp_id'] or c['ptnr2_label_comp_id'])
            fk=frozenset([k1,k2]); anyc.setdefault(fk,set()).add(c['conn_type_id'])
            if (c['conn_type_id'] or '').lower().startswith('covale'):
                rec.setdefault(fk,[]).append({k1:c['ptnr1_label_atom_id'],k2:c['ptnr2_label_atom_id'],'d':c['pdbx_dist_value']})
        st=gemmi.make_structure_from_block(b); st.setup_entities(); st.remove_hydrogens()
        # nonpolymer components present + heavy atoms
        mpro_ents=ENT.get(e,{})
        ent_of={}
        for en in st.entities:
            for sub in en.subchains: ent_of[sub]=en
        del doc,b
        m=st[0]
        try: ns=gemmi.NeighborSearch(m,st.cell,5).populate(); res['symmetry_checked']=True
        except RuntimeError:
            st.cell=gemmi.UnitCell(); ns=gemmi.NeighborSearch(m,st.cell,5).populate(); res['symmetry_checked']=False
        crystal=st.cell.is_crystal()
        comps={}
        polychains={}
        for ch in m:
            for r in ch:
                en=ent_of.get(r.subchain)
                if r.entity_type==gemmi.EntityType.NonPolymer:
                    comps.setdefault(r.name,0); comps[r.name]=max(comps[r.name],sum(1 for a in r if not a.is_hydrogen()))
        res['nonpoly_comps']=comps
        def is_mpro_res(r):
            en=ent_of.get(r.subchain); return en is not None and en.name in mpro_ents and r.entity_type==gemmi.EntityType.Polymer
        # heavy atoms per partner polymer subchain (non-Mpro)
        sub_heavy={}
        for ch in m:
            for r in ch:
                if r.entity_type==gemmi.EntityType.Polymer and not is_mpro_res(r):
                    sub_heavy[r.subchain]=sub_heavy.get(r.subchain,0)+sum(1 for a in r if not a.is_hydrogen())
        lig_ok={}; out=[]; n_mpro_chains=0; cys145=[]
        for ch in m:
            seen_mpro=False
            for r in ch:
                if not is_mpro_res(r): continue
                seen_mpro=True
                p145=mpro_ents[ent_of[r.subchain].name]
                pos=None if r.label_seq is None else r.label_seq-(p145-145)
                targets=[]
                if r.name=='CYS': targets=['SG']
                elif r.name=='SER': targets=['OG']
                elif r.name=='HIS' and pos==41: targets=['NE2','ND1']
                if pos==145: cys145.append(f'{ch.name}{r.seqid.num}{r.name}')
                if not targets: continue
                nk=rk(ch.name,r.seqid.num,r.seqid.icode,r.name)
                best={}
                for a in r:
                    if a.name not in targets or a.occ<=0: continue
                    for mk in ns.find_atoms(a.pos,'\0',radius=R):
                        cra=mk.to_cra(m); lr=cra.residue; la=cra.atom
                        if lr.name in WATER or la.is_hydrogen() or la.occ<=0: continue
                        if lr.entity_type==gemmi.EntityType.NonPolymer:
                            lk=(cra.chain.name,lr.seqid.num,lr.seqid.icode,lr.name)
                            if lk not in lig_ok:
                                hv=[x for x in lr if not x.is_hydrogen()]
                                lig_ok[lk]=(len(hv)>1 and not all(x.element.is_metal for x in hv), len(hv))
                            if not lig_ok[lk][0] or la.element.is_metal: continue
                            ptype='nonpolymer'; nheavy=lig_ok[lk][1]
                        elif lr.entity_type==gemmi.EntityType.Polymer and not is_mpro_res(lr):
                            lk=(cra.chain.name,lr.seqid.num,lr.seqid.icode,lr.name); ptype='polymer_nonMpro'; nheavy=sub_heavy.get(lr.subchain)
                        else: continue
                        if a.has_altloc() and la.has_altloc() and a.altloc!=la.altloc: continue
                        if mk.image_idx==0 or not crystal: d=a.pos.dist(la.pos); img=0
                        else: d=a.pos.dist(st.cell.find_nearest_pbc_position(a.pos,la.pos,mk.image_idx)); img=mk.image_idx
                        if d>R: continue
                        key=(lk,a.name)
                        if key not in best or d<best[key][0]:
                            best[key]=(d,la.name,la.element.name,a.altloc,la.altloc,la.occ,a.occ,img,ptype,nheavy,lr.subchain)
                for (lk,an),(d,lan,lel,calt,lalt,locc,cocc,img,ptype,nheavy,lsub) in best.items():
                    lkey=rk(lk[0],lk[1],lk[2],lk[3]); fk=frozenset([nk,lkey]); recs=rec.get(fk,[])
                    oth=[k for f2 in rec for k in f2 if nk in f2 and k!=nk and k!=lkey]
                    out.append({'nuc_chain':ch.name,'nuc_seq':r.seqid.num,'nuc_ins':r.seqid.icode.strip(),'nuc_res':r.name,'nuc_atom':an,'nsp5_pos':pos,
                        'nuc_alt':calt.strip() if calt!='\0' else '','partner_type':ptype,'lig_comp':lk[3],'lig_chain':lk[0],'lig_seq':lk[1],'lig_ins':lk[2].strip(),
                        'lig_label_asym':lsub,'lig_atom':lan,'lig_el':lel,'lig_alt':lalt.strip() if lalt!='\0' else '','lig_occ':round(locc,2),'nuc_occ':round(cocc,2),
                        'lig_heavy_atoms':nheavy,'dist':round(d,3),'sym_image':img,'recorded':bool(recs),
                        'rec_atom_exact':any(x.get(nk)==an and x.get(lkey)==lan for x in recs),'rec_dist':';'.join(str(x['d']) for x in recs),
                        'conn_types_pair':';'.join(sorted(t for t in anyc.get(fk,set()) if t)),
                        'nuc_cov_other':';'.join(f'{k[0]}:{k[3]}{k[1]}' for k in oth)})
            if seen_mpro: n_mpro_chains+=1
        res['n_mpro_chains']=n_mpro_chains; res['pos145_residues']=';'.join(cys145)
        # every covale row touching an Mpro Cys145 / Ser145
        c145=set()
        for ch in m:
            for r in ch:
                if is_mpro_res(r) and r.label_seq is not None and r.label_seq-(mpro_ents[ent_of[r.subchain].name]-145)==145: c145.add(rk(ch.name,r.seqid.num,r.seqid.icode,r.name))
        res['covale_rows_pos145']=[ {k:c[k] for k in ('id','ptnr1_auth_asym_id','ptnr1_auth_seq_id','ptnr1_auth_comp_id','ptnr1_label_atom_id','ptnr2_auth_asym_id','ptnr2_auth_seq_id','ptnr2_auth_comp_id','ptnr2_label_atom_id','pdbx_dist_value')}
            for c in cl if (c['conn_type_id'] or '').lower().startswith('covale') and (rk(c['ptnr1_auth_asym_id'],c['ptnr1_auth_seq_id'],c['pdbx_ptnr1_PDB_ins_code'],c['ptnr1_auth_comp_id']) in c145 or rk(c['ptnr2_auth_asym_id'],c['ptnr2_auth_seq_id'],c['pdbx_ptnr2_PDB_ins_code'],c['ptnr2_auth_comp_id']) in c145)]
        res['n_covale_rows_total']=sum(1 for c in cl if (c['conn_type_id'] or '').lower().startswith('covale'))
        res['pairs']=out; res['status']='ok'
    except Exception as ex:
        import traceback; res['status']='error'; res['error']=traceback.format_exc()[-400:]
    return res
if __name__=='__main__':
    ids=open(f'{D}/out/mpro_entries.txt').read().split()
    if len(sys.argv)>1: ids=sys.argv[1:]
    outp=f'{D}/out/scan.jsonl' if len(sys.argv)==1 else '/dev/stdout'
    with open(outp,'w') as f, mp.Pool(8) as pool:
        for i,r in enumerate(pool.imap_unordered(scan,ids,chunksize=4)):
            f.write(json.dumps(r)+'\n')
            if i%250==0: print(i,time.strftime('%H:%M:%S'),file=sys.stderr,flush=True)
