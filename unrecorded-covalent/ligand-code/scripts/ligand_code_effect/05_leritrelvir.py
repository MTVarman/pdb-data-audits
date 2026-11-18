"""Leritrelvir case: per entry and per ligand copy, CCD code, Cys145 SG-C7 distance (nearest image), geometry at
warhead carbon C7 (C7-O6 length, sum of the three heavy-atom bond angles at C7: ~360 deg = planar ketone,
~328 deg = tetrahedral), whether a covale row exists, close-contact listing, has_protein_modification,
modification_feature rows, dates, site, authors."""
import gemmi, glob, json, math, pandas as pd
O='/workspace/ligand_code_effect'
CC=json.load(open(f'{O}/out/ccd_parsed.json'))
def v(x):
    s=gemmi.cif.as_string(x) if x is not None else None; return None if s in (None,'?','.','') else s
def ang(a,b,c):
    return math.degrees(gemmi.calculate_angle(a,b,c))
rows=[]
for p in sorted(glob.glob(f'{O}/data/cif/*.cif.gz')):
    doc=gemmi.cif.read(p); b=doc.sole_block(); e=b.name.upper()
    g=lambda t:v(b.find_value(t))
    st=gemmi.make_structure_from_block(b); st.setup_entities(); st.remove_hydrogens(); m=st[0]
    ns=gemmi.NeighborSearch(m,st.cell,5).populate()
    covs=[]; t=b.find('_struct_conn.',['conn_type_id','ptnr1_label_comp_id','ptnr1_label_atom_id','ptnr1_auth_asym_id','ptnr2_label_comp_id','ptnr2_label_atom_id','ptnr2_auth_asym_id','pdbx_dist_value'])
    for r in t: covs.append([v(r[i]) for i in range(8)])
    cc=b.find('_pdbx_validate_close_contact.',['auth_comp_id_1','auth_atom_id_1','auth_comp_id_2','auth_atom_id_2','dist'])
    ccl=[[v(r[i]) for i in range(5)] for r in cc]
    mf=len(b.find('_pdbx_modification_feature.',['type'])) if b.find('_pdbx_modification_feature.',['type']) else 0
    revs=[v(x) for x in b.find_values('_pdbx_audit_revision_history.revision_date')]
    au=[v(x) for x in b.find_values('_audit_author.name')]
    for ch in m:
        for res in ch:
            if res.name not in ('7ON','A1EN0'): continue
            at={a.name:a for a in res}
            if 'C7' not in at: continue
            c7=at['C7'].pos
            best=None
            for mk in ns.find_atoms(c7,'\0',radius=4.5):
                cra=mk.to_cra(m)
                if cra.residue.name=='CYS' and cra.atom.name=='SG' and cra.residue.seqid.num==145:
                    dd=st.cell.find_nearest_image(c7,cra.atom.pos).dist()
                    if best is None or dd<best[0]: best=(dd,cra.chain.name)
            nbrs=[n for n in ('C6','C8','O6') if n in at]
            # find heavy neighbours of C7 by distance (<1.7 A within residue)
            hn=[a for a in res if a.name!='C7' and a.pos.dist(c7)<1.75]
            angs=[ang(hn[i].pos,c7,hn[j].pos) for i in range(len(hn)) for j in range(i+1,len(hn))]
            rec=any(x[0] and x[0].startswith('covale') and ((x[1]==res.name and x[2]=='C7') or (x[4]==res.name and x[5]=='C7')) for x in covs)
            recany=any(x[0] and x[0].startswith('covale') and res.name in (x[1],x[4]) for x in covs)
            inclose=any(((x[0]==res.name and x[1]=='C7' and x[3]=='SG') or (x[2]==res.name and x[3]=='C7' and x[1]=='SG')) for x in ccl)
            r=CC[res.name]
            rows.append(dict(entry=e,lig_comp=res.name,lig_chain=ch.name,lig_seq=res.seqid.num,sg145_c7_dist=round(best[0],3) if best else None,
                c7_o6=round(at['C7'].pos.dist(at['O6'].pos),3) if 'O6' in at else None,
                c7_n_heavy_nb=len(hn),c7_angle_sum=round(sum(angs),1) if len(hn)==3 else None,
                covale_row_c7=rec,covale_row_any=recany,sg_c7_in_close_contact=inclose,
                has_protein_modification=g('_pdbx_entry_details.has_protein_modification'),n_modification_feature=mf,
                process_site=g('_pdbx_database_status.process_site'),deposit_site=g('_pdbx_database_status.deposit_site'),
                deposit_date=g('_pdbx_database_status.recvd_initial_deposition_date'),initial_release=min(revs),last_revision=max(revs),
                last_author=au[-1],resolution=g('_refine.ls_d_res_high'),title=(g('_struct.title') or '')[:110],
                ccd_formula=r['formula'],ccd_n_atoms=r['n_atoms'],ccd_n_heavy=r['n_heavy'],ccd_leaving_atoms=';'.join(r['leaving_atoms']) or 'none',
                ccd_initial_date=r['pdbx_initial_date'],ccd_release_date=r['ccd_release_date'],ccd_create_site=r['create_site'],ccd_defining_entry=r['pdbx_model_coordinates_db_code'],
                ccd_type=r['type'],ccd_pdbx_type=r['pdbx_type'],ccd_pcm=r['pdbx_pcm'] or 'blank',ccd_pcm_link=';'.join(x[1] for x in r['pcm_rows'] if x[1]) or 'none',
                ccd_name=r['name'][:60]))
d=pd.DataFrame(rows); d.to_csv(f'{O}/out/leritrelvir_copies.tsv',sep='\t',index=False)
pd.set_option('display.width',250)
print(d[['entry','lig_comp','lig_chain','sg145_c7_dist','c7_o6','c7_angle_sum','covale_row_c7','covale_row_any','sg_c7_in_close_contact','has_protein_modification','n_modification_feature','deposit_date','initial_release','last_author']].to_string())
