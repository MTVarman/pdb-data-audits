import csv,gzip,re,collections,gemmi,xml.etree.ElementTree as ET
rows=list(csv.DictReader(open('archive_removed_811_sites.csv')))
COF=set("""GDP GTP ATP ADP AMP ANP GNP GSP GCP AGS ACP APC CTP UTP TTP DTP CMP UMP NAD NAI NAP NDP NADP FAD FMN HEM HEC HEA HEB COA ACO SAM SAH PLP PMP TPP TDP BTN B12 CNC SF4 FES F3S CLA BCL CYC PEB PUB BLA BLR RET MGD PQQ""".split())
ADD=set("EDO GOL PEG PG4 PGE 1PE DMS ACT FMT ACE MPD BME DTT TRS MES HEPES EPE IMD CIT SO4 PO4 NO3 FLC TAR".split())
def elem(cif_block):pass
out=[]
cache={}
def load(e):
    if e in cache: return cache[e]
    doc=gemmi.cif.read(f'cif/{e}.cif.gz'); b=doc.sole_block()
    st=gemmi.make_structure_from_block(b); st.setup_entities()
    desc={r[0]:r[1] for r in b.find('_entity.',['id','pdbx_description'])}
    unp={}
    t=b.find('_struct_ref.',['entity_id','db_name','pdbx_db_accession'])
    for r in t:
        if r[1]=='UNP': unp.setdefault(r[0],r[2])
    cct={}
    for r in b.find('_chem_comp.',['id','type','name']): cct[r[0]]=(r[1],r[2])
    meth=b.find_value('_exptl.method') or ''
    res=b.find_value('_refine.ls_d_res_high') or b.find_value('_em_3d_reconstruction.resolution') or ''
    val={}
    try:
        t2=ET.parse(gzip.open(f'val/{e}.xml.gz')).getroot()
        for m in t2.iter('ModelledSubgroup'):
            if m.get('model','1')!='1': continue
            val[(m.get('chain'),m.get('resnum'),m.get('resname'))]=(m.get('rscc'),m.get('rsr'))
    except Exception as ex: pass
    cache[e]=(st,desc,unp,cct,meth.strip("'"),res,val); return cache[e]
for x in rows:
    e=x['entry'].lower(); st,desc,unp,cct,meth,resol,val=load(e)
    m=st[0]; ch=x['prot_chain']
    pr,pa=x['prot'].split(':'); lr,la=x['lig'].split(':')
    prn,prnum=re.match(r'([A-Z0-9]+?)(-?\d+)$',pr).groups()
    P=L=None; pres=lres=pchain=lchain=None
    for c in m:
        for r in c:
            key=r.name+str(r.seqid.num)
            if c.name==ch and key==pr and P is None:
                a=r.find_atom(pa,'*')
                if a: P,pres,pchain=a,r,c
            if key==lr:
                a=r.find_atom(la,'*')
                if a and (L is None or (P is not None and a.pos.dist(P.pos)<L.pos.dist(P.pos))): L,lres,lchain=a,r,c
    if P is None or L is None:
        out.append(dict(note='atom_not_found',**x)); continue
    d=P.pos.dist(L.pos); pe,le=P.element.name,L.element.name
    rc=gemmi.Element(pe).covalent_r+gemmi.Element(le).covalent_r
    lname=lres.name; ltype,lfull=cct.get(lname,('',''))
    het=pe in('S','O','N','Se') ; 
    if pe=='C' and le=='C': chem='C-C_inspect'
    elif het and le in('C','B','P','Si','S','Se') and d<=rc+0.4: chem='plausible_bond'
    elif het and le in('O','N') : chem='heteroatom_pair_not_bond'
    else: chem='other'
    lt=ltype.lower()
    if lname in COF or 'linking' in lt or 'nucleotide' in lfull.lower(): lcat='nucleotide_or_cofactor'
    elif 'saccharide' in lt: lcat='sugar'
    elif lname in ADD: lcat='additive'
    else: lcat='ligand'
    pv=val.get((pchain.name,str(pres.seqid.num),pres.name),(None,None)); lv=val.get((lchain.name,str(lres.seqid.num),lres.name),(None,None))
    def ok(v):
        try: return float(v[0])>=0.8 and (v[1] is None or float(v[1])<=0.3)
        except: return None
    po,lo=ok(pv),ok(lv)
    dens='supported' if po and lo else ('unknown' if po is None or lo is None else 'not_supported')
    ent=pchain.get_subchain(pres.subchain).entity_id if False else None
    eid=None
    for en in st.entities:
        if pres.subchain in en.subchains: eid=en.name
    out.append(dict(entry=x['entry'],chain=ch,residue=pres.name,resnum=pres.seqid.num,prot_atom=pa,prot_el=pe,
      ligand=lname,lig_atom=la,lig_el=le,dist=round(d,3),rcov_sum=round(rc,2),chem=chem,lig_category=lcat,lig_name=lfull[:60],
      protein=desc.get(eid,'')[:60],uniprot=unp.get(eid,''),method=meth,resolution=resol,
      prot_rscc=pv[0],prot_rsr=pv[1],lig_rscc=lv[0],lig_rsr=lv[1],density=dens,removed_on=x['latest_revision']))
keys=list(out[0].keys())
with open('removed_811_classified.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=keys,extrasaction='ignore'); w.writeheader(); w.writerows([o for o in out if 'chem' in o])
bad=[o for o in out if 'chem' not in o]; print('not found',len(bad))
ok_=[o for o in out if 'chem' in o]
C=collections.Counter
print(C(o['chem'] for o in ok_)); print(C(o['lig_category'] for o in ok_)); print(C(o['density'] for o in ok_))
core=[o for o in ok_ if o['chem']=='plausible_bond' and o['lig_category']=='ligand']
print('plausible+ligand',len(core), C(o['residue'] for o in core).most_common())
print('… density supported', C(o['residue'] for o in core if o['density']=='supported').most_common())
