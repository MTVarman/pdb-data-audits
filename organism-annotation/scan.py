import json, gzip, sys, os, re, csv, time, urllib.request, urllib.parse
from concurrent.futures import ThreadPoolExecutor
import gemmi

ENTRIES=json.load(open("entries.json"))
GROUP=sys.argv[1]            # dog|cat|horse|human_sample
OUT=f"scan_{GROUP}.tsv"
os.makedirs("cif",exist_ok=True); os.makedirs("uniprot",exist_ok=True)

# --- SIFTS: pdb -> chain -> accession
sifts={}
with gzip.open("pdb_chain_uniprot.tsv.gz","rt") as f:
    for line in f:
        if line.startswith("#") or line.startswith("PDB"): continue
        p=line.rstrip("\n").split("\t")
        sifts.setdefault(p[0].upper(),{}).setdefault(p[1],set()).add(p[2])

# --- UniProt fetch (cached): reviewed flag, canonical seq, organism taxid
def _jdump(obj,fn):
    import threading, uuid
    tmp=fn+f".{os.getpid()}.{threading.get_ident()}.{uuid.uuid4().hex[:6]}.tmp"
    try:
        json.dump(obj,open(tmp,"w")); os.replace(tmp,fn)
    except FileNotFoundError:
        pass

def uniprot(acc):
    fn=f"uniprot/{acc}.json"
    if os.path.exists(fn):
        try: return json.load(open(fn))
        except Exception: pass
    url=f"https://rest.uniprot.org/uniprotkb/{acc}.json?fields=accession,reviewed,sequence,organism_id,gene_primary"
    for attempt in range(3):
        try:
            d=json.load(urllib.request.urlopen(url,timeout=60))
            rec={"acc":acc,"reviewed":d.get("entryType","").startswith("UniProtKB reviewed"),
                 "seq":d.get("sequence",{}).get("value",""),"taxid":d.get("organism",{}).get("taxonId"),
                 "gene":(d.get("genes") or [{}])[0].get("geneName",{}).get("value")}
            _jdump(rec,fn); return rec
        except Exception as e:
            if "404" in str(e):
                rec={"acc":acc,"reviewed":None,"seq":"","taxid":None,"gene":None,"missing":True}
                _jdump(rec,fn); return rec
            time.sleep(2)
    return {"acc":acc,"reviewed":None,"seq":"","taxid":None,"gene":None,"error":True}


def _get(url):
    for attempt in range(3):
        try: return json.load(urllib.request.urlopen(url,timeout=60))
        except Exception as e:
            if "404" in str(e): return None
            time.sleep(2)
    return None

def reviewed_for(acc, taxid=None, identity=0.9):
    """reviewed Swiss-Prot member of the unreviewed acc's UniRef90 cluster (fallback UniRef50); prefer same lineage; cached"""
    fn=f"uniprot/rev_{acc}.json"
    if os.path.exists(fn):
        try: return json.load(open(fn))
        except Exception: pass
    best=None
    for ident in (0.9,0.5):
        d=_get("https://rest.uniprot.org/uniref/search?query="+urllib.parse.quote(f"uniprot_id:{acc} AND identity:{ident}")+"&fields=id,members&format=json&size=1")
        res=(d or {}).get("results",[])
        if not res: continue
        rep=res[0]["id"].split("_",1)[1]
        members=[rep]+[m for m in res[0].get("members",[]) if not m.startswith("UPI") and m!=rep]
        for i in range(0,len(members),100):
            chunk=members[i:i+100]
            q="("+" OR ".join(f"accession:{m}" for m in chunk)+") AND reviewed:true"
            r=_get("https://rest.uniprot.org/uniprotkb/search?query="+urllib.parse.quote(q)+"&fields=accession,organism_id,lineage_ids&format=json&size=100")
            for e in (r or {}).get("results",[]):
                a=e["primaryAccession"]; oid=e.get("organism",{}).get("taxonId")
                lin=[x.get("taxonId") for x in e.get("organism",{}).get("lineage",[])] if isinstance(e.get("organism",{}).get("lineage"),list) else []
                same = taxid and (oid==taxid or taxid in lin)
                if same: best=a; break
                if best is None: best=a
            if best and taxid and same: break
        if best: break
    _jdump(best,fn); return best

def fetch_cif(pdb):
    fn=f"cif/{pdb}.cif.gz"
    if os.path.exists(fn): return fn
    url=f"https://files.rcsb.org/download/{pdb}.cif.gz"
    for attempt in range(3):
        try:
            urllib.request.urlretrieve(url,fn); return fn
        except Exception: time.sleep(2)
    return None

def align_subs(entity_seq, ref_seq):
    """Global alignment of entity (canonical one-letter) vs reference; return substitutions list [(ref_pos,ref_aa,ent_aa)], aligned length, identity"""
    import gemmi
    res=gemmi.align_string_sequences(list(entity_seq), list(ref_seq), [], gemmi.AlignmentScoring("b"))
    subs=[]; aligned=0; ident=0
    i=j=0
    for c in res.cigar_str():
        pass
    # walk cigar
    cig=res.cigar_str()
    for n,op in re.findall(r"(\d+)([MID])",cig):
        n=int(n)
        if op=="M":
            for k in range(n):
                a=entity_seq[i]; b=ref_seq[j]
                aligned+=1
                if a==b: ident+=1
                elif a!="X" and b!="X": subs.append((j+1,b,a))
                i+=1; j+=1
        elif op=="I": i+=n
        else: j+=n
    return subs, aligned, ident

def scan(pdb):
    pdb=pdb.upper()
    fn=fetch_cif(pdb)
    if not fn: return [{"pdb":pdb,"error":"fetch"}]
    try:
        doc=gemmi.cif.read(fn); blk=doc.sole_block()
    except Exception as e:
        return [{"pdb":pdb,"error":"parse"}]
    # entity info
    ent_mut={}; ent_desc={}
    t=blk.find("_entity.",["id","pdbx_mutation","pdbx_description"])
    for r in t: ent_mut[r[0]]=gemmi.cif.as_string(r[1]); ent_desc[r[0]]=gemmi.cif.as_string(r[2])
    ent_seq={}; ent_chains={}
    t=blk.find("_entity_poly.",["entity_id","type","pdbx_seq_one_letter_code_can","pdbx_strand_id"])
    for r in t:
        if "polypeptide" not in gemmi.cif.as_string(r[1]): continue
        ent_seq[r[0]]=gemmi.cif.as_string(r[2]).replace("\n","").replace(" ","").upper()
        ent_chains[r[0]]=gemmi.cif.as_string(r[3]).split(",")
    # organism per entity
    ent_tax={}
    for cat in ("_entity_src_gen.","_entity_src_nat.","_pdbx_entity_src_syn."):
        key={"_entity_src_gen.":"pdbx_gene_src_ncbi_taxonomy_id","_entity_src_nat.":"pdbx_ncbi_taxonomy_id","_pdbx_entity_src_syn.":"ncbi_taxonomy_id"}[cat]
        t=blk.find(cat,["entity_id",key])
        for r in t: ent_tax.setdefault(r[0],gemmi.cif.as_string(r[1]))
    # depositor struct_ref
    dep_ref={}
    t=blk.find("_struct_ref.",["id","entity_id","db_name","pdbx_db_accession"])
    for r in t:
        if gemmi.cif.as_string(r[2]).upper()=="UNP": dep_ref.setdefault(r[1],gemmi.cif.as_string(r[3]))
    # struct_ref_seq_dif
    dif={}
    t=blk.find("_struct_ref_seq_dif.",["pdbx_pdb_id_code","mon_id","db_mon_id","pdbx_seq_db_accession_code","details","pdbx_auth_seq_num"])
    for r in t:
        d=gemmi.cif.as_string(r[4]).upper()
        dif.setdefault("all",[]).append((gemmi.cif.as_string(r[1]),gemmi.cif.as_string(r[2]),d))
    dif_subs=[x for x in dif.get("all",[]) if x[0]!=x[1] and x[0] not in("?",".") and x[1] not in("?",".") and not any(k in x[2] for k in("TAG","CLONING","EXPRESSION","LINKER","INSERTION","DELETION","INITIATING","LEADER"))]
    rows=[]
    for eid,seq in ent_seq.items():
        chains=ent_chains.get(eid,[])
        accs=set()
        for ch in chains: accs|=sifts.get(pdb,{}).get(ch,set())
        if not accs: 
            rows.append({"pdb":pdb,"entity":eid,"taxid":ent_tax.get(eid),"sifts_acc":"","dep_acc":dep_ref.get(eid,""),"reviewed":"","status":"no_sifts","declared_mut":ent_mut.get(eid,""),"n_dif":len(dif_subs),"len":len(seq)})
            continue
        acc=sorted(accs)[0]
        up=uniprot(acc)
        if not up["seq"]:
            rows.append({"pdb":pdb,"entity":eid,"taxid":ent_tax.get(eid),"sifts_acc":acc,"dep_acc":dep_ref.get(eid,""),"reviewed":up.get("reviewed"),"status":"uniprot_missing","declared_mut":ent_mut.get(eid,""),"n_dif":len(dif_subs),"len":len(seq)}); continue
        subs,aligned,ident=align_subs(seq,up["seq"])
        declared=ent_mut.get(eid,"") or ""
        decl_set=set(re.findall(r"[A-Z]\d+[A-Z]",declared.upper().replace(" ","")))
        undeclared=[f"{b}{p}{a}" for p,b,a in subs if f"{b}{p}{a}" not in decl_set]
        label="WT" if declared.strip() in("","?",".") else "MUT"
        rev_acc=""; n_rev=""; rev_und=""
        if up["reviewed"] is False:
            ra=reviewed_for(acc, up["taxid"])
            if ra:
                rv=uniprot(ra)
                if rv["seq"]:
                    s2,a2,i2=align_subs(seq,rv["seq"])
                    und2=[f"{b}{p}{a}" for p,b,a in s2 if f"{b}{p}{a}" not in decl_set]
                    rev_acc=ra; n_rev=len(und2); rev_und=";".join(und2[:20])
        rows.append({"rev_acc":rev_acc,"n_undeclared_vs_reviewed":n_rev,"undeclared_vs_reviewed":rev_und,"pdb":pdb,"entity":eid,"taxid":ent_tax.get(eid),"sifts_acc":acc,"dep_acc":dep_ref.get(eid,""),"reviewed":up["reviewed"],"acc_taxid":up["taxid"],"gene":up["gene"],
                     "status":"ok","label":label,"declared_mut":declared,"n_subs_vs_ref":len(subs),"n_undeclared":len(undeclared),"undeclared":";".join(undeclared[:20]),
                     "aligned":aligned,"identity":round(ident/aligned,4) if aligned else 0,"n_dif":len(dif_subs),"len":len(seq)})
    return rows

ids=ENTRIES[GROUP]
done=set()
if os.path.exists(OUT):
    for r in csv.DictReader(open(OUT),delimiter="\t"): done.add(r["pdb"])
todo=[i for i in ids if i.upper() not in done]
print(GROUP,len(ids),"todo",len(todo),file=sys.stderr)
fields=["pdb","entity","taxid","sifts_acc","dep_acc","reviewed","acc_taxid","gene","status","label","declared_mut","n_subs_vs_ref","n_undeclared","undeclared","aligned","identity","rev_acc","n_undeclared_vs_reviewed","undeclared_vs_reviewed","n_dif","len","error"]
newfile=not os.path.exists(OUT)
with open(OUT,"a") as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",extrasaction="ignore")
    if newfile: w.writeheader()
    with ThreadPoolExecutor(int(os.environ.get("THREADS","8"))) as ex:
        for n,rows in enumerate(ex.map(scan,todo)):
            for r in rows: w.writerow(r)
            if n%50==0: f.flush(); print(n,file=sys.stderr)
print("done",GROUP,file=sys.stderr)
