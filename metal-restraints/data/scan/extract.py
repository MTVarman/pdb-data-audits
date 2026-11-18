import gzip,json,re,sys,os,time
from concurrent.futures import ProcessPoolExecutor
IDX=os.path.expanduser('~/mnt/morphogenesis/results/2026-10-31_whole-pdb-index/index')
OUT=os.path.expanduser('~/metalc')
METALS=set("ZN CA MG NA MN K FE CU CD HG NI CO SR CS PT BA TL PB SM AU RB YB LI MO LU CR OS GD TB LA AG HO GA CE W RU RE PR IR EU AL V PD U".split())
STD=set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL A C G U DA DC DG DT HOH MSE".split())
row=re.compile(r'\["MetalC",[^\]]*\]')
hdr={k:re.compile(r'"%s":("?[^,"]*"?)'%k) for k in ('method','resolution','deposited')}
def el_of(name):
    n=re.sub(r'[^A-Z]','',name.upper())
    return n
def side(res,atom):
    a=el_of(atom)
    if res not in STD:
        if a in METALS: return a
        if a[:2] in METALS and len(a)>=2 and a[:2] not in ('CA',) : return a[:2]
    return None
def lig_el(res,atom):
    a=atom.upper()
    if res in ('CL','BR','IOD','F'): return 'X'
    if a.startswith('SE'): return 'SE'
    return a[0] if a[0] in 'NOSCP' else '?'
def do(n):
    p=f'{IDX}/shard_{n:02d}.jsonl.gz'; o=f'{OUT}/shard_{n:02d}.tsv'; k=0; ne=0
    with gzip.open(p,'rt') as f, open(o+'.tmp','w') as w:
        for line in f:
            ne+=1
            if '"MetalC"' not in line: continue
            pid=line[10:14]
            head=line[:3000]
            m=hdr['method'].search(head); meth=(m.group(1).strip('"') if m else '')
            r=hdr['resolution'].search(head); res=(r.group(1).strip('"') if r else '')
            d=hdr['deposited'].search(head); yr=(d.group(1).strip('"')[:4] if d else '')
            mt='E' if meth.startswith('ELECTRON MICRO') else ('X' if meth.startswith('X-RAY') else 'O')
            for mm in row.finditer(line):
                try: x=json.loads(mm.group(0))
                except Exception: continue
                if len(x)<15: continue
                r1,a1,r2,a2,dist=x[2],x[5],x[8],x[11],x[14]
                m1=side(r1,a1); m2=side(r2,a2)
                if m1 and not m2: M,LR,LA=m1,r2,a2; MR=r1; LS=x[9]; LC=x[7]; MS=x[3]; MC=x[1]
                elif m2 and not m1: M,LR,LA=m2,r1,a1; MR=r2; LS=x[3]; LC=x[1]; MS=x[9]; MC=x[7]
                else: continue
                if not isinstance(dist,(int,float)): continue
                w.write(f"{pid}\t{mt}\t{res}\t{yr}\t{M}\t{MR}\t{lig_el(LR,LA)}\t{LR}\t{LA}\t{dist}\t{LC}\t{LS}\t{MC}\t{MS}\n"); k+=1
    os.replace(o+'.tmp',o); return n,ne,k
if __name__=='__main__':
    t=time.time(); ns=[int(a) for a in sys.argv[1:]]
    with ProcessPoolExecutor(3) as ex:
        for n,ne,k in ex.map(do,ns): print('shard',n,'entries',ne,'metalc rows',k,flush=True)
    print('sec',round(time.time()-t))
