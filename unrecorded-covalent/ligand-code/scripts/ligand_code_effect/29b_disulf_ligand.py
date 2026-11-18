"""Count _struct_conn disulf rows that join a Cys SG to a non-polymer ligand atom (ligand S-Cys S pairs the
distance-table note counts as recorded, while scripts 12-29 count only covale rows). Writes out/step6/29b_disulf_ligand.tsv."""
import gemmi, json, sys, os, gzip, re
RX=re.compile(r'^loop_\s*\n(?=(#|loop_|data_|\Z))',re.M)   # empty loops in header files, as in script 12
from multiprocessing import Pool
H='/workspace/ligand_code_effect/archive/headers'
files=[l.strip() for l in open(sys.argv[1])]
NP={}
for l in open('/workspace/ligand_code_effect/archive/parsed.jsonl'):
    r=json.loads(l)
    if 'np' in r: NP[r['e']]=set(r['np'])
def f(fn):
    e=fn.split('.')[0]; out=[]
    try: b=gemmi.cif.read_string(RX.sub('',gzip.open(os.path.join(H,fn),'rt').read())).sole_block()
    except Exception: return [['ERR',e]]
    t=b.find('_struct_conn.',['conn_type_id','ptnr1_label_comp_id','ptnr1_label_atom_id','ptnr1_auth_asym_id','ptnr1_auth_seq_id',
                              'ptnr2_label_comp_id','ptnr2_label_atom_id','ptnr2_auth_asym_id','ptnr2_auth_seq_id','pdbx_dist_value'])
    np_=NP.get(e,set()); nd=[0]
    for x in t:
        if x[0]!='disulf': continue
        nd[0]+=1
        for a,bb in ((1,5),(5,1)):
            if x[a]=='CYS' and x[a+1]=='SG' and x[bb] in np_:
                out.append([e,x[a+2],x[a+3],x[bb],x[bb+1],x[bb+2],x[bb+3],x[9]])
    out.append(['N',e,nd[0]]); return out
with Pool(8) as p: R=[y for z in p.map(f,files,chunksize=200) for y in z]
nd=sum(r[2] for r in R if r[0]=='N'); R=[r for r in R if r[0]!='N']; print('disulf rows read',nd)
print('parse errors',sum(r[0]=='ERR' for r in R)); R=[r for r in R if r[0]!='ERR']
with open('out/step6/29b_disulf_ligand.tsv','w') as o:
    o.write('entry\tnuc_chain\tnuc_seq\tlig_comp\tlig_atom\tlig_chain\tlig_seq\tdist\n')
    for r in R: o.write('\t'.join(map(str,r))+'\n')
print('disulf rows joining Cys SG to a non-polymer ligand:',len(R),'entries',len({r[0] for r in R}),'codes',len({r[3] for r in R}))
