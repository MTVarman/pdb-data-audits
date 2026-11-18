"""Fetch every polymer entity sequence for the union of candidate entries (RCSB GraphQL) and keep entities whose
sequence aligns to SARS-CoV-2 nsp5 (P0DTD1 3264-3569) better than to SARS-CoV-1 nsp5 (P0C6X7 3241-3546),
with >= 150 aligned nsp5 residues, >= 90% identity over the aligned region, and nsp5 position 145 covered.
Writes out/mpro_entities.tsv and out/mpro_entries.txt"""
import json, urllib.request, time
from Bio import Align
from Bio.Align import substitution_matrices
D='/workspace/mpro_covalent_census'
d=json.load(open(f'{D}/data/search_raw.json'))
ents=set(d['A'])|set(d['B']); entries=sorted({x.split('_')[0] for x in ents}|set(d['C']))
print('union candidate entries',len(entries), 'A',len({x.split('_')[0] for x in d['A']}),'B',len({x.split('_')[0] for x in d['B']}),'C',len(d['C']))
def fetch(u, data=None, hdr=None):
    for i in range(6):
        try:
            r=urllib.request.Request(u,data=data,headers=hdr or {})
            with urllib.request.urlopen(r,timeout=300) as f: return f.read()
        except Exception as e: err=e; time.sleep(3*(i+1))
    raise err
s1=json.loads(fetch('https://rest.uniprot.org/uniprotkb/P0C6X7.json'))['sequence']['value']
SARS2=d['nsp5']
# locate SARS-CoV-1 nsp5 by its start 'SGFRKMAF'
i=s1.find('SGFRKMAF'); SARS1=s1[i:i+306]; assert SARS1.endswith('VTFQ'),SARS1[-10:]
Q='''query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id polymer_entities{rcsb_id entity_poly{pdbx_seq_one_letter_code_can rcsb_entity_polymer_type pdbx_strand_id}
 rcsb_entity_source_organism{ncbi_taxonomy_id} rcsb_polymer_entity{pdbx_description}}}}'''
rows=[]
for k in range(0,len(entries),200):
    b=entries[k:k+200]
    r=json.loads(fetch('https://data.rcsb.org/graphql',json.dumps({'query':Q,'variables':{'ids':b}}).encode(),{'Content-Type':'application/json'}))
    for e in r['data']['entries'] or []:
        for pe in e['polymer_entities'] or []:
            ep=pe['entity_poly'] or {}
            rows.append((e['rcsb_id'],pe['rcsb_id'],ep.get('rcsb_entity_polymer_type'),ep.get('pdbx_strand_id'),(ep.get('pdbx_seq_one_letter_code_can') or '').replace('\n',''),
                ';'.join(str(o['ncbi_taxonomy_id']) for o in (pe['rcsb_entity_source_organism'] or [])),(pe['rcsb_polymer_entity'] or {}).get('pdbx_description')))
    print(k,len(rows),flush=True)
al=Align.PairwiseAligner(); al.mode='local'; al.substitution_matrix=substitution_matrices.load('BLOSUM62'); al.open_gap_score=-10; al.extend_gap_score=-1
def stat(seq,ref):
    a=al.align(ref,seq)[0]; ident=0; n=0; p145=None
    for (rs,re_),(ss,se) in zip(*a.aligned):
        for j in range(re_-rs):
            n+=1; ident+= ref[rs+j]==seq[ss+j]
            if rs+j==144: p145=(ss+j+1, seq[ss+j])
    return n,ident,p145
out=open(f'{D}/out/mpro_entities.tsv','w'); out.write('entry\tentity\tstrands\tlength\taligned\tident_sars2\tident_sars1\tseq_pos145\tres145\ttaxids\tdescription\n')
keep=set(); n_prot=0
for e,pe,typ,strands,seq,tax,desc in rows:
    if typ!='Protein' or len(seq)<100: continue
    n_prot+=1
    seq=''.join(c if c in 'ARNDCQEGHILKMFPSTWYVBZX' else 'X' for c in seq.upper())
    n2,i2,p=stat(seq,SARS2)
    if n2<150 or i2/n2<0.80: continue
    n1,i1,_=stat(seq,SARS1)
    if not (i2/n2>=0.90 and i2>=i1 and p): 
        out.write(f'#REJECT\t{e}\t{pe}\t{len(seq)}\t{n2}\t{i2/n2:.3f}\t{i1/max(n1,1):.3f}\t{p}\t{tax}\t{desc}\n'); continue
    keep.add(e)
    out.write(f'{e}\t{pe}\t{strands}\t{len(seq)}\t{n2}\t{i2/n2:.4f}\t{i1/n1:.4f}\t{p[0]}\t{p[1]}\t{tax}\t{desc}\n')
out.close()
open(f'{D}/out/mpro_entries.txt','w').write('\n'.join(sorted(keep))+'\n')
print('protein entities checked',n_prot,'Mpro entries',len(keep))
srcA={x.split('_')[0] for x in d['A']}; srcB={x.split('_')[0] for x in d['B']}; srcC=set(d['C'])
print('kept found by A',len(keep&srcA),'B',len(keep&srcB),'C',len(keep&srcC),'only A',len(keep-srcB-srcC),'only B',len(keep-srcA-srcC),'only C',len(keep-srcA-srcB))
