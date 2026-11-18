import collections,statistics as st,csv,json
sites=collections.defaultdict(list); meta={}
for line in open('all.tsv'):
    pid,mt,res,yr,M,MR,X,LR,LA,d,LC,LS,MC,MS=line.rstrip('\n').split('\t'); d=float(d)
    if d>3.5 or d<1.2: continue
    meta[pid]=(mt,res,yr); sites[(pid,M,X,LR,LS,LA,MR,MS)].append(d)
S=[(k[0],k[1],k[2],k[3],k[6],st.median(v)) for k,v in sites.items()]   # pid,M,X,LR,MR,d
print('unique sites',len(S))
def pooled(M,X,T,h=0.005,label=''):
    for mt in ('E','X'):
        v=[(p,d) for (p,m,x,lr,mr,d) in S if m==M and (X is None or x==X) and meta[p][0]==mt]
        inw=[(p,d) for p,d in v if abs(d-T)<=h]; fl=sum(1 for p,d in v if 0.015<=abs(d-T)<=0.035)
        exp=fl*(2*h)/0.040; ents=collections.Counter(p for p,d in inw)
        print(f"  {M}-{X or 'any':3s} @ {T:.3f} {mt}: sites {len(v):7d}  in window {len(inw):6d}  expected {exp:8.1f}  excess {len(inw)-exp:8.1f}  entries w/ >=1 {len(ents):5d}, >=2 {sum(1 for c in ents.values() if c>=2):4d}, >=3 {sum(1 for c in ents.values() if c>=3):4d}")
for M,X,T in [('MG','N',1.897),('MG','N',2.303),('ZN','N',2.301),('ZN','O',2.301),('NI',None,1.973),('NI',None,2.301),('CD',None,2.301),('HG',None,2.301),('CA','N',2.25),('CA','N',2.301),('NA','N',2.08),('NA','N',2.301),('K','N',2.444),('MN','N',2.301),('FE','N',2.301),('CU','N',2.301),('CO','N',2.301),('BA','O',2.75),('SR','O',2.40),('MG','O',2.301),('MG','O',2.18),('MG','O',1.95),('ZN','N',2.25),('MG','N',2.25)]:
    pooled(M,X,T)
# union of clearly-wrong pinned entries from pinned.tsv
R=list(csv.DictReader(open('pinned.tsv'),delimiter='\t'))
def sel(metal,els,lo,hi): return {(r['pdb'],r['method'],r['year']) for r in R if r['metal']==metal and (els is None or r['lig_el'] in els) and lo<=float(r['center'])<=hi and float(r['sd_in_window'])<=0.012}
G={'Mg-N 1.895':sel('MG','N',1.885,1.905),'Mg-N 2.30':sel('MG','N',2.291,2.315),'Zn-N 2.30':sel('ZN','N',2.291,2.315),'Ni 1.972':sel('NI',None,1.962,1.984),'Cd 2.30':sel('CD',None,2.291,2.315),'Ba-O 2.75':sel('BA','O',2.74,2.76),'Na-O 2.30':sel('NA','O',2.291,2.315)}
U=set().union(*G.values())
for k,v in G.items(): print(k,'entries',len(v),'EM',sum(1 for e in v if e[1]=='E'),'Xray',sum(1 for e in v if e[1]=='X'),'other',sum(1 for e in v if e[1]=='O'))
print('union entries',len({e[0] for e in U}),'by year',sorted(collections.Counter(e[2] for e in U).items()))
# Mg-N 1.895: which residues carry the Mg
mgp={e[0] for e in G['Mg-N 1.895']}
c=collections.Counter(); 
for (p,m,x,lr,mr,d) in S:
    if p in mgp and m=='MG' and x=='N' and abs(d-1.897)<=0.012: c[mr]+=1
print('Mg-N 1.895 sites by metal residue:',c.most_common(8)); ce=collections.defaultdict(set)
for (p,m,x,lr,mr,d) in S:
    if p in mgp and m=='MG' and x=='N' and abs(d-1.897)<=0.012: ce[mr].add(p)
print(' entries by metal residue:',{k:len(v) for k,v in ce.items()})
json.dump({k:sorted(v) for k,v in G.items()},open('fallback_groups.json','w'))
