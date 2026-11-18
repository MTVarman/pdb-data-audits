import collections,json,statistics as st,sys,os
QM=json.loads('{"CA|N": 2.25, "CA|S": 2.712, "CA|SE": 2.807, "CA|X": 2.578, "CA|C": 2.519, "MG|N": 1.895, "MG|S": 2.318, "MG|SE": 2.413, "MG|X": 2.213, "MG|C": 2.105, "NA|N": 2.08, "NA|S": 2.514, "NA|X": 2.397, "NA|C": 2.328, "K|N": 2.444, "K|S": 2.923, "K|X": 2.775, "K|C": 2.761, "NI|N": 1.972, "NI|O": 1.972, "NI|S": 1.972, "NI|X": 1.972, "SR|N": 2.493, "SR|O": 2.4, "SR|S": 2.831, "CS|N": 2.75, "CS|O": 2.75, "CS|S": 2.75, "PT|N": 2.75, "PT|O": 2.75, "PT|S": 2.75, "BA|N": 2.75, "BA|O": 2.75, "BA|S": 2.75, "TL|N": 2.75, "TL|O": 2.75, "TL|S": 2.75, "PB|N": 2.75, "PB|O": 2.75, "PB|S": 2.75, "AU|N": 2.75, "AU|O": 2.75, "AU|S": 2.75, "RB|N": 2.565, "RB|O": 2.45, "RB|S": 2.903, "LI|N": 1.749, "LI|O": 1.59, "LI|S": 2.189, "MO|N": 2.75, "MO|O": 2.75, "MO|S": 2.75, "CR|N": 2.146, "CR|O": 2.146, "CR|S": 2.146, "OS|N": 2.11, "OS|O": 2.11, "LA|N": 2.75, "LA|O": 2.75, "AG|N": 2.75, "AG|O": 2.75, "AG|S": 2.75, "GA|N": 1.816, "GA|O": 1.781, "GA|S": 2.206, "W|N": 2.75, "W|O": 2.75, "W|S": 2.75, "RU|N": 2.75, "RU|O": 2.75, "RU|S": 2.75, "RE|N": 2.75, "RE|O": 2.75, "IR|N": 2.75, "IR|O": 2.75, "AL|N": 1.77, "AL|O": 1.693, "AL|S": 2.195, "V|N": 2.103, "V|O": 2.103, "V|S": 2.103, "PD|N": 2.75, "PD|O": 2.75, "PD|S": 2.75}')
T1=set('NA MG K CA'.split()); T2=set('MN FE CO CU ZN'.split())
def fallback(M,X):
    if M in T1: return X!='O'
    if M in T2: return X not in ('O','N','S')
    return True
W=0.03; NMIN=int(sys.argv[1]) if len(sys.argv)>1 else 5; FMIN=float(sys.argv[2]) if len(sys.argv)>2 else 0.7
sites=collections.defaultdict(list); meta={}
for line in open('all.tsv'):
    pid,mt,res,yr,M,MR,X,LR,LA,d,LC,LS,MC,MS=line.rstrip('\n').split('\t'); d=float(d)
    if d>3.5 or d<1.2: continue
    meta[pid]=(mt,res,yr)
    sites[(pid,M,X,LR,LS,LA,MR,MS)].append(d)
groups=collections.defaultdict(list)
for (pid,M,X,LR,LS,LA,MR,MS),v in sites.items(): groups[(pid,M,X)].append(st.median(v))
tested=collections.Counter(); pinned=[]
for (pid,M,X),v in groups.items():
    n=len(v)
    if n<NMIN: continue
    mt=meta[pid][0]; tested[(M,X,mt)]+=1
    v.sort(); best=(0,0,0); j=0
    for i in range(n):
        while v[i]-v[j]>W: j+=1
        if i-j+1>best[0]: best=(i-j+1,j,i)
    c,j,i=best
    if c>=NMIN and c/n>=FMIN:
        w=v[j:i+1]; cen=st.median(w); sd=st.pstdev(w)
        fb=fallback(M,X); q=QM.get(f'{M}|{X}')
        cls='other'
        if fb and q is not None and abs(cen-q)<=0.012: cls='fallback_qm'
        elif fb and abs(cen-2.303)<=0.012: cls='fallback_2.30'
        elif not fb: cls='table_pair'
        pinned.append((pid,mt,meta[pid][1],meta[pid][2],M,X,n,c,round(cen,3),round(sd,4),cls))
out=os.path.expanduser('~/metalc/pinned.tsv')
with open(out,'w') as f:
    f.write('pdb\tmethod\tres\tyear\tmetal\tlig_el\tn_sites\tn_in_window\tcenter\tsd_in_window\tclass\n')
    for r in sorted(pinned): f.write('\t'.join(map(str,r))+'\n')
print('groups tested (>=%d sites):'%NMIN,sum(tested.values()),'| pinned groups:',len(pinned),'| entries with any pinned group:',len({r[0] for r in pinned}))
for cls in ('fallback_qm','fallback_2.30','other','table_pair'):
    r=[x for x in pinned if x[10]==cls]
    print(f' {cls:14s} groups {len(r):5d} entries {len({x[0] for x in r}):5d}  EM {len({x[0] for x in r if x[1]=="E"}):4d}  Xray {len({x[0] for x in r if x[1]=="X"}):4d}')
print('\nfallback hits by metal|ligand element (entries, EM/Xray):')
for cls in ('fallback_qm','fallback_2.30'):
    c=collections.defaultdict(set)
    for x in pinned:
        if x[10]==cls: c[(x[4],x[5])].add((x[0],x[1]))
    print(' ',cls,sorted([(k[0]+'-'+k[1],len(v),sum(1 for e in v if e[1]=="E"),sum(1 for e in v if e[1]=="X")) for k,v in c.items()],key=lambda t:-t[1])[:25])
print('\ntop pinned values, fallback pairs only (metal-el, center, entries):')
c=collections.defaultdict(set)
for x in pinned:
    if fallback(x[4],x[5]): c[(x[4],x[5],round(x[8],2))].add(x[0])
print(' ',sorted([(k[0]+'-'+k[1],k[2],len(v)) for k,v in c.items()],key=lambda t:-t[2])[:40])
print('\ntop pinned values, table pairs (metal-el, center, entries):')
c=collections.defaultdict(set)
for x in pinned:
    if not fallback(x[4],x[5]): c[(x[4],x[5],round(x[8],2))].add(x[0])
print(' ',sorted([(k[0]+'-'+k[1],k[2],len(v)) for k,v in c.items()],key=lambda t:-t[2])[:30])
print('\nfallback hits by year:',sorted(collections.Counter(y for (p,y) in {(x[0],x[3]) for x in pinned if x[10].startswith('fallback')}).items()))
