import collections,json
site=collections.defaultdict(list); meta={}
for line in open('all.tsv'):
    pid,mt,res,yr,M,MR,X,LR,LA,d,LC,LS,MC,MS=line.rstrip('\n').split('\t')
    if M!='ZN' or MR!='ZN': continue
    d=float(d)
    if d>3.0: continue
    meta[pid]=mt; site[(pid,MC,MS)].append((LR,LA,d))
print('zinc sites',len(site))
def w(v,t,h=0.004): return sum(1 for x in v if abs(x-t)<=h)
for mt,name in (('E','cryo-EM'),('X','X-ray')):
    print('\n==',name)
    cls=collections.defaultdict(lambda:{'N':[],'S':[],'n':0})
    for (pid,mc,ms),L in site.items():
        if meta[pid]!=mt: continue
        nC=sum(1 for lr,la,d in L if la=='SG'); nH=len({(lr,la) for lr,la,d in L if la in('ND1','NE2')}) and sum(1 for lr,la,d in L if la in('ND1','NE2'))
        k=(nC,nH) if nC+nH<=4 else ('>4',)
        c=cls[k]; c['n']+=1
        c['N']+=[d for lr,la,d in L if la in('ND1','NE2')]; c['S']+=[d for lr,la,d in L if la=='SG']
    print(' class(nCys,nHis)  sites |  Zn-N: n  @2.300 @2.002 @2.040 @2.074 | Zn-S: n  @2.330 @2.318 @2.306 @2.298 @2.300')
    for k in [(4,0),(3,1),(2,2),(1,3),(0,4),(3,0),(2,1),(1,2),(0,3),(2,0),(1,1),(0,2),(1,0),(0,1)]:
        c=cls.get(k)
        if not c: continue
        N,S=c['N'],c['S']
        print(f"  {str(k):8s} {c['n']:6d} | {len(N):6d}  {w(N,2.300):5d} {w(N,2.002):5d} {w(N,2.040):5d} {w(N,2.074):5d} | {len(S):6d}  {w(S,2.330):5d} {w(S,2.318):5d} {w(S,2.306):5d} {w(S,2.298):5d} {w(S,2.300):5d}")
