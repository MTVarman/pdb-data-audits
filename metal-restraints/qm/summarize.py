import json,numpy as np
from pyscf import gto
from pyscf.dispersion import dftd3
XC={'b':'b3lyp','c':'cam-b3lyp'}; BS={'s':'def2-SVP','t':'def2-TZVP'}
def load(fn,pig):
    R=[]
    for l in open(fn).read().strip().split('\n'):
        p=l.split(); st=[x.split('/') for x in p[6:]]
        R.append(dict(pig=pig,xc=XC[p[0]],basis=BS[p[1]],tda=p[2]=='1',d=float(p[3]),e=float(p[4]),qy=float(st[0][0]),fy=float(st[0][1]),qx=float(st[1][0]),fx=float(st[1][1]),s3=float(st[2][0]),t1=st[0][2],t2=st[1][2]))
    return R
R=load('rows_chla.txt','Chl a')+load('rows_bchla.txt','BChl a')
# state identity check: S1 must be HOMO->LUMO dominated (00), S2 HOMO-1->LUMO (10)
assert all(r['t1'].startswith('00') and r['t2'].startswith('10') for r in R)
# D3(BJ) for both pigments
def d3curve(modelfn,xc,ds):
    m=json.load(open(modelfn)); X0=np.array(m['xyz']); sym=m['sym']; mg=m['mg']; n3=m['n3']; im=m['im']
    u=X0[n3]-X0[mg]; d0=np.linalg.norm(u); u/=d0; out={}
    for d in ds:
        X=X0.copy(); X[im]+=u*(d-d0)
        mol=gto.M(atom=[(a,x) for a,x in zip(sym,X)],basis='sto-3g',verbose=0)
        out[d]=float(dftd3.DFTD3Dispersion(mol,xc=xc,version='d3bj').get_dispersion()['energy'])
    return out
D3={}
for pig,fn in (('Chl a','chla_imidazole_model.json'),('BChl a','bchla_imidazole_model.json')):
    for xc in ('b3lyp','cam-b3lyp'):
        D3[(pig,xc)]=d3curve(fn,xc,[2.5,2.35,2.2,2.05,1.9,1.76,1.6])
H=627.509; CM=8065.544
out=[]
for r in R:
    ref=[q for q in R if q['pig']==r['pig'] and q['xc']==r['xc'] and q['basis']==r['basis'] and q['tda']==r['tda'] and q['d']==2.2][0]
    r['dQy_meV']=round((r['qy']-ref['qy'])*1000,2); r['dQy_cm']=round((r['qy']-ref['qy'])*CM,1)
    r['dQx_meV']=round((r['qx']-ref['qx'])*1000,2); r['dQx_cm']=round((r['qx']-ref['qx'])*CM,1)
    r['strain']=round((r['e']-ref['e'])*H,2)
    dd=D3[(r['pig'],r['xc'])]; r['strain_d3']=round((r['e']-ref['e']+dd[r['d']]-dd[2.2])*H,2)
json.dump(R,open('results_all.json','w'),indent=1)
def show(pig,xc,basis,tda):
    rows=sorted([r for r in R if r['pig']==pig and r['xc']==xc and r['basis']==basis and r['tda']==tda],key=lambda r:-r['d'])
    if not rows: return
    print(f"\n{pig} | {xc} | {basis} | {'TDA' if tda else 'full TDDFT'}")
    print(' d(A)   Qy(eV)  f     dQy(meV) dQy(cm-1)  Qx(eV)  dQx(meV) dQx(cm-1)  strain  strain+D3 (kcal/mol)')
    for r in rows: print(f" {r['d']:.2f}  {r['qy']:.5f} {r['fy']:.3f}  {r['dQy_meV']:+7.2f} {r['dQy_cm']:+8.1f}   {r['qx']:.5f} {r['dQx_meV']:+7.2f} {r['dQx_cm']:+8.1f}   {r['strain']:+7.2f} {r['strain_d3']:+7.2f}")
for pig in ('Chl a','BChl a'):
    for xc in ('b3lyp','cam-b3lyp'):
        for basis in ('def2-SVP','def2-TZVP'):
            for tda in (True,False): show(pig,xc,basis,tda)
print('\n2.20 -> 1.90 summary')
for r in R:
    if r['d']==1.9: print(f" {r['pig']:7s} {r['xc']:10s} {r['basis']:9s} {'TDA ' if r['tda'] else 'TDDFT'}  dQy {r['dQy_meV']:+6.2f} meV ({r['dQy_cm']:+6.1f} cm-1)  dQx {r['dQx_meV']:+7.2f} meV ({r['dQx_cm']:+7.1f} cm-1)  strain {r['strain']:5.2f} / D3 {r['strain_d3']:5.2f}")
