# Rigid Mg-N(imidazole) scan, TDA / TDDFT excitation energies on GPU (gpu4pyscf 1.9.1, cupy-cuda12x, pyscf 2.15.0).
# Executed version (with outputs) is the Colab notebook chl_mg_ligand_squish_tddft_2026-11-11.ipynb.
# Setup on Colab: pip uninstall -y cupy-cuda13x; pip install pyscf gpu4pyscf-cuda12x
import json,time,sys,numpy as np
from pyscf import gto
from gpu4pyscf.dft import rks
m=json.load(open(sys.argv[1]))   # chla_imidazole_model.json or bchla_imidazole_model.json
X0=np.array(m['xyz']); sym=m['sym']; mg=m['mg']; n3=m['n3']; im=m['im']
u=X0[n3]-X0[mg]; d0=float(np.linalg.norm(u)); u=u/d0
def geom(d):
    X=X0.copy(); X[im]+=u*(d-d0); return [(a,tuple(float(v) for v in x)) for a,x in zip(sym,X)]
def tonp(a): return a.get() if hasattr(a,'get') else np.asarray(a)
DM={}
def run(d,xc,basis='def2-svp',nst=3,tda=True):
    t0=time.time(); mol=gto.M(atom=geom(d),basis=basis,verbose=0)
    mf=rks.RKS(mol,xc=xc).density_fit(); mf.grids.level=3; mf.conv_tol=1e-9
    e=float(mf.kernel(dm0=DM.get((xc,basis)))); DM[(xc,basis)]=mf.make_rdm1(); t1=time.time()
    td=mf.TDA() if tda else mf.TDDFT(); td.nstates=nst; td.kernel()
    es=[float(x)*27.211386 for x in tonp(td.e)]; f=[float(x) for x in tonp(td.oscillator_strength())]
    nocc=mol.nelectron//2; st=[]
    for k in range(nst):
        w=2*tonp(td.xy[k][0])**2; idx=np.dstack(np.unravel_index(np.argsort(-w.ravel())[:2],w.shape))[0]
        st.append([round(es[k],5),round(f[k],4)]+[f'H-{nocc-1-i}>L+{a}:{w[i,a]:.2f}' for i,a in idx])
    return dict(d=d,xc=xc,basis=basis,tda=tda,nbf=int(mol.nao),conv=bool(mf.converged),e=e,states=st,t=[round(t1-t0),round(time.time()-t1)])
if __name__=='__main__':
    for xc in ('b3lyp','camb3lyp'):
        for d in (2.20,1.90,2.35,2.05,1.76,2.50,1.60):
            print(json.dumps(run(d,xc)),flush=True)
