import gemmi,numpy as np,json
from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from ase.optimize import LBFGS
from tblite.ase import TBLite
b=gemmi.cif.read('CLA.cif').sole_block()
t=b.find('_chem_comp_atom.',['atom_id','type_symbol','pdbx_model_Cartn_x_ideal','pdbx_model_Cartn_y_ideal','pdbx_model_Cartn_z_ideal','model_Cartn_x','model_Cartn_y','model_Cartn_z'])
at={}
for r in t:
    xyz=[r[2],r[3],r[4]] if '?' not in (r[2],r[3],r[4]) else [r[5],r[6],r[7]]
    at[r[0].strip('"')]=(r[1].capitalize(),np.array([float(v) for v in xyz]))
bonds=[(r[0].strip('"'),r[1].strip('"')) for r in b.find('_chem_comp_bond.',['atom_id_1','atom_id_2'])]
nb={k:set() for k in at}
for a,c in bonds: nb[a].add(c); nb[c].add(a)
kill=set(['C%d'%i for i in range(1,21)]+['CAA','CBA','CGA','O1A','O2A','CGD','O1D','O2D','CED','CBC'])
killH={h for k in kill for h in nb[k] if at[h][0]=='H'}
caps=[('C2A','CAA'),('CBD','CGD'),('CAC','CBC')]
names=[];sym=[];xyz=[]
for k,(e,p) in at.items():
    if k in kill or k in killH: continue
    names.append(k);sym.append(e);xyz.append(p)
for keep,gone in caps:
    v=at[gone][1]-at[keep][1]; names.append('Hcap_'+keep);sym.append('H');xyz.append(at[keep][1]+1.09*v/np.linalg.norm(v))
xyz=np.array(xyz); ix={n:i for i,n in enumerate(names)}
N4=xyz[[ix[k] for k in('NA','NB','NC','ND')]]; cen=N4.mean(0); n=np.linalg.svd(N4-cen)[2][2]
# imidazole
im=Chem.AddHs(Chem.MolFromSmiles('c1c[nH]cn1')); AllChem.EmbedMolecule(im,randomSeed=7); AllChem.MMFFOptimizeMolecule(im)
P=im.GetConformer().GetPositions(); el=[a.GetSymbol() for a in im.GetAtoms()]
n3=[a.GetIdx() for a in im.GetAtoms() if a.GetSymbol()=='N' and all(x.GetSymbol()!='H' for x in a.GetNeighbors())][0]
ring=[a.GetIdx() for a in im.GetAtoms() if a.GetSymbol()!='H']; u=P[n3]-P[ring].mean(0); u/=np.linalg.norm(u)
def rot(a,bv):
    v=np.cross(a,bv); c=np.dot(a,bv); s=np.linalg.norm(v)
    if s<1e-8: return np.eye(3)
    K=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]]); return np.eye(3)+K+K@K*((1-c)/s**2)
Rm=rot(u,-n); P=(Rm@(P-P[n3]).T).T + xyz[ix['MG']]+n*2.2
natp=len(sym); sym+=el; xyz=np.vstack([xyz,P]); names+=['IM_'+e+str(i) for i,e in enumerate(el)]
print('atoms',len(sym),'heavy',sum(s!='H' for s in sym),'formula check electrons',sum(gemmi.Element(s).atomic_number for s in sym))
A=Atoms(sym,positions=xyz); A.calc=TBLite(method='GFN2-xTB',verbosity=0)
LBFGS(A,logfile='opt.log').run(fmax=0.05,steps=600)
X=A.get_positions(); N4=X[[ix[k] for k in('NA','NB','NC','ND')]]; cen=N4.mean(0); n=np.linalg.svd(N4-cen)[2][2]
iN=natp+n3; 
if np.dot(X[iN]-cen,n)<0: n=-n
print('relaxed: Mg-N(im) %.3f  Mg out of ring %.3f  ligand height %.3f  Mg-N(ring) %.3f'%(np.linalg.norm(X[iN]-X[ix['MG']]),np.dot(X[ix['MG']]-cen,n),np.dot(X[iN]-cen,n),np.mean(np.linalg.norm(N4-X[ix['MG']],axis=1))))
json.dump(dict(sym=sym,xyz=X.tolist(),mg=ix['MG'],n3=iN,im=list(range(natp,len(sym))),names=names),open('model.json','w'))
A.write('model_relaxed.xyz')
