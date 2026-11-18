"""Parse each CCD definition: type, pdbx_type, formula, atom counts, leaving atoms, dates, creating site,
PCM annotation (linking atom), InChIKey, and an RDKit molecule built from the CCD atom/bond tables.
Keys for grouping the same drug across codes:
  inchikey14  - InChIKey connectivity block (identical connectivity incl. H count)
  skeleton    - canonical SMILES of heavy atoms only, all bonds set to single, no charge/stereo/aromaticity,
                leaving atoms removed. Same skeleton = same heavy-atom graph, so a pre-reaction form
                (ketone, nitrile) and a post-reaction form (hemithioketal, thioimidate) of one drug collapse."""
import gemmi, glob, json, os
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
O='/workspace/ligand_code_effect'
def v(x):
    if x is None: return None
    s=gemmi.cif.as_string(x); return None if s in ('?','.','') else s
BO={'SING':1,'DOUB':2,'TRIP':3,'AROM':1.5}
out={}
for p in sorted(glob.glob(f'{O}/data/ccd/*.cif')):
    b=gemmi.cif.read(p).sole_block(); cid=b.name
    g=lambda t:v(b.find_value(t)) if b.find_value(t) is not None else None
    r={k:g('_chem_comp.'+k) for k in ('name','type','pdbx_type','formula','pdbx_initial_date','pdbx_modified_date',
       'pdbx_processing_site','pdbx_pcm','pdbx_model_coordinates_db_code','pdbx_release_status','pdbx_replaced_by','pdbx_replaces','pdbx_ambiguous_flag')}
    at=b.find('_chem_comp_atom.',['atom_id','type_symbol','charge','pdbx_leaving_atom_flag','pdbx_aromatic_flag'])
    atoms=[(v(a[0]),(v(a[1]) or 'X').capitalize(),int(v(a[2]) or 0) if (v(a[2]) or '0').lstrip('-').isdigit() else 0,v(a[3]) or 'N',v(a[4]) or 'N') for a in at]
    bt=b.find('_chem_comp_bond.',['atom_id_1','atom_id_2','value_order','pdbx_aromatic_flag'])
    bonds=[(v(x[0]),v(x[1]),(v(x[2]) or 'SING').upper(),v(x[3]) or 'N') for x in bt]
    r['n_atoms']=len(atoms); r['n_heavy']=sum(a[1]!='H' for a in atoms)
    r['leaving_atoms']=[a[0] for a in atoms if a[3]=='Y']; r['leaving_heavy']=[a[0] for a in atoms if a[3]=='Y' and a[1]!='H']
    pc=b.find('_pdbx_chem_comp_pcm.',['modified_residue_id','comp_id_linking_atom','modified_residue_id_linking_atom','category'])
    r['pcm_rows']=[[v(x[i]) for i in range(4)] for x in pc]
    au=b.find('_pdbx_chem_comp_audit.',['action_type','date','processing_site'])
    r['audit']=[[v(x[0]),v(x[1]),v(x[2])] for x in au]
    cr=[a for a in r['audit'] if a[0]=='Create component']; r['create_site']=cr[0][2] if cr else r['pdbx_processing_site']
    rel=[a[1] for a in r['audit'] if a[0]=='Initial release']; r['ccd_release_date']=rel[0] if rel else None
    r['pcm_modify_dates']=[a[1] for a in r['audit'] if a[0]=='Modify PCM']
    de=b.find('_pdbx_chem_comp_descriptor.',['type','program','descriptor'])
    r['inchikey']=next((v(x[2]) for x in de if v(x[0])=='InChIKey'),None)
    r['smiles']=next((v(x[2]) for x in de if v(x[0])=='SMILES_CANONICAL' and 'CACTVS' in (v(x[1]) or '')),None) or next((v(x[2]) for x in de if v(x[0]) and v(x[0]).startswith('SMILES')),None)
    # neighbour table (heavy + H) for warhead classification
    nb={a[0]:[] for a in atoms}; el={a[0]:a[1] for a in atoms}; arom={a[0]:a[4]=='Y' for a in atoms}
    for a1,a2,o,ar in bonds:
        if a1 in nb and a2 in nb: nb[a1].append((a2,o,ar)); nb[a2].append((a1,o,ar))
    r['nb']={k:[(n,o) for n,o,_ in x] for k,x in nb.items()}; r['el']=el; r['arom']={k:x for k,x in arom.items() if x}
    # skeleton key
    try:
        m=Chem.RWMol(); idx={}
        lv=set(r['leaving_heavy'])
        for a in atoms:
            if a[1]=='H' or a[0] in lv: continue
            at_=Chem.Atom(a[1] if a[1]!='X' else 'C'); at_.SetNoImplicit(True); idx[a[0]]=m.AddAtom(at_)
        for a1,a2,o,ar in bonds:
            if a1 in idx and a2 in idx: m.AddBond(idx[a1],idx[a2],Chem.BondType.SINGLE)
        mm=m.GetMol(); mm.UpdatePropertyCache(strict=False)
        frags=Chem.GetMolFrags(mm,asMols=True,sanitizeFrags=False)
        big=max(frags,key=lambda x:x.GetNumAtoms())
        r['skeleton']=Chem.MolToSmiles(big,canonical=True); r['skeleton_n']=big.GetNumAtoms()
    except Exception as e: r['skeleton']=None; r['skeleton_err']=str(e)
    out[cid]=r
import sys
json.dump(out,open(sys.argv[1] if len(sys.argv)>1 else f'{O}/out/ccd_parsed.json','w'))
print(len(out), sum(1 for r in out.values() if r.get('skeleton')), sum(1 for r in out.values() if r['pdbx_pcm']=='Y'), sum(1 for r in out.values() if r['leaving_heavy']))
