import sys,io,re
import iotbx.pdb, mmtbx.model
fn,mode=sys.argv[1],sys.argv[2]
inp=iotbx.pdb.input(fn); log=io.StringIO()
model=mmtbx.model.manager(model_input=inp,log=log)
p=mmtbx.model.manager.get_default_pdb_interpretation_params()
if mode=='on': p.pdb_interpretation.automatic_linking.link_metals=True
model.process(pdb_interpretation_params=p,make_restraints=True)
mg=[a for a in model.get_hierarchy().atoms() if a.name.strip()=='MG'][0]
geo=model.restraints_as_geo(); open(fn.rsplit('.',1)[0]+'_'+fn.rsplit('.',1)[1]+f'_{mode}.geo','w').write(geo)
L=geo.split('\n'); hit=[i for i,l in enumerate(L) if l.startswith('bond') and i+1<len(L) and ('NE2 HIS' in l+L[i+1]) and ('MG   CLA' in l+L[i+1])]
print(f"input={fn}  link_metals={'True' if mode=='on' else 'Auto (default)'}  element read as {mg.element!r}")
if not hit: print('   no Mg-NE2 bond restraint')
for i in hit: print('\n'.join('   '+x for x in L[i:i+4]))
