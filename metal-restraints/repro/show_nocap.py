# Same as show.py, but with the 2023-03-07 element-uppercasing step (commit 8575b15fd) switched off,
# to show what that one step changes. Everything else is cctbx-base 2025.11 at default settings.
import sys,io
import iotbx.pdb, mmtbx.model
import mmtbx.monomer_library.pdb_interpretation as pi
pi.build_all_chain_proxies._check_for_capitalised_element = lambda self: None
fn=sys.argv[1]
inp=iotbx.pdb.input(fn); log=io.StringIO()
model=mmtbx.model.manager(model_input=inp,log=log)
p=mmtbx.model.manager.get_default_pdb_interpretation_params()
model.process(pdb_interpretation_params=p,make_restraints=True)
mg=[a for a in model.get_hierarchy().atoms() if a.name.strip()=='MG'][0]
geo=model.restraints_as_geo(); open(fn.rsplit('.',1)[0]+'_'+fn.rsplit('.',1)[1]+'_nocap.geo','w').write(geo)
L=geo.split('\n'); hit=[i for i,l in enumerate(L) if l.startswith('bond') and i+1<len(L) and ('NE2 HIS' in l+L[i+1]) and ('MG   CLA' in l+L[i+1])]
print(f"input={fn}  all defaults, uppercasing step disabled  element read as {mg.element!r}")
if not hit: print('   no Mg-NE2 bond restraint')
for i in hit: print('\n'.join('   '+x for x in L[i:i+4]))
