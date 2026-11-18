import sys,io
import iotbx.pdb, mmtbx.model
fn,mode=sys.argv[1],sys.argv[2]
inp=iotbx.pdb.input(fn); log=io.StringIO()
model=mmtbx.model.manager(model_input=inp,log=log)
p=mmtbx.model.manager.get_default_pdb_interpretation_params()
if mode=='on': p.pdb_interpretation.automatic_linking.link_metals=True
if mode=='nomcl': p.pdb_interpretation.restraints_library.mcl=False
model.process(pdb_interpretation_params=p,make_restraints=True)
geo=model.restraints_as_geo(); open(fn.rsplit('.',1)[0]+f'_{mode}.geo','w').write(geo)
L=geo.split('\n'); hit=[i for i,l in enumerate(L) if l.startswith('bond') and i+1<len(L) and 'ZN ' in (l+L[i+1])]
print(f"input={fn}  settings={ {'auto':'all defaults','on':'link_metals=True','nomcl':'mcl=False'}[mode] }")
if not hit: print('   no Zn bond restraint')
for i in hit: print('\n'.join('   '+x for x in L[i:i+4]))
