import sys,io
import iotbx.pdb, mmtbx.model
from mmtbx.monomer_library import linking_utils
for fn in sys.argv[1:]:
    model=mmtbx.model.manager(model_input=iotbx.pdb.input(fn),log=io.StringIO())
    model.process(pdb_interpretation_params=mmtbx.model.manager.get_default_pdb_interpretation_params(),make_restraints=True)
    for a in model.get_hierarchy().atoms():
        if a.name.strip() in('MG','ZN') or (a.name.strip()=='NE2' and a.parent().resname=='HIS' and a.parent().parent().resseq.strip() in('52','116')):
            c=linking_utils.get_classes(a)
            print(fn, a.quote(), 'important_only=%r'%c.important_only, {k:getattr(c,k) for k in dir(c) if not k.startswith('_') and isinstance(getattr(c,k),bool) and getattr(c,k)})
