import sys,io,traceback
import iotbx.pdb, mmtbx.model
fn=sys.argv[1]
inp=iotbx.pdb.input(fn)
model=mmtbx.model.manager(model_input=inp,log=io.StringIO())
mg=[a for a in model.get_hierarchy().atoms() if a.name.strip()=='MG'][0]
state={'el':mg.element,'done':False}
def tr(frame,event,arg):
    if state['done']: return None
    if event in('line',):
        if mg.element!=state['el']:
            print('element changed %r -> %r'%(state['el'],mg.element))
            # the change happened on the previous line of this frame or in a callee that just returned
            st=traceback.extract_stack(frame)
            for s in st[-6:]: print('   %s:%d %s | %s'%(s.filename.split('site-packages/')[-1],s.lineno,s.name,s.line))
            state['done']=True; sys.settrace(None); return None
    return tr
sys.settrace(tr)
model.process(pdb_interpretation_params=mmtbx.model.manager.get_default_pdb_interpretation_params(),make_restraints=True)
sys.settrace(None)
