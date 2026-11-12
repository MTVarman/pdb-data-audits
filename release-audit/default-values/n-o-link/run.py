import sys, io
import iotbx.pdb, mmtbx.model
for fn in sys.argv[1:]:
    inp = iotbx.pdb.input(fn); log = io.StringIO()
    model = mmtbx.model.manager(model_input=inp, log=log)
    p = mmtbx.model.manager.get_default_pdb_interpretation_params()
    model.process(pdb_interpretation_params=p, make_restraints=True)
    L = model.restraints_as_geo().split('\n')
    hit = [i for i, l in enumerate(L) if l.startswith('bond') and i + 1 < len(L) and 'NZ  LYS' in l + L[i+1] and 'ATP' in l + L[i+1]]
    print(f'input={fn}')
    if not hit: print('   no Lys NZ - ATP bond restraint')
    for i in hit: print('\n'.join('   ' + x for x in L[i:i+4]))
    for l in log.getvalue().split('\n'):
        if 'link' in l.lower() and ('ATP' in l or 'LYS' in l or 'Number' in l): print('   log:', l.strip()[:120])
