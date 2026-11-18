# Run on Google Colab (CPU). Part 1 builds the inputs; part 2 installs each cctbx-base build and prints the restraints.
# The two fragments are the files frag_MG.pdb and frag_ZN_his3.pdb in this folder (uploaded to /content).
import subprocess, os, sys, time
os.chdir('/content')
SETUP = '''
set -e
cd /content
[ -x bin/micromamba ] || (curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba)
[ -d mon_lib ] || git clone -q https://github.com/rlabduke/mon_lib
(cd mon_lib && git checkout -q c382734a31f6d673cf20683712f5315133406303)
[ -d geostd ] || git clone -q --filter=blob:none --sparse https://github.com/phenix-project/geostd
(cd geostd && git sparse-checkout set a c h l t list rna_dna && git checkout -q f4d21cc5f4dbf84f55cc972fb6b287354300ab2d)
rm -rf merged; mkdir -p merged/a merged/c merged/h merged/l merged/t merged/z merged/list
cp geostd/a/data_ALA.cif geostd/a/data_ASN.cif merged/a/
cp geostd/c/data_CLA.cif merged/c/; cp geostd/h/data_HIS.cif merged/h/; cp geostd/l/data_LEU.cif merged/l/; cp geostd/t/data_TRP.cif merged/t/
cp mon_lib/z/ZN.cif merged/z/; cp mon_lib/ener_lib.cif merged/; cp geostd/ener_lib.cif merged/geostd_ener_lib.cif
cp mon_lib/list/mon_lib_list.cif merged/list/; cp geostd/list/geostd_list.cif merged/list/
'''
print(subprocess.run(['bash', '-c', SETUP], capture_output=True, text=True).stderr[-2000:])
TEST = r'''
import sys, io
import iotbx.pdb, mmtbx.model
fn, metal = sys.argv[1], sys.argv[2]
log = io.StringIO()
how = ''
try:
    model = mmtbx.model.manager(model_input=iotbx.pdb.input(fn), log=log)
    model.process(make_restraints=True)
    how = 'process()'
except Exception as e1:
    try:
        model = mmtbx.model.manager(model_input=iotbx.pdb.input(fn), log=log, build_grm=True)
        how = 'build_grm'
    except Exception as e2:
        model = mmtbx.model.manager(model_input=iotbx.pdb.input(fn), log=log)
        model.process_input_model(make_restraints=True)
        how = 'process_input_model'
try:
    geo = model.restraints_as_geo()
except Exception as e3:
    out = io.StringIO()
    model.get_restraints_manager().geometry.show_sorted(sites_cart=model.get_sites_cart(), site_labels=[a.id_str() for a in model.get_hierarchy().atoms()], f=out)
    geo = out.getvalue()
el = [a.element for a in model.get_hierarchy().atoms() if a.name.strip() == metal]
L = geo.split('\n')
hits = []
for i, l in enumerate(L):
    if l.lstrip().startswith('bond') and i + 4 < len(L):
        two = l + L[i+1]
        if metal == 'MG' and 'NE2 HIS' in two and 'MG   CLA' in two: hits.append(i)
        if metal == 'ZN' and 'HIS' in two and 'ZN    ZN' in two: hits.append(i)
print('RES how=%s element_after=%r n_restraints=%d' % (how, el, len(hits)))
for i in hits:
    print('RES   ' + L[i].strip()[:40] + ' | ' + L[i+3].strip())
'''
open('test.py', 'w').write(TEST)
VERS = [('2020.8','3.8'),('2020.9','3.8'),('2020.10','3.8'),('2020.11','3.8'),('2020.12','3.8'),('2021.1','3.8'),
        ('2021.2','3.9'),('2022.1','3.9'),('2023.10','3.10'),('2024.5','3.10'),('2024.6','3.10')]
env = dict(os.environ, MMTBX_CCP4_MONOMER_LIB='/content/merged')
for v, py in VERS:
    p = '/content/env_' + v
    t0 = time.time()
    if not os.path.exists(p + '/bin/python'):
        r = subprocess.run(['./bin/micromamba','create','-y','-q','-p',p,'-c','conda-forge','cctbx-base=' + v,'python=' + py], capture_output=True, text=True)
        if r.returncode != 0:
            print('RES VERSION', v, 'INSTALL FAILED', r.stderr[-400:].replace('\n', ' ')); continue
    print('RES VERSION cctbx-base', v, 'python', py, 'installed in %.0fs' % (time.time() - t0)); sys.stdout.flush()
    for fn, metal in (('frag_MG.pdb', 'MG'), ('frag_ZN_his3.pdb', 'ZN')):
        r = subprocess.run(['./bin/micromamba','run','-p',p,'python','test.py',fn,metal], capture_output=True, text=True, env=env)
        got = [l for l in r.stdout.split('\n') if l.startswith('RES')]
        print('RES  input', fn)
        print('\n'.join(got) if got else 'RES   FAILED: ' + r.stderr[-500:].replace('\n', ' | '))
        sys.stdout.flush()
print('ALL_DONE')
