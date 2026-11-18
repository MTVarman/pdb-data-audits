"""Fetch CCD files for every comp that appears in a nucleophile covale row or a nucleophile close contact anywhere
in archive/parsed.jsonl (adds to data/ccd)."""
import json, os, urllib.request, concurrent.futures as cf, time
O='/workspace/ligand_code_effect'
comps=set()
import sys
for l in open(sys.argv[1] if len(sys.argv)>1 else f'{O}/archive/parsed.jsonl'):
    r=json.loads(l)
    for x in r.get('cov',[])+r.get('cc',[]): comps.add(x[5])
need=sorted(c for c in comps if not os.path.exists(f'{O}/data/ccd/{c}.cif'))
def get(c):
    for u in (f'https://files.rcsb.org/ligands/download/{c}.cif',f'https://files.wwpdb.org/pub/pdb/refdata/chem_comp/{c[-1]}/{c}/{c}.cif'):
        for t in range(3):
            try: open(f'{O}/data/ccd/{c}.cif','wb').write(urllib.request.urlopen(u,timeout=60).read()); return c,'ok'
            except Exception as e: err=str(e)[:60]; time.sleep(1)
    return c,'FAIL '+err
with cf.ThreadPoolExecutor(16) as ex: res=list(ex.map(get,need))
with open(f'{O}/out/ccd_download_log.tsv','a') as f: f.write(f'# archive step {time.strftime("%Y-%m-%d %H:%M %Z")}\n'+''.join(f'{c}\t{s}\n' for c,s in res))
print(len(comps),len(need),[r for r in res if r[1].startswith('FAIL')][:20])
