"""Download CCD definition files (files.rcsb.org/ligands/download/XXX.cif) for every ligand comp in contacts_merged."""
import pandas as pd, os, urllib.request, concurrent.futures as cf, hashlib, time
O='/workspace/ligand_code_effect'
comps=sorted(set(pd.read_csv(f'{O}/out/contacts_merged.tsv',sep='\t',dtype=str).lig_comp))
def get(c):
    p=f'{O}/data/ccd/{c}.cif'
    if os.path.exists(p) and os.path.getsize(p)>100: return c,'cached'
    for u in (f'https://files.rcsb.org/ligands/download/{c}.cif',f'https://files.wwpdb.org/pub/pdb/refdata/chem_comp/{c[-1]}/{c}/{c}.cif'):
        for t in range(3):
            try:
                b=urllib.request.urlopen(u,timeout=60).read(); open(p,'wb').write(b); return c,'ok'
            except Exception as e: err=str(e); time.sleep(1)
    return c,'FAIL '+err
with cf.ThreadPoolExecutor(16) as ex: res=list(ex.map(get,comps))
open(f'{O}/out/ccd_download_log.tsv','w').write(f'# fetched {time.strftime("%Y-%m-%d %H:%M %Z")}\n'+''.join(f'{c}\t{s}\n' for c,s in res))
print(sum(s!='FAIL' and not s.startswith('FAIL') for c,s in res), [r for r in res if r[1].startswith('FAIL')])
