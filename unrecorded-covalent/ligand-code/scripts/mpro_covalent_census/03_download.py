"""Download each Mpro entry's mmCIF (gz) from files.rcsb.org; log URL, UTC time, bytes, SHA-256."""
import os, time, hashlib, urllib.request
from concurrent.futures import ThreadPoolExecutor
D='/workspace/mpro_covalent_census'; CIF=f'{D}/data/cif'; LOG=f'{D}/out/download_log.tsv'
ids=open(f'{D}/out/mpro_entries.txt').read().split()
def get(e):
    p=f'{CIF}/{e}.cif.gz'; url=f'https://files.rcsb.org/download/{e}.cif.gz'
    for i in range(5):
        try:
            t=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
            with urllib.request.urlopen(url,timeout=120) as r: b=r.read()
            open(p,'wb').write(b); return f'{e}\t{url}\t{t}\t{len(b)}\t{hashlib.sha256(b).hexdigest()}'
        except Exception as ex: err=repr(ex); time.sleep(2*(i+1))
    return f'{e}\t{url}\tFAILED\t0\t{err}'
with open(LOG,'w') as f, ThreadPoolExecutor(16) as ex:
    f.write('entry\turl\tretrieved_utc\tbytes\tsha256\n')
    for l in ex.map(get,ids): f.write(l+'\n')
