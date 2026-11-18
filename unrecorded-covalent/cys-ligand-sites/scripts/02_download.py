"""Download every candidate mmCIF (gz) from files.rcsb.org once; log URL, UTC time, bytes, SHA-256."""
import os, sys, time, hashlib, urllib.request
from concurrent.futures import ThreadPoolExecutor
D = os.path.join(os.path.dirname(__file__), '..'); CIF = f'{D}/cif'; LOG = f'{D}/out/download_log.tsv'
ids = [l.split()[0] for l in open(f'{D}/out/candidates.tsv').read().split('\n')[1:] if l.strip()]
done = set()
if os.path.exists(LOG):
    done = {l.split('\t')[0] for l in open(LOG) if l.strip() and not l.startswith('entry')}
def get(e):
    p = f'{CIF}/{e}.cif.gz'
    if e in done and os.path.exists(p): return None
    url = f'https://files.rcsb.org/download/{e}.cif.gz'
    for i in range(5):
        try:
            t = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
            with urllib.request.urlopen(url, timeout=120) as r: b = r.read()
            open(p + '.tmp', 'wb').write(b); os.replace(p + '.tmp', p)
            return f'{e}\t{url}\t{t}\t{len(b)}\t{hashlib.sha256(b).hexdigest()}'
        except Exception as ex:
            err = repr(ex); time.sleep(2 * (i + 1))
    return f'{e}\t{url}\tFAILED\t0\t{err}'
new = not os.path.exists(LOG)
with open(LOG, 'a') as f, ThreadPoolExecutor(int(sys.argv[1]) if len(sys.argv)>1 else 16) as ex:
    if new: f.write('entry\turl\tretrieved_utc\tbytes\tsha256\n')
    n = 0
    for line in ex.map(get, ids):
        if line: f.write(line + '\n'); f.flush()
        n += 1
        if n % 1000 == 0: print(n, time.strftime('%H:%M:%S'), flush=True)
print('done', n)
