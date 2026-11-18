"""Download header-only mmCIF (no atom_site; includes _struct_conn, _pdbx_validate_close_contact, status) for every
experimental PDB entry with >=1 non-polymer entity (RCSB search rcsb_entry_info.nonpolymer_entity_count>0).
Stored gzip-compressed as served (Accept-Encoding: gzip)."""
import urllib.request, os, sys, time, concurrent.futures as cf, gzip
O='/workspace/ligand_code_effect/archive'
ids=[l.strip() for l in open(f'{O}/entry_ids_nonpolymer.txt') if l.strip()]
if len(sys.argv)>1: ids=ids[:int(sys.argv[1])]
def get(e):
    p=f'{O}/headers/{e}.cif.gz'
    if os.path.exists(p) and os.path.getsize(p)>200: return e,'cached'
    for t in range(4):
        try:
            r=urllib.request.urlopen(urllib.request.Request(f'https://files.rcsb.org/header/{e}.cif',headers={'Accept-Encoding':'gzip'}),timeout=60)
            b=r.read()
            if r.headers.get('Content-Encoding')!='gzip': b=gzip.compress(b)
            open(p+'.tmp','wb').write(b); os.replace(p+'.tmp',p); return e,'ok'
        except Exception as ex: err=str(ex)[:80]; time.sleep(2*(t+1))
    return e,'FAIL '+err
t0=time.time(); n=0; fails=[]
with cf.ThreadPoolExecutor(48) as ex:
    for e,s in ex.map(get,ids):
        n+=1
        if s.startswith('FAIL'): fails.append((e,s))
        if n%5000==0: print(n,round(time.time()-t0),len(fails),flush=True)
open(f'{O}/header_fetch_fails.tsv','w').write(''.join(f'{e}\t{s}\n' for e,s in fails))
print('done',n,len(fails),round(time.time()-t0))
