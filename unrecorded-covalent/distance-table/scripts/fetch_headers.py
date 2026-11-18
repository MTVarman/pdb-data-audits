"""Read _struct_conn and the two close-contact validation categories from the header file of each entry.
Usage: python3 -I fetch_headers.py ids.txt out.jsonl [workers]
One JSON line per entry: url, retrieval time (UTC), size and SHA-256 of the file as served, and the rows.
Files are parsed in memory and not kept.  Resumable: entries already in out.jsonl are skipped."""
import gzip, hashlib, json, re, sys, threading, time, urllib.request, queue, os
import gemmi
ids_path, out_path = sys.argv[1], sys.argv[2]
workers = int(sys.argv[3]) if len(sys.argv) > 3 else 16
ids = [l.strip().upper() for l in open(ids_path) if l.strip()]
done = set()
if os.path.exists(out_path):
    for line in open(out_path):
        try: done.add(json.loads(line)['entry'])
        except Exception: pass
todo = [e for e in ids if e not in done]
print(f'{len(ids)} entries, {len(done)} done, {len(todo)} to fetch', flush=True)
CONN = ['conn_type_id', 'ptnr1_auth_asym_id', 'ptnr1_auth_comp_id', 'ptnr1_auth_seq_id', 'pdbx_ptnr1_PDB_ins_code', 'ptnr1_label_atom_id', 'pdbx_ptnr1_label_alt_id', 'ptnr1_symmetry',
        'ptnr2_auth_asym_id', 'ptnr2_auth_comp_id', 'ptnr2_auth_seq_id', 'pdbx_ptnr2_PDB_ins_code', 'ptnr2_label_atom_id', 'pdbx_ptnr2_label_alt_id', 'ptnr2_symmetry', 'pdbx_dist_value']
CLOSE = ['PDB_model_num', 'auth_asym_id_1', 'auth_comp_id_1', 'auth_seq_id_1', 'PDB_ins_code_1', 'label_alt_id_1', 'auth_atom_id_1',
         'auth_asym_id_2', 'auth_comp_id_2', 'auth_seq_id_2', 'PDB_ins_code_2', 'label_alt_id_2', 'auth_atom_id_2', 'dist']
SYMM = CLOSE + ['site_symmetry_1', 'site_symmetry_2']
def table(block, cat, cols):
    t = block.find(cat, [cols[0]] + ['?' + c for c in cols[1:]])
    out = []
    for row in t:
        vals = []
        for i in range(len(cols)):
            v = row[i] if row.has(i) else ''
            v = gemmi.cif.as_string(v) if v not in ('?', '.') else ''
            vals.append(v)
        out.append(vals)
    return out
def one(e):
    url = f'https://files.rcsb.org/header/{e}.cif'
    req = urllib.request.Request(url, headers={'User-Agent': 'pdb-struct-audit/1.0', 'Accept-Encoding': 'gzip'})
    last = None; raw = None
    for a in range(5):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                raw = r.read()
                if r.headers.get('Content-Encoding', '') == 'gzip': raw = gzip.decompress(raw)
            break
        except Exception as ex:
            last = repr(ex)[:200]; raw = None; time.sleep(2 + 3 * a)
    if raw is None: return {'entry': e, 'url': url, 'status': 'failed', 'error': last}
    try:
        text = raw.decode('utf8', 'replace')
        # header files end with loop_ stubs where the coordinate loops were cut; drop the empty stubs
        text = re.sub(r'(?m)^loop_[ \t]*\r?\n(?=#[ \t]*\r?$)', '', text)
        block = gemmi.cif.read_string(text).sole_block()
        return {'entry': e, 'url': url, 'retrieved_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'bytes': len(raw),
                'sha256': hashlib.sha256(raw).hexdigest(), 'status': 'ok',
                'conn': table(block, '_struct_conn.', CONN), 'close': table(block, '_pdbx_validate_close_contact.', CLOSE),
                'symm': table(block, '_pdbx_validate_symm_contact.', SYMM)}
    except Exception as ex:
        return {'entry': e, 'url': url, 'status': 'parse_failed', 'bytes': len(raw), 'error': repr(ex)[:200]}
q = queue.Queue(); lock = threading.Lock()
for e in todo: q.put(e)
n = [0]; t0 = time.time()
out = open(out_path, 'a')
def work():
    while True:
        try: e = q.get_nowait()
        except queue.Empty: return
        rec = one(e)
        with lock:
            out.write(json.dumps(rec, separators=(',', ':')) + '\n'); out.flush(); n[0] += 1
            if n[0] % 500 == 0: print(f'  {n[0]}/{len(todo)}  {time.time() - t0:.0f}s', flush=True)
ts = [threading.Thread(target=work, daemon=True) for _ in range(workers)]
for t in ts: t.start()
for t in ts: t.join()
out.close()
print('done', n[0], f'{time.time() - t0:.0f}s', flush=True)
