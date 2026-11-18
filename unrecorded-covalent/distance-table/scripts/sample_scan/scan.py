"""Past-the-line scan for non-cysteine protein nucleophiles.

Usage: python3 -I scan.py ids.txt out_dir [threads]

Per entry, writes one JSON line into out_dir/shard_NN.jsonl. Resumable.
Kept per entry: method, resolution, title, keywords, deposition date, every
_struct_conn row, the PDB's own close-contact table, and every contact between a
protein side-chain nucleophile atom and a ligand C/N/O/S/P/B/Se atom within
2.6 A (model 1, asymmetric unit, image 0).

Nucleophile list: SER OG, THR OG1, TYR OH, LYS NZ, HIS ND1/NE2, ASP OD1/OD2,
GLU OE1/OE2, SEC SE, plus CYS SG as an internal control.
"""
import gzip, hashlib, json, os, sys, time, threading, urllib.request, zlib, tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
import gemmi

NUC = {('CYS', 'SG'), ('SER', 'OG'), ('THR', 'OG1'), ('TYR', 'OH'), ('LYS', 'NZ'),
       ('HIS', 'ND1'), ('HIS', 'NE2'), ('ASP', 'OD1'), ('ASP', 'OD2'),
       ('GLU', 'OE1'), ('GLU', 'OE2'), ('SEC', 'SE')}
LIG_EL = {'C', 'N', 'O', 'S', 'P', 'B', 'Se'}
RADIUS = 2.6
WATER = {'HOH', 'DOD', 'WAT', 'H2O'}
NSHARD = 16
HUGE_BYTES = 20_000_000
TMP = tempfile.mkdtemp(prefix='scan_')
HEAVY = threading.Semaphore(3)


def clean(v):
    if v in ('?', '.'):
        return None
    if len(v) >= 2 and v[0] == v[-1] and v[0] in '\'"':
        v = v[1:-1]
    if v.startswith(';'):
        v = v.strip(';').strip()
    return v


def table(block, prefix, tags):
    t = block.find_mmcif_category(prefix)
    if not len(t):
        return []
    have = [x[len(prefix):] for x in t.tags]
    idx = [have.index(x) if x in have else None for x in tags]
    return [[clean(r[i]) if i is not None else None for i in idx] for r in t]


def fetch(pdbid, tries=4):
    url = f'https://files.rcsb.org/download/{pdbid.upper()}.cif.gz'
    last = None
    for i in range(tries):
        try:
            t = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
            with urllib.request.urlopen(url, timeout=180) as r:
                return url, t, r.read()
        except Exception as e:  # noqa
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f'download failed: {last}')


def key(chain, res, atom):
    return (chain, res.seqid.num, res.seqid.icode, res.name, atom)


def process(pdbid):
    url, t, raw = fetch(pdbid)
    if len(raw) > HUGE_BYTES:
        return {'pdbid': pdbid, 'status': 'too_big', 'bytes': len(raw)}
    if len(raw) > 3_000_000:
        with HEAVY:
            return _process(pdbid, url, t, raw)
    return _process(pdbid, url, t, raw)


def _process(pdbid, url, t, raw):
    rec = {'pdbid': pdbid, 'url': url, 'retrieved_utc': t, 'bytes': len(raw),
           'sha256': hashlib.sha256(raw).hexdigest()}
    tmp = os.path.join(TMP, f'{pdbid}.cif.gz')
    with open(tmp, 'wb') as f:
        f.write(raw)
    del raw
    try:
        doc = gemmi.cif.read(tmp)
    finally:
        os.remove(tmp)
    block = doc.sole_block()
    st = gemmi.make_structure_from_block(block)
    st.setup_entities()
    methods = [clean(x) for x in block.find_values('_exptl.method')]
    rec['method'] = methods[0] if methods else ''
    rec['resolution'] = round(st.resolution, 3) if st.resolution else None
    rec['r_free'] = clean(block.find_value('_refine.ls_R_factor_R_free') or '?')
    rec['title'] = clean(block.find_value('_struct.title') or '?')
    rec['keywords'] = clean(block.find_value('_struct_keywords.pdbx_keywords') or '?')
    rec['deposited'] = clean(block.find_value('_pdbx_database_status.recvd_initial_deposition_date') or '?')
    rec['close'] = table(block, '_pdbx_validate_close_contact.',
                         ['PDB_model_num', 'auth_asym_id_1', 'auth_comp_id_1', 'auth_seq_id_1', 'auth_atom_id_1',
                          'auth_asym_id_2', 'auth_comp_id_2', 'auth_seq_id_2', 'auth_atom_id_2', 'dist'])
    conn = []
    for c in st.connections:
        a, b = c.partner1, c.partner2
        conn.append([c.type.name, a.chain_name, a.res_id.name, a.res_id.seqid.num, a.atom_name,
                     b.chain_name, b.res_id.name, b.res_id.seqid.num, b.atom_name,
                     c.asu.name, round(c.reported_distance, 3) if c.reported_distance else None])
    rec['conn'] = conn
    del block, doc
    if len(st) == 0:
        rec.update(status='no_model', contacts=[])
        return rec
    model = st[0]
    try:
        ns = gemmi.NeighborSearch(model, st.cell, 5).populate()
    except RuntimeError:
        ns = gemmi.NeighborSearch(model, gemmi.UnitCell(), 5).populate()
        rec['ns_no_symmetry'] = True
    pair, respair = {}, {}
    for c in st.connections:
        a, b = c.partner1, c.partner2
        ka = (a.chain_name, a.res_id.seqid.num, a.res_id.seqid.icode, a.res_id.name, a.atom_name)
        kb = (b.chain_name, b.res_id.seqid.num, b.res_id.seqid.icode, b.res_id.name, b.atom_name)
        for x, y in ((ka, kb), (kb, ka)):
            pair.setdefault((x, y), set()).add(c.type.name)
            respair.setdefault((x[:4], y[:4]), set()).add(f'{c.type.name}:{x[4]}-{y[4]}')
    contacts = []
    for ch in model:
        for res in ch:
            et = res.entity_type
            if res.name in WATER or et not in (gemmi.EntityType.NonPolymer, gemmi.EntityType.Branched):
                continue
            if not any(a.element.name == 'C' for a in res):
                continue  # organic ligands only; bare metals and salts excluded
            nheavy = sum(1 for a in res if a.element.name not in ('H', 'D'))
            for la in res:
                if la.element.name not in LIG_EL:
                    continue
                for mark in ns.find_atoms(la.pos, '\0', radius=RADIUS):
                    if mark.image_idx != 0:
                        continue
                    cra = mark.to_cra(model)
                    pres, pa = cra.residue, cra.atom
                    if pres.entity_type != gemmi.EntityType.Polymer or (pres.name, pa.name) not in NUC:
                        continue
                    d = la.pos.dist(pa.pos)
                    if d >= RADIUS:
                        continue
                    kp = key(cra.chain.name, pres, pa.name)
                    kl = key(ch.name, res, la.name)
                    contacts.append({
                        'prot': f'{cra.chain.name}:{pres.name}{pres.seqid.num}{pres.seqid.icode.strip()}:{pa.name}',
                        'prot_el': pa.element.name,
                        'prot_altloc': pa.altloc.strip('\0'), 'prot_occ': round(pa.occ, 2),
                        'prot_b': round(pa.b_iso, 1),
                        'lig': f'{res.name}_{ch.name}_{res.seqid.num}{res.seqid.icode.strip()}',
                        'lig_entity': 'branched' if et == gemmi.EntityType.Branched else 'nonpolymer',
                        'lig_nheavy': nheavy,
                        'lig_atom': la.name, 'lig_el': la.element.name,
                        'lig_altloc': la.altloc.strip('\0'), 'lig_occ': round(la.occ, 2),
                        'lig_b': round(la.b_iso, 1),
                        'dist': round(d, 3),
                        'conn_pair': sorted(pair.get((kl, kp), ())),
                        'conn_residues': sorted(respair.get((kl[:4], kp[:4]), ()))})
    rec.update(status='ok', contacts=contacts)
    return rec


def shard(pdbid):
    return zlib.crc32(pdbid.encode()) % NSHARD


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    ids = [l.strip().upper() for l in open(args[0]) if l.strip()]
    od = args[1]
    os.makedirs(od, exist_ok=True)
    threads = int(args[2]) if len(args) > 2 else 24
    done = set()
    for i in range(NSHARD):
        p = os.path.join(od, f'shard_{i:02d}.jsonl')
        if os.path.exists(p):
            good = []
            with open(p) as f:
                for line in f:
                    try:
                        done.add(json.loads(line)['pdbid'])
                        good.append(line)
                    except Exception:
                        pass
            with open(p, 'w') as f:
                f.writelines(good)
    todo = [e for e in ids if e not in done]
    print(f'{len(ids)} ids, {len(done)} done, {len(todo)} to do', flush=True)
    if not todo:
        print('ALL DONE', flush=True)
        return
    fhs = [open(os.path.join(od, f'shard_{i:02d}.jsonl'), 'a') for i in range(NSHARD)]
    lock = threading.Lock()
    n = 0
    t0 = time.time()
    with ThreadPoolExecutor(threads) as ex:
        it = iter(todo)
        futs = {}

        def submit_more():
            while len(futs) < threads * 3:
                try:
                    e = next(it)
                except StopIteration:
                    return
                futs[ex.submit(process, e)] = e
        submit_more()
        while futs:
            for f in as_completed(list(futs)):
                e = futs.pop(f)
                try:
                    r = f.result()
                except Exception as err:
                    r = {'pdbid': e, 'status': 'error', 'error': str(err)[:300]}
                with lock:
                    fh = fhs[shard(e)]
                    fh.write(json.dumps(r, separators=(',', ':')) + '\n')
                    fh.flush()
                    n += 1
                    if n % 500 == 0:
                        el = time.time() - t0
                        print(f'{n}/{len(todo)} {el:.0f}s ({n/el:.2f}/s) eta {(len(todo)-n)/(n/el)/60:.1f}min', flush=True)
                submit_more()
                break
    for fh in fhs:
        fh.close()
    print('done', n, f'{time.time()-t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
