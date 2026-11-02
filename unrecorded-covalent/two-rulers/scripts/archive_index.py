"""Whole-PDB index: download every entry once, keep what later questions need.

Usage: python3 -I archive_index.py ids.txt out_dir [threads]
Writes out_dir/shard_NN.jsonl (32 shards, one JSON line per entry). Resumable.

Per entry:
  provenance  url, retrieval time, bytes, SHA-256
  meta        method, resolution, R-free, title, keywords, deposition date,
              revision history (+ which categories each revision touched)
  entities    type, description, pdbx_mutation, sequence, chains, source organism/taxid
  refs        _struct_ref (UniProt etc. + the reference sequence segment),
              _struct_ref_seq (alignment ranges), _struct_ref_seq_dif (every declared
              difference: engineered mutation / conflict / expression tag ...)
  unobs       unobserved polymer residues as label_seq ranges per chain
  mod_res     _pdbx_struct_mod_residue
  conn        every _struct_conn row (type, both partners, symmetry, distance)
  close       the PDB's own close-contact table (_pdbx_validate_close_contact)
  ligands     inventory of non-polymer and branched residues (atoms, occupancy, altlocs, B)
  contacts    ligand C/B/P/S/Se within 2.5 A of a protein nucleophile (same as
              archive_sweep.py; branched entities included and flagged)
  pp_contacts protein side-chain nucleophile atom within 2.2 A of a polymer heavy atom
              of a residue >= 2 apart or on another chain (unrecorded cross-link scan)
"""
import gzip, hashlib, json, os, sys, time, threading, urllib.request, zlib
from concurrent.futures import ThreadPoolExecutor, as_completed
import gemmi

NUC = {('CYS', 'SG'), ('SER', 'OG'), ('THR', 'OG1'), ('TYR', 'OH'), ('LYS', 'NZ'),
       ('HIS', 'ND1'), ('HIS', 'NE2'), ('ASP', 'OD1'), ('ASP', 'OD2'),
       ('GLU', 'OE1'), ('GLU', 'OE2'), ('MET', 'SD'), ('SEC', 'SE'),
       ('ARG', 'NH1'), ('ARG', 'NH2'), ('ARG', 'NE'), ('ASN', 'ND2'), ('GLN', 'NE2'),
       ('TRP', 'NE1')}
LIG_EL = {'C', 'B', 'P', 'S', 'Se'}
RADIUS, PP_RADIUS = 2.5, 2.2
WATER = {'HOH', 'DOD', 'WAT', 'H2O'}
NSHARD = 32
BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}

CATS = {
    'entity': ('_entity.', ['id', 'type', 'pdbx_description', 'pdbx_mutation', 'pdbx_fragment', 'pdbx_ec', 'pdbx_number_of_molecules']),
    'entity_poly': ('_entity_poly.', ['entity_id', 'type', 'pdbx_seq_one_letter_code_can', 'pdbx_strand_id', 'nstd_monomer']),
    'src_gen': ('_entity_src_gen.', ['entity_id', 'pdbx_gene_src_scientific_name', 'pdbx_gene_src_ncbi_taxonomy_id', 'pdbx_gene_src_gene', 'pdbx_host_org_scientific_name']),
    'src_nat': ('_entity_src_nat.', ['entity_id', 'pdbx_organism_scientific', 'pdbx_ncbi_taxonomy_id']),
    'src_syn': ('_pdbx_entity_src_syn.', ['entity_id', 'organism_scientific', 'ncbi_taxonomy_id']),
    'struct_ref': ('_struct_ref.', ['id', 'db_name', 'db_code', 'pdbx_db_accession', 'pdbx_db_isoform', 'entity_id', 'pdbx_seq_one_letter_code', 'pdbx_align_begin']),
    'struct_ref_seq': ('_struct_ref_seq.', ['align_id', 'ref_id', 'pdbx_strand_id', 'seq_align_beg', 'seq_align_end', 'pdbx_db_accession', 'db_align_beg', 'db_align_end', 'pdbx_auth_seq_align_beg', 'pdbx_auth_seq_align_end']),
    'seq_dif': ('_struct_ref_seq_dif.', ['align_id', 'pdbx_pdb_strand_id', 'mon_id', 'seq_num', 'db_mon_id', 'pdbx_seq_db_seq_num', 'pdbx_seq_db_accession_code', 'details', 'pdbx_auth_seq_num']),
    'mod_res': ('_pdbx_struct_mod_residue.', ['auth_asym_id', 'auth_comp_id', 'auth_seq_id', 'parent_comp_id', 'details']),
    'close': ('_pdbx_validate_close_contact.', ['PDB_model_num', 'auth_asym_id_1', 'auth_comp_id_1', 'auth_seq_id_1', 'auth_atom_id_1', 'label_alt_id_1',
                                                 'auth_asym_id_2', 'auth_comp_id_2', 'auth_seq_id_2', 'auth_atom_id_2', 'label_alt_id_2', 'dist']),
    'revisions': ('_pdbx_audit_revision_history.', ['ordinal', 'major_revision', 'minor_revision', 'revision_date']),
    'revision_categories': ('_pdbx_audit_revision_category.', ['revision_ordinal', 'category']),
    'revision_groups': ('_pdbx_audit_revision_group.', ['revision_ordinal', 'group']),
}

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

def unobs_ranges(block):
    rows = table(block, '_pdbx_unobs_or_zero_occ_residues.', ['PDB_model_num', 'polymer_flag', 'occupancy_flag', 'auth_asym_id', 'label_seq_id'])
    per = {}
    for m, pf, of, ch, s in rows:
        if m != '1' or pf != 'Y' or of != '1':
            continue
        if s is None or not s.lstrip('-').isdigit():
            continue
        per.setdefault(ch, []).append(int(s))
    out = {}
    for ch, v in per.items():
        v = sorted(set(v)); rng = []
        for x in v:
            if rng and x == rng[-1][1] + 1:
                rng[-1][1] = x
            else:
                rng.append([x, x])
        out[ch] = rng
    return out

import tempfile
HEAVY = threading.Semaphore(2)   # gz 2-15 MB: at most two at a time
HUGE_BYTES = 15_000_000          # larger entries are deferred to a final single-file pass
BIG_MODE = '--big' in sys.argv
TMP = tempfile.mkdtemp(prefix='idx_')

def process(pdbid):
    url, t, raw = fetch(pdbid)
    if len(raw) > HUGE_BYTES and not BIG_MODE:
        return {'pdbid': pdbid, 'status': 'deferred_big', 'bytes': len(raw)}
    if len(raw) > 2_000_000:
        with HEAVY:
            return _process(pdbid, url, t, raw)
    return _process(pdbid, url, t, raw)

def _process(pdbid, url, t, raw):
    rec = {'pdbid': pdbid, 'url': url, 'retrieved_utc': t, 'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}
    tmp = os.path.join(TMP, f'{pdbid}.cif.gz')
    with open(tmp, 'wb') as f:
        f.write(raw)
    del raw
    try:
        doc = gemmi.cif.read(tmp)   # gemmi reads the gz directly, no Python-side copy
    finally:
        os.remove(tmp)
    block = doc.sole_block()
    st = gemmi.make_structure_from_block(block)
    st.setup_entities()
    methods = [clean(x) for x in block.find_values('_exptl.method')]
    rec['method'] = methods[0] if methods else ''
    rec['methods'] = methods
    rec['resolution'] = round(st.resolution, 3) if st.resolution else None
    rec['r_free'] = clean(block.find_value('_refine.ls_R_factor_R_free') or '?')
    rec['title'] = clean(block.find_value('_struct.title') or '?')
    rec['keywords'] = clean(block.find_value('_struct_keywords.pdbx_keywords') or '?')
    rec['deposited'] = clean(block.find_value('_pdbx_database_status.recvd_initial_deposition_date') or '?')
    for name, (prefix, tags) in CATS.items():
        rec[name] = table(block, prefix, tags)
    rec['unobs'] = unobs_ranges(block)
    conn = []
    for c in st.connections:
        a, b = c.partner1, c.partner2
        conn.append([c.type.name, a.chain_name, a.res_id.name, a.res_id.seqid.num, a.res_id.seqid.icode.strip(), a.atom_name, a.altloc.strip('\0'),
                     b.chain_name, b.res_id.name, b.res_id.seqid.num, b.res_id.seqid.icode.strip(), b.atom_name, b.altloc.strip('\0'),
                     c.asu.name, round(c.reported_distance, 3) if c.reported_distance else None])
    rec['conn'] = conn
    del block, doc  # tables are extracted; free the parsed text before the neighbour search
    if len(st) == 0:
        rec.update(status='no_model', contacts=[], pp_contacts=[], ligands=[]); return rec
    model = st[0]
    try:
        ns = gemmi.NeighborSearch(model, st.cell, 5).populate()
    except RuntimeError:
        # non-standard crystal frame (SCALEn): search without symmetry; only image 0 is used anyway
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
    contacts, ligands, pp = [], [], []
    for ch in model:
        for res in ch:
            et = res.entity_type
            if et == gemmi.EntityType.Polymer:
                # protein-protein unrecorded cross-link scan, from side-chain nucleophiles
                for pa in res:
                    if (res.name, pa.name) not in NUC:
                        continue
                    for mark in ns.find_atoms(pa.pos, '\0', radius=PP_RADIUS):
                        if mark.image_idx != 0:
                            continue
                        cra = mark.to_cra(model)
                        r2, a2 = cra.residue, cra.atom
                        if r2.entity_type != gemmi.EntityType.Polymer or a2.element.name in ('H', 'D'):
                            continue
                        if cra.chain.name == ch.name and abs(r2.seqid.num - res.seqid.num) < 2 and r2.seqid.num != res.seqid.num:
                            continue
                        if cra.chain.name == ch.name and r2.seqid.num == res.seqid.num and r2.seqid.icode == res.seqid.icode:
                            continue
                        k1, k2 = key(ch.name, res, pa.name), key(cra.chain.name, r2, a2.name)
                        if k1 > k2 and (r2.name, a2.name) in NUC:
                            continue  # count nucleophile-nucleophile pairs once
                        d = pa.pos.dist(a2.pos)
                        if d >= PP_RADIUS:
                            continue
                        pp.append([f'{ch.name}:{res.name}{res.seqid.num}{res.seqid.icode.strip()}:{pa.name}', pa.altloc.strip('\0'), round(pa.occ, 2),
                                   f'{cra.chain.name}:{r2.name}{r2.seqid.num}{r2.seqid.icode.strip()}:{a2.name}', a2.altloc.strip('\0'), round(a2.occ, 2),
                                   a2.name in BACKBONE, round(d, 3), sorted(pair.get((k1, k2), ()))])
                continue
            if res.name in WATER or et not in (gemmi.EntityType.NonPolymer, gemmi.EntityType.Branched):
                continue
            heavy = [a for a in res if a.element.name not in ('H', 'D')]
            if not heavy:
                continue
            ligands.append([res.name, ch.name, f'{res.seqid.num}{res.seqid.icode.strip()}', 'branched' if et == gemmi.EntityType.Branched else 'nonpolymer',
                            len(heavy), round(sum(a.occ for a in heavy) / len(heavy), 2), sorted({a.altloc for a in heavy if a.altloc != '\0'}),
                            round(sum(a.b_iso for a in heavy) / len(heavy), 1)])
            if not any(a.element.name == 'C' for a in res):
                continue
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
                    kp, kl = key(cra.chain.name, pres, pa.name), key(ch.name, res, la.name)
                    contacts.append({
                        'prot': f'{cra.chain.name}:{pres.name}{pres.seqid.num}{pres.seqid.icode.strip()}:{pa.name}',
                        'prot_altloc': pa.altloc.strip('\0'), 'prot_occ': round(pa.occ, 2),
                        'lig': f'{res.name}_{ch.name}_{res.seqid.num}{res.seqid.icode.strip()}',
                        'lig_entity': 'branched' if et == gemmi.EntityType.Branched else 'nonpolymer',
                        'lig_atom': la.name, 'lig_el': la.element.name,
                        'lig_altloc': la.altloc.strip('\0'), 'lig_occ': round(la.occ, 2), 'dist': round(d, 3),
                        'conn_pair': sorted(pair.get((kl, kp), ())),
                        'conn_residues': sorted(respair.get((kl[:4], kp[:4]), ()))})
    rec.update(status='ok', contacts=contacts, pp_contacts=pp, ligands=ligands)
    return rec

def shard(pdbid):
    return zlib.crc32(pdbid.encode()) % NSHARD

def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    maxn = next((int(a.split('=')[1]) for a in sys.argv if a.startswith('--max=')), None)
    ids = [l.strip() for l in open(args[0]) if l.strip()]
    od = args[1]; os.makedirs(od, exist_ok=True)
    threads = int(args[2]) if len(args) > 2 else 20
    dpath = os.path.join(od, 'deferred_big.txt')
    deferred = set(l.strip() for l in open(dpath)) if os.path.exists(dpath) else set()
    done = set()
    for i in range(NSHARD):
        p = os.path.join(od, f'shard_{i:02d}.jsonl')
        if os.path.exists(p):
            good = []
            with open(p) as f:
                for line in f:
                    try:
                        done.add(json.loads(line)['pdbid']); good.append(line)
                    except Exception:
                        pass  # truncated last line from a killed run
            with open(p, 'w') as f:
                f.writelines(good)
    if BIG_MODE:
        todo = [e for e in ids if e in deferred and e not in done]
    else:
        todo = [e for e in ids if e not in done and e not in deferred]
    if not todo:
        print('ALL DONE', len(done), 'done,', len(deferred - done), 'deferred big left', flush=True)
        return
    if maxn:
        todo = todo[:maxn]
    print(f'{len(ids)} ids, {len(done)} done, {len(deferred)} deferred, {len(todo)} to do this run', flush=True)
    dfh = open(dpath, 'a')
    fhs = [open(os.path.join(od, f'shard_{i:02d}.jsonl'), 'a') for i in range(NSHARD)]
    lock = threading.Lock(); n = 0; t0 = time.time()
    with ThreadPoolExecutor(threads) as ex:
        it = iter(todo); futs = {}
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
                    if r.get('status') == 'deferred_big':
                        dfh.write(e + '\n'); dfh.flush(); n += 1
                        submit_more(); break
                    fh = fhs[shard(e)]
                    fh.write(json.dumps(r, separators=(',', ':')) + '\n'); fh.flush()
                    n += 1
                    if n % 2000 == 0:
                        el = time.time() - t0
                        print(f'{n}/{len(todo)} {el:.0f}s ({n/el:.1f}/s)', flush=True)
                submit_more()
                break
    for fh in fhs:
        fh.close()
    dfh.close()
    print('done', n, f'{time.time()-t0:.0f}s', flush=True)

if __name__ == '__main__':
    main()
