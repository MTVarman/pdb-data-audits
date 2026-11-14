"""Cys SG -> non-polymer ligand heavy-atom contact scan with _struct_conn match.

Conventions follow Varman's anchor-bond scan (antenna note, scan.py) and two-rulers index:
  model 1 only; hydrogens removed; occupancy > 0 on both atoms; alternate conformers
  must match or be blank; symmetry-aware distance via gemmi NeighborSearch + nearest PBC image.
Ligand = residue of a non-polymer entity, not water, more than one heavy atom (drops ions/metals),
and not composed only of metal atoms. Every pair (Cys residue, ligand residue instance) with any
SG-ligand heavy-atom distance <= 2.6 A is written, keeping the closest atom pair for that ligand
instance (over alternate conformers).
Record match (residue level, the primary call): a _struct_conn row with conn_type_id 'covale*'
joining this Cys (auth chain, auth seq, ins code) to this ligand instance (auth chain, auth seq,
ins code, comp), in either partner order, any symmetry operator. Atom-level and looser flags kept too.

usage: 03_scan.py <entries.txt> <out.jsonl> [workers]
"""
import sys, os, json, gzip, time
import multiprocessing as mp
import gemmi
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'); CIF = f'{D}/cif'
R = 2.6
WATER = {'HOH', 'DOD', 'WAT', 'H2O'}
MAXBYTES = int(os.environ.get('MAXBYTES', 25_000_000))

def v(x):
    if x is None: return None
    x = gemmi.cif.as_string(x)
    return None if x in ('?', '.', '') else x

def first(b, tag):
    x = b.find_value(tag)
    return v(x) if x is not None else None

def meta(b):
    m = {'process_site': first(b, '_pdbx_database_status.process_site'),
         'deposit_date': first(b, '_pdbx_database_status.recvd_initial_deposition_date'),
         'method': None, 'resolution': None}
    meths = [v(x) for x in b.find_values('_exptl.method')]
    m['method'] = ';'.join(x for x in meths if x) or None
    for tag in ('_refine.ls_d_res_high', '_em_3d_reconstruction.resolution', '_reflns.d_resolution_high'):
        x = first(b, tag)
        if x:
            try: m['resolution'] = float(x); break
            except ValueError: pass
    revs = [v(x) for x in b.find_values('_pdbx_audit_revision_history.revision_date')]
    revs = [x for x in revs if x]
    m['initial_release'] = min(revs) if revs else None
    m['last_revision'] = max(revs) if revs else None
    return m

def conns(b):
    cols = ['conn_type_id', 'ptnr1_auth_asym_id', 'ptnr1_auth_seq_id', 'pdbx_ptnr1_PDB_ins_code', 'ptnr1_label_comp_id', 'ptnr1_label_atom_id', 'pdbx_ptnr1_label_alt_id',
            'ptnr2_auth_asym_id', 'ptnr2_auth_seq_id', 'pdbx_ptnr2_PDB_ins_code', 'ptnr2_label_comp_id', 'ptnr2_label_atom_id', 'pdbx_ptnr2_label_alt_id',
            'ptnr1_symmetry', 'ptnr2_symmetry', 'pdbx_dist_value', 'ptnr1_auth_comp_id', 'ptnr2_auth_comp_id']
    cat = b.get_mmcif_category('_struct_conn.')   # dict tag -> list (None for ?/.)
    out = []
    if not cat: return out
    n = len(next(iter(cat.values())))
    for i in range(n):
        d = {}
        for c in cols:
            x = cat.get(c)
            x = x[i] if x is not None else None
            d[c] = None if x in (None, False, True, '?', '.', '') else str(x)
        out.append(d)
    return out

def rkey(chain, seq, ins):
    return (chain, str(seq), (ins or '').strip())

def scan_entry(e):
    p = f'{CIF}/{e}.cif.gz'
    res = {'entry': e}
    try:
        if not os.path.exists(p): res['status'] = 'no_file'; return res
        if os.path.getsize(p) > MAXBYTES: res['status'] = 'skipped_too_large'; return res
        doc = gemmi.cif.read(p); b = doc.sole_block()
        res.update(meta(b))
        cl = conns(b)
        cov = [c for c in cl if (c['conn_type_id'] or '').lower().startswith('covale')]
        pairs_rec = {}   # frozenset of residue keys -> list of (atom1, atom2, dist, types)
        res_cov = {}     # residue key -> set of partner (key, comp)
        anytype = {}
        for c in cl:
            k1 = rkey(c['ptnr1_auth_asym_id'], c['ptnr1_auth_seq_id'], c['pdbx_ptnr1_PDB_ins_code']) + (c['ptnr1_auth_comp_id'] or c['ptnr1_label_comp_id'],)
            k2 = rkey(c['ptnr2_auth_asym_id'], c['ptnr2_auth_seq_id'], c['pdbx_ptnr2_PDB_ins_code']) + (c['ptnr2_auth_comp_id'] or c['ptnr2_label_comp_id'],)
            fk = frozenset([k1, k2])
            anytype.setdefault(fk, set()).add(c['conn_type_id'])
            if (c['conn_type_id'] or '').lower().startswith('covale'):
                pairs_rec.setdefault(fk, []).append({k1: c['ptnr1_label_atom_id'], k2: c['ptnr2_label_atom_id'], 'd': c['pdbx_dist_value'], 'sym': (c['ptnr1_symmetry'], c['ptnr2_symmetry'])})
                res_cov.setdefault(k1, set()).add(k2); res_cov.setdefault(k2, set()).add(k1)
        st = gemmi.make_structure_from_block(b); st.setup_entities(); st.remove_hydrogens()
        del doc, b
        if len(st) == 0: res['status'] = 'no_model'; return res
        m = st[0]
        res['n_atoms'] = m.count_atom_sites()
        res['symmetry_checked'] = True
        try:
            ns = gemmi.NeighborSearch(m, st.cell, 5).populate()
        except RuntimeError:          # non-standard crystal frame: fall back to asymmetric unit only
            st.cell = gemmi.UnitCell(); res['symmetry_checked'] = False
            ns = gemmi.NeighborSearch(m, st.cell, 5).populate()
        crystal = st.cell.is_crystal()
        lig_ok = {}
        out = []
        for ch in m:
            for r in ch:
                if r.name != 'CYS': continue
                best = {}
                for a in r:
                    if a.name != 'SG' or a.occ <= 0: continue
                    for mk in ns.find_atoms(a.pos, '\0', radius=R):
                        cra = mk.to_cra(m)
                        lr = cra.residue
                        if lr.entity_type != gemmi.EntityType.NonPolymer or lr.name in WATER: continue
                        lk = (cra.chain.name, lr.seqid.num, lr.seqid.icode, lr.name)
                        if lk not in lig_ok:
                            heavy = [x for x in lr if not x.is_hydrogen()]
                            lig_ok[lk] = len(heavy) > 1 and not all(x.element.is_metal for x in heavy)
                        if not lig_ok[lk]: continue
                        la = cra.atom
                        if la.occ <= 0 or la.is_hydrogen() or la.element.is_metal: continue
                        ec = 'C' if la.element.name == 'C' else ('S' if la.element.name == 'S' else 'X')
                        if a.has_altloc() and la.has_altloc() and a.altloc != la.altloc: continue
                        if mk.image_idx == 0 or not crystal:
                            d = a.pos.dist(la.pos); img = 0
                        else:
                            d = a.pos.dist(st.cell.find_nearest_pbc_position(a.pos, la.pos, mk.image_idx)); img = mk.image_idx
                        if d > R: continue
                        if (lk, ec) not in best or d < best[(lk, ec)][0]:
                            best[(lk, ec)] = (d, la.name, la.element.name, a.altloc, la.altloc, la.occ, a.occ, img, la.b_iso, a.b_iso)
                if not best: continue
                ck = rkey(ch.name, r.seqid.num, r.seqid.icode) + ('CYS',)
                nearest = {}
                for (lk, ec), val in best.items():
                    if ec not in nearest or val[0] < best[(nearest[ec], ec)][0]: nearest[ec] = lk
                for (lk, ec), (d, lan, lel, calt, lalt, locc, cocc, img, lb, cb) in best.items():
                    lkey = rkey(lk[0], lk[1], lk[2]) + (lk[3],)
                    fk = frozenset([ck, lkey])
                    recs = pairs_rec.get(fk, [])
                    atom_exact = any(x.get(ck) == 'SG' and x.get(lkey) == lan for x in recs)
                    sg_any = any(x.get(ck) == 'SG' for x in recs)
                    cys_to_comp = any(p[3] == lk[3] for p in res_cov.get(ck, ()))         # Cys bonded to same comp in any copy
                    lig_other = [p for p in res_cov.get(lkey, ()) if p != ck and p[3] != 'HOH']
                    lig_other_prot = [p for p in lig_other if p[3] in AA]
                    cys_other = [p for p in res_cov.get(ck, ()) if p != lkey]
                    out.append({'cys_chain': ch.name, 'cys_seq': r.seqid.num, 'cys_ins': r.seqid.icode.strip(), 'cys_alt': calt.strip() if calt != '\0' else '',
                                'lig_comp': lk[3], 'lig_chain': lk[0], 'lig_seq': lk[1], 'lig_ins': lk[2].strip(), 'lig_atom': lan, 'lig_el': lel,
                                'lig_alt': lalt.strip() if lalt != '\0' else '', 'lig_occ': round(locc, 2), 'cys_occ': round(cocc, 2),
                                'lig_b': round(lb, 1), 'sg_b': round(cb, 1),
                                'dist': round(d, 3), 'sym_image': img, 'el_class': ec, 'nearest_lig_for_cys': lk == nearest[ec],
                                'recorded': bool(recs), 'rec_atom_exact': atom_exact, 'rec_sg_any': sg_any,
                                'rec_dist': ';'.join(str(x['d']) for x in recs),
                                'conn_types_pair': ';'.join(sorted(t for t in anytype.get(fk, set()) if t)),
                                'cys_cov_same_comp_any_copy': cys_to_comp,
                                'lig_cov_other_protein_res': ';'.join(f'{p[0]}:{p[3]}{p[1]}' for p in lig_other_prot),
                                'cys_cov_other': ';'.join(f'{p[0]}:{p[3]}{p[1]}' for p in cys_other)})
        res['pairs'] = out; res['status'] = 'ok'
    except Exception as ex:
        res['status'] = 'error'; res['error'] = repr(ex)[:300]
    return res

AA = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','SEC','MSE'}

if __name__ == '__main__':
    lst, outp = sys.argv[1], sys.argv[2]; nw = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    ids = [x.strip().split()[0] for x in open(lst) if x.strip() and not x.startswith('entry')]
    done = set()
    if os.path.exists(outp):
        for l in open(outp):
            try: done.add(json.loads(l)['entry'])
            except Exception: pass
    todo = [e for e in ids if e not in done]
    print('todo', len(todo), 'done', len(done), flush=True)
    with open(outp, 'a') as f, mp.Pool(nw, maxtasksperchild=200) as pool:
        for i, r in enumerate(pool.imap_unordered(scan_entry, todo, chunksize=4)):
            f.write(json.dumps(r) + '\n')
            if i % 500 == 0: f.flush(); print(i, time.strftime('%H:%M:%S'), flush=True)
    print('finished', flush=True)
