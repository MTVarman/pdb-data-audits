"""Add to each archive residue pair what the entry's header file says about it:
which _struct_conn record types join the two residues, which join the closest atom pair itself, and whether
that atom pair is listed in _pdbx_validate_close_contact (same asymmetric unit) or _pdbx_validate_symm_contact.
Residues are matched on chain, residue name and residue number, atoms on atom name as well.
Usage: python3 -I add_header_flags.py residue_pairs_in.csv.gz headers.jsonl residue_pairs_out.csv.gz header_files_out.csv.gz"""
import csv, gzip, json, sys, collections
src, hdr, dst, files_out = sys.argv[1:5]
H = {}
with gzip.open(files_out, 'wt', newline='') as f:
    w = csv.writer(f)
    w.writerow(['entry', 'url', 'retrieved_utc', 'bytes', 'sha256', 'status', 'n_struct_conn', 'n_close_contact', 'n_symm_contact'])
    seen = set()
    for line in open(hdr):
        r = json.loads(line)
        if r['entry'] in seen: continue
        seen.add(r['entry'])
        w.writerow([r['entry'], r['url'], r.get('retrieved_utc', ''), r.get('bytes', ''), r.get('sha256', ''), r['status'],
                    len(r.get('conn', [])), len(r.get('close', [])), len(r.get('symm', []))])
        if r['status'] != 'ok': continue
        conn = collections.defaultdict(set); conn_atoms = collections.defaultdict(set)
        for c in r['conn']:
            a = (c[1], c[2], c[3]); b = (c[8], c[9], c[10])
            conn[frozenset([a, b])].add(c[0].lower())
            conn_atoms[frozenset([a + (c[5],), b + (c[12],)])].add(c[0].lower())
        close, close_res, symm = {}, set(), {}
        for c in r['close']:
            a = (c[1], c[2], c[3], c[6]); b = (c[7], c[8], c[9], c[12])
            close[frozenset([a, b])] = c[13]; close_res.add(frozenset([a[:3], b[:3]]))
        for c in r['symm']:
            a = (c[1], c[2], c[3], c[6]); b = (c[7], c[8], c[9], c[12])
            symm[frozenset([a, b])] = c[13]
        H[r['entry']] = (conn, close, close_res, symm, conn_atoms)
rows = list(csv.DictReader(gzip.open(src, 'rt', newline='')))
extra = ['header_read', 'conn_types', 'conn_types_same_atoms', 'disulf', 'close_listed', 'close_dist', 'close_listed_other_atoms', 'symm_listed']
cols = [c for c in rows[0].keys() if c not in extra] + extra
n = collections.Counter()
with gzip.open(dst, 'wt', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
    for r in rows:
        h = H.get(r['entry'])
        if h is None:
            r.update(header_read=0, conn_types='', conn_types_same_atoms='', disulf='', close_listed='', close_dist='', close_listed_other_atoms='', symm_listed='')
            n['no_header'] += 1
        else:
            conn, close, close_res, symm, conn_atoms = h
            p = (r['prot_chain'], r['prot_res'], r['prot_seq']); q = (r['lig_chain'], r['lig_comp'], r['lig_seq'])
            types = sorted(t for t in conn.get(frozenset([p, q]), ()))
            k = frozenset([p + (r['prot_atom'],), q + (r['lig_atom'],)])
            cl = k in close
            r.update(header_read=1, conn_types='|'.join(types), conn_types_same_atoms='|'.join(sorted(conn_atoms.get(k, ()))), disulf=int('disulf' in types), close_listed=int(cl), close_dist=close.get(k, ''),
                     close_listed_other_atoms=int((not cl) and frozenset([p, q]) in close_res), symm_listed=int(k in symm))
            n['with_header'] += 1
            cov = any(t.startswith('covale') for t in types)
            n[('scan_recorded', r['recorded'], 'header_covale', cov)] += 1
        w.writerow(r)
for k, v in sorted(n.items(), key=str): print(k, v)
