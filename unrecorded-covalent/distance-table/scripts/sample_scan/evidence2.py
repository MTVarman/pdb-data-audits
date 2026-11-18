"""Stage 2 of the evidence layer: wwPDB validation values for the window rows.

Usage: python3 -I evidence2.py contacts_all.csv work_dir out_csv
Each entry's validation report is downloaded, parsed, used, and discarded before
the next one, so only one report is in memory at a time.
Density thresholds, fixed before any report was read: both partners RSCC >= 0.80
and RSR <= 0.30 and occupancy > 0.
"""
import collections, csv, gzip, json, os, re, sys, time, urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

WINDOW_LO, WINDOW_HI = -0.45, 0.45


def get(url, path, tries=3, timeout=120):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return open(path, 'rb').read()
    last = None
    for i in range(tries):
        try:
            data = urllib.request.urlopen(url, timeout=timeout).read()
            open(path, 'wb').write(data)
            return data
        except Exception as e:  # noqa
            last = e
            time.sleep(1 + i)
    raise RuntimeError(str(last))


def site(subs, row):
    pch, pres, patom = row['protein_atom'].split(':')
    m = re.match(r'([A-Z]+)(-?\d+)(\w?)', pres)
    pname, pnum = m.group(1), m.group(2)
    lname, lch, lnum = row['ligand'].rsplit('_', 2)
    lnum = re.match(r'-?\d+', lnum).group(0)

    def pick(ch, nm_num, nm):
        return [s for s in subs if s.get('model', '1') == '1' and s.get('chain') == ch
                and s.get('resnum') == nm_num and s.get('resname') == nm]
    pc, lc = pick(pch, pnum, pname), pick(lch, lnum, lname)

    def worst(c):
        return min((s for s in c if s.get('rscc')), key=lambda s: float(s.get('rscc')),
                   default=c[0] if c else None)
    pg, lg = worst(pc), worst(lc)
    cids = {cl.get('cid') for s in pc for cl in s.findall('clash') if cl.get('atom') == patom}
    clash = ''
    for s in lc:
        for cl in s.findall('clash'):
            if cl.get('cid') in cids and cl.get('atom') == row['ligand_atom']:
                clash = cl.get('clashmag')
    allcids = {cl.get('cid') for s in pc for cl in s.findall('clash')}
    loose = any(cl.get('cid') in allcids for s in lc for cl in s.findall('clash'))

    def f(g, k):
        return (g.get(k) or '') if g is not None else ''
    out = {'prot_rscc': f(pg, 'rscc'), 'prot_rsr': f(pg, 'rsr'),
           'lig_rscc': f(lg, 'rscc'), 'lig_rsr': f(lg, 'rsr'),
           'clash_listed': bool(clash), 'clashmag': clash,
           'clash_elsewhere_same_residue_pair': loose,
           'in_report': bool(pc) and bool(lc)}
    try:
        def ok(g):
            return float(g.get('rscc')) >= 0.80 and float(g.get('rsr')) <= 0.30
        if pg is not None and lg is not None and pg.get('rscc') and lg.get('rscc'):
            okp, okl = ok(pg), ok(lg)
            occ_ok = float(row['prot_occ']) > 0 and float(row['lig_occ']) > 0
            out['density_class'] = ('density_supported' if (okp and okl and occ_ok) else
                                    'ligand_uncertain' if okp else
                                    'protein_uncertain' if okl else 'unsupported')
        else:
            out['density_class'] = 'no_density_data'
    except Exception:
        out['density_class'] = 'no_density_data'
    return out


BLANK = {'prot_rscc': '', 'prot_rsr': '', 'lig_rscc': '', 'lig_rsr': '', 'clash_listed': '',
         'clashmag': '', 'clash_elsewhere_same_residue_pair': '', 'in_report': False,
         'density_class': 'no_report'}


def do_entry(args):
    pdbid, rows, wd = args
    e = pdbid.lower()
    p = os.path.join(wd, f'{e}_validation.xml.gz')
    try:
        data = get(f'https://files.rcsb.org/pub/pdb/validation_reports/{e[1:3]}/{e}/{e}_validation.xml.gz', p)
        subs = ET.fromstring(gzip.decompress(data)).findall('ModelledSubgroup')
    except Exception as ex:  # noqa
        return [{**r, **BLANK} for r in rows], str(ex)[:120]
    out = []
    for r in rows:
        try:
            out.append({**r, **site(subs, r)})
        except Exception:
            out.append({**r, **BLANK})
    del subs
    return out, None


def main():
    src, wd, dst = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(wd, exist_ok=True)
    rows = []
    with open(src) as f:
        for r in csv.DictReader(f):
            off = float(r['offset_from_line_A'])
            if WINDOW_LO <= off < WINDOW_HI:
                rows.append(r)
    by = collections.defaultdict(list)
    for r in rows:
        by[r['pdbid']].append(r)
    print(f'window rows {len(rows)}, entries {len(by)}', flush=True)
    tasks = [(k, v, wd) for k, v in sorted(by.items())]
    done, errs, n = [], {}, 0
    t0 = time.time()
    with ThreadPoolExecutor(12) as ex:
        for res, err in ex.map(do_entry, tasks):
            done += res
            n += 1
            if err:
                errs[res[0]['pdbid']] = err
            if n % 400 == 0:
                print(f'{n}/{len(tasks)} entries {time.time()-t0:.0f}s', flush=True)
    done.sort(key=lambda x: (x['pair'], float(x['dist_A'])))
    with open(dst, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(done[0].keys()))
        w.writeheader()
        w.writerows(done)
    print('written', len(done), 'rows;', 'report failures', len(errs), flush=True)
    print('in_report', sum(1 for x in done if x['in_report']),
          'density_supported', sum(1 for x in done if x['density_class'] == 'density_supported'), flush=True)
    json.dump({'window_rows': len(done), 'entries': len(tasks), 'report_failures': errs,
               'thresholds': {'rscc_min': 0.80, 'rsr_max': 0.30, 'occ_min': 0.0}},
              open(dst.replace('.csv', '_run.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
