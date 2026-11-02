"""Evidence layer for the whole-PDB sweep.

Usage: python3 -I archive_evidence.py archive.jsonl triage_contacts.csv work_dir out_prefix

1. Probe-line test (pre-registered before reading: clash-listed share ~0% below
   probe's bond line, ~100% above it). All Cys SG - ligand C pairs, no bond record,
   full occupancy, X-ray, 1.70-2.50 A; read each entry's wwPDB validation report.
2. Candidates: Cys SG - ligand C/S, 2.0-2.5 A, full occupancy, not class 1/3.
   For each: density (RSCC/RSR both partners), MolProbity clash at the pair, entry
   title, and whether the same ligand component is recorded as covalently bonded
   (_struct_conn covale) to a protein residue in any entry of the sweep.
"""
import csv, gzip, json, os, re, sys, urllib.request, collections
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

def get(url, path):
    if os.path.exists(path):
        return open(path, 'rb').read()
    for i in range(3):
        try:
            data = urllib.request.urlopen(url, timeout=120).read()
            open(path, 'wb').write(data)
            return data
        except Exception as e:  # noqa
            last = e
    raise RuntimeError(str(last))

def val_report(pdbid, wd):
    e = pdbid.lower()
    p = os.path.join(wd, f'{e}_validation.xml.gz')
    try:
        data = get(f'https://files.rcsb.org/pub/pdb/validation_reports/{e[1:3]}/{e}/{e}_validation.xml.gz', p)
        return ET.fromstring(gzip.decompress(data)).findall('ModelledSubgroup')
    except Exception:
        return None

def title(pdbid, wd):
    p = os.path.join(wd, f'{pdbid.lower()}.json')
    try:
        d = json.loads(get(f'https://data.rcsb.org/rest/v1/core/entry/{pdbid.upper()}', p))
        kw = (d.get('struct_keywords') or {}).get('pdbx_keywords', '') or ''
        return d.get('struct', {}).get('title', ''), kw
    except Exception:
        return '', ''

def site(subs, row):
    pch, pres, patom = row['protein_atom'].split(':')
    m = re.match(r'([A-Z]+)(-?\d+)(\w?)', pres)
    pname, pnum = m.group(1), m.group(2)
    lname, lch, lnum = row['ligand'].rsplit('_', 2)
    lnum = re.match(r'-?\d+', lnum).group(0)
    pick = lambda ch, num, nm: [s for s in subs if s.get('model', '1') == '1' and s.get('chain') == ch
                                and s.get('resnum') == num and s.get('resname') == nm]
    pc, lc = pick(pch, pnum, pname), pick(lch, lnum, lname)
    worst = lambda c: min((s for s in c if s.get('rscc')), key=lambda s: float(s.get('rscc')), default=c[0] if c else None)
    pg, lg = worst(pc), worst(lc)
    cids = {cl.get('cid') for s in pc for cl in s.findall('clash') if cl.get('atom') == patom}
    clash = ''
    for s in lc:
        for cl in s.findall('clash'):
            if cl.get('cid') in cids and cl.get('atom') == row['ligand_atom']:
                clash = cl.get('clashmag')
    f = lambda g, k: (g.get(k) or '') if g is not None else ''
    out = {'prot_rscc': f(pg, 'rscc'), 'prot_rsr': f(pg, 'rsr'), 'lig_rscc': f(lg, 'rscc'), 'lig_rsr': f(lg, 'rsr'),
           'clash_listed': bool(clash), 'clashmag': clash, 'in_report': bool(pc) and bool(lc)}
    try:
        ok = lambda g: float(g.get('rscc')) >= 0.8 and float(g.get('rsr')) <= 0.3
        if pg is not None and lg is not None and pg.get('rscc') and lg.get('rscc'):
            out['density_class'] = {(True, True): 'density_supported', (True, False): 'ligand_uncertain',
                                    (False, True): 'protein_uncertain', (False, False): 'unsupported'}[(ok(pg), ok(lg))]
        else:
            out['density_class'] = 'no_density_data'
    except Exception:
        out['density_class'] = 'no_density_data'
    return out

ADDITIVE = set('BME DTT DTU DTV DTD SGM MPO TCE EDO ME2 1PE MRD'.split())

def main():
    raw, tri, wd, pre = sys.argv[1:5]
    os.makedirs(wd, exist_ok=True)
    # covalent-elsewhere index: component -> entries where it has a covale record to a nucleophile
    cov_comp = collections.defaultdict(set)
    recorded = []   # recorded Cys SG - ligand C covale bonds, X-ray, full occupancy, 1.95-2.50 A
    for line in open(raw):
        r = json.loads(line)
        for c in r.get('contacts', []):
            if any(t == 'Covale' for t in c['conn_pair']):
                cov_comp[c['lig'].split('_')[0]].add(r['pdbid'])
                if (c['prot'].endswith(':SG') and ':CYS' in c['prot'] and c['lig_el'] == 'C' and 1.95 <= c['dist'] < 2.50
                        and 'X-RAY' in (r.get('method') or '') and c['prot_occ'] >= 0.99 and c['lig_occ'] >= 0.99
                        and not (c['prot_altloc'] and c['lig_altloc'] and c['prot_altloc'] != c['lig_altloc'])):
                    recorded.append({'pdbid': r['pdbid'], 'protein_atom': c['prot'], 'ligand': c['lig'], 'ligand_atom': c['lig_atom'],
                                     'dist_A': c['dist'], 'resolution': r.get('resolution'), 'method': r.get('method')})
    rows = list(csv.DictReader(open(tri)))
    for r in rows:
        r['dist_A'] = float(r['dist_A']); r['resolution'] = float(r['resolution']) if r['resolution'] not in ('', 'None') else None
        r['res_linked'] = bool(re.search('Covale|Disulf', r['residue_level_conn'] or ''))
        if r['ligand'].split('_')[0] in ADDITIVE and not r['family']:
            r['family'] = 'additive_reducing_agent'
    cys_all = [r for r in rows if r['residue'] == 'CYS' and r['protein_atom'].endswith(':SG') and r['pair'] in ('C-S', 'S-S')
               and r['occupancy_class'] == 'full']
    probe_set = [r for r in cys_all if r['pair'] == 'C-S' and 'X-RAY' in r['method'] and 1.70 <= r['dist_A'] < 2.50]
    cand = [r for r in cys_all if 2.0 <= r['dist_A'] < 2.5 and r['family'] == '' and not r['res_linked']]
    fam_counts = collections.Counter(r['family'] or ('residue_already_bonded' if r['res_linked'] else 'candidate')
                                     for r in cys_all if 2.0 <= r['dist_A'] < 2.5)
    print('Cys SG 2.0-2.5 A full-occupancy rows by family:', dict(fam_counts))
    need = sorted({r['pdbid'] for r in probe_set} | {r['pdbid'] for r in cand} | {r['pdbid'] for r in recorded})
    print('recorded C-S bonds 1.95-2.50 A (X-ray, full occ):', len(recorded), 'rows', len({r['pdbid'] for r in recorded}), 'entries')
    print('probe-line set', len(probe_set), 'rows', len({r['pdbid'] for r in probe_set}), 'entries;',
          'candidates', len(cand), 'rows', len({r['pdbid'] for r in cand}), 'entries; reports needed', len(need), flush=True)
    with ThreadPoolExecutor(16) as ex:
        reps = dict(zip(need, ex.map(lambda e: val_report(e, wd), need)))
    with ThreadPoolExecutor(16) as ex:
        tids = sorted({r['pdbid'] for r in cand})
        titles = dict(zip(tids, ex.map(lambda e: title(e, wd), tids)))
    # 1. probe-line test
    pl = []
    for r in probe_set:
        subs = reps.get(r['pdbid'])
        if not subs:
            continue
        s = site(subs, r)
        if not s['in_report']:
            continue
        pl.append({**r, **s})
    bins = collections.OrderedDict()
    for lo in [x / 100 for x in range(170, 250, 5)]:
        b = [x for x in pl if lo <= x['dist_A'] < lo + 0.05]
        bins[f'{lo:.2f}-{lo+0.05:.2f}'] = (len(b), sum(1 for x in b if x['clash_listed']))
    print('\nPROBE LINE (C-S, line 2.01 A): bin  n  clash-listed')
    for k, (n, c) in bins.items():
        print(f'  {k}  {n:5d}  {c:5d}  {100*c/n if n else 0:5.1f}%')
    below = [x for x in pl if x['dist_A'] <= 2.01]; above = [x for x in pl if x['dist_A'] > 2.01]
    print(f'  <=2.01: {sum(x["clash_listed"] for x in below)}/{len(below)} clash-listed;  >2.01: {sum(x["clash_listed"] for x in above)}/{len(above)}')
    with open(f'{pre}_probe_line_test.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(pl[0].keys())); w.writeheader(); w.writerows(pl)
    rl = []
    for r in recorded:
        subs = reps.get(r['pdbid'])
        if not subs:
            continue
        s2 = site(subs, r)
        if s2['in_report']:
            rl.append({**r, **s2})
    print('\nRECORDED bonds (struct_conn covale present), same distance bins:')
    for lo in [x / 100 for x in range(195, 250, 5)]:
        b = [x for x in rl if lo <= x['dist_A'] < lo + 0.05]
        print(f'  {lo:.2f}-{lo+0.05:.2f}  {len(b):5d}  {sum(1 for x in b if x["clash_listed"]):5d}')
    print(f'  recorded >2.01: {sum(x["clash_listed"] for x in rl if x["dist_A"] > 2.01)}/{sum(1 for x in rl if x["dist_A"] > 2.01)} clash-listed;'
          f'  recorded <=2.01: {sum(x["clash_listed"] for x in rl if x["dist_A"] <= 2.01)}/{sum(1 for x in rl if x["dist_A"] <= 2.01)}')
    if rl:
        with open(f'{pre}_probe_line_recorded.csv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rl[0].keys())); w.writeheader(); w.writerows(rl)
    # 2. candidates
    out = []
    for r in cand:
        subs = reps.get(r['pdbid'])
        s = site(subs, r) if subs else {'density_class': 'no_report', 'clash_listed': '', 'clashmag': '', 'in_report': False,
                                        'prot_rscc': '', 'prot_rsr': '', 'lig_rscc': '', 'lig_rsr': ''}
        comp = r['ligand'].split('_')[0]
        others = sorted(cov_comp.get(comp, set()) - {r['pdbid']})
        t, kw = titles.get(r['pdbid'], ('', ''))
        out.append({**r, **s, 'same_ligand_recorded_covalent_in': ' '.join(others[:10]) + (f' (+{len(others)-10})' if len(others) > 10 else ''),
                    'n_entries_same_ligand_covalent': len(others),
                    'same_ligand_covalent_in_this_entry': r['pdbid'] in cov_comp.get(comp, set()),
                    'title_says_covalent': bool(re.search(r'covalent|irreversib|adduct|acyl-enzyme|thioester|warhead', t, re.I)),
                    'title': t, 'keywords': kw})
    out.sort(key=lambda x: x['dist_A'])
    with open(f'{pre}_cys_candidates.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
    xr = [x for x in out if 'X-RAY' in x['method'] and x['resolution'] and x['resolution'] <= 2.5]
    print('\nCANDIDATES (Cys SG, 2.0-2.5 A, full occupancy, not class 1/3):', len(out), 'rows', len({x['pdbid'] for x in out}), 'entries')
    print(' by method', collections.Counter(x['method'] for x in out).most_common(5))
    print(' X-ray <=2.5 A:', len(xr), 'rows', len({x['pdbid'] for x in xr}), 'entries')
    print('  density', collections.Counter(x['density_class'] for x in xr))
    ds = [x for x in xr if x['density_class'] == 'density_supported']
    print('  density_supported:', len(ds), 'rows', len({x['pdbid'] for x in ds}), 'entries')
    ev = lambda x: x['n_entries_same_ligand_covalent'] > 0 or x['same_ligand_covalent_in_this_entry'] or x['title_says_covalent']
    print('  ...with outside evidence (same ligand recorded covalent elsewhere / in this entry, or title):',
          sum(1 for x in ds if ev(x)), 'rows', len({x['pdbid'] for x in ds if ev(x)}), 'entries')
    print('  clash-listed among X-ray candidates in report:', sum(1 for x in xr if x['clash_listed'] is True), '/', sum(1 for x in xr if x['in_report']))

if __name__ == '__main__':
    main()
