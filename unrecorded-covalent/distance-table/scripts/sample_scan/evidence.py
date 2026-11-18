"""Contact table + wwPDB evidence for the evidence-sample scan.

Usage: python3 -I evidence.py scan_dir work_dir out_prefix

Builds one row per (entry, ligand instance, protein nucleophile atom) keeping the
closest ligand atom, then attaches:
  - probe's bond line for the element pair (covalent radii sum + 0.2 A) and the
    band (line .. line+0.40)
  - _struct_conn status for that atom pair, and for the residue pair
  - CCD-derived ligand-atom tier (reactive centre / sp3-with-heteroatom / neither)
  - for rows in a wide distance window: wwPDB validation-report RSCC, RSR,
    occupancy for both partners and whether the pair is a clash outlier
Density thresholds, fixed before reading any report: both partners RSCC >= 0.80
and RSR <= 0.30 and occupancy > 0.
"""
import collections, csv, glob, gzip, json, os, re, sys, time, urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
import gemmi

COV = {'C': 0.77, 'N': 0.70, 'O': 0.66, 'S': 1.04, 'B': 0.86, 'P': 1.10, 'Se': 1.17, 'F': 0.58}
FUDGE = 0.20
BAND = 0.40
WINDOW_LO, WINDOW_HI = -0.45, 0.45   # offset from line, rows that get a validation report
VAL = {'C': 4, 'B': 3, 'P': 5, 'S': 2, 'Se': 2}
ORD = {'SING': 1, 'DOUB': 2, 'TRIP': 3, 'AROM': 1.5}

COF = set('''PEB PUB PCB CYC PVB BLA BLV LBV BLR RBV
HEC HEM HEB HAS HEO
FAD FMN FDA FAE
PLP PMP P1T PLR PDP PXG PLG PPG
BTN BTI LPA LPB
RET
PLM MYR STE OLA LAU DAO PAM MYA
FAR GER GRG FPP GPP
PNS PPS 4PS
GSH GTT GSF GDS
ACE FOR NAD NAP NAI TPP TDP MQ7 MQ8 UQ1 UQ2'''.split())
NUCLEOTIDE_LIKE = set('''ATP ADP AMP ANP ACP APC AGS ADX GTP GDP GMP GNP GCP GSP GAV
CTP CDP CMP UTP UDP UMP TTP TDP TMP DTP DCP DGT DUT ITP IDP
NAD NAP NDP NAI FAD FMN COA ACO SAM SAH'''.split())
ADDITIVE = set('BME DTT DTU DTV DTD SGM MPO TCE EDO ME2 1PE MRD PEG PGE GOL ACT ACY EPE TRS MES'.split())
HBOND_PAIRS = {('O', 'O'), ('O', 'N'), ('N', 'O'), ('N', 'N')}


def get(url, path, tries=3, timeout=120):
    if os.path.exists(path):
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


def val_report(pdbid, wd):
    e = pdbid.lower()
    p = os.path.join(wd, f'{e}_validation.xml.gz')
    try:
        data = get(f'https://files.rcsb.org/pub/pdb/validation_reports/{e[1:3]}/{e}/{e}_validation.xml.gz', p)
        return ET.fromstring(gzip.decompress(data)).findall('ModelledSubgroup')
    except Exception:
        return None


def site(subs, row):
    """RSCC / RSR / occupancy for both partners and clash at the flagged atom pair."""
    pch, pres, patom = row['protein_atom'].split(':')
    m = re.match(r'([A-Z]+)(-?\d+)(\w?)', pres)
    pname, pnum = m.group(1), m.group(2)
    lname, lch, lnum = row['ligand'].rsplit('_', 2)
    lnum = re.match(r'-?\d+', lnum).group(0)

    def pick(ch, num, nm):
        return [s for s in subs if s.get('model', '1') == '1' and s.get('chain') == ch
                and s.get('resnum') == num and s.get('resname') == nm]
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
    # any clash at all between the two residues (looser match, reported separately)
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
            occ_ok = row['prot_occ'] > 0 and row['lig_occ'] > 0
            out['density_class'] = ('density_supported' if (okp and okl and occ_ok) else
                                    'ligand_uncertain' if okp else
                                    'protein_uncertain' if okl else 'unsupported')
        else:
            out['density_class'] = 'no_density_data'
    except Exception:
        out['density_class'] = 'no_density_data'
    return out


def ccd_fetch(comp, d):
    p = os.path.join(d, f'{comp}.cif')
    if not os.path.exists(p):
        try:
            data = urllib.request.urlopen(f'https://files.rcsb.org/ligands/download/{comp}.cif', timeout=60).read()
            open(p, 'wb').write(data)
        except Exception:
            return None
    return p


def ccd_atoms(path):
    b = gemmi.cif.read_file(path).sole_block()
    atoms = {r[0].strip('"'): (r[1].capitalize() if len(r[1]) > 1 else r[1], r[2])
             for r in b.find('_chem_comp_atom.', ['atom_id', 'type_symbol', 'pdbx_aromatic_flag'])}
    nb = collections.defaultdict(list)
    for r in b.find('_chem_comp_bond.', ['atom_id_1', 'atom_id_2', 'value_order', 'pdbx_aromatic_flag']):
        a1, a2 = r[0].strip('"'), r[1].strip('"')
        nb[a1].append((a2, r[2].upper(), r[3]))
        nb[a2].append((a1, r[2].upper(), r[3]))
    name = (b.find_value('_chem_comp.name') or '').strip('"; \n')
    ctype = (b.find_value('_chem_comp.type') or '').strip('"; \n')
    return atoms, nb, name, ctype


def classify(atoms, nb, atom):
    if atom not in atoms:
        return '?', 'atom not in CCD'
    el, arom = atoms[atom]
    bonds = nb.get(atom, [])
    val = sum(ORD.get(o, 1) for _, o, _ in bonds)
    els = [atoms.get(o, ('?', ''))[0] for o, _, _ in bonds]
    nH = sum(1 for e in els if e in ('H', 'D'))
    heavy = [e for e in els if e not in ('H', 'D')]
    orders = [o for _, o, _ in bonds]
    if el in ('B', 'P', 'Se'):
        return 'A', el
    if el == 'S':
        if nH >= 1 or val < VAL['S']:
            return 'A', 'thiol / open S'
        return ('C', 'thioether / ring S') if val <= 2 else ('B', 'oxidised S (sulfonyl etc.)')
    if el == 'O':
        if nH >= 1 or val < 2:
            return 'B', 'hydroxyl / open O'
        return 'C', 'ether / carbonyl O'
    if el == 'N':
        if val < 3:
            return 'A', 'open-valence N'
        return 'C', 'amine / amide N'
    if el != 'C':
        return 'C', el
    if val < 4:
        return 'A', 'open valence C'
    if arom == 'Y' or any(f == 'Y' for _, _, f in bonds):
        return 'C', 'aromatic C'
    if 'TRIP' in orders:
        return 'A', 'C#N / C#C'
    if 'DOUB' in orders:
        dbl = [atoms.get(o, ('?', ''))[0] for o, ord_, _ in bonds if ord_ == 'DOUB']
        return 'A', 'C=' + ('O' if 'O' in dbl else 'N' if 'N' in dbl else 'S' if 'S' in dbl else 'C')
    if any(e in ('O', 'N', 'S', 'F', 'Cl', 'CL', 'Br', 'BR', 'I', 'Se', 'SE', 'B', 'P') for e in heavy):
        return 'B', 'sp3 C with heteroatom'

    def has_hetero_double(a):
        return any(o == 'DOUB' and atoms.get(x, ('?', ''))[0] in ('O', 'N', 'S') for x, o, _ in nb.get(a, []))
    n1 = [x for x, _, _ in bonds if atoms.get(x, ('?', ''))[0] not in ('H', 'D')]
    n2 = {y for x in n1 for y, _, _ in nb.get(x, []) if y != atom and atoms.get(y, ('?', ''))[0] not in ('H', 'D')}
    if any(has_hetero_double(x) for x in n1):
        return 'B', 'sp3 C alpha to C=O/S=O/C=N'
    if any(has_hetero_double(y) for y in n2):
        return 'B', 'sp3 C beta to C=O/S=O/C=N'
    if nH >= 3:
        return 'C', 'methyl'
    return 'C', 'sp3 C, C/H only'


def main():
    scan_dir, wd, pre = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(wd, exist_ok=True)
    ccd_dir = os.path.join(wd, 'ccd')
    os.makedirs(ccd_dir, exist_ok=True)
    arms = json.load(open('/home/claude/work/lane3/arms.json'))['arms']
    arm_of = {}
    for a, lst in arms.items():
        for e in lst:
            arm_of.setdefault(e, []).append(a[0])

    rows, meta, st = {}, {}, collections.Counter()
    cov_comp = collections.defaultdict(set)      # component -> entries where it has a covale to a nucleophile
    for p in sorted(glob.glob(os.path.join(scan_dir, 'shard_*.jsonl'))):
        for line in open(p):
            r = json.loads(line)
            st[r.get('status')] += 1
            meta[r['pdbid']] = {'method': r.get('method', ''), 'resolution': r.get('resolution'),
                                'title': r.get('title') or '', 'keywords': r.get('keywords') or '',
                                'deposited': r.get('deposited') or '', 'r_free': r.get('r_free') or ''}
            for c in r.get('contacts', []):
                if any(t == 'Covale' for t in c['conn_pair']):
                    cov_comp[c['lig'].split('_')[0]].add(r['pdbid'])
                k = (r['pdbid'], c['lig'], c['prot'])
                if k not in rows or c['dist'] < rows[k]['dist']:
                    rows[k] = {**c, 'pdbid': r['pdbid']}
    print('entries parsed', sum(st.values()), dict(st), flush=True)
    print('raw contact rows (closest atom per ligand-instance/nucleophile-atom):', len(rows), flush=True)

    comps = sorted({k[1].split('_')[0] for k in rows})
    print('components to fetch from the CCD:', len(comps), flush=True)
    with ThreadPoolExecutor(16) as ex:
        paths = dict(zip(comps, ex.map(lambda c: ccd_fetch(c, ccd_dir), comps)))
    cache = {}

    out = []
    for (pdbid, lig, prot), c in rows.items():
        comp = lig.split('_')[0]
        if comp not in cache:
            try:
                cache[comp] = ccd_atoms(paths[comp]) if paths.get(comp) else None
            except Exception:
                cache[comp] = None
        if cache[comp]:
            tier, kind = classify(cache[comp][0], cache[comp][1], c['lig_atom'])
            cname, ctype = cache[comp][2], cache[comp][3]
        else:
            tier, kind, cname, ctype = '?', 'CCD unavailable', '', ''
        pres = prot.split(':')[1][:3]
        pel = c['prot_el']
        lel = c['lig_el']
        line = round(COV.get(lel, 0) + COV.get(pel, 0) + FUDGE, 3)
        d = c['dist']
        off = round(d - line, 3)
        fam = ''
        if comp in COF:
            fam = 'class1_cofactor_ptm'
        elif comp in ADDITIVE:
            fam = 'additive_cryoprotectant'
        elif (comp in NUCLEOTIDE_LIKE or lel == 'P') and pres in ('ASP', 'GLU', 'LYS', 'SER', 'THR', 'HIS', 'TYR'):
            fam = 'class3_nucleotide_phosphate'
        pa, la = c['prot_altloc'], c['lig_altloc']
        occ_class = ('altloc_crossing' if (pa and la and pa != la) else
                     'full' if c['prot_occ'] >= 0.99 and c['lig_occ'] >= 0.99 else 'partial')
        recorded = any(t == 'Covale' for t in c['conn_pair'])
        out.append({
            'pdbid': pdbid, 'arm': ''.join(sorted(arm_of.get(pdbid, []))),
            'ligand': lig, 'ligand_comp': comp, 'ligand_name': cname[:90], 'ligand_type': ctype[:40],
            'ligand_nheavy': c['lig_nheavy'], 'lig_entity': c['lig_entity'],
            'protein_atom': prot, 'residue': pres, 'ligand_atom': c['lig_atom'],
            'pair': f'{lel}-{pel}', 'lig_el': lel, 'prot_el': pel,
            'dist_A': d, 'probe_line_A': line, 'offset_from_line_A': off,
            'in_band': line <= d < line + BAND,
            'past_line': d >= line,
            'recorded_covale': recorded,
            'conn_pair': ' '.join(c['conn_pair']),
            'residue_level_conn': ' '.join(c['conn_residues']),
            'residue_already_bonded': bool(re.search('Covale|Disulf|Metalc', ' '.join(c['conn_residues']))) and not recorded,
            'occupancy_class': occ_class, 'prot_altloc': pa, 'lig_altloc': la,
            'prot_occ': c['prot_occ'], 'lig_occ': c['lig_occ'],
            'prot_b': c['prot_b'], 'lig_b': c['lig_b'],
            'ligand_atom_tier': tier, 'ligand_atom_kind': kind, 'family': fam,
            'hbond_capable_pair': (lel, pel) in HBOND_PAIRS,
            'method': meta[pdbid]['method'], 'resolution': meta[pdbid]['resolution'],
            'deposited': meta[pdbid]['deposited'],
            'title': meta[pdbid]['title'][:160], 'keywords': meta[pdbid]['keywords'][:80],
            'n_entries_same_comp_recorded_covalent': len(cov_comp.get(comp, set()) - {pdbid}),
            'same_comp_recorded_covalent_in_this_entry': pdbid in cov_comp.get(comp, set()),
        })
    out.sort(key=lambda x: (x['pair'], x['dist_A'], x['pdbid']))
    with open(f'{pre}_contacts_all.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print('contact rows written:', len(out), flush=True)

    # validation reports for rows inside the wide window around the line
    win = [x for x in out if WINDOW_LO <= x['offset_from_line_A'] < WINDOW_HI and x['residue'] != 'CYS']
    cys = [x for x in out if WINDOW_LO <= x['offset_from_line_A'] < WINDOW_HI and x['residue'] == 'CYS']
    need = sorted({x['pdbid'] for x in win} | {x['pdbid'] for x in cys})
    print(f'rows in window (non-Cys) {len(win)}, Cys control {len(cys)}, reports needed {len(need)}', flush=True)
    with ThreadPoolExecutor(16) as ex:
        reps = dict(zip(need, ex.map(lambda e: val_report(e, wd), need)))
    miss = [e for e in need if not reps.get(e)]
    print('validation reports retrieved', len(need) - len(miss), 'missing', len(miss), flush=True)
    ev = []
    for x in win + cys:
        subs = reps.get(x['pdbid'])
        s = (site(subs, x) if subs else
             {'prot_rscc': '', 'prot_rsr': '', 'lig_rscc': '', 'lig_rsr': '', 'clash_listed': '',
              'clashmag': '', 'clash_elsewhere_same_residue_pair': '', 'in_report': False,
              'density_class': 'no_report'})
        ev.append({**x, **s})
    ev.sort(key=lambda x: (x['pair'], x['dist_A']))
    with open(f'{pre}_window_evidence.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(ev[0].keys()))
        w.writeheader()
        w.writerows(ev)
    print('window evidence rows written:', len(ev), flush=True)
    json.dump({'scan_status': dict(st), 'rows_all': len(out), 'rows_window': len(ev),
               'reports_missing': miss, 'thresholds': {'rscc_min': 0.80, 'rsr_max': 0.30, 'occ_min': 0.0},
               'fudge': FUDGE, 'band': BAND},
              open(f'{pre}_run.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
