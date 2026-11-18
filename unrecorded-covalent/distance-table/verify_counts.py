#!/usr/bin/env python3
"""Recompute the numbers quoted in the note from the shipped tables and compare each with the value in the text.
Standard library only.  Run from this folder:  python3 verify_counts.py
Prints PASS or FAIL per check and a total.  Takes about a minute (it reads the 968,000-row atom-pair table).

recorded       a _struct_conn row of type covale or disulf joins the pair
               (archive scan: the two residues; evidence sample: the two atoms)
close contact  the entry file's own _pdbx_validate_close_contact table has a row naming the two atoms
distances      as stored in the tables (0.001 A), except the pairs in boundary_exact_distances.csv, whose rounded
               distance lies exactly on a boundary; for those the unrounded distance is used"""
import collections, csv, gzip, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, 'data')


def rows(name):
    path = os.path.join(D, name)
    fh = gzip.open(path, 'rt', newline='') if name.endswith('.gz') else open(path, newline='')
    with fh:
        yield from csv.DictReader(fh)


RESULTS = []


def check(label, got, want):
    ok = got == want
    RESULTS.append(ok)
    print(('PASS  ' if ok else 'FAIL  ') + label + ('' if ok else f'   computed {got!r}, note says {want!r}'))


def pct(k, n):
    return round(100.0 * k / n, 1)


# ------------------------------------------------------------------------------------------------ the rules
TABLE1 = {  # pair: (window lower, window upper, validator line, between) as printed in Table 1
    'C-S': ('1.55', '2.00', '2.01', '0.01'), 'C-O': ('1.10', '1.55', '1.63', '0.08'), 'C-N': ('1.20', '1.60', '1.67', '0.07'),
    'B-O': ('1.20', '1.65', '1.72', '0.07'), 'B-N': ('1.20', '1.65', '1.76', '0.11'), 'S-S': ('1.85', '2.20', '2.28', '0.08'),
    'O-P': ('1.40', '1.75', '1.96', '0.21'), 'N-P': ('1.54', '1.77', '2.00', '0.23'), 'O-S': ('1.35', '1.80', '1.90', '0.10'),
    'N-S': ('1.45', '1.90', '1.94', '0.04'), 'C-Se': ('1.75', '2.15', '2.14', '-0.01')}
maxit = collections.defaultdict(set)
releases = {}
for r in rows('maxit_bond_table_by_release.csv'):
    a, b = r['pair'].split('_')
    maxit['-'.join([a, b.capitalize()])].add((r['lower_limit_A'], r['upper_limit_A']))
    releases[r['maxit_release']] = r['release_date']
check('MAXIT releases tabulated: 11, from 8.120 (2008-02-12) to 11.400 (2025-10-16)',
      (len(releases), releases.get('8.120'), releases.get('11.400')), (11, '2008-02-12', '2025-10-16'))
# probe radii from the shipped rows of atomprops.h
cov, vdw = {}, {}
for line in open(os.path.join(HERE, 'sources', 'probe_atomprops_rows.txt')):
    m = re.match(r'\s*\{atom\w+,\s*\d+,\s*"(\w+)",\s*([\w.]+),\s*[\w.]+,\s*[\d.]+,\s*([\d.]+),', line)
    if m:
        cov[m.group(1)] = float(m.group(3))
        vdw[m.group(1)] = 1.70 if m.group(2) == 'ATOMC_EXPLICIT_VDW' else float(m.group(2))
for pair, (lo, hi, line, gap) in TABLE1.items():
    check(f'Table 1 {pair}: one window in all 11 releases, {lo}-{hi}', maxit[pair], {(lo, hi)})
    a, b = pair.split('-')
    pl = cov[a] + cov[b] + 0.2
    check(f'Table 1 {pair}: validator line {line}', f'{pl:.2f}', line)
    check(f'Table 1 {pair}: between the two {gap}', f'{pl - float(hi):.2f}', gap)
check('probe: boron van der Waals radius 0.20', vdw['B'], 0.20)
check('probe: oxygen van der Waals radius 1.40', vdw['O'], 1.40)
check('boron-oxygen pair touches at 1.60 and reaches a 0.4 A overlap only at 1.20 or less',
      (round(vdw['B'] + vdw['O'], 2), round(vdw['B'] + vdw['O'] - 0.4, 2)), (1.60, 1.20))

WIN = {p: (float(v[0]), float(v[1])) for p, v in TABLE1.items()}
LINE = {p: float(v[2]) for p, v in TABLE1.items()}
CUT = 2.195    # MAXIT: a pair at 2.195 A or more is not listed as a close contact


def zone(d, pair):
    lo, hi = WIN[pair]
    if d < lo:
        return 'below'
    if d <= hi:
        return 'inside'
    if d <= LINE[pair]:
        return 'between'
    if d <= LINE[pair] + 0.40:
        return 'past'
    return 'far'


def pairname(e1, e2):
    return '-'.join(sorted([e1.capitalize(), e2.capitalize()]))


# ------------------------------------------------------------------------------------------------ exact distances on a boundary
BXA, BXS = {}, {}
for b in rows('boundary_exact_distances.csv'):
    k = (b['entry'], b['protein_chain'], b['protein_residue'], b['protein_seq'], b['protein_atom'],
         b['ligand_chain'], b['ligand_comp'], b['ligand_seq'], b['ligand_atom'])
    exact = float(b['dist_exact'])
    if b['source'] == 'archive':
        BXA[k] = exact
    else:
        BXS[k + (b['protein_altloc'], b['ligand_altloc'])] = exact
    if abs(exact - float(b['dist_rounded'])) > 0.00051:
        print('      boundary distance does not round to the stored value:', b)
check('re-measured on a boundary: 63 archive pairs, 14 sample contacts', (len(BXA), len(BXS)), (63, 14))

# ------------------------------------------------------------------------------------------------ archive scan
ns = collections.Counter(r['status'] for r in rows('archive_not_scanned.csv'))
check('archive: 161,735 entries scanned of 161,748; 8 parse errors, 5 skipped for size',
      (sum(1 for _ in rows('archive_per_entry.csv.gz')), sum(1 for _ in rows('archive_per_entry.csv.gz')) + sum(ns.values()), ns['parse error'], ns['skipped for size']),
      (161735, 161748, 8, 5))
hf = list(rows('archive_header_files.csv.gz'))
check('header files read for all 64,856 entries with a pair', (len(hf), sum(r['status'] == 'ok' for r in hf)), (64856, 64856))
A = []
scan_vs_header = collections.Counter()
for r in rows('archive_residue_pairs.csv.gz'):
    try:
        res = float(r['resolution'])
    except ValueError:
        res = None
    try:
        dep = int(r['dep_year'])
    except ValueError:
        dep = None
    k = (r['entry'], r['prot_chain'], r['prot_res'], r['prot_seq'], r['prot_atom'], r['lig_chain'], r['lig_comp'], r['lig_seq'], r['lig_atom'])
    types = r['conn_types'].split('|')
    scan_vs_header[(r['recorded'] == '1') == any(t.startswith('covale') for t in types)] += 1
    A.append({'entry': r['entry'], 'res': r['prot_res'], 'atom': r['prot_atom'], 'comp': r['lig_comp'], 'latom': r['lig_atom'],
              'pair': pairname(r['prot_el'], r['lig_el']), 'd': BXA.get(k, float(r['dist'])), 'dep': dep,
              'rec': r['recorded'] == '1' or r['disulf'] == '1',
              'rec_atoms': any(t.startswith('covale') or t == 'disulf' for t in r['conn_types_same_atoms'].split('|')),
              'x': r['method1'].startswith('X-RAY') and res is not None and res <= 2.5,
              'close': r['close_listed'] == '1', 'symm': r['symm_listed'] == '1'})
check('archive: residue pairs', len(A), 226753)
check('archive: entries with a pair', len({r['entry'] for r in A}), 64856)
check('the scan and the header file agree on a covale record for every pair', scan_vs_header[True], 226753)
X = [r for r in A if r['x']]
check('archive, X-ray at 2.5 A or better: 118,373 residue pairs in 46,742 entries', (len(X), len({r['entry'] for r in X})), (118373, 46742))

GROUPS = [('Cys S to ligand C', ('CYS',), 'C-S'), ('Ser, Thr, Tyr O to C', ('SER', 'THR', 'TYR'), 'C-O'), ('Asp, Glu O to C', ('ASP', 'GLU'), 'C-O'),
          ('Lys N to C', ('LYS',), 'C-N'), ('His N to C', ('HIS',), 'C-N'), ('Cys S to ligand S', ('CYS',), 'S-S'), ('Ser, Thr, Tyr O to P', ('SER', 'THR', 'TYR'), 'O-P')]
TABLE2 = {  # recorded/all: inside, between, past (2014 on); between (before 2014)
    'Cys S to ligand C': ((6796, 7096), (2, 39), (31, 682), (4, 45)), 'Ser, Thr, Tyr O to C': ((860, 895), (5, 23), (5, 52), (11, 14)),
    'Asp, Glu O to C': ((386, 436), (0, 17), (0, 64), (38, 78)), 'Lys N to C': ((1965, 2010), (7, 18), (31, 54), (16, 19)),
    'His N to C': ((191, 237), (0, 13), (9, 26), (18, 21)), 'Cys S to ligand S': ((390, 417), (0, 14), (0, 22), (4, 43)),
    'Ser, Thr, Tyr O to P': ((226, 234), (1, 3), (0, 1), (4, 6))}
tot = collections.Counter()
for label, residues, pair in GROUPS:
    g = [r for r in X if r['res'] in residues and r['pair'] == pair and r['dep'] is not None]

    def cell(z, late):
        s = [r for r in g if zone(r['d'], pair) == z and (r['dep'] >= 2014) == late]
        return (sum(r['rec'] for r in s), len(s))
    got = (cell('inside', True), cell('between', True), cell('past', True), cell('between', False))
    check(f'Table 2 {label}', got, TABLE2[label])
    for name, (k, n) in zip(('in', 'bt', 'pa', 'bt_old'), got):
        tot[name + '_k'] += k
        tot[name + '_n'] += n
check('Table 2 all seven, inside: 10,814 / 11,325 (95.5%)', (tot['in_k'], tot['in_n'], pct(tot['in_k'], tot['in_n'])), (10814, 11325, 95.5))
check('Table 2 all seven, between: 15 / 127 (11.8%)', (tot['bt_k'], tot['bt_n'], pct(tot['bt_k'], tot['bt_n'])), (15, 127, 11.8))
check('Table 2 all seven, past: 76 / 901 (8.4%)', (tot['pa_k'], tot['pa_n'], pct(tot['pa_k'], tot['pa_n'])), (76, 901, 8.4))
check('Table 2 all seven, between, before 2014: 95 / 226 (42.0%)', (tot['bt_old_k'], tot['bt_old_n'], pct(tot['bt_old_k'], tot['bt_old_n'])), (95, 226, 42.0))
check('inside the window with no record since 2014: 4.5%', round(100 - pct(tot['in_k'], tot['in_n']), 1), 4.5)
G7 = [r for r in X if r['dep'] is not None and r['dep'] >= 2014 and any(r['res'] in g[1] and r['pair'] == g[2] for g in GROUPS)]
bey = [r for r in G7 if r['rec'] and zone(r['d'], r['pair']) in ('between', 'past')]
cc = collections.Counter(r['comp'] for r in bey)
check('recorded beyond the window since 2014: 91 in 50 entries; PLP 37, FAD 8, PEB 7',
      (len(bey), len({r['entry'] for r in bey}), cc['PLP'], cc['FAD'], cc['PEB']), (91, 50, 37, 8, 7))
s = [r for r in X if r['res'] == 'LYS' and r['pair'] == 'C-N' and r['dep'] is not None and r['dep'] >= 2014 and 1.60 < r['d'] <= 2.07 and r['rec']]
check('recorded lysine contacts between 1.60 and 2.07 A since 2014: 38, of which PLP 37', (len(s), sum(r['comp'] == 'PLP' for r in s)), (38, 37))
atoms = tuple(sum(r['rec_atoms'] for r in G7 if zone(r['d'], r['pair']) == z) for z in ('inside', 'between', 'past'))
check('Table 2 since 2014 counting only records on the two atoms: 10,809 / 15 / 73', atoms, (10809, 15, 73))
old = [r for r in X if r['res'] == 'CYS' and r['pair'] == 'C-S' and r['dep'] is not None and r['dep'] <= 2004]
lo_ = [r for r in old if 1.90 <= r['d'] <= 2.00]
hi_ = [r for r in old if 2.00 < r['d'] <= 2.10]
check('cysteine deposited up to 2004: 108 of 111 recorded at 1.90-2.00 A, 5 of 60 at 2.00-2.10 A',
      (sum(r['rec'] for r in lo_), len(lo_), sum(r['rec'] for r in hi_), len(hi_)), (108, 111, 5, 60))

BONDPAIRS = ('C-S', 'C-O', 'C-N', 'S-S', 'O-P', 'O-S', 'N-S', 'N-P')


def grouped(data):
    g = collections.defaultdict(list)
    for r in data:
        if r['pair'] in BONDPAIRS:
            g[(r['entry'], r['res'], r['atom'], r['comp'], r['latom'], r['pair'])].append(r)
    return g


mixed = collections.Counter()
for key, v in grouped(X).items():
    rec = [r['d'] for r in v if r['rec']]
    un = [r['d'] for r in v if not r['rec']]
    if rec and un:
        lo, hi = WIN[key[5]]
        sep = all(lo <= d <= hi for d in rec) and not any(lo <= d <= hi for d in un)
        late = v[0]['dep'] is not None and v[0]['dep'] >= 2014
        mixed['n'] += 1
        mixed['sep'] += sep
        mixed[key[5] + ' n'] += 1
        mixed[key[5] + ' sep'] += sep
        mixed[('late' if late else 'early') + ' n'] += 1
        mixed[('late' if late else 'early') + ' sep'] += sep
check('inside one entry: 250 mixed groups, 230 separated exactly by the window', (mixed['n'], mixed['sep']), (250, 230))
for p, want in (('C-S', (173, 182)), ('C-N', (26, 28)), ('C-O', (12, 14)), ('S-S', (15, 21)), ('O-P', (4, 4)), ('O-S', (0, 1))):
    check(f'inside one entry, {p}: {want[0]} of {want[1]}', (mixed[p + ' sep'], mixed[p + ' n']), want)
check('inside one entry: 122 of 134 before 2014, 108 of 116 since', (mixed['early sep'], mixed['early n'], mixed['late sep'], mixed['late n']), (122, 134, 108, 116))

strad = collections.Counter()
examples = {}
for key, v in grouped(A).items():
    hi = WIN[key[5]][1]
    below = [r for r in v if hi - 0.03 < r['d'] <= hi]
    above = [r for r in v if hi < r['d'] <= hi + 0.03]
    if below and above:
        b = max(below, key=lambda r: r['d'])
        a = min(above, key=lambda r: r['d'])
        kind = 'split' if b['rec'] and not a['rec'] else 'both' if b['rec'] and a['rec'] else 'neither' if not b['rec'] and not a['rec'] else 'reversed'
        strad[kind] += 1
        strad['n'] += 1
        if kind == 'split':
            strad['split_longer_listed'] += a['close']
            if not a['close']:
                examples['unlisted'] = (key[0], key[5], a['d'])
        examples[(key[0], key[3])] = (b['d'], b['rec'], a['d'], a['rec'])
check('at the limit: 91 groups; 67 split, 17 both recorded, 7 neither, 0 reversed',
      (strad['n'], strad['split'], strad['both'], strad['neither'], strad['reversed']), (91, 67, 17, 7, 0))
check('example 6QVM heme c: 1.999 recorded, 2.004 not', examples[('6QVM', 'HEC')], (1.999, True, 2.004, False))
check('example 4XLO flavin-histidine: 1.594 recorded, 1.604 not', examples[('4XLO', 'FAD')], (1.594, True, 1.604, False))
check('example 8ERJ lysine adduct: 1.597 recorded, 1.608 not', examples[('8ERJ', 'WQF')], (1.597, True, 1.608, False))
check('split groups: the longer copy is a close contact in 66 of 67; the other is S-S at 2.23 A',
      (strad['split_longer_listed'], examples.get('unlisted')), (66, ('2WOW', 'S-S', 2.23)))

short = [r for r in A if r['res'] == 'CYS' and r['pair'] == 'C-S' and r['d'] < 1.55]
inside = [r for r in A if r['res'] == 'CYS' and r['pair'] == 'C-S' and 1.55 <= r['d'] <= 2.00]
xs = [r for r in short if r['x']]
xi = [r for r in inside if r['x']]
check('below the window, X-ray: 24 of 64 Cys-C contacts with a record on the two atoms; inside 11,406 of 11,795',
      (sum(r['rec_atoms'] for r in xs), len(xs), sum(r['rec_atoms'] for r in xi), len(xi)), (24, 64, 11406, 11795))
check('below the window, all methods: 2,493 in 152 entries, 1,625 in 7EZX; 37 on the two atoms, 561 per residue pair',
      (len(short), len({r['entry'] for r in short}), sum(r['entry'] == '7EZX' for r in short), sum(r['rec_atoms'] for r in short), sum(r['rec'] for r in short)),
      (2493, 152, 1625, 37, 561))


def between_unrecorded(data):
    s = [r for r in data if r['pair'] in BONDPAIRS and zone(r['d'], r['pair']) == 'between' and not r['rec']]
    return (len(s), len({r['entry'] for r in s}))


check('archive, between the lines, no record: 246 pairs in 187 X-ray entries', between_unrecorded(X), (246, 187))
check('archive, between the lines, no record: 367 pairs in 258 entries, all methods', between_unrecorded(A), (367, 258))

# ------------------------------------------------------------------------------------------------ close-contact tables, archive
B = [r for r in A if r['pair'] in BONDPAIRS]
ra = [r for r in B if r['rec_atoms']]
check('Table 3 archive: a record joins the two atoms, 37,357 pairs in 10,139 entries, 0 listed',
      (len(ra), len({r['entry'] for r in ra}), sum(r['close'] for r in ra)), (37357, 10139, 0))
un = [r for r in B if not r['rec']]
lo = [r for r in un if r['d'] < CUT]
hi = [r for r in un if r['d'] >= CUT]
check('Table 3 archive: no record, under 2.195 A: 8,276 pairs, 8,137 listed in 1,608 entries',
      (len(lo), sum(r['close'] for r in lo), len({r['entry'] for r in lo if r['close']})), (8276, 8137, 1608))
nl = [r for r in lo if not r['close']]
check('of the other 139: 42 in the symmetry-contact table', sum(r['symm'] for r in nl), 42)
ul = list(rows('unlisted_pairs_check.csv'))
check('and 97 with both atoms below full occupancy', (len(ul), sum(r['reason'] == 'both atoms below full occupancy' for r in ul), len(nl) - sum(r['symm'] for r in nl)), (97, 97, 97))
check('Table 3 archive: no record, 2.195 A or more: 5,319 pairs, 0 listed', (len(hi), sum(r['close'] for r in hi)), (5319, 0))
st1 = [r for r in un if 2.190 <= r['d'] < 2.195]
st2 = [r for r in un if 2.195 <= r['d'] < 2.200]
check('archive: 52 of 53 listed at 2.190-2.195 A, 0 of 50 at 2.195-2.200 A', (sum(r['close'] for r in st1), len(st1), sum(r['close'] for r in st2), len(st2)), (52, 53, 0, 50))
xl = [r for r in X if r['pair'] in BONDPAIRS and not r['rec'] and r['d'] < CUT]
check('X-ray at 2.5 A or better: 2,310 unrecorded pairs under the cut-off: 91 below, 780 inside, 193 between, 1,152 past, 94 further',
      (len(xl),) + tuple(sum(zone(r['d'], r['pair']) == z for r in xl) for z in ('below', 'inside', 'between', 'past', 'far')), (2310, 91, 780, 193, 1152, 94))

# ------------------------------------------------------------------------------------------------ 4DL1
e = [(float(r['distance_A']), r['recorded'] == 'True', r['listed_as_close_contact'] == 'True', r['close_contact_row']) for r in rows('esters_4DL1.csv')]
check('4DL1: sixteen esters, eight recorded', (len(e), sum(x[1] for x in e)), (16, 8))
check('4DL1: longest recorded 1.54999, shortest unrecorded 1.55026', (max(x[0] for x in e if x[1]), min(x[0] for x in e if not x[1])), (1.54999, 1.55026))
check('4DL1: all between 1.531 and 1.562', (round(min(x[0] for x in e), 3), round(max(x[0] for x in e), 3)), (1.531, 1.562))
check('4DL1: the eight unrecorded are close-contact rows 1 to 8; no recorded one is listed',
      (sorted(int(x[3]) for x in e if not x[1] and x[2]), sum(x[2] for x in e if x[1])), ([1, 2, 3, 4, 5, 6, 7, 8], 0))

# ------------------------------------------------------------------------------------------------ evidence sample
sel = json.load(open(os.path.join(D, 'sample_selection.json')))
check('sample: arms 11,426 / 6,477 / 34,880 of a pool of 139,519; union 48,072',
      (sel['counts']['A_text'], sel['counts']['B_warhead'], sel['counts']['C_baseline'], sel['baseline_pool_size'], sel['union']), (11426, 6477, 34880, 139519, 48072))


def full(r):
    return r['occupancy_class'] == 'full' and float(r['prot_occ'] or 0) > 0 and float(r['lig_occ'] or 0) > 0


CC = collections.defaultdict(set)
ncc = 0
for c in rows('sample_close_contact_rows.csv.gz'):
    CC[c['pdbid']].add(frozenset([(c['chain_1'], c['comp_1'], c['seq_1'], c['atom_1']), (c['chain_2'], c['comp_2'], c['seq_2'], c['atom_2'])]))
    ncc += 1
check('sample: 91,505 close-contact rows in 10,839 entries', (ncc, sum(1 for v in CC.values() if v)), (91505, 10839))


def prep(r):
    ch, resx, atom = r['protein_atom'].split(':')
    seq = resx[len(r['residue']):]                     # residue number, with any insertion code
    num = re.match(r'-?\d+', seq).group(0)             # the close-contact rows carry the number without it
    comp, lch, lseq = r['ligand'].rsplit('_', 2)
    lseq = re.match(r'-?\d+', lseq).group(0)
    r['cp'] = pairname(r['lig_el'], r['prot_el'])
    r['d'] = BXS.get((r['pdbid'], ch, r['residue'], seq, atom, lch, comp, lseq, r['ligand_atom'], r['prot_altloc'], r['lig_altloc']), float(r['dist_A']))
    r['rec'] = r['conn_pair'] in ('Covale', 'Disulf')
    r['hb'] = r['hbond_capable_pair'] == 'True'
    r['close'] = frozenset([(ch, r['residue'], num, atom), (lch, comp, lseq, r['ligand_atom'])]) in CC[r['pdbid']]
    return r


S = [prep(r) for r in rows('sample_contacts.csv.gz')]
check('sample: 55,913 contacts in 20,425 entries, 46,653 at full occupancy', (len(S), len({r['pdbid'] for r in S}), sum(1 for r in S if full(r))), (55913, 20425, 46653))
SB = [r for r in S if full(r) and r['cp'] in WIN and not r['hb']]
sr = [r for r in SB if r['rec']]
check('Table 3 sample: a record joins the two atoms, 8,364 in 3,779 entries, 0 listed', (len(sr), len({r['pdbid'] for r in sr}), sum(r['close'] for r in sr)), (8364, 3779, 0))
su = [r for r in SB if not r['rec'] and r['d'] < CUT]
check('Table 3 sample: no record, under 2.195 A, 990 in 479 entries, 990 listed', (len(su), len({r['pdbid'] for r in su}), sum(r['close'] for r in su)), (990, 479, 990))
so = [r for r in SB if not r['rec'] and r['d'] >= CUT]
check('Table 3 sample: no record, 2.195 A or more, 812, 0 listed', (len(so), sum(r['close'] for r in so)), (812, 0))
s1 = [r for r in SB if not r['rec'] and 2.190 <= r['d'] < 2.195]
s2 = [r for r in SB if not r['rec'] and 2.195 <= r['d'] < 2.200]
check('sample: 7 of 7 listed at 2.190-2.195 A, 0 of 6 at 2.195-2.200 A', (sum(r['close'] for r in s1), len(s1), sum(r['close'] for r in s2), len(s2)), (7, 7, 0, 6))

E = []
for r in rows('sample_window_evidence.csv.gz'):
    prep(r)
    r['clash'] = r['clash_listed'] == 'True'
    r['inrep'] = r['in_report'] == 'True'
    r['dens'] = r['density_class'] == 'density_supported'
    E.append(r)
check('sample: validation evidence read for 11,792 contacts', len(E), 11792)
EF = [r for r in E if full(r) and r['cp'] in WIN and not r['hb']]
zb = [r for r in EF if zone(r['d'], r['cp']) == 'between']
zbu = [r for r in zb if not r['rec']]
check('Table 4 all pairs: 76 recorded; 141 unrecorded in 95 entries, all in a report, 125 density-supported, 0 clash, 117 close contacts',
      (sum(r['rec'] for r in zb), len(zbu), len({r['pdbid'] for r in zbu}), sum(r['inrep'] for r in zbu), sum(r['dens'] for r in zbu), sum(r['clash'] for r in zbu), sum(r['close'] for r in zbu)),
      (76, 141, 95, 141, 125, 0, 117))
TABLE4 = [('Cys, C-S', ('CYS',), 'C-S', (5, 27, 23, 25, 0, 27)), ('Cys, S-S', ('CYS',), 'S-S', (4, 24, 18, 15, 0, 0)), ('Ser and Thr, C-O', ('SER', 'THR'), 'C-O', (10, 7, 7, 6, 0, 7)),
          ('Asp and Glu, C-O', ('ASP', 'GLU'), 'C-O', (5, 45, 24, 44, 0, 45)), ('Lys and His, C-N', ('LYS', 'HIS'), 'C-N', (17, 11, 8, 10, 0, 11)),
          ('Ser and Thr, B-O', ('SER', 'THR'), 'B-O', (27, 22, 11, 22, 0, 22))]
named = 0
for label, residues, cp, want in TABLE4:
    s = [r for r in zbu if r['residue'] in residues and r['cp'] == cp]
    sr_ = [r for r in zb if r['rec'] and r['residue'] in residues and r['cp'] == cp]
    named += len(s)
    check(f'Table 4 {label}', (len(sr_), len(s), len({r['pdbid'] for r in s}), sum(r['dens'] for r in s), sum(r['clash'] for r in s), sum(r['close'] for r in s)), want)
check('Table 4: five unrecorded contacts on pairs not listed', len(zbu) - named, 5)
check('Table 4: every unrecorded contact missing from the close-contact column is S-S at 2.195 A or more',
      all(r['cp'] == 'S-S' and r['d'] >= CUT for r in zbu if not r['close']), True)

inrep = [r for r in EF if r['inrep']]
ins = [r for r in inrep if not r['rec'] and r['d'] <= LINE[r['cp']]]
past = [r for r in inrep if not r['rec'] and zone(r['d'], r['cp']) == 'past']
past_rec = [r for r in inrep if r['rec'] and zone(r['d'], r['cp']) == 'past']


def frac(s):
    return (sum(r['clash'] for r in s), len(s), sum(r['close'] for r in s))


check('Table 5: unrecorded up to the validator line: clash 0 / 500, close contact 476', frac(ins), (0, 500, 476))
check('Table 5: past the line, cysteine: clash 460 / 470, close contact 309', frac([r for r in past if r['residue'] == 'CYS']), (460, 470, 309))
check('Table 5: past the line, other residues, pairs without boron: clash 85 / 86, close contact 80',
      frac([r for r in past if r['residue'] != 'CYS' and 'B' not in r['cp'].split('-')]), (85, 86, 80))
sb = [r for r in past if r['residue'] == 'SER' and r['cp'] == 'B-O']
check('Table 5: past the line, serine to boron: clash 0 / 68, close contact 68', frac(sb), (0, 68, 68))
check('Table 5: recorded contacts past the line: clash 0 / 136, close contact 0', frac(past_rec), (0, 136, 0))
check('Table 5: every unrecorded contact missing from the close-contact column is at 2.195 A or more',
      all(r['d'] >= CUT for r in ins + past if not r['close']), True)
check('serine-boron past the line: 21 entries, 67 density-supported, deposited 2009 to 2025',
      (len({r['pdbid'] for r in sb}), sum(r['dens'] for r in sb), min(r['deposited'][:4] for r in sb), max(r['deposited'][:4] for r in sb)), (21, 67, '2009', '2025'))
check('no other residue has an unrecorded boron contact past the line', len([r for r in past if r['residue'] != 'SER' and 'B' in r['cp'].split('-')]), 0)

bo = [r for r in S if full(r) and r['cp'] == 'B-O' and r['residue'] in ('SER', 'THR') and int(r['deposited'][:4]) >= 2014]
bz = {z: [r for r in bo if zone(r['d'], 'B-O') == z] for z in ('inside', 'between', 'past')}
check('serine and threonine boron since 2014: inside 232 / 238, between 16 / 37, past 5 / 72',
      tuple((sum(r['rec'] for r in bz[z]), len(bz[z])) for z in ('inside', 'between', 'past')), ((232, 238), (16, 37), (5, 72)))
check('14 of the 16 recorded between the lines are in 6N80', sum(1 for r in bz['between'] if r['rec'] and r['pdbid'] == '6N80'), 14)

bn = [r for r in rows('boron_neighbour_clashes.csv') if r['zone'] == 'past_probe_line_to_2.12']
nuc = [int(r['n_records_nucleophile_O_or_H_to_same_boronate']) for r in bn]
els = collections.Counter()
for r in bn:
    for p in r['clash_partners'].split('; '):
        if re.match(r'^(OG1?|HG1?)~', p):
            els[re.search(r'\((\w+)\)', p).group(1)] += 1
check('serine-boron past the line: 68 contacts, none with a clash record naming the boron', (len(bn), sum(int(r['clash_records_naming_the_boron']) for r in bn)), (68, 0))
check('in 57 the serine O or H is listed against another atom of the boronate: 97 records, 63 to O, 34 to C', (sum(1 for n in nuc if n), sum(nuc), els['O'], els['C']), (57, 97, 63, 34))
check('all 97 name an atom bonded to the boron', sum(int(r['n_of_those_naming_an_atom_bonded_to_the_boron']) for r in bn), 97)
bl = list(rows('boron_ligand_clash_records.csv'))
check('211 entries, 536 boron-containing ligand residues, 561 clash records, none naming a boron',
      (len({r['pdbid'] for r in bl}), len(bl), sum(int(r['clash_records']) for r in bl), sum(int(r['clash_records_naming_boron']) for r in bl)), (211, 536, 561, 0))
check('the same 211 entries as the boron contact table', {r['pdbid'] for r in bl} == {r['pdbid'] for r in rows('boron_clash_check.csv')}, True)

tri = list(rows('past_line_triage.csv'))
c = collections.Counter(r['triage'] for r in tri)
check('other contacts past the line: 221 in 134 entries', (len(tri), len({r['pdbid'] for r in tri})), (221, 134))
ALL_LINE = dict(LINE, **{'O-O': 1.52, 'N-O': 1.56, 'N-N': 1.60})   # covalent radii plus 0.2 A, for the pairs that can also be hydrogen bonds
sel_tri = {(r['pdbid'], r['protein_atom'], r['ligand'], r['ligand_atom']) for r in E
           if full(r) and r['cp'] in ALL_LINE and not r['rec'] and r['residue'] != 'CYS' and ALL_LINE[r['cp']] < r['d'] <= ALL_LINE[r['cp']] + 0.40}
check('the triage table is exactly the unrecorded non-cysteine contacts from the line to 0.40 A beyond it',
      sel_tri == {(r['pdbid'], r['protein_atom'], r['ligand'], r['ligand_atom']) for r in tri}, True)
check('triage: 86 candidate bonds (68 serine-boron), 42 cofactor, 67 O/N pairs, 12 bond or near-attack, 14 other',
      (c['candidate bond'], sum(1 for r in tri if r['triage'] == 'candidate bond' and r['pair'] == 'B-O'), c['cofactor with known covalent attachment'], c['O-O, N-O or N-N pair'],
       c['bond or near-attack complex'], c['unreactive ligand atom'] + c['additive'] + c['phosphoryl or nucleotide']), (86, 68, 42, 67, 12, 14))

# vildagliptin, 6B1E
v = {r['protein_atom']: r for r in E if r['pdbid'] == '6B1E' and r['ligand_comp'] == 'LF7' and r['pair'] == 'C-O'}
check('6B1E chain B: 1.422 A, recorded', (v['B:SER630:OG']['d'], v['B:SER630:OG']['rec']), (1.422, True))
check('6B1E chain A: 1.569 A, not recorded, a close contact, not a clash',
      (v['A:SER630:OG']['d'], v['A:SER630:OG']['rec'], v['A:SER630:OG']['close'], v['A:SER630:OG']['clash']), (1.569, False, True, False))
check('6B1E: both at full occupancy, ligand RSCC 0.96',
      ([full(x) for x in v.values()], sorted(round(float(x['lig_rscc']), 2) for x in v.values())), ([True, True], [0.96, 0.96]))

# vaborbactam
vb = list(rows('vaborbactam_entries.csv'))
rec_ = [r for r in vb if r['recorded_covale'] == 'True']
unr = [r for r in vb if r['recorded_covale'] != 'True']
check('vaborbactam: 13 entries, 22 serine-boron contacts', (len({r['pdbid'] for r in vb}), len(vb)), (13, 22))
check('vaborbactam: the 17 at 1.622 A or less are recorded and none is a close contact',
      (len(rec_), max(float(r['distance_A']) for r in rec_), any(r['listed_as_close_contact'] == 'True' for r in rec_)), (17, 1.6217, False))
check('vaborbactam: the 5 longer ones are close contacts and none is recorded',
      (len(unr), min(float(r['distance_A']) for r in unr), all(r['listed_as_close_contact'] == 'True' for r in unr)), (5, 1.6544, True))
check('vaborbactam: 8BTW 1.683, 8BV4 1.685, 8TMT 1.878',
      sorted((r['pdbid'], round(float(r['distance_A']), 3)) for r in unr if r['pdbid'] != '6V7H'), [('8BTW', 1.683), ('8BV4', 1.685), ('8TMT', 1.878)])
h = sorted((r['ligand_conformer'], round(float(r['distance_A']), 3), r['recorded_covale']) for r in vb if r['pdbid'] == '6V7H')
check('6V7H: conformer A recorded at 1.611 and 1.622, conformer B unrecorded at 1.654 and 1.663',
      h, [('A', 1.611, 'True'), ('A', 1.622, 'True'), ('B', 1.654, 'False'), ('B', 1.663, 'False')])

# ------------------------------------------------------------------------------------------------ the two scans against each other
L3 = {}
ents = set()
for r in S:
    ch, resx, atom = r['protein_atom'].split(':')
    comp, lch, lseq = r['ligand'].rsplit('_', 2)
    L3[(r['pdbid'], ch, r['residue'], resx[len(r['residue']):], atom, lch, comp, lseq, r['ligand_atom'])] = (
        float(r['dist_A']), r['recorded_covale'] == 'True' or 'Covale' in (r['residue_level_conn'] or ''))
    ents.add(r['pdbid'])
both = agree = close = 0
for r in rows('archive_atom_pairs.csv.gz'):
    if r['entry'] not in ents:
        continue
    k = (r['entry'], r['prot_chain'], r['prot_res'], r['prot_seq'], r['prot_atom'], r['lig_chain'], r['lig_comp'], r['lig_seq'], r['lig_atom'])
    if k in L3:
        both += 1
        agree += (r['recorded'] == '1') == L3[k][1]
        close += abs(L3[k][0] - float(r['dist'])) <= 0.0011
check('the two scans: 49,538 atom pairs in both, covale status agrees for all, distance to 0.001 A for 49,347', (both, agree, close), (49538, 49538, 49347))

print(f'\n{sum(RESULTS)} of {len(RESULTS)} checks pass')
sys.exit(0 if all(RESULTS) else 1)
