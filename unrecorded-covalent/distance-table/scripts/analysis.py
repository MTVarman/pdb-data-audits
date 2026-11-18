"""All numbers quoted in the note, computed from the shipped tables.
Usage: python3 -I analysis.py data_dir out_dir
data_dir holds: archive_residue_pairs.csv.gz, sample_contacts.csv.gz, sample_window_evidence.csv.gz,
sample_close_contact_rows.csv.gz, boundary_exact_distances.csv, boron_neighbour_clashes.csv, boron_clash_check.csv,
sample_selection.json

Distances are the scan's, rounded to 0.001 A, except for the pairs whose rounded distance lies exactly on a boundary:
those are taken unrounded from boundary_exact_distances.csv.

Definitions used throughout:
  recorded      a _struct_conn row of type covale or disulf joins the pair
                (archive scan: the two residues; evidence sample: the two atoms)
  close contact the entry file's own _pdbx_validate_close_contact table has a row naming the two atoms
  clash         the entry's wwPDB validation report lists the two atoms as a clash"""
import csv, json, sys, os, re, collections
import gzip
def op(path):
    """Open a table that may be gzip-compressed."""
    return gzip.open(path, 'rt', newline='') if path.endswith('.gz') else open(path, newline='')
D, O = sys.argv[1], sys.argv[2]
os.makedirs(O, exist_ok=True)
# --- the rules -------------------------------------------------------------------------------------
# MAXIT filterlib BondUtil_global.h, _allowed_bonds[]: (lower, upper); identical v8.120 (2008) .. v11.400 (2025)
MAXIT = {'C-S': (1.55, 2.00), 'C-O': (1.10, 1.55), 'C-N': (1.20, 1.60), 'B-O': (1.20, 1.65), 'B-N': (1.20, 1.65),
         'S-S': (1.85, 2.20), 'O-P': (1.40, 1.75), 'O-S': (1.35, 1.80), 'N-S': (1.45, 1.90), 'N-P': (1.54, 1.77),
         'C-Se': (1.75, 2.15), 'O-O': (1.30, 1.60), 'N-O': (1.18, 1.50), 'N-N': (1.15, 1.55)}
# probe atomprops.h covalent radii; bonded if d <= r1 + r2 + 0.2 (COVRADFUDGE)
COV = {'C': 0.77, 'N': 0.70, 'O': 0.66, 'S': 1.04, 'P': 1.10, 'B': 0.86, 'Se': 1.17}
# MAXIT rolodist-util xtal-util.C, contact_filter_check(): a pair at 2.195 A or more is not listed as a close contact
CLOSE_CUT = 2.195
BONDPAIRS = ('C-S', 'C-O', 'C-N', 'S-S', 'O-P', 'O-S', 'N-S', 'N-P')      # archive scan (boron not scanned there)
def canon(a, b):
    a, b = sorted([a, b]); return f'{a}-{b}'
def probe_line(pair):
    a, b = pair.split('-'); return round(COV[a] + COV[b] + 0.2, 2)
def zone(d, pair):
    lo, hi = MAXIT[pair]; pl = probe_line(pair)
    if d < lo: return 'below'
    if d <= hi: return 'inside'
    if d <= pl: return 'between'
    if d <= pl + 0.40: return 'past'
    return 'far'
N = collections.OrderedDict()
def put(k, v): N[k] = v
def dump(name, rows_):
    with open(f'{O}/{name}', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows_[0])); w.writeheader(); w.writerows(rows_)
rules = []
for p in ['C-S', 'C-O', 'C-N', 'B-O', 'B-N', 'S-S', 'O-P', 'N-P', 'O-S', 'N-S', 'C-Se']:
    lo, hi = MAXIT[p]; pl = probe_line(p)
    rules.append({'pair': p, 'maxit_lower': f'{lo:.2f}', 'maxit_upper': f'{hi:.2f}', 'probe_line': f'{pl:.2f}', 'gap': f'{pl - hi:+.2f}'})
    put(f'rule.{p}.upper', hi); put(f'rule.{p}.probe', pl)
dump('table1_rules.csv', rules)

# --- exact distances for pairs on a boundary -----------------------------------------------------
BXA, BXS = {}, {}
for b in csv.DictReader(open(f'{D}/boundary_exact_distances.csv', newline='')):
    k = (b['entry'], b['protein_chain'], b['protein_residue'], b['protein_seq'], b['protein_atom'], b['ligand_chain'], b['ligand_comp'], b['ligand_seq'], b['ligand_atom'])
    if b['source'] == 'archive': BXA[k] = float(b['dist_exact'])
    else: BXS[k + (b['protein_altloc'], b['ligand_altloc'])] = float(b['dist_exact'])
put('boundary.archive_pairs', len(BXA)); put('boundary.sample_contacts', len(BXS))
# --- archive-wide residue pairs --------------------------------------------------------------------
A = []
for r in csv.DictReader(op(f'{D}/archive_residue_pairs.csv.gz')):
    r['d'] = BXA.get((r['entry'], r['prot_chain'], r['prot_res'], r['prot_seq'], r['prot_atom'], r['lig_chain'], r['lig_comp'], r['lig_seq'], r['lig_atom']), float(r['dist']))
    r['rec'] = r['recorded'] == '1' or r['disulf'] == '1'
    r['hdr'] = r['header_read'] == '1'; r['close'] = r['close_listed'] == '1'; r['symm'] = r['symm_listed'] == '1'
    r['close_other'] = r['close_listed_other_atoms'] == '1'
    r['rec_atoms'] = any(t.startswith('covale') or t == 'disulf' for t in r['conn_types_same_atoms'].split('|'))
    r['pair'] = canon(r['prot_el'].capitalize(), r['lig_el'].capitalize())
    try: r['res_A'] = float(r['resolution'])
    except Exception: r['res_A'] = None
    r['xray25'] = r['method1'].startswith('X-RAY') and r['res_A'] is not None and r['res_A'] <= 2.5
    try: r['dep'] = int(r['dep_year'])
    except Exception: r['dep'] = None
    A.append(r)
put('archive.residue_pairs', len(A)); put('archive.entries', len({r['entry'] for r in A}))
put('archive.header_read.pairs', sum(r['hdr'] for r in A)); put('archive.header_read.entries', len({r['entry'] for r in A if r['hdr']}))
put('archive.recorded_as_disulf_only', sum(1 for r in A if r['recorded'] != '1' and r['disulf'] == '1'))
X = [r for r in A if r['xray25']]
put('archive.xray25.residue_pairs', len(X)); put('archive.xray25.entries', len({r['entry'] for r in X}))
GROUPS = [('cys_CS', 'Cys SG to C', ('CYS',), 'C-S'), ('hydroxyl_CO', 'Ser/Thr/Tyr O to C', ('SER', 'THR', 'TYR'), 'C-O'),
          ('carboxylate_CO', 'Asp/Glu O to C', ('ASP', 'GLU'), 'C-O'), ('lys_CN', 'Lys NZ to C', ('LYS',), 'C-N'),
          ('his_CN', 'His N to C', ('HIS',), 'C-N'), ('cys_SS', 'Cys SG to S', ('CYS',), 'S-S'), ('ser_OP', 'Ser/Thr/Tyr O to P', ('SER', 'THR', 'TYR'), 'O-P')]
t2 = []
for gid, label, residues, pair in GROUPS:
    for era in ('all', '2014on', 'pre2014'):
        g = [r for r in X if r['prot_res'] in residues and r['pair'] == pair and
             (era == 'all' or (r['dep'] is not None and ((r['dep'] >= 2014) == (era == '2014on'))))]
        row = {'group': gid, 'label': label, 'pair': pair, 'era': era}
        for z in ('below', 'inside', 'between', 'past'):
            s = [r for r in g if zone(r['d'], pair) == z]
            ents = collections.defaultdict(list)
            for r in s: ents[r['entry']].append(r['rec'])
            k = sum(r['rec'] for r in s)
            row[f'{z}_recorded'] = k; row[f'{z}_n'] = len(s)
            row[f'{z}_entries'] = len(ents); row[f'{z}_entries_none_recorded'] = sum(1 for v in ents.values() if not any(v))
            put(f'A.{gid}.{era}.{z}', [k, len(s), len(ents), row[f'{z}_entries_none_recorded']])
        t2.append(row)
dump('table2_archive_zones.csv', t2)
for era in ('all', '2014on', 'pre2014'):
    for z in ('inside', 'between', 'past'):
        k = sum(N[f'A.{g[0]}.{era}.{z}'][0] for g in GROUPS); n = sum(N[f'A.{g[0]}.{era}.{z}'][1] for g in GROUPS)
        put(f'A.total.{era}.{z}', [k, n, round(100.0 * k / n, 1)])
# the same totals counting a record only when it names the two atoms of the contact
for era in ('2014on', 'pre2014'):
    for z in ('below', 'inside', 'between', 'past'):
        g = [r for r in X if r['dep'] is not None and ((r['dep'] >= 2014) == (era == '2014on')) and
             any(r['prot_res'] in G[2] and r['pair'] == G[3] for G in GROUPS) and zone(r['d'], r['pair']) == z]
        put(f'A.total_atoms.{era}.{z}', [sum(r['rec_atoms'] for r in g), sum(r['rec'] for r in g), len(g)])
s = [r for r in X if r['prot_res'] == 'LYS' and r['pair'] == 'C-N' and r['dep'] is not None and r['dep'] >= 2014 and 1.60 < r['d'] <= 2.07 and r['rec']]
put('lys.recorded_beyond_window_2014on', [len(s), sum(r['lig_comp'] == 'PLP' for r in s)])
# what the records beyond the window are, since 2014 (the seven groups, between and past zones)
bey = [r for r in X if r['dep'] is not None and r['dep'] >= 2014 and r['rec'] and
       any(r['prot_res'] in g[2] and r['pair'] == g[3] and zone(r['d'], g[3]) in ('between', 'past') for g in GROUPS)]
put('beyond.recorded_2014on', [len(bey), len({r['entry'] for r in bey})])
put('beyond.recorded_2014on.components', collections.Counter(r['lig_comp'] for r in bey).most_common(8))
# the 2.00 A switch for cysteine by deposition period: recorded / all at 1.90-2.00 and at 2.00-2.10 A
for lab, a, b in (('to2004', 0, 2004), ('2005-2013', 2005, 2013), ('2014on', 2014, 9999)):
    g = [r for r in X if r['prot_res'] == 'CYS' and r['pair'] == 'C-S' and r['dep'] is not None and a <= r['dep'] <= b]
    lo = [r for r in g if 1.90 <= r['d'] <= 2.00]; hi = [r for r in g if 2.00 < r['d'] <= 2.10]
    put(f'switch.cys_CS.{lab}', [sum(r['rec'] for r in lo), len(lo), sum(r['rec'] for r in hi), len(hi)])

# within-entry groups: same entry, residue type, nucleophile atom, ligand component, ligand atom
def groups(rows_):
    g = collections.defaultdict(list)
    for r in rows_:
        if r['pair'] in BONDPAIRS:
            g[(r['entry'], r['prot_res'], r['prot_atom'], r['lig_comp'], r['lig_atom'], r['pair'])].append(r)
    return g
mixed = []; cnt = collections.Counter()
for k, v in groups(X).items():
    rec = [r['d'] for r in v if r['rec']]; un = [r['d'] for r in v if not r['rec']]
    if not rec or not un: continue
    lo, hi = MAXIT[k[5]]
    ins = lambda d: lo <= d <= hi
    sep = all(ins(d) for d in rec) and all(not ins(d) for d in un)
    cnt['mixed'] += 1; cnt['separated'] += sep; cnt[f'mixed.{k[5]}'] += 1; cnt[f'separated.{k[5]}'] += sep
    mixed.append({'entry': k[0], 'residue': k[1], 'nucleophile_atom': k[2], 'ligand_comp': k[3], 'ligand_atom': k[4], 'pair': k[5],
                  'n_recorded': len(rec), 'n_unrecorded': len(un), 'recorded_min': min(rec), 'recorded_max': max(rec),
                  'unrecorded_min': min(un), 'unrecorded_max': max(un), 'window_separates_exactly': sep})
put('within.mixed_groups', cnt['mixed']); put('within.separated_exactly', cnt['separated'])
dep_of = {}
for r in X: dep_of[r['entry']] = r['dep']
for era in ('pre2014', '2014on'):
    g = [m for m in mixed if dep_of.get(m['entry']) is not None and ((dep_of[m['entry']] >= 2014) == (era == '2014on'))]
    put(f'within.{era}', [sum(m['window_separates_exactly'] for m in g), len(g)])
for p in BONDPAIRS:
    if cnt[f'mixed.{p}']: put(f'within.{p}', [cnt[f'separated.{p}'], cnt[f'mixed.{p}']])
dump('within_entry_mixed_groups.csv', sorted(mixed, key=lambda m: (m['pair'], m['entry'])))
# copies of the same bond within 0.03 A either side of the upper limit, all methods
st = []; sc = collections.Counter()
for k, v in groups(A).items():
    hi = MAXIT[k[5]][1]
    below = [r for r in v if hi - 0.03 < r['d'] <= hi]; above = [r for r in v if hi < r['d'] <= hi + 0.03]
    if not below or not above: continue
    b = max(below, key=lambda r: r['d']); a = min(above, key=lambda r: r['d'])
    kind = ('split_as_predicted' if b['rec'] and not a['rec'] else 'both_recorded' if b['rec'] and a['rec'] else
            'neither_recorded' if not b['rec'] and not a['rec'] else 'reversed')
    sc[kind] += 1; sc['all'] += 1; sc[f'{kind}.{k[5]}'] += 1; sc[f'all.{k[5]}'] += 1
    st.append({'entry': k[0], 'residue': k[1], 'nucleophile_atom': k[2], 'ligand_comp': k[3], 'ligand_atom': k[4], 'pair': k[5],
               'upper_limit': hi, 'below_dist': b['d'], 'below_recorded': b['rec'], 'below_site': b['prot_chain'] + b['prot_seq'],
               'above_dist': a['d'], 'above_recorded': a['rec'], 'above_site': a['prot_chain'] + a['prot_seq'], 'kind': kind,
               'above_listed_as_close_contact': a['close'] if a['hdr'] else '',
               'method': b['method1'], 'resolution': b['resolution'], 'dep_year': b['dep_year']})
for kd in ('all', 'split_as_predicted', 'both_recorded', 'neither_recorded', 'reversed'): put(f'straddle.{kd}', sc[kd])
for p in BONDPAIRS:
    if sc[f'all.{p}']: put(f'straddle.{p}', [sc[f'split_as_predicted.{p}'], sc[f'both_recorded.{p}'], sc[f'neither_recorded.{p}'], sc[f'reversed.{p}'], sc[f'all.{p}']])
sp = [m for m in st if m['kind'] == 'split_as_predicted' and m['above_listed_as_close_contact'] != '']
put('straddle.split.above_listed_as_close_contact', [sum(1 for m in sp if m['above_listed_as_close_contact']), len(sp)])
dump('straddling_copies.csv', sorted(st, key=lambda m: (m['pair'], m['above_dist'] - m['below_dist'])))
# between the lines, no record: how many, archive-wide
def between_unrec(rows_):
    s = [r for r in rows_ if r['pair'] in BONDPAIRS and zone(r['d'], r['pair']) == 'between' and not r['rec']]
    return [len(s), len({r['entry'] for r in s})]
put('A.between_unrecorded.xray25', between_unrec(X)); put('A.between_unrecorded.all_methods', between_unrec(A))
# lower limit, cysteine
for lab, rows_ in (('xray25', X), ('all', A)):
    g = [r for r in rows_ if r['prot_res'] == 'CYS' and r['pair'] == 'C-S' and r['d'] < 1.55]
    put(f'lower.cys_CS.{lab}', [sum(r['rec'] for r in g), len(g), len({r['entry'] for r in g}), sum(r['entry'] == '7EZX' for r in g)])
    put(f'lower.cys_CS.{lab}.7EZX', [sum(r['rec'] for r in g if r['entry'] == '7EZX'), sum(r['entry'] == '7EZX' for r in g)])
    put(f'lower.cys_CS.{lab}.record_on_these_atoms', [sum(r['rec_atoms'] for r in g), len(g)])
    put(f'lower.cys_CS.{lab}.other_entries', [sum(r['rec'] for r in g if r['entry'] != '7EZX'), sum(r['entry'] != '7EZX' for r in g)])
    g = [r for r in rows_ if r['prot_res'] == 'CYS' and r['pair'] == 'C-S' and 1.55 <= r['d'] <= 2.00]
    put(f'lower.cys_CS_inside.{lab}', [sum(r['rec'] for r in g), len(g)])
    put(f'lower.cys_CS_inside.{lab}.record_on_these_atoms', [sum(r['rec_atoms'] for r in g), len(g)])

# the archive's close-contact tables: where the unrecorded pairs are listed
def ledger(rows_, tag):
    b = [r for r in rows_ if r['hdr'] and r['pair'] in BONDPAIRS]
    ra = [r for r in b if r['rec_atoms']]                       # a record on the two atoms themselves
    put(f'{tag}.record_on_these_atoms', [sum(r['close'] for r in ra), len(ra), len({r['entry'] for r in ra})])
    ro = [r for r in b if r['rec'] and not r['rec_atoms']]      # the residues are joined, but by other atoms
    put(f'{tag}.record_on_other_atoms', [sum(r['close'] for r in ro), len(ro), len({r['entry'] for r in ro})])
    un = [r for r in b if not r['rec']]
    lo = [r for r in un if r['d'] < CLOSE_CUT]; hi = [r for r in un if r['d'] >= CLOSE_CUT]
    listed = [r for r in lo if r['close']]
    put(f'{tag}.unrecorded_under_cut', {'pairs': len(lo), 'entries': len({r['entry'] for r in lo}), 'close_contact': len(listed),
                                        'close_contact_entries': len({r['entry'] for r in listed}),
                                        'pct': round(100.0 * len(listed) / len(lo), 1) if lo else None,
                                        'symmetry_contact_table': sum(r['symm'] and not r['close'] for r in lo),
                                        'close_contact_on_other_atoms': sum(r['close_other'] and not r['symm'] and not r['close'] for r in lo)})
    put(f'{tag}.unrecorded_at_or_over_cut', [sum(r['close'] for r in hi), len(hi), sum(r['close'] for r in hi if r['d'] >= 2.2)])
    for z in ('below', 'inside', 'between', 'past', 'far'):
        s = [r for r in lo if zone(r['d'], r['pair']) == z]
        put(f'{tag}.unrecorded_under_cut.{z}', [sum(r['close'] for r in s), len(s)])
    for a, c in ((2.15, 2.19), (2.19, 2.195), (2.195, 2.20), (2.20, 2.25)):
        s = [r for r in un if a <= r['d'] < c]
        put(f'{tag}.step.{a}-{c}', [sum(r['close'] for r in s), len(s)])
ledger(A, 'A.close.all_methods'); ledger(X, 'A.close.xray25')
ledger([r for r in X if r['dep'] is not None and r['dep'] >= 2014], 'A.close.xray25.2014on')
ledger([r for r in X if r['dep'] is not None and r['dep'] < 2014], 'A.close.xray25.pre2014')

# --- evidence sample -------------------------------------------------------------------------------
S = list(csv.DictReader(op(f'{D}/sample_contacts.csv.gz')))
sel = json.load(open(f'{D}/sample_selection.json'))
put('sample.rows', len(S)); put('sample.entries_with_contact', len({r['pdbid'] for r in S}))
def full(r): return r['occupancy_class'] == 'full' and float(r['prot_occ'] or 0) > 0 and float(r['lig_occ'] or 0) > 0
# the files' own close-contact rows, keyed by entry and atom pair
CC = collections.defaultdict(set); ncc = 0
for c in csv.DictReader(op(f'{D}/sample_close_contact_rows.csv.gz')):
    a = (c['chain_1'], c['comp_1'], c['seq_1'], c['atom_1']); b = (c['chain_2'], c['comp_2'], c['seq_2'], c['atom_2'])
    CC[c['pdbid']].add(frozenset([a, b])); ncc += 1
put('sample.close_contact_rows', ncc)
def atoms(r):
    pch, pres, patom = r['protein_atom'].split(':')
    m = re.match(r'([A-Z0-9]{3})(-?\d+)', pres)
    lname, lch, lnum = r['ligand'].rsplit('_', 2)
    return frozenset([(pch, m.group(1), m.group(2), patom), (lch, lname, re.match(r'-?\d+', lnum).group(0), r['ligand_atom'])])
def skey(r):
    pch, pres, patom = r['protein_atom'].split(':')
    comp, lch, lseq = r['ligand'].rsplit('_', 2)
    return (r['pdbid'], pch, r['residue'], pres[len(r['residue']):], patom, lch, comp, re.match(r'-?\d+', lseq).group(0), r['ligand_atom'], r['prot_altloc'], r['lig_altloc'])
def prep(r):
    r['d'] = BXS.get(skey(r), float(r['dist_A'])); r['rec'] = r['conn_pair'] in ('Covale', 'Disulf')
    r['cp'] = canon(r['lig_el'], r['prot_el']); r['hb'] = r['hbond_capable_pair'] == 'True'
    r['close'] = atoms(r) in CC[r['pdbid']]
    r['dep'] = int(r['deposited'][:4]) if r['deposited'][:4].isdigit() else None
for r in S: prep(r)
SF = [r for r in S if full(r)]
put('sample.full_rows', len(SF)); put('sample.full_entries', len({r['pdbid'] for r in SF}))
put('sample.recorded_as_disulf', sum(1 for r in S if r['conn_pair'] == 'Disulf'))
# the close-contact table against the record, every full-occupancy contact on a pair in the link table
B = [r for r in SF if r['cp'] in MAXIT and not r['hb']]
rec = [r for r in B if r['rec']]
put('S.close.recorded', [sum(r['close'] for r in rec), len(rec), len({r['pdbid'] for r in rec})])
SREC = rec
un = [r for r in B if not r['rec']]
lo = [r for r in un if r['d'] < CLOSE_CUT]; hi = [r for r in un if r['d'] >= CLOSE_CUT]
put('S.close.unrecorded_under_cut', [sum(r['close'] for r in lo), len(lo), len({r['pdbid'] for r in lo})])
put('S.close.unrecorded_at_or_over_cut', [sum(r['close'] for r in hi), len(hi)])
t3 = []
for z in ('below', 'inside', 'between', 'past', 'far'):
    a = [r for r in rec if zone(r['d'], r['cp']) == z]; u = [r for r in un if zone(r['d'], r['cp']) == z]
    ul = [r for r in u if r['d'] < CLOSE_CUT]
    put(f'S.close.zone.{z}', [sum(r['close'] for r in a), len(a), sum(r['close'] for r in u), len(u), sum(r['close'] for r in ul), len(ul)])
    t3.append({'zone': z, 'recorded_listed': sum(r['close'] for r in a), 'recorded_n': len(a), 'unrecorded_listed': sum(r['close'] for r in u), 'unrecorded_n': len(u),
               'unrecorded_under_2.195_listed': sum(r['close'] for r in ul), 'unrecorded_under_2.195_n': len(ul)})
dump('close_contact_by_zone_sample.csv', t3)
AC = N['A.close.all_methods.record_on_these_atoms']; AU = N['A.close.all_methods.unrecorded_under_cut']; AO = N['A.close.all_methods.unrecorded_at_or_over_cut']
SU = N['S.close.unrecorded_under_cut']; SO = N['S.close.unrecorded_at_or_over_cut']
dump('table3_close_contact_listing.csv', [
    {'set': 'archive scan, all methods', 'pairs': 'a record joins the two atoms', 'n': AC[1], 'listed_as_close_contact': AC[0], 'notes': ''},
    {'set': 'archive scan, all methods', 'pairs': 'no record, under 2.195 A', 'n': AU['pairs'], 'listed_as_close_contact': AU['close_contact'],
     'notes': f"{AU['symmetry_contact_table']} in the symmetry-contact table; the rest have both atoms below full occupancy (data/unlisted_pairs_check.csv)"},
    {'set': 'archive scan, all methods', 'pairs': 'no record, 2.195 A or more', 'n': AO[1], 'listed_as_close_contact': AO[0], 'notes': ''},
    {'set': 'evidence sample, full occupancy', 'pairs': 'a record joins the two atoms', 'n': len(SREC), 'listed_as_close_contact': sum(r['close'] for r in SREC), 'notes': ''},
    {'set': 'evidence sample, full occupancy', 'pairs': 'no record, under 2.195 A', 'n': SU[1], 'listed_as_close_contact': SU[0], 'notes': ''},
    {'set': 'evidence sample, full occupancy', 'pairs': 'no record, 2.195 A or more', 'n': SO[1], 'listed_as_close_contact': SO[0], 'notes': ''}])
for a, c in ((2.15, 2.19), (2.19, 2.195), (2.195, 2.20), (2.20, 2.25)):
    s = [r for r in un if a <= r['d'] < c]
    put(f'S.close.step.{a}-{c}', [sum(r['close'] for r in s), len(s)])
hbp = [r for r in SF if r['cp'] in MAXIT and r['hb'] and not r['rec']]
put('S.close.hbond_pairs_unrecorded_under_cut', [sum(r['close'] for r in hbp if r['d'] < CLOSE_CUT), sum(1 for r in hbp if r['d'] < CLOSE_CUT)])
put('S.close.hbond_pairs_unrecorded_at_or_over_cut', [sum(r['close'] for r in hbp if r['d'] >= CLOSE_CUT), sum(1 for r in hbp if r['d'] >= CLOSE_CUT)])
# boron, serine and threonine
bo = [r for r in SF if r['cp'] == 'B-O' and r['residue'] in ('SER', 'THR')]
put('boron.contacts', len(bo)); put('boron.entries', len({r['pdbid'] for r in bo}))
tb = []
for era in ('all', '2014on', 'pre2014'):
    g = [r for r in bo if era == 'all' or (r['dep'] is not None and ((r['dep'] >= 2014) == (era == '2014on')))]
    row = {'era': era}
    for z in ('inside', 'between', 'past'):
        s = [r for r in g if zone(r['d'], 'B-O') == z]
        ents = collections.defaultdict(list)
        for r in s: ents[r['pdbid']].append(r['rec'])
        k = sum(r['rec'] for r in s)
        row[f'{z}_recorded'] = k; row[f'{z}_n'] = len(s); row[f'{z}_entries'] = len(ents)
        row[f'{z}_entries_none_recorded'] = sum(1 for v in ents.values() if not any(v))
        put(f'boron.{era}.{z}', [k, len(s), len(ents), row[f'{z}_entries_none_recorded']])
    tb.append(row)
dump('boron_zones.csv', tb)
bb = [r for r in bo if zone(r['d'], 'B-O') == 'between' and r['rec'] and r['dep'] is not None and r['dep'] >= 2014]
put('boron.2014on.between.recorded_by_entry', collections.Counter(r['pdbid'] for r in bb).most_common(5))
# rows with validation evidence
E = list(csv.DictReader(op(f'{D}/sample_window_evidence.csv.gz')))
for r in E:
    prep(r); r['clash'] = r['clash_listed'] == 'True'; r['inrep'] = r['in_report'] == 'True'
    r['dens'] = r['density_class'] == 'density_supported'
EF = [r for r in E if full(r)]
put('evidence.rows', len(E)); put('evidence.entries', len({r['pdbid'] for r in E}))
# between the annotation limit and the validator's line
zb = [r for r in EF if r['cp'] in MAXIT and not r['hb'] and zone(r['d'], r['cp']) == 'between']
zbu = [r for r in zb if not r['rec']]
put('between.contacts', len(zb)); put('between.recorded', sum(r['rec'] for r in zb))
put('between.unrecorded', len(zbu)); put('between.unrecorded.entries', len({r['pdbid'] for r in zbu}))
put('between.unrecorded.in_report', sum(r['inrep'] for r in zbu)); put('between.unrecorded.density_supported', sum(r['dens'] for r in zbu))
put('between.unrecorded.clash', sum(r['clash'] for r in zbu)); put('between.unrecorded.close_contact', sum(r['close'] for r in zbu))
put('between.unrecorded.under_cut.close_contact', [sum(r['close'] for r in zbu if r['d'] < CLOSE_CUT), sum(1 for r in zbu if r['d'] < CLOSE_CUT)])
t4 = []; named = 0
for lab, residues, cp in (('Cys, C-S', ('CYS',), 'C-S'), ('Cys, S-S', ('CYS',), 'S-S'), ('Ser and Thr, C-O', ('SER', 'THR'), 'C-O'), ('Asp and Glu, C-O', ('ASP', 'GLU'), 'C-O'),
                          ('Lys and His, C-N', ('LYS', 'HIS'), 'C-N'), ('Ser and Thr, B-O', ('SER', 'THR'), 'B-O')):
    s = [r for r in zbu if r['residue'] in residues and r['cp'] == cp]; named += len(s)
    sr = [r for r in zb if r['rec'] and r['residue'] in residues and r['cp'] == cp]
    row = {'residue_and_pair': lab, 'band': f"{MAXIT[cp][1]:.2f}-{probe_line(cp):.2f}", 'recorded_in_band': len(sr), 'unrecorded': len(s), 'entries': len({r['pdbid'] for r in s}),
           'density_supported': sum(r['dens'] for r in s), 'clash': sum(r['clash'] for r in s), 'close_contact': sum(r['close'] for r in s)}
    t4.append(row); put(f'between.{lab}', [row['recorded_in_band'], row['unrecorded'], row['entries'], row['density_supported'], row['clash'], row['close_contact']])
put('between.other_pairs', len(zbu) - named)
t4.append({'residue_and_pair': 'All pairs', 'band': '', 'recorded_in_band': sum(r['rec'] for r in zb), 'unrecorded': len(zbu), 'entries': len({r['pdbid'] for r in zbu}),
           'density_supported': sum(r['dens'] for r in zbu), 'clash': sum(r['clash'] for r in zbu), 'close_contact': sum(r['close'] for r in zbu)})
dump('table4_between_the_lines.csv', t4)
# past the validator's line, to +0.40, no record: clash in the report, close contact in the file
def flip(rows_):
    s = [r for r in rows_ if r['inrep']]
    return [sum(r['clash'] for r in s), len(s), sum(r['close'] for r in s), sum(r['close'] for r in s if r['d'] < CLOSE_CUT), sum(1 for r in s if r['d'] < CLOSE_CUT)]
pc = [r for r in EF if r['cp'] in MAXIT and not r['rec'] and zone(r['d'], r['cp']) == 'past']
t5 = []
for lab, sub in (('inside the validator line, all pairs', [r for r in EF if r['cp'] in MAXIT and not r['hb'] and not r['rec'] and r['d'] <= probe_line(r['cp'])]),
                 ('past the line: cysteine', [r for r in pc if r['residue'] == 'CYS']),
                 ('past the line: other residues, pairs without boron', [r for r in pc if r['residue'] != 'CYS' and not r['hb'] and 'B' not in r['cp'].split('-')]),
                 ('past the line: serine to boron', [r for r in pc if r['residue'] == 'SER' and r['cp'] == 'B-O'])):
    v = flip(sub); put(f'clash.{lab}', v)
    t5.append({'unrecorded_contacts': lab, 'clash': v[0], 'n_in_report': v[1], 'close_contact': v[2], 'close_contact_under_2.195': v[3], 'n_under_2.195': v[4]})
dump('table5_clash_listing.csv', t5)
pr = [r for r in EF if r['cp'] in MAXIT and not r['hb'] and r['rec'] and zone(r['d'], r['cp']) == 'past' and r['inrep']]
put('clash.past the line: recorded', [sum(r['clash'] for r in pr), len(pr)])
t5.append({'unrecorded_contacts': 'recorded contacts past the line (for comparison)', 'clash': sum(r['clash'] for r in pr), 'n_in_report': len(pr),
           'close_contact': sum(r['close'] for r in pr), 'close_contact_under_2.195': sum(r['close'] for r in pr if r['d'] < CLOSE_CUT), 'n_under_2.195': sum(1 for r in pr if r['d'] < CLOSE_CUT)})
dump('table5_clash_listing.csv', t5)
# the triage set (past_line_triage.csv) is defined from the validator's line alone: line < d <= line + 0.40, residues other than cysteine.
# (For O-O the validator's line, 1.52 A, lies below the link table's upper limit, 1.60 A, so the window-first zone() would call three of them inside.)
tri = [r for r in EF if r['cp'] in MAXIT and not r['rec'] and r['residue'] != 'CYS' and probe_line(r['cp']) < r['d'] <= probe_line(r['cp']) + 0.40]
put('triage.selection', [len(tri), len({r['pdbid'] for r in tri})])
T = list(csv.DictReader(open(f'{D}/past_line_triage.csv', newline='')))
tk = {(r['pdbid'], r['protein_atom'], r['ligand'], r['ligand_atom']) for r in T}
put('triage.table_rows_matching_selection', [len(tk & {(r['pdbid'], r['protein_atom'], r['ligand'], r['ligand_atom']) for r in tri}), len(T)])
put('triage.classes', collections.Counter(r['triage'] for r in T).most_common())
put('triage.candidate_serine_boron', sum(1 for r in T if r['triage'] == 'candidate bond' and r['pair'] == 'B-O' and ':SER' in r['protein_atom']))
put('clash.past.boron_other_residues', len([r for r in pc if r['residue'] != 'SER' and 'B' in r['cp'].split('-') and r['inrep']]))
sb = [r for r in pc if r['cp'] == 'B-O' and r['residue'] == 'SER']
put('boron.past.unrecorded', len(sb)); put('boron.past.entries', len({r['pdbid'] for r in sb})); put('boron.past.density_supported', sum(r['dens'] for r in sb))
put('boron.past.deposited', [min(r['deposited'][:4] for r in sb), max(r['deposited'][:4] for r in sb)])
# neighbour clashes
BN = list(csv.DictReader(open(f'{D}/boron_neighbour_clashes.csv')))
for z, lab in (('past_probe_line_to_2.12', 'past'), ('between_limit_and_probe_line', 'between')):
    s = [r for r in BN if r['zone'] == z]
    hit = [[p for p in r['clash_partners'].split('; ') if re.match(r'^(OG1?|HG1?)~', p)] for r in s]
    put(f'boron.{lab}.neigh', [len(s), sum(1 for h in hit if h), sum(len(h) for h in hit), sum(int(r['clash_records_naming_the_boron']) for r in s)])
    el = collections.Counter(re.search(r'\((\w+)\)', p).group(1) for h in hit for p in h)
    put(f'boron.{lab}.neigh_elements', dict(sorted(el.items())))
bc = list(csv.DictReader(open(f'{D}/boron_clash_check.csv')))
put('boron.table.entries', len({r['pdbid'] for r in bc}))
json.dump(N, open(f'{O}/numbers.json', 'w'), indent=1)
for k, v in N.items(): print(k, '=', v)
