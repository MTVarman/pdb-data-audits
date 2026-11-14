#!/usr/bin/env python3
"""Recompute every headline number in README.md from the shipped tables and print PASS/FAIL.

Record: Varman, M. T. (2026). Unrecorded covalent bonds to cysteine concentrate by wwPDB processing
site: a 24,144-entry scan of _struct_conn records, with two matched cases. v1.0.

Usage: python3 verify_counts.py   (keep it at the top level of the unzipped archive, next to out/)
Standard library only.
"""
import csv, os, sys
from math import lgamma, exp
HERE = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(HERE, 'out')
csv.field_size_limit(10**8)
def rows(name):
    p = os.path.join(O, name)
    if not os.path.exists(p): sys.exit(f'missing {p}: keep verify_counts.py next to out/')
    with open(p, newline='') as fh: return list(csv.DictReader(fh, delimiter='\t'))
res = []
def check(label, got, exp_):
    ok = got == exp_; res.append(ok)
    print(f'{"PASS" if ok else "FAIL"}  {label}: got {got}, expected {exp_}')
def B(x): return x == 'True'
def fr(xs): return (sum(xs), len(xs))
def _lp(a, b, c, d):
    return (lgamma(a+b+1)+lgamma(c+d+1)+lgamma(a+c+1)+lgamma(b+d+1)-lgamma(a+b+c+d+1)-lgamma(a+1)-lgamma(b+1)-lgamma(c+1)-lgamma(d+1))
def fisher(a, b, c, d):
    r1, c1, n = a+b, a+c, a+b+c+d; p0 = _lp(a, b, c, d); p = 0.0
    for x in range(max(0, r1+c1-n), min(r1, c1)+1):
        lp = _lp(x, r1-x, c1-x, n-r1-c1+x)
        if lp <= p0 + 1e-9: p += exp(lp)
    return (a*d)/(b*c), min(p, 1.0)
def sci(p): return f'{p:.0e}'

print('-- sample --')
cand = rows('candidates.tsv'); ent = rows('entries_scanned.tsv')
check('candidate entries', len(cand), 24144)
check('arm A "covalent cysteine"', sum('A' in r['arms'] for r in cand), 12563)
check('arm B "covalent inhibitor"', sum('B' in r['arms'] for r in cand), 16281)
check('entries scanned with status ok', sum(r['status'] == 'ok' for r in ent), 24144)
check('entries with any SG contact <= 2.6 A', sum(int(r['n_pairs']) > 0 for r in ent), 3909)
check('entries scanned without symmetry (3LDH)', [r['entry'] for r in ent if r['symmetry_checked'] == 'False'], ['3LDH'])

A = rows('contacts_all.tsv')
for r in A:
    r['d'] = float(r['dist']); r['rec'] = B(r['recorded']); r['site'] = r['process_site'] or '(blank)'
    try: r['yr'] = int(r['deposit_year'])
    except ValueError: r['yr'] = None
base = [r for r in A if B(r['nearest_lig_for_cys']) and 1.6 <= r['d'] <= 2.4]
def pop(el, cats): return [r for r in base if r['el_class'] == el and r['category'] in cats]
def bands(x):
    return (fr([r['rec'] for r in x if r['d'] < 2.0]), fr([r['rec'] for r in x if 2.0 <= r['d'] < 2.1]), fr([r['rec'] for r in x if r['d'] >= 2.1]))
DL = pop('C', {'other_ligand', 'named_drug'})

print('\n-- main table: drug-like ligands, carbon partner --')
check('drug-like contacts', len(DL), 4249)
check('drug-like entries', len({r['entry'] for r in DL}), 2516)
b = bands(DL)
check('recorded 1.60-<2.00', b[0], (3897, 4094)); check('recorded 2.00-<2.10', b[1], (27, 69)); check('recorded 2.10-2.40', b[2], (26, 86))
x = fr([r['rec'] for r in DL if 1.95 <= r['d'] < 2.0]); y = fr([r['rec'] for r in DL if 2.0 <= r['d'] < 2.05])
check('1.95-2.00 A recorded', x, (50, 62)); check('2.00-2.05 A recorded', y, (20, 44))
check('step Fisher p (1.95-2.00 vs 2.00-2.05)', sci(fisher(x[0], x[1]-x[0], y[0], y[1]-y[0])[1]), '3e-04')
lg = fr([r['rec'] for r in DL if r['d'] >= 2.0]); check('recorded >= 2.0 A', lg, (53, 155))
old = fr([r['rec'] for r in DL if r['d'] >= 2.0 and r['yr'] is not None and r['yr'] < 2010])
new = fr([r['rec'] for r in DL if r['d'] >= 2.0 and r['yr'] is not None and r['yr'] >= 2020])
check('>= 2.0 A recorded, deposited before 2010', old, (23, 40)); check('>= 2.0 A recorded, deposited 2020-2026', new, (16, 66))
check('unrecorded and bonded elsewhere', sum(B(r['bonded_elsewhere']) for r in DL), 14)
# sensitivity: max change in any band percentage
def pct(t): return 100 * t[0] / t[1]
full = [pct(t) for t in b]; worst = 0.0
for keep in (lambda r: not B(r['bonded_elsewhere']), lambda r: 'X-RAY' in r['method'],
             lambda r: float(r['lig_occ']) >= 1 and float(r['cys_occ']) >= 1, lambda r: r['cys_alt'] == '' and r['lig_alt'] == '',
             lambda r: r['sym_image'] == '0'):
    s = bands([r for r in DL if keep(r)]); worst = max(worst, max(abs(pct(s[i]) - full[i]) for i in range(3)))
check('largest sensitivity shift <= 4.0 points', round(worst, 1) <= 4.0, True)

print('\n-- processing site --')
lo = [r for r in DL if r['d'] < 2.0]; hi = [r for r in DL if r['d'] >= 2.0]
exp_lo = {'RCSB': (2173, 2225), 'PDBE': (1039, 1068), 'PDBC': (233, 251), 'PDBJ': (429, 526), 'BNL': (23, 24)}
exp_hi = {'RCSB': (40, 101), 'PDBE': (10, 30), 'PDBC': (1, 2), 'PDBJ': (2, 22)}
for s, v in exp_lo.items(): check(f'{s} < 2.0 A recorded', fr([r['rec'] for r in lo if r['site'] == s]), v)
for s, v in exp_hi.items(): check(f'{s} 2.0-2.4 A recorded', fr([r['rec'] for r in hi if r['site'] == s]), v)
miss = [r for r in lo if not r['rec']]
check('unrecorded < 2.0 A (all sites)', len(miss), 197); check('of which PDBj', sum(r['site'] == 'PDBJ' for r in miss), 97)
check('PDBj share of < 2.0 A contacts (%)', round(100 * sum(r['site'] == 'PDBJ' for r in lo) / len(lo)), 13)
pj = [r for r in lo if r['site'] == 'PDBJ']; ot = [r for r in lo if r['site'] != 'PDBJ']
check('other sites < 2.0 A recorded', fr([r['rec'] for r in ot]), (3468, 3568))
orr, p = fisher(len(pj) - sum(r['rec'] for r in pj), sum(r['rec'] for r in pj), len(ot) - sum(r['rec'] for r in ot), sum(r['rec'] for r in ot))
check('contact-level odds ratio PDBj vs others', round(orr, 1), 7.8); check('contact-level Fisher p', sci(p), '1e-37')
E = {}
for r in lo: E.setdefault((r['entry'], r['site']), []).append(r['rec'])
ej = [not all(v) for k, v in E.items() if k[1] == 'PDBJ']; eo = [not all(v) for k, v in E.items() if k[1] != 'PDBJ']
check('PDBj entries with >= 1 miss', (sum(ej), len(ej)), (47, 324)); check('other-site entries with >= 1 miss', (sum(eo), len(eo)), (58, 2118))
exp_ent = {'RCSB': (28, 1325), 'PDBE': (18, 631), 'PDBC': (11, 140), 'BNL': (1, 22)}
for s, v in exp_ent.items():
    z = [not all(vv) for k, vv in E.items() if k[1] == s]; check(f'{s} entries with >= 1 miss', (sum(z), len(z)), v)
orr, p = fisher(sum(ej), len(ej) - sum(ej), sum(eo), len(eo) - sum(eo))
check('entry-level odds ratio', round(orr, 1), 6.0); check('entry-level Fisher p', sci(p), '4e-16')
check('PDBj < 2.0 A recorded, deposited 2020-2022', fr([r['rec'] for r in pj if r['yr'] and 2020 <= r['yr'] <= 2022]), (153, 194))
check('PDBj < 2.0 A recorded, deposited 2023-2026', fr([r['rec'] for r in pj if r['yr'] and r['yr'] >= 2023]), (122, 162))

print('\n-- nirmatrelvir (4WI) and Cys12 --')
w = [r for r in lo if r['lig_comp'] == '4WI']
for s, v in {'RCSB': (49, 49), 'PDBC': (23, 23), 'PDBE': (2, 2), 'PDBJ': (7, 30)}.items():
    check(f'4WI < 2.0 A recorded at {s}', fr([r['rec'] for r in w if r['site'] == s]), v)
c12 = [r for r in lo if r['cys_seq'] == '12']
check('Cys12 < 2.0 A recorded at PDBj', fr([r['rec'] for r in c12 if r['site'] == 'PDBJ']), (5, 15))
check('Cys12 < 2.0 A recorded elsewhere', fr([r['rec'] for r in c12 if r['site'] != 'PDBJ']), (163, 166))

print('\n-- named drugs (KRAS G12C, BTK, EGFR/HER: all under 2.0 A and recorded) --')
ND = [r for r in DL if r['category'] == 'named_drug']
for cl in ('KRAS G12C', 'BTK', 'EGFR/HER'):
    z = [r for r in ND if r['drug_class'] == cl and r['d'] < 2.0]
    check(f'{cl} contacts < 2.0 A recorded', all(r['rec'] for r in z) and len(z) > 0, True)
check('only named-drug contact >= 2.0 A outside Mpro', [(r['entry'], r['lig_comp'], r['dist'], r['rec']) for r in ND if r['d'] >= 2.0 and r['drug_class'] != 'SARS-CoV-2 Mpro'], [('7JXM', 'YY3', '2.1', False)])

print('\n-- controls in the same sample --')
H = bands(pop('C', {'heme_porphyrin'})); check('heme/porphyrin < 2.0 A recorded', H[0], (3019, 3019))
check('heme/porphyrin >= 2.0 A recorded', (H[1][0] + H[2][0], H[1][1] + H[2][1]), (0, 291))
bil = [r for r in pop('C', {'bilin'}) if r['d'] < 2.0 and not r['rec']]
check('bilin < 2.0 A misses', len(bil), 1454); check('bilin < 2.0 A misses at PDBj', sum(r['site'] == 'PDBJ' for r in bil), 1423)
F = bands(pop('C', {'flavin'})); check('flavin < 2.0 A recorded', F[0], (114, 120))
check('flavin >= 2.0 A recorded', (F[1][0] + F[2][0], F[1][1] + F[2][1]), (2, 11))
S = [r for r in base if r['el_class'] == 'S']
for lo_, hi_, v in ((1.95, 2.0, (58, 75)), (2.0, 2.05, (204, 249)), (2.05, 2.1, (102, 124)), (2.1, 2.2, (62, 93)), (2.2, 2.4000001, (17, 146))):
    check(f'ligand S {lo_:.2f}-{min(hi_, 2.4):.2f} A recorded', fr([r['rec'] for r in S if lo_ <= r['d'] < hi_]), v)

print('\n-- overlap with strict-cut-491 --')
U = rows('unrecorded_under_2.0A_druglike.tsv')
check('unrecorded < 2.0 A rows / entries', (len(U), len({r['entry'] for r in U})), (197, 105))
check('already in strict-cut-491', sum(r.get('in_strict_cut_491') == 'True' for r in U), 131)

print('\n-- example rows quoted in README --')
def find(e, ch, comp):
    return [(r['dist'], r['rec']) for r in DL if r['entry'] == e and r['cys_chain'] == ch and r['lig_comp'] == comp]
for e, ch, comp, d, rec in [('6OIM', 'A', 'MOV', '1.805', True), ('5P9J', 'A', '8E8', '1.85', True), ('6JX4', 'A', 'YY3', '1.811', True),
                            ('8FD9', 'A', 'XQQ', '1.683', True), ('4G5J', 'A', '0WN', '1.841', True), ('7YCC', 'A', 'IQC', '1.818', False),
                            ('6JZ0', 'A', 'CKO', '1.932', False), ('7VJL', 'B', '7IF', '1.981', True), ('7VJL', 'A', '7IF', '2.0', False),
                            ('9PZF', 'A', 'A1CQ2', '1.943', True), ('9PZF', 'B', 'A1CQ2', '2.003', False), ('7MAT', 'B', 'H37', '1.766', True),
                            ('7MAT', 'A', 'H37', '2.006', False), ('8TWY', 'A', 'RNX', '2.0', True), ('8ANF', 'A', 'N1R', '2.0', True)]:
    check(f'{e} {ch} {comp}', find(e, ch, comp), [(d, rec)])
check('7YCC IQC unrecorded chains', sum(1 for r in DL if r['entry'] == '7YCC' and r['lig_comp'] == 'IQC' and not r['rec'] and r['d'] < 2), 3)
z = [r['d'] for r in DL if r['entry'] == '8X6R' and r['lig_comp'] == 'Y9D' and not r['rec']]
check('8X6R Y9D unrecorded chains and range', (len(z), min(z), max(z)), (6, 1.784, 1.81))
z = [r['d'] for r in DL if r['entry'] in ('7WCT', '7WCW', '7WCX') and r['lig_comp'] == '90F' and not r['rec']]
check('7WCT/7WCW/7WCX 90F unrecorded, range', (len(z), min(z), max(z)), (3, 1.627, 1.651))
z = [r['d'] for r in lo if r['lig_comp'] == '4WI' and r['site'] == 'PDBJ' and not r['rec']]
check('4WI PDBj misses distance range', (len(z), min(z), max(z)), (23, 1.761, 1.846))



print('\n-- deposition year of the unrecorded contacts --')
check('PDBj unrecorded, deposited 2021-2023', fr([bool(r['yr'] and 2021 <= r['yr'] <= 2023) for r in lo if r['site'] == 'PDBJ' and not r['rec']]), (76, 97))
check('RCSB unrecorded, deposited 2006-2013', fr([bool(r['yr'] and 2006 <= r['yr'] <= 2013) for r in lo if r['site'] == 'RCSB' and not r['rec']]), (35, 52))

print('\n-- per depositing group (out/depositor_groups.tsv) --')
G = rows('depositor_groups.tsv'); gk = {r['entry']: r for r in G}
check('entries with a group key', len(G), 2442)
check('group table matches contact table', all(int(gk[e]['n_unrecorded_lt2']) == sum(1 for r in lo if r['entry'] == e and not r['rec']) for e in list(gk)[::25]), True)
for key, exp_, tot_ in (('citation_group', {'RCSB': (23, 546), 'PDBE': (15, 275), 'PDBC': (7, 60), 'PDBJ': (30, 182), 'BNL': (1, 14)}, ((30, 182), (46, 895), '3.6', '1e-06')),
                        ('last_author_group', {'RCSB': (23, 337), 'PDBE': (13, 179), 'PDBC': (5, 29), 'PDBJ': (25, 127), 'BNL': (1, 15)}, ((25, 127), (42, 560), '3.0', '1e-04'))):
    tot, hit = {}, {}
    for r in G:
        tot.setdefault(r['process_site'], set()).add(r[key])
        if int(r['n_unrecorded_lt2']) > 0: hit.setdefault(r['process_site'], set()).add(r[key])
    for s_, v in exp_.items(): check(f'{key} {s_} groups with a miss', (len(hit.get(s_, ())), len(tot[s_])), v)
    a = len(hit['PDBJ']); b_ = len(tot['PDBJ']) - a
    c = sum(len(v) for k, v in hit.items() if k != 'PDBJ'); d = sum(len(v) for k, v in tot.items() if k != 'PDBJ') - c
    orr, pv = fisher(a, b_, c, d)
    check(f'{key} PDBj vs others', ((a, a + b_), (c, c + d), f'{orr:.1f}', sci(pv)), tot_)

print('\n-- matched case 1: nirmatrelvir on SARS-CoV-2 Mpro Cys145 (out/nirmatrelvir_targets.tsv) --')
T = {r['entry']: r for r in rows('nirmatrelvir_targets.tsv')}
N = [r for r in lo if r['lig_comp'] == '4WI']
check('nirmatrelvir contacts / entries under 2.0 A', (len(N), len({r['entry'] for r in N})), (104, 64))
check('distinct source organisms', len({T[r['entry']]['source_organism'] for r in N}), 9)
S2 = [r for r in N if T[r['entry']]['source_organism'] == 'Severe acute respiratory syndrome coronavirus 2' and r['cys_seq'] == '145']
check('SARS-CoV-2 Cys145 contacts / entries', (len(S2), len({r['entry'] for r in S2})), (80, 56))
for s_, v in (('RCSB', (44, 44, 29, 14)), ('PDBC', (10, 10, 10, 2)), ('PDBE', (2, 2, 2, 1)), ('PDBJ', (6, 24, 15, 7))):
    g = [r for r in S2 if r['site'] == s_]
    check(f'SARS-CoV-2 Cys145 {s_} recorded / contacts / entries / papers',
          (sum(r['rec'] for r in g), len(g), len({r['entry'] for r in g}), len({gk[r['entry']]['citation_group'] for r in g})), v)
g = [r for r in S2 if r['site'] != 'PDBJ']
check('SARS-CoV-2 Cys145 outside PDBj: recorded / contacts / entries / papers',
      (sum(r['rec'] for r in g), len(g), len({r['entry'] for r in g}), len({gk[r['entry']]['citation_group'] for r in g})), (56, 56, 41, 17))
g = [r for r in S2 if r['site'] == 'PDBJ']
check('SARS-CoV-2 Cys145 PDBj papers: recorded-only / unrecorded',
      (len({gk[r['entry']]['citation_group'] for r in g if r['rec']}), len({gk[r['entry']]['citation_group'] for r in g if not r['rec']})), (3, 4))
z = [r['d'] for r in g if not r['rec']]
check('SARS-CoV-2 Cys145 PDBj unrecorded: n, range', (len(z), min(z), max(z)), (18, 1.761, 1.846))
for s_, v in (('RCSB', (49, 49)), ('PDBC', (23, 23)), ('PDBE', (2, 2)), ('PDBJ', (7, 30))):
    check(f'nirmatrelvir, all targets, {s_}', fr([r['rec'] for r in N if r['site'] == s_]), v)
g = [r for r in N if r['site'] == 'PDBJ' and not r['rec']]
check('nirmatrelvir PDBj unrecorded: contacts / entries / papers',
      (len(g), len({r['entry'] for r in g}), len({gk[r['entry']]['citation_group'] for r in g})), (23, 13, 5))
check('UQ0: 8I01/8I05/8I08 share a citation, 8BEO does not',
      (len({gk[e]['citation_group'] for e in ('8I01', '8I05', '8I08')}), gk['8BEO']['citation_group'] == gk['8I01']['citation_group']), (1, False))

print('\n-- composition of the PDBj excess --')
pj_mis = [r for r in lo if r['site'] == 'PDBJ' and not r['rec']]
check('PDBj unrecorded < 2.0 A', len(pj_mis), 97)
check('PDBj unrecorded: distinct components / entries',
      (len({r['lig_comp'] for r in pj_mis}), len({r['entry'] for r in pj_mis})), (31, 47))
check('PDBj unrecorded 4WI', sum(r['lig_comp'] == '4WI' for r in pj_mis), 23)
check('PDBj unrecorded UQ0', sum(r['lig_comp'] == 'UQ0' for r in pj_mis), 23)
ex = [r for r in lo if r['lig_comp'] not in ('4WI', 'UQ0')]
check('excl 4WI+UQ0, PDBj recorded < 2.0 A', fr([r['rec'] for r in ex if r['site'] == 'PDBJ']), (422, 473))
check('excl 4WI+UQ0, other sites recorded < 2.0 A', fr([r['rec'] for r in ex if r['site'] != 'PDBJ']), (3388, 3488))

print('\n-- matched case 2: UQ0 on glyoxylate carboligase --')
for e, site, n_, k_ in (('8BEO', 'PDBE', 6, 6), ('8I01', 'PDBJ', 6, 0), ('8I05', 'PDBJ', 6, 0), ('8I08', 'PDBJ', 11, 0)):
    g = [r for r in lo if r['entry'] == e and r['lig_comp'] == 'UQ0']
    check(f'{e} UQ0 contacts / recorded / site', (len(g), sum(r['rec'] for r in g), g[0]['site'] if g else None), (n_, k_, site))
z = [r['d'] for r in lo if r['entry'] == '8BEO' and r['lig_comp'] == 'UQ0']
check('8BEO UQ0 distance range', (round(min(z), 3), round(max(z), 3)), (1.664, 1.718))
z = [r['d'] for r in lo if r['entry'] in ('8I01', '8I05', '8I08') and r['lig_comp'] == 'UQ0']
check('8I01/8I05/8I08 UQ0 distance range', (round(min(z), 3), round(max(z), 3)), (1.765, 1.776))

print('\n-- cctbx generic C-S link value (1.766 +- 0.010 A) --')
ON = lambda r: abs(r['d'] - 1.766) <= 0.010
check('all < 2.0 A on value', (sum(ON(r) for r in lo), len(lo)), (938, 4094))
check('unrecorded < 2.0 A on value', (sum(ON(r) for r in lo if not r['rec']), sum(1 for r in lo if not r['rec'])), (65, 197))
check('recorded < 2.0 A on value', (sum(ON(r) for r in lo if r['rec']), sum(1 for r in lo if r['rec'])), (873, 3897))
for s_, v in (('RCSB', (2, 52)), ('PDBE', (4, 29)), ('PDBC', (3, 18)), ('PDBJ', (55, 97)), ('BNL', (1, 1))):
    g = [r for r in lo if r['site'] == s_ and not r['rec']]
    check(f'{s_} unrecorded on value', (sum(ON(r) for r in g), len(g)), v)
g = [r for r in lo if r['lig_comp'] == '4WI' and r['site'] == 'PDBJ' and not r['rec']]
check('4WI PDBj unrecorded on value', (sum(ON(r) for r in g), len(g)), (22, 23))
g = [r for r in lo if r['lig_comp'] == '4WI' and r['rec']]
check('4WI recorded on value', (sum(ON(r) for r in g), len(g)), (51, 81))
g = [r for r in lo if r['lig_comp'] == 'UQ0' and r['site'] == 'PDBJ' and not r['rec']]
check('UQ0 PDBj unrecorded on value', (sum(ON(r) for r in g), len(g)), (21, 23))
g = [r for r in lo if r['lig_comp'] == 'UQ0' and r['rec']]
check('UQ0 8BEO recorded on value', (sum(ON(r) for r in g), len(g)), (0, 6))
be = {}
for r in lo:
    if not r['rec']: be.setdefault(r['entry'], []).append(r)
check('entries whose every unrecorded contact is on value',
      (sum(1 for v in be.values() if all(ON(r) for r in v)), len(be)), (30, 105))

n = sum(res); print(f'\n{n}/{len(res)} checks passed'); sys.exit(0 if n == len(res) else 1)
