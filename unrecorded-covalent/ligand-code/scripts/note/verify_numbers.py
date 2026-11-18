"""Recompute the numbers quoted in the addendum from the bundled tables.
Run from anywhere:  python3 scripts/note/verify_numbers.py
Needs pandas and scipy. Prints one line per check and a final tally."""
import os, re, sys
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

R = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'data'))
A, L, MP, D = (os.path.join(R, x) for x in ('archive', 'ligand_code_effect', 'mpro_census', 'downstream'))
checks = []
def ck(label, got, want):
    ok = (got == want) if not isinstance(want, float) else abs(got - want) < 1e-9
    checks.append(ok); print(('PASS ' if ok else 'FAIL ') + f'{label}: got {got!r}, expected {want!r}')
def txt(p): return open(os.path.join(L, p)).read()

# ---- archive-wide (header files) ----
M = pd.read_csv(f'{A}/entries_meta.tsv', sep='\t', dtype=str, keep_default_na=False)
ck('entries parsed', len(M), 199326)
ck('entries with close-contact table', int((M.has_cc_table == 'True').sum()), 94342)
B = pd.read_csv(f'{A}/bonds_all.tsv', sep='\t', dtype={'entry': str}, low_memory=False)
ck('nucleophile-ligand bonds (any distance)', len(B), 45450)
ck('bonds < 2.0 A', int(B.lt2.sum()), 43488)
ck('unrecorded bonds < 2.0 A', int((B.lt2 & ~B.recorded).sum()), 5594)
lt = B[B.lt2]
for nuc, tot, un in [('CYS', 31396, 4362), ('SER', 4316, 516), ('LYS', 6636, 484), ('HIS', 1140, 232)]:
    s = lt[lt.nuc_res == nuc]; ck(f'{nuc} bonds < 2.0 A (total, unrecorded)', (len(s), int((~s.recorded).sum())), (tot, un))
for ps, tot, un in [('PDBJ', 17580, 3950), ('RCSB', 14308, 874), ('PDBE', 9633, 481), ('PDBC', 997, 268), ('BNL', 939, 21)]:
    s = lt[lt.ps == ps]; ck(f'{ps} bonds < 2.0 A (total, unrecorded)', (len(s), int((~s.recorded).sum())), (tot, un))
pe = lt[lt.lig_comp == 'PEB']
ck('PEB bonds < 2.0 A, unrecorded, entries', (len(pe), int((~pe.recorded).sum()), pe.entry.nunique()), (11022, 2859, 42))
U = set(lt[~lt.recorded].entry); Rr = set(B[B.recorded].entry); only = sorted(U - Rr)
ck('entries with >=1 unrecorded bond < 2.0 A', len(U), 1057)
ck('... of which no recorded nucleophile-ligand covale row', len(only), 832)
h = M.set_index('e').loc[only].hpm.value_counts().to_dict()
ck('has_protein_modification in those 832 (Y, N, absent)', (h.get('Y', 0), h.get('N', 0), h.get('', 0)), (338, 239, 255))
E = pd.read_csv(f'{L}/step3/entry_level_bonded.tsv', sep='\t', dtype={'entry': str}, low_memory=False)
ck('(code, entry) units, misses, codes', (len(E), int(E.miss.sum()), E.lig_comp.nunique()), (12199, 1099, 3864))
dl = E[E.cat == 'drug-like']; ck('drug-like units, misses', (len(dl), int(dl.miss.sum())), (6583, 571))
for site, n, m in [('PDBJ', 935, 199), ('RCSB', 3325, 194), ('PDBE', 1987, 123), ('PDBC', 185, 42), ('BNL', 151, 13)]:
    s = dl[dl.ps == site]; ck(f'drug-like units at {site} (units, misses)', (len(s), int(s.miss.sum())), (n, m))
S = pd.read_csv(f'{L}/step3/attribution_sensitivity.tsv', sep='\t').set_index('subset')
ck('born-unbonded OR, all drug-like', (S.loc['drug-like', 'unbonded_OR'], S.loc['drug-like', 'unbonded_CI']), (2.98, '1.85-4.80'))
ck('born-unbonded OR, specific + bond-like', (S.loc['drug-like, specific + bond-like', 'unbonded_OR'], S.loc['drug-like, specific + bond-like', 'unbonded_CI'], S.loc['drug-like, specific + bond-like', 'unbonded_p']), (1.47, '0.67-3.21', 0.33))
ck('born-unbonded later uses / misses (strict)', (int(S.loc['drug-like, specific + bond-like', 'unbonded_later_units']), int(S.loc['drug-like, specific + bond-like', 'unbonded_later_misses'])), (258, 22))
t20 = txt('step3/20_stdout.txt')
ck('first-bonded-use OR (strict)', bool(re.search(r'first_miss OR 14\.23 \(4\.28-47\.38\) p 1\.5e-05', t20)), True)
ck('different-lab OR (strict)', bool(re.search(r'different-lab later uses only: first_miss OR 8\.71 \(2\.27-33\.39\) p 0\.0016', t20)), True)
t18 = txt('step3/18_stdout.txt')
ck('lab+day same-drug pairs 507/486/15/5', bool(re.search(r"'pairs': 507, 'concordant_both_rec': np.int64\(486\), 'concordant_both_miss': np.int64\(15\), 'discordant_full': np.int64\(5\)", t18)), True)
ck('fixed / never / flip-flop among 29', "{'never recorded': 14, 'fixed': 11, 'flip-flop': 4}" in t18, True)

# ---- step 4 (strict set) ----
t22 = txt('step4/22_stdout.txt')
ck('strict set size', 'strict units 5355 misses 334 codes 3320 labs 1463' in t22, True)
for s in ['15/40 (38%) vs 86/1289 (7%) OR 8.5 (4.3-16.5), Fisher p=4.8e-08',
          '4/19 (21%) vs 58/693 (8%) OR 3.2 (1.1-9.3), Fisher p=0.075',
          '7/26 (27%) vs 36/982 (4%) OR 10.0 (4.0-24.6), Fisher p=5.6e-05',
          '4/15 (27%) vs 13/363 (4%) OR 10.2 (3.0-34.3), Fisher p=0.0028',
          'same code: 22/49 (45%)  labs=23', 'other code, same drug (InChIKey-14 or skeleton): 0/5 (0%)  labs=4',
          'other code, different drug: 68/599 (11%)  labs=72', 'baseline (later units of labs with earlier units, none missed): 87/2578 (3%)',
          'c_first_miss: OR 10.22 (2.74-38.13) p=0.00054', 'c_first_miss: OR 6.98 (1.75-27.78)', 'c_first_miss: OR 7.28 (1.87-28.37)',
          'c_first_miss: OR 12.53 (3.59-43.70) p=7.3e-05', 'lab_prior_miss_other: OR 2.94 (1.22-7.09) p=0.016',
          'lab_prior_miss_other: OR 0.75 (0.20-2.82) p=0.67', 'pdbj: OR 5.71 (3.01-10.85)']:
    ck('step4: ' + s, s in t22, True)
t23 = txt('step4/23_stdout.txt')
for s in ['carried (first missed, later missed)           12', 'escaped (first missed, later all recorded)      9',
          'first recorded, later all recorded            406', 'pcm: carried 4/12 vs always-recorded 402/406  Fisher p=1.1e-11']:
    ck('step4 CCD: ' + s.split('  ')[0], s in t23, True)
t24 = txt('step4/24_stdout.txt')
ck('PCM bulk date 2024-09-27 count', "'2024-09-27': 3454" in t24, True)
ck('post-update misses without/with PCM', ('code had no PCM: 1/28 missed' in t24) and ('code had PCM: 6/271 missed' in t24), True)
ck('first-missed codes lacking PCM', "['59H', '80I', 'CPR', 'HUR', 'MKZ', 'PBG', 'TSU', 'XVE']" in t24, True)
ck('pre-reaction carried vs always-recorded', 'carried 9 / 12 vs always-recorded 220 / 406 Fisher p 0.24' in t24, True)

# ---- step 5 ----
t27 = txt('step5/27_stdout.txt')
for s in ['c_first_miss: OR 13.72 (4.52-41.72) p=3.9e-06', 'pdbj: OR 6.61 (3.25-13.45) p=1.9e-07',
          'struct_conn touched in a later revision: recorded 1673 / 5021  missed 136 / 334',
          'deposited after the first entry was released: 46 / 55 | missed among after: 18 / 46  among before: 7 / 9']:
    ck('step5: ' + s[:40], s in t27, True)
t28 = txt('step5/28_stdout.txt')
ck('leaving flag vs CCD (1481, 3559, 1)', bool(re.search(r'none\s+3559\s+1\s*\n\s*one\s+0\s+1481', t28)), True)
ck('atom pair differs: 108 units, 45 codes', 'atom pair differs from first recorded: 108 units in 45 codes' in t28, True)
FM = pd.read_csv(f'{L}/step5/first_missed_later_units.tsv', sep='\t', dtype={'entry': str})
ck('first-missed codes with later use', FM.code.nunique(), 21)
g = FM.groupby('code').miss.sum(); ck('carried codes (>=1 later miss)', int((g > 0).sum()), 12)
for code, ents in [('MKZ', ['7RP3', '7RP4']), ('XVE', ['4Q24']), ('59H', ['8VEQ'])]:
    ck(f'{code} later uses all missed', sorted(FM[FM.code == code].entry), ents)
    ck(f'{code} later misses', int(FM[FM.code == code].miss.sum()), len(ents))
def bond(e, c): s = B[(B.entry == e) & (B.lig_comp == c)]; return sorted(zip(s.nuc_res, s.nuc_seq.astype(int), s.dist.round(2), s.recorded))
ck('6T5V MKZ', bond('6T5V', 'MKZ'), [('CYS', 12, 1.82, False)])
ck('7RP3 MKZ', bond('7RP3', 'MKZ'), [('CYS', 12, 1.84, False)])
ck('7RP4 MKZ', bond('7RP4', 'MKZ'), [('CYS', 12, 1.81, False), ('CYS', 12, 1.82, False)])
ck('2XYP XVE', bond('2XYP', 'XVE'), [('CYS', 163, 1.79, False)])
ck('4Q24 XVE', bond('4Q24', 'XVE'), [('CYS', 37, 1.77, False)])
ck('5DF7 59H', bond('5DF7', '59H'), [('SER', 294, 1.39, False), ('SER', 294, 1.39, False)])
ck('8VEQ 59H', bond('8VEQ', '59H'), [('SER', 310, 1.37, False)])

# ---- Cys-sample feature models (step 1) ----
LM = pd.read_csv(f'{L}/logistic_models.tsv', sep='\t')
def lm(model, term): r = LM[(LM.model == model) & (LM.term == term)].iloc[0]; return (r.OR, r.ci_lo, r.ci_hi)
ck('pre-reaction adjusted OR', lm('site+year+unreacted', 'unreacted'), (0.95, 0.44, 2.06))
ck('any leaving atom adjusted OR', lm('site+year+leaving', 'leaving'), (0.88, 0.37, 2.12))
ck('leaving atom on attachment atom OR', lm('site+year+leaving_at_attach', 'leaving_at_attach'), (0.14, 0.04, 0.52))
ck('code created for this entry OR', lm('site+year+newcode', 'newcode'), (0.36, 0.2, 0.64))
ck('PDBj OR site+year', lm('site+year', 'C(site, Treatment("RCSB"))[T.PDBJ]'), (12.97, 6.35, 26.5))
ch = pd.read_csv(f'{L}/code_history_logistic.tsv', sep='\t')
ck('born-unlinked OR (Cys sample) present', bool(((ch.OR.round(2) == 5.89) | (ch.OR.round(1) == 5.9)).any()), True)

# ---- leritrelvir ----
C = pd.read_csv(f'{L}/leritrelvir_copies.tsv', sep='\t', dtype={'entry': str})
x = C[C.deposit_date == '2026-02-11']; a, b = x[x.lig_comp == '7ON'], x[x.lig_comp == 'A1EN0']
ck('Xiong batch entries 7ON / A1EN0', (a.entry.nunique(), b.entry.nunique()), (11, 5))
ck('Xiong batch copies recorded 7ON', (int(a.covale_row_c7.sum()), len(a)), (22, 22))
ck('Xiong batch copies recorded A1EN0', (int(b.covale_row_c7.sum()), len(b)), (0, 10))
ck('A1EN0 copies in close-contact table', int(b.sg_c7_in_close_contact.sum()), 10)
ck('7ON copies in close-contact table', int(a.sg_c7_in_close_contact.sum()), 0)
ck('7ON SG-C7 range', (round(a.sg145_c7_dist.min(), 3), round(a.sg145_c7_dist.max(), 3)), (1.849, 1.92))
ck('A1EN0 SG-C7 range', (round(b.sg145_c7_dist.min(), 3), round(b.sg145_c7_dist.max(), 3)), (1.815, 1.932))
ck('7ON C7 angle-sum range', (round(a.c7_angle_sum.min(), 1), round(a.c7_angle_sum.max(), 1)), (347.3, 356.7))
ck('A1EN0 C7 angle-sum range', (round(b.c7_angle_sum.min(), 1), round(b.c7_angle_sum.max(), 1)), (348.5, 358.8))
ck('Mann-Whitney distance p', round(mannwhitneyu(a.sg145_c7_dist, b.sg145_c7_dist).pvalue, 2), 0.82)
ck('Mann-Whitney angle p', round(mannwhitneyu(a.c7_angle_sum, b.c7_angle_sum).pvalue, 2), 0.70)
ck('Fisher per copy', f'{fisher_exact([[22, 0], [0, 10]])[1]:.1e}', '1.6e-08')
ck('Fisher per entry', f'{fisher_exact([[11, 0], [0, 5]])[1]:.1e}', '2.3e-04')
ck('has_protein_modification by code (Xiong batch)', (set(a.has_protein_modification), set(b.has_protein_modification)), ({'Y'}, {'N'}))

# ---- Mpro census ----
AF = pd.read_csv(f'{MP}/annotation_fields.tsv', sep='\t'); l2 = AF[AF.band == '<2.0']
ck('Mpro bonds < 2.0 A recorded / total', (int(l2.recorded.sum()), len(l2)), (1037, 1109))
ck('Mpro unrecorded in close-contact table', int((~l2.recorded & l2.flagged_close_contact).sum()), 72)
ck('Mpro recorded in close-contact table', int((l2.recorded & l2.flagged_close_contact).sum()), 0)
for site, rec, tot in [('PDBJ', 112, 168), ('PDBC', 81, 93), ('PDBE', 156, 159), ('RCSB', 688, 689)]:
    s = l2[l2.process_site == site]; ck(f'Mpro {site} recorded / total', (int(s.recorded.sum()), len(s)), (rec, tot))
PF = pd.read_csv(f'{MP}/protein_modification_flag.tsv', sep='\t', keep_default_na=False)
ck('Mpro covalent entries', len(PF), 807)
ck('Mpro entries with recorded bond: HPM Y', int(((PF.all_unrec == False) & (PF.has_protein_modification == 'Y')).sum()), 754)
nr = PF[PF.all_unrec == True]
ck('Mpro entries with no recorded bond (N, absent, mod rows)', (int((nr.has_protein_modification == 'N').sum()), int((nr.has_protein_modification == '').sum()), int(nr.n_modification_feature.astype(int).sum())), (28, 25, 0))
PD = pd.read_csv(f'{MP}/per_drug_by_site.tsv', sep='\t')
w = PD[PD.comp == '4WI'].set_index('site')
ck('nirmatrelvir Cys145 recorded outside PDBj', (int(w.loc[['RCSB', 'PDBE', 'PDBC'], 'recorded_lt2'].sum()), int(w.loc[['RCSB', 'PDBE', 'PDBC'], 'contacts_lt2'].sum())), (56, 56))
ck('nirmatrelvir Cys145 recorded at PDBj', (int(w.loc['PDBJ', 'recorded_lt2']), int(w.loc['PDBJ', 'contacts_lt2'])), (6, 24))
N = pd.read_csv(f'{D}/cys_scan_named_drug_contacts.tsv', sep='\t'); n4 = N[(N.lig_comp == '4WI') & (N.dist < 2.0)]
o = n4[n4.process_site != 'PDBJ']; j = n4[n4.process_site == 'PDBJ']
ck('nirmatrelvir, Cys-scan sample (all Cys targets): outside PDBj', (int(o.recorded.sum()), len(o)), (74, 74))
ck('nirmatrelvir, Cys-scan sample: PDBj', (int(j.recorded.sum()), len(j)), (7, 30))
ck('Mpro PDBj contact Fisher', f'{fisher_exact([[56, 112], [16, 925]])[0]:.1f}', '28.9')

# ---- downstream ----
BS = pd.read_csv(f'{D}/all_benchmarks_headline_with_without.tsv', sep='\t', dtype=str)
dd = pd.to_numeric(BS.delta_pp, errors='coerce'); ms = pd.to_numeric(BS.maxshift_pp, errors='coerce')
ck('largest observed headline change (pp)', float(dd.abs().max()), 0.33)
ck('largest possible headline change (pp)', float(ms.max()), 0.99)
t = open(f'{D}/plinder_summary_numbers.txt').read()
ck('PLINDER: 0/180 unrecorded flagged', 'in PLINDER: 180; flagged is_covalent via this Cys: 0/180' in t, True)
ck('PLINDER: recorded flagged 3716/3754', 'flagged is_covalent via this Cys: 3716/3754' in t, True)

# ---- v1.1: contacts against the MAXIT windows (Varman distance-table note v1.1, doi:10.5281/zenodo.24295951, Table 1) ----
W = pd.read_csv(f'{L}/step6/contacts_windowed.tsv', sep='\t', dtype={'entry': str}, low_memory=False)
ck('windowed contacts total / v1.0 rows / 2.0-2.2 A extension / Cys SG-ligand O,N extension', (len(W), int((W.src == 'bonds_all').sum()), int((W.src == 'cc_2.0-2.2').sum()), int((W.src == 'cc_SG_O/N').sum())), (47776, 45450, 2201, 125))
ck('Cys SG-ligand O,N extension: all unrecorded, all Cys SG, ligand O or N', (int(W[W.src == 'cc_SG_O/N'].recorded.sum()), set(W[W.src == 'cc_SG_O/N'].nuc_atom), set(W[W.src == 'cc_SG_O/N'].lig_el)), (0, {'SG'}, {'N', 'O'}))
T1 = {('C', 'S'): (1.55, 2.00), ('C', 'O'): (1.10, 1.55), ('C', 'N'): (1.20, 1.60), ('B', 'O'): (1.20, 1.65), ('B', 'N'): (1.20, 1.65), ('S', 'S'): (1.85, 2.20),
      ('O', 'P'): (1.40, 1.75), ('N', 'P'): (1.54, 1.77), ('O', 'S'): (1.35, 1.80), ('N', 'S'): (1.45, 1.90), ('C', 'Se'): (1.75, 2.15)}   # Table 1 as published
WIN = {}
for (x, y), w in T1.items(): WIN[(x, y)] = w; WIN[(y, x)] = w
def wclass(r):
    w = WIN.get((r.nuc_el, r.lig_el))
    if w is None: return 'no window'
    return 'below' if r.dist < w[0] else ('above' if r.dist > w[1] else 'inside')
ck('window class recomputed from Table 1 limits (unordered pairs) matches the table', bool((W.apply(wclass, axis=1) == W.window).all()), True)
ck('no contact with a Se nucleophile, so C-Se cannot occur', int((W.nuc_el == 'Se').sum()), 0)
I = W[W.window == 'inside']; pc = lambda u, n: round(100 * u / n, 1)
nu = lambda x: int((~x.recorded).sum())
ck('in-window contacts, unrecorded, %', (len(I), nu(I), pc(nu(I), len(I))), (40899, 2971, 7.3))
ck('in-window unrecorded: entries, codes', (I[~I.recorded].entry.nunique(), I[~I.recorded].lig_comp.nunique()), (703, 368))
I = I.assign(t1=['-'.join(k) if k in T1 else '-'.join(k[::-1]) for k in zip(I.lig_el, I.nuc_el)])
for pr, n, u, pct in [('C-S', 28706, 2187, 7.6), ('C-N', 5919, 270, 4.6), ('S-S', 2095, 210, 10.0), ('C-O', 1963, 107, 5.5), ('O-P', 1081, 38, 3.5), ('B-O', 746, 95, 12.7),
                      ('O-S', 177, 30, 16.9), ('N-S', 121, 26, 21.5), ('N-P', 54, 8, 14.8), ('B-N', 37, 0, 0.0)]:
    s = I[I.t1 == pr]; ck(f'in-window {pr} (contacts, unrecorded, %)', (len(s), nu(s), pc(nu(s), len(s))), (n, u, pct))
ck('four smaller windows (O-S, N-S, N-P, B-N): contacts', len(I[I.t1.isin(['O-S', 'N-S', 'N-P', 'B-N'])]), 389)
for pr, n, u in [('S-O', 121, 7), ('O-S', 56, 23), ('S-N', 56, 8), ('N-S', 65, 18)]:
    s = I[I.pair == pr]; ck(f'in-window ligand-first {pr} (contacts, unrecorded)', (len(s), nu(s)), (n, u))
nb = I[I.cat != 'bilin']; ck('in-window without bilins (contacts, unrecorded, %)', (len(nb), nu(nb), pc(nu(nb), len(nb))), (28851, 1316, 4.6))
ck('in-window non-bilin unrecorded at PDBj', int((~nb.recorded & (nb.ps == 'PDBJ')).sum()), 433)
ck('distance-table note: 10,814 of 11,325 recorded -> 511 unrecorded, 4.5%', (11325 - 10814, round(100 * (11325 - 10814) / 11325, 1)), (511, 4.5))
def onb(r): w = WIN[(r.nuc_el, r.lig_el)]; return abs(r.dist - w[0]) < 1e-9 or abs(r.dist - w[1]) < 1e-9
ob = I[I.apply(onb, axis=1)]; ck('in-window contacts on a bound; C-S at 2.00 A', (len(ob), int(((ob.pair == 'C-S') & (ob.dist == 2.0)).sum())), (128, 90))
ex = I.drop(ob.index); ck('bounds exclusive (contacts, unrecorded, %)', (len(ex), nu(ex), pc(nu(ex), len(ex))), (40771, 2856, 7.0))
ck('in-window contacts at >= 2.0 A (all, unrecorded)', (int((I.dist >= 2.0).sum()), nu(I[I.dist >= 2.0])), (1786, 250))
lt2w = W[(W.dist < 2.0) & (W.src == 'bonds_all')]
ck('v1.0 bonds < 2.0 A and unrecorded (from windowed table)', (len(lt2w), nu(lt2w)), (43488, 5594))
for wc, n, u in [('inside', 39072, 2680), ('below', 2384, 2278), ('above', 963, 622), ('no window', 1069, 14)]:
    s = lt2w[lt2w.window == wc]; ck(f'bonds < 2.0 A, {wc} (total, unrecorded)', (len(s), nu(s)), (n, u))
nw = lt2w[lt2w.window == 'no window']
s = nw[nw.lig_el.isin(['O', 'N']) & nw.nuc_el.isin(['O', 'N'])]; ck('no-window: ligand O/N on Ser/Lys/His (contacts, unrecorded)', (len(s), nu(s)), (1003, 0))
ck('no-window unrecorded are all As or Se ligands', (nu(nw), set(nw[~nw.recorded].lig_el)), (14, {'As', 'Se'}))
wo = W[W.window.isin(['below', 'above'])]; ck('windowed pairs outside window (total, recorded, %)', (len(wo), int(wo.recorded.sum()), pc(int(wo.recorded.sum()), len(wo))), (5523, 648, 11.7))
for site, n, u, pct in [('PDBJ', 14779, 1965, 13.3), ('PDBC', 907, 166, 18.3), ('RCSB', 14440, 550, 3.8), ('PDBE', 9696, 266, 2.7)]:
    s = I[I.ps == site]; ck(f'in-window at {site} (contacts, unrecorded, %)', (len(s), nu(s), pc(nu(s), len(s))), (n, u, pct))
for cat, n, u, pct in [('bilin', 12048, 1655, 13.7), ('lipid/acyl', 956, 141, 14.7), ('drug-like', 12245, 735, 6.0), ('heme', 8831, 65, 0.7)]:
    s = I[I.cat == cat]; ck(f'in-window {cat} (contacts, unrecorded, %)', (len(s), nu(s), pc(nu(s), len(s))), (n, u, pct))
ck('in-window PEB unrecorded', nu(I[I.lig_comp == 'PEB']), 1200)
DS = pd.read_csv(f'{L}/step6/29b_disulf_ligand.tsv', sep='\t', dtype=str)
kk = set(zip(W.entry, W.nuc_chain.astype(str), W.nuc_seq.astype(str), W.lig_comp, W.lig_chain.astype(str), W.lig_seq.astype(str)))
inn = sum((r.entry, r.nuc_chain, r.nuc_seq, r.lig_comp, r.lig_chain, r.lig_seq) in kk for r in DS.itertuples())
ck('ligand-Cys disulf rows: count, residue pairs, all in S-S window, absent from contact table', (len(DS), len(DS.drop_duplicates(['entry', 'nuc_chain', 'nuc_seq', 'lig_comp', 'lig_chain', 'lig_seq'])), bool(DS.dist.astype(float).between(1.85, 2.20).all()), len(DS) - inn), (34, 34, True, 29))
xb = W[W.entry.isin(['23MB', '23MD', '23MF', '23MG', '23MI']) & (W.lig_comp == 'A1EN0')]
ck('A1EN0 2026 batch: contacts, inside C-S window, recorded, header distance range', (len(xb), int((xb.window == 'inside').sum()), int(xb.recorded.sum()), round(xb.dist.min(), 2), round(xb.dist.max(), 2)), (10, 10, 0, 1.81, 1.93))
x7 = W[W.entry.isin(['23LW', '23LX', '23LY', '23LZ', '23MA', '23MC', '23ME', '23MH', '23MJ', '23MK', '23ML']) & (W.lig_comp == '7ON')]
ck('7ON 2026 batch: contacts, inside, recorded', (len(x7), int((x7.window == 'inside').sum()), int(x7.recorded.sum())), (22, 22, 22))
ck('A1EN0 coordinate copies inside the C-S window', int(b.sg145_c7_dist.between(1.55, 2.00).sum()), 10)
ck('7ON coordinate copies inside the C-S window', int(a.sg145_c7_dist.between(1.55, 2.00).sum()), 22)
h9 = C[C.entry == '9MCH']; ck('9MCH SG-C7 range', (round(h9.sg145_c7_dist.min(), 2), round(h9.sg145_c7_dist.max(), 2)), (2.74, 2.81))
ck('7ON / A1EN0 heavy atoms', (set(C[C.lig_comp == '7ON'].ccd_n_heavy), set(C[C.lig_comp == 'A1EN0'].ccd_n_heavy)), ({45}, {45}))
# Mpro census against the C-S window
ck('Mpro Cys145 bonds < 2.0 A inside 1.55-2.00 A (all, unrecorded)', (int(l2.dist.between(1.55, 2.0).sum()), int(l2[~l2.recorded].dist.between(1.55, 2.0).sum())), (1109, 72))
u72 = l2[~l2.recorded]; ck('Mpro unrecorded distance range', (round(u72.dist.min(), 3), round(u72.dist.max(), 3)), (1.696, 1.99))
ck('Mpro unrecorded entries', u72.entry.nunique(), 45)
YS = pd.read_csv(f'{MP}/by_deposit_year_site_entries.tsv', sep='\t', header=[0, 1], index_col=0)
yy = YS[YS.index.astype(str).str.isdigit()]; yy.index = yy.index.astype(int); late = yy[yy.index >= 2024]
ck('Mpro PDBj covalent entries deposited 2024+: with a miss / total', (int(late[('sum', 'PDBJ')].sum()), int(late[('count', 'PDBJ')].sum())), (0, 14))
UV = pd.read_csv(f'{MP}/unrecorded_version_history.tsv', sep='\t', dtype={'entry': str})
ck('Mpro unrecorded contacts (<2.4 A) never recorded in an archived version', (len(UV), int((UV.covale_in_any_earlier_version == False).sum())), (89, 89))
pj = UV[(UV.process_site == 'PDBJ') & (UV.band == '<2.0')]; dd = M.set_index('e').loc[sorted(set(pj.entry)), 'dep']
ck('Mpro PDBj misses: deposition months', (dd.min()[:7], dd.max()[:7]), ('2020-12', '2023-02'))

# ---- v1.1: carry-forward on in-window contacts (script 30) ----
t30 = txt('step6/30_stdout.txt')
def blk(a, b=None, t=None): x = (t or t30).split(a, 1)[1]; return x.split(b, 1)[0] if b else x
ck('in-window contacts < 2.0 A in no bonded unit', 'in-window contacts < 2.0 A 39113 unrecorded 2721 ; not in any bonded unit (left out): 31 unrecorded 31' in t30, True)
ck('in-window units / misses / codes', 'units with >=1 in-window contact < 2.0 A 10777 misses 552 codes 3725' in t30, True)
ck('in-window strict set', 'in-window strict units 5611 misses 299 codes 3475 labs 1527' in t30, True)
ck('overlap with v1.0 strict set', 'overlap with the original strict set: 5291 units' in t30, True)
vo = blk('########## VALIDATION', '########## IN-WINDOW')
for s_ in ['[orig] step-3 GLM first_miss: OR 14.23 (4.28-47.38)', '15/40 (38%) vs 86/1289 (7%) OR 8.5 (4.3-16.5)', 'c_first_miss: OR 10.22 (2.74-38.13)', 'c_first_miss: OR 12.53 (3.59-43.70)']:
    ck('validation reproduces v1.0: ' + s_[:40], s_ in vo, True)
w3 = blk('########## IN-WINDOW', '[win] == strict')
for s_ in ['later uses 1813; first missed 25/55 (45%) codes=19; first recorded 77/1758 (4%)', 'step-3 GLM first_miss: OR 16.65 (5.11-54.29) p=3.1e-06',
           'other-lab later uses only: first_miss OR 10.46 (2.64-41.40)']:
    ck('in-window step 3: ' + s_[:45], s_ in w3, True)
ws = blk('[win] == strict', '[win] == drop HUR/80I')
for s_ in ['first-missed codes with a later use 19', '15/41 (37%) vs 73/1369 (5%) OR 10.3 (5.3-20.1)', '4/18 (22%) vs 50/738 (7%) OR 4.2 (1.4-12.7)',
           'same code: 20/43 (47%)', 'other code, same drug: 0/5 (0%) labs=4', 'other code, different drug: 67/583', 'baseline labs (earlier units, none missed): 82/2773',
           'c_first_miss: OR 10.59 (2.77-40.51)', 'lab_prior_miss_other: OR 3.50 (1.32-9.26)', 'pdbj: OR 7.29 (3.67-14.46)',
           '59H 1/1, 80I 7/7, ALD 0/12, CPR 2/2, EOF 0/1, FFQ 0/4, FL3 0/2, FNO 0/1, H7N 0/1, HDS 1/3, HUR 8/8, ISP 1/5, MKZ 2/2, MXC 0/1, MXF 0/1, O2F 0/1, TSU 1/1, XVE 1/1, YDD 1/1']:
    ck('in-window strict: ' + s_[:45], s_ in ws, True)
fm = dict((c, tuple(map(int, f.split('/')))) for c, f in re.findall(r'(\w+) (\d+/\d+)', ws.split('(code: later missed/later): ')[1].split('\n')[0]))
ck('first-missed codes: carried (>=1 later miss), escaped', (sorted(c for c, v in fm.items() if v[0] > 0), sorted(c for c, v in fm.items() if v[0] == 0)),
   (['59H', '80I', 'CPR', 'HDS', 'HUR', 'ISP', 'MKZ', 'TSU', 'XVE', 'YDD'], ['ALD', 'EOF', 'FFQ', 'FL3', 'FNO', 'H7N', 'MXC', 'MXF', 'O2F']))
ck('later misses: total, from HUR and 80I', (sum(v[0] for v in fm.values()), fm['HUR'][0] + fm['80I'][0]), (25, 15))
wh = blk('[win] == drop HUR/80I', '[win] == drop Li J. lab')
for s_ in ['7/33 (21%) vs 73/1369 (5%) OR 5.0 (2.1-11.6)', 'c_first_miss: OR 6.74 (1.59-28.51) p=0.0095']: ck('in-window without HUR/80I: ' + s_[:40], s_ in wh, True)
wl = blk('[win] == drop Li J. lab', '[win] == drop PDBj'); ck('in-window without Li J.: joint OR', 'c_first_miss: OR 7.79 (1.78-34.03) p=0.0063' in wl, True)
wp = blk('[win] == drop PDBj', '[sensitivity win+bondlike]')
for s_ in ['first-missed codes with a later use 11', '6/25 (24%) vs 27/1044 (3%) OR 12.3 (4.7-32.4)', 'c_first_miss: OR 15.19 (3.82-60.32)',
           'lab_prior_miss_other: OR 0.62 (0.11-3.59)', '59H 1/1, ALD 0/9, CPR 2/2, FFQ 0/4, H7N 0/1, HDS 1/3, ISP 1/5, MKZ 2/2, O2F 0/1, TSU 1/1, XVE 1/1']:
    ck('in-window without PDBj: ' + s_[:40], s_ in wp, True)
fp = dict((c, int(f.split('/')[0])) for c, f in re.findall(r'(\w+) (\d+/\d+)', wp.split('(code: later missed/later): ')[1].split('\n')[0]))
ck('without PDBj: later misses, codes carrying them', (sum(fp.values()), sum(v > 0 for v in fp.values())), (9, 7))
wb = blk('[sensitivity win+bondlike]', '########## ATTRIBUTION')
for s_ in ['units 5291 misses 282', 'step-3 GLM first_miss: OR 17.82 (5.12-61.96)']: ck('window + bond-like: ' + s_, s_ in wb, True)
ck('window + bond-like, without HUR/80I', 'c_first_miss: OR 7.65 (1.61-36.39)' in blk('[win+bondlike] == drop HUR/80I', '[win+bondlike] == drop Li J. lab'), True)
WR = pd.read_csv(f'{L}/step6/window_results.tsv', sep='\t')
sw = WR[(WR['set'] == 'win') & (WR.subset == 'strict') & (WR.design == 'lab after first miss: other code, same drug')].stat.iloc[0]
ck('in-window code switches (same five units as v1.0)', sw, '1FM->WNT 8YP4 miss=0;J3U->3KI 7FEA miss=0;NOH->7K8 5M4Z miss=0;O2F->MER 8RWS miss=0;O2F->MER 9HPW miss=0')
AT = pd.read_csv(f'{L}/step6/window_attribution.tsv', sep='\t').set_index('set')
r = AT.loc['in-window strict set']
ck('attribution, strict', (int(r.misses), int(r.first_day_misses), int(r.later_misses_first_recorded), int(r.later_misses_first_missed), int(r.later_units_first_missed), float(r.expected_if_first_recorded), float(r.excess), float(r.share_later_first_missed), float(r.share_excess)),
   (299, 197, 77, 25, 55, 4.2, 20.8, 0.084, 0.07))
ck('attribution, strict: misses at PDBj/PDBc (first day, later first recorded)', (int(r['first day: misses at PDBj or PDBc']), int(r['later, first recorded: misses at PDBj or PDBc'])), (85, 50))
r = AT.loc['all in-window units']
ck('attribution, all in-window units', (int(r.units), int(r.misses), int(r.later_misses_first_missed), float(r.excess), float(r.share_later_first_missed), float(r.share_excess)), (10777, 552, 61, 54.1, 0.111, 0.098))
ck('contact level: in-window unrecorded < 2.0 A in units; in later uses of first-missed codes', 'in-window unrecorded contacts < 2.0 A 2690; in later uses of first-missed codes 86' in t30, True)
SW = pd.read_csv(f'{L}/step6/window_strict_units.tsv', sep='\t', dtype={'entry': str}, low_memory=False)
def hist(c): g = SW[SW.lig_comp == c].sort_values('dep'); return [(e, int(m), str(d)[:4]) for e, m, d in zip(g.entry, g.miss, g.dep)]
ck('TSU in-window strict history', hist('TSU'), [('2CHA', 1, '1975'), ('1EST', 1, '1976')])
ck('TSU uses: laboratories', SW[SW.lig_comp == 'TSU'].lab.nunique(), 2)
ck('HDS in-window strict history', hist('HDS'), [('1LPO', 1, '1995'), ('1LPP', 1, '1995'), ('1EXW', 0, '2000'), ('1QZ3', 0, '2003')])
ck('HDS first two entries: one laboratory', SW[SW.entry.isin(['1LPO', '1LPP']) & (SW.lig_comp == 'HDS')].lab.nunique(), 1)
ck('ISP in-window strict history', hist('ISP'), [('1NTP', 1, '1987'), ('1TRY', 0, '1994'), ('1GMH', 0, '1994'), ('1TRN', 0, '1995'), ('1O6E', 0, '2002'), ('1XLU', 1, '2004')])
ck('TSU, HDS contacts are S-O (ligand S on Ser OG); ISP is P-O', (set(I[I.lig_comp == 'TSU'].pair), set(I[(I.lig_comp == 'HDS') & I.entry.isin(['1LPO', '1EXW', '1QZ3'])].pair), set(I[I.lig_comp == 'ISP'].pair), set(I[I.lig_comp.isin(['TSU', 'ISP']) | ((I.lig_comp == 'HDS') & (I.pair == 'S-O'))].nuc_res)), ({'S-O'}, {'S-O'}, {'P-O'}, {'SER'}))
import numpy as np, statsmodels.api as sm, statsmodels.formula.api as smf
Dx = SW.sort_values(['lig_comp', 'dep', 'rel', 'entry']).copy(); f1 = Dx.groupby('lig_comp').head(1).set_index('lig_comp')
Dx['first_miss'] = Dx.lig_comp.map(f1.miss); Dx['first_dep'] = Dx.lig_comp.map(f1.dep); Dx['first_lab'] = Dx.lig_comp.map(f1.au_last)
xx = Dx[Dx.dep > Dx.first_dep].copy(); xx['same_lab'] = xx.au_last == xx.first_lab; xx['code'] = xx.lig_comp.astype('category').cat.codes
sl_ = [v for v in sorted(xx.site.unique()) if v != 'RCSB']   # site factor with RCSB as reference, coded by hand (C is the leritrelvir table here)
for v in sl_: xx['s_' + v] = (xx.site == v).astype(int)
xx['same_lab'] = xx.same_lab.astype(int)
mm = smf.glm('miss ~ first_miss + same_lab + ' + ' + '.join('s_' + v for v in sl_) + ' + yearc', xx, family=sm.families.Binomial()).fit(cov_type='cluster', cov_kwds={'groups': xx.code})
ck('step-3 in-window OR recomputed from window_strict_units.tsv (1 decimal)', (len(xx), round(float(np.exp(mm.params['first_miss'])), 1), tuple(round(float(v), 1) for v in np.exp(mm.conf_int().loc['first_miss']))), (1813, 16.7, (5.1, 54.3)))
ck('v1.0 bond-like joint model PDBj term (Table 3)', 'pdbj: OR 5.71 (3.01-10.85)' in blk('[orig] == strict', '[orig] == drop HUR/80I'), True)
# four-window sensitivity (data/ligand_code_effect/step6/four_window)
t4 = txt('step6/four_window/30_stdout.txt')
ck('four windows: without HUR/80I joint OR', 'c_first_miss: OR 5.37 (0.78-36.73) p=0.087' in blk('[win] == drop HUR/80I', '[win] == drop Li J. lab', t4), True)
ck('four windows: without Li J. joint OR', 'c_first_miss: OR 6.54 (0.91-46.96) p=0.062' in blk('[win] == drop Li J. lab', '[win] == drop PDBj', t4), True)

# ---- v1.1: retained numbers not checked in v1.0 ----
ck('Mpro PDBj contact Fisher p', f'{fisher_exact([[56, 112], [16, 925]])[1]:.1e}', '1.1e-35')
def cp(e): s = C[C.entry == e].sort_values('sg145_c7_dist'); return (sorted(s.sg145_c7_dist.round(3)), int(s.covale_row_c7.sum()))
ck('23MC 7ON distances, recorded', cp('23MC'), ([1.879, 1.885], 2)); ck('23MD A1EN0 distances, recorded', cp('23MD'), ([1.862, 1.898], 0))
g7 = C[C.lig_comp == '7ON'].iloc[0]; gA = C[C.lig_comp == 'A1EN0'].iloc[0]
ck('7ON CCD created / site / defining entry', (g7.ccd_initial_date, g7.ccd_create_site, g7.ccd_defining_entry), ('2023-03-01', 'PDBC', '8IGN'))
ck('A1EN0 CCD created / site / defining entry', (gA.ccd_initial_date, gA.ccd_create_site, gA.ccd_defining_entry), ('2025-03-24', 'PDBC', '9MCH'))
t27r = t27.split('miss ~ c_first_miss + C(c_first_era) + C(ref')[1].split('miss ~')[0]
ps_ = [float(m) for m in re.findall(r'C\(ref, Treatment\("PHENIX"\)\)\[T\.\w+\]: OR [\d.]+ \([\d.]+-[\d.]+\) p=([\d.e-]+)', t27r)]
ck('refinement-program terms: n, min p', (len(ps_), min(ps_)), (4, 0.097))
for code, n, labs in [('ALD', 12, 4), ('FFQ', 4, 3)]:
    s_ = FM[FM.code == code]; ck(f'{code} escaped: later units recorded, labs', (int((s_.miss == 0).sum()), len(s_), s_.lab.nunique()), (n, n, labs))
s_ = FM[FM.code == 'ALD']; ck('ALD later uses: sites, first year', (sorted(set(s_.site)), str(s_.dep.min())[:4]), (['PDBE', 'PDBJ'], '2020'))
SL = pd.read_csv(f'{D}/slipped_into_noncovalent_or_test_sets.tsv', sep='\t', dtype={'entry': str}); sl = SL[SL.band == '<2.00']
ck("Runs N' Poses entries with unrecorded Cys bond < 2.0 A (all PDBj)", (sorted(set(sl[sl['where'].str.startswith('Runs')].entry)), set(sl[sl['where'].str.startswith('Runs')].site)), (['7EI4', '7VVT', '7YCC', '7YCE', '8X6R'], {'PDBJ'}))
ck('LP-PDBBind test split entries < 2.0 A, labelled False', (sorted(set(sl[sl['where'].str.startswith('LP')].entry)), set(sl[sl['where'].str.startswith('LP')].label_in_set)), (['2PRE', '3O0U', '3WZU', '5Z1D', '6JZ0'], {'LP covalent=False'}))
ck('DiffDock time-split test entries', sorted(set(SL[SL['where'].str.startswith('DiffDock')].entry)), ['6I9A', '6S9W', '6S9X'])
CS = pd.read_csv(f'{D}/cross_dataset_summary.tsv', sep='\t', dtype=str)
ck('BioLiP2 holds 194 of the 197', CS[CS.dataset == 'BioLiP2'].iloc[0, 5].split(' / ')[0], '194')

ck('CCD definitions downloaded (ok rows in log)', sum(1 for l in open(f'{L}/ccd_download_log.tsv') if l.rstrip('\n').endswith('\tok')), 4907)
ck('Mpro entries scanned', 'Mpro entries scanned 1957 ' in open(f'{MP}/summary_numbers.txt').read(), True)

print(f'\n{sum(checks)}/{len(checks)} checks passed')
sys.exit(0 if all(checks) else 1)
