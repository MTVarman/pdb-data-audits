"""Figures for the addendum, drawn from the bundled tables.  python3 scripts/note/make_figures.py
Writes figures/fig1..fig6 (fig1 and fig4 use the in-window tables in data/ligand_code_effect/step6) (.png at 300 dpi and .pdf). Figure 2 is the Mpro census figure, copied as produced
by scripts/mpro_covalent_census/13_figure.py."""
import os, re, shutil
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

H = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
R = f'{H}/data'; F = f'{H}/figures'; os.makedirs(F, exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.titlesize': 8.5, 'axes.titleweight': 'bold',
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.titlelocation': 'left'})
REC, MISS, GREY, DARK, PDBJ = '#2a7f62', '#c0392b', '#8c8c8c', '#2c3e50', '#d4882b'
def save(fig, name):
    fig.savefig(f'{F}/{name}.png', dpi=300, bbox_inches='tight'); fig.savefig(f'{F}/{name}.pdf', bbox_inches='tight'); plt.close(fig)
def ci(s):
    m = re.search(r'OR ([\d.]+) \(([\d.]+)-([\d.]+)\)', s); return tuple(float(v) for v in m.groups())

# ---------- Figure 1: contacts against the MAXIT windows (Varman distance-table note v1.1, Table 1) ----------
W = pd.read_csv(f'{R}/ligand_code_effect/step6/contacts_windowed.tsv', sep='\t', dtype={'entry': str}, low_memory=False)
# Table 1 pairs; the contact table writes the ligand element first, so unordered pairs collect both orders
T1 = [('C–S', ['C-S'], (1.55, 2.00), 'Cys SG – ligand C'), ('C–N', ['C-N'], (1.20, 1.60), 'Lys/His N – ligand C'), ('C–O', ['C-O'], (1.10, 1.55), 'Ser OG – ligand C'),
      ('S–S', ['S-S'], (1.85, 2.20), 'Cys SG – ligand S'), ('O–P', ['P-O'], (1.40, 1.75), 'Ser OG – ligand P'),
      ('B–O', ['B-O'], (1.20, 1.65), 'Ser OG – ligand B'), ('O–S', ['S-O', 'O-S'], (1.35, 1.80), 'Ser OG–S / Cys SG–O'), ('N–S', ['S-N', 'N-S'], (1.45, 1.90), 'Lys/His N–S / SG–N'),
      ('N–P', ['P-N'], (1.54, 1.77), 'Lys/His N – ligand P'), ('B–N', ['B-N'], (1.20, 1.65), 'Lys/His N – ligand B')]
fig = plt.figure(figsize=(7.2, 7.0)); gs = fig.add_gridspec(3, 10, hspace=1.25, wspace=1.9, height_ratios=[1, 1, 1.15])
bins = np.arange(1.0, 2.2201, 0.02)
for j, (p, keys, (lo, hi), lab) in enumerate(T1):
    a = fig.add_subplot(gs[j // 5, 2 * (j % 5):2 * (j % 5) + 2]); s = W[W.pair.isin(keys)]
    a.axvspan(lo, hi, color='#e8eef5', zorder=0)
    a.hist([s[s.recorded].dist.clip(1.0, 2.22), s[~s.recorded].dist.clip(1.0, 2.22)], bins=bins, stacked=True, color=[REC, MISS], lw=0)
    a.set_yscale('log'); a.set_xlim(1.0, 2.22); a.set_xticks([1.2, 1.6, 2.0]); a.tick_params(labelsize=6)
    ins = s[s.window == 'inside']; u = int((~ins.recorded).sum())
    a.set_title(f'{"ABCDEFGHIJ"[j]}  {p} ({lo:.2f}–{hi:.2f})\n{lab}\nin window: {u:,} / {len(ins):,}', fontsize=5.6)
    a.set_xlabel('distance (Å)', fontsize=6)
    if j % 5 == 0: a.set_ylabel('contacts (log)', fontsize=6)
a = fig.add_subplot(gs[2, :5])
I = W[W.window == 'inside']; sites = ['RCSB', 'PDBE', 'PDBJ', 'PDBC']; names = {'RCSB': 'RCSB PDB', 'PDBE': 'PDBe', 'PDBJ': 'PDBj', 'PDBC': 'PDBc'}
for i, sn in enumerate(sites):
    s = I[I.ps == sn]; o = s[s.cat != 'bilin']; b = s[s.cat == 'bilin']
    uo, ub = int((~o.recorded).sum()), int((~b.recorded).sum()); n = len(s)
    a.barh(i, uo / n * 100, color=PDBJ if sn in ('PDBJ', 'PDBC') else GREY); a.barh(i, ub / n * 100, left=uo / n * 100, color='#e4b98a' if sn in ('PDBJ', 'PDBC') else '#cfcfcf')
    a.text((uo + ub) / n * 100 + 0.4, i, f'{uo + ub:,} / {n:,}', va='center', fontsize=6.3)
a.set_yticks(range(4)); a.set_yticklabels([names[x] for x in sites]); a.invert_yaxis(); a.set_xlim(0, 24)
a.set_xlabel('% of in-window contacts with no covale row'); a.set_title('K  In-window contacts by processing site')
a.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=GREY), plt.Rectangle((0, 0), 1, 1, color='#cfcfcf')], labels=['non-bilin ligands', 'bilins'], fontsize=5.8, frameon=False, loc='upper right')
a = fig.add_subplot(gs[2, 5:])
cats = [('bilin', 'bilins'), ('drug-like', 'drug-like'), ('lipid/acyl', 'lipid / acyl'), ('heme', 'hemes'), ('cofactor', 'cofactors'), ('flavin', 'flavins')]
for i, (c, lab) in enumerate(cats):
    s = I[I.cat == c]; u = int((~s.recorded).sum()); n = len(s)
    a.barh(i, u / n * 100, color=MISS); a.text(u / n * 100 + 0.4, i, f'{u:,} / {n:,}', va='center', fontsize=6.3)
a.set_yticks(range(len(cats))); a.set_yticklabels([c[1] for c in cats]); a.invert_yaxis(); a.set_xlim(0, 24)
a.set_xlabel('% of in-window contacts with no covale row'); a.set_title('L  In-window contacts by ligand class')
fig.text(0, -0.035, 'Header mmCIFs of 199,326 entries (16 Nov 2026). One contact per (nucleophile residue, ligand residue): recorded = covale row (pdbx_dist_value);\n'
         'unrecorded = closest pair listed only in _pdbx_validate_close_contact (≤ 2.2 Å). Shaded: MAXIT 11.400 _allowed_bonds window (Varman distance-table note v1.1, Table 1).\n'
         'C–Se (1.75–2.15 Å) cannot occur with these nucleophiles. Unrecorded contacts with a ligand O or N are collected for Cys SG only (O–S, N–S panels).', fontsize=6.3, color='#555')
save(fig, 'fig1_window')

# ---------- Figure 2: Mpro census (copied) ----------
src = f'{R}/mpro_census/mpro_covalent_census.png'
if os.path.exists(src): shutil.copy(src, f'{F}/fig2_mpro_census.png')

# ---------- Figure 3: leritrelvir ----------
C = pd.read_csv(f'{R}/ligand_code_effect/leritrelvir_copies.tsv', sep='\t', dtype={'entry': str})
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={'width_ratios': [1.2, 1], 'wspace': 0.95})
mk = {'2026-02-11': 'o', '2025-03-17': 's', '2023-02-21': 'D'}
for (code, dep), g in C.groupby(['lig_comp', 'deposit_date']):
    col = REC if code == '7ON' else PDBJ
    ax[0].scatter(g.sg145_c7_dist, g.c7_angle_sum, marker=mk[dep], s=22, facecolors=col if dep == '2026-02-11' else 'none', edgecolors=col, linewidths=0.9)
ax[0].axvspan(1.55, 2.00, color='#e8eef5', zorder=0); ax[0].text(1.47, 0.02, 'MAXIT C–S window\n1.55–2.00 Å (shaded)', transform=ax[0].get_xaxis_transform(), ha='left', fontsize=5.3, color='#555'); ax[0].set_xlabel('Cys145 SG – ligand C7 distance (Å)'); ax[0].set_ylabel('sum of bond angles at C7 (°)\n360 = planar ketone')
ax[0].set_title('A  Leritrelvir copies: geometry by CCD code and batch')
h = [Line2D([], [], marker='o', ls='', color=REC, label='7ON (bound form)'), Line2D([], [], marker='o', ls='', color=PDBJ, label='A1EN0 (ketone form)'),
     Line2D([], [], marker='o', ls='', mfc='k', mec='k', label='filled: 11 Feb 2026 batch (16 entries)'), Line2D([], [], marker='s', ls='', mfc='none', mec='k', label='open square: 17 Mar 2025 batch'),
     Line2D([], [], marker='D', ls='', mfc='none', mec='k', label='open diamond: 8IGN (2023)')]
ax[0].set_xlim(1.45, 3.3); ax[0].legend(handles=h, fontsize=5.5, frameon=False, loc='center right', markerscale=0.7)
x = C[C.deposit_date == '2026-02-11']; a, b = x[x.lig_comp == '7ON'], x[x.lig_comp == 'A1EN0']
cats = ['covale row present', 'pair listed in\nclose-contact table', 'has_protein_modification = Y']
va = [a.covale_row_c7.mean(), a.sg_c7_in_close_contact.mean(), (a.drop_duplicates('entry').has_protein_modification == 'Y').mean()]
vb = [b.covale_row_c7.mean(), b.sg_c7_in_close_contact.mean(), (b.drop_duplicates('entry').has_protein_modification == 'Y').mean()]
la = [f'{int(a.covale_row_c7.sum())}/{len(a)}', f'{int(a.sg_c7_in_close_contact.sum())}/{len(a)}', f'{int((a.drop_duplicates("entry").has_protein_modification=="Y").sum())}/{a.entry.nunique()}']
lb = [f'{int(b.covale_row_c7.sum())}/{len(b)}', f'{int(b.sg_c7_in_close_contact.sum())}/{len(b)}', f'{int((b.drop_duplicates("entry").has_protein_modification=="Y").sum())}/{b.entry.nunique()}']
yy = np.arange(3)
ax[1].barh(yy - 0.18, np.array(va) * 100, 0.34, color=REC, label='7ON (11 entries, 22 copies)'); ax[1].barh(yy + 0.18, np.array(vb) * 100, 0.34, color=PDBJ, label='A1EN0 (5 entries, 10 copies)')
for i in range(3):
    ax[1].text(va[i] * 100 + 2, i - 0.18, la[i], va='center', fontsize=6.5); ax[1].text(vb[i] * 100 + 2, i + 0.18, lb[i], va='center', fontsize=6.5)
ax[1].set_yticks(yy); ax[1].set_yticklabels(cats); ax[1].invert_yaxis(); ax[1].set_xlim(0, 125); ax[1].set_xticks([0, 50, 100])
ax[1].set_xlabel('% of copies (rows 1–2) or entries (row 3)'); ax[1].set_title('B  11 Feb 2026 batch: what the file records'); ax[1].legend(fontsize=6, frameon=False, loc='upper left', bbox_to_anchor=(0, -0.2))
save(fig, 'fig3_leritrelvir')

# ---------- Figure 4: carry-forward, contacts inside the MAXIT windows ----------
t30 = open(f'{R}/ligand_code_effect/step6/30_stdout.txt').read()
tw = t30.split('[win] == strict')[1].split('[win] == drop')[0]
def frac(lbl):
    m = re.search(re.escape(lbl) + r'\s*:\s*(\d+)/(\d+)', tw); return int(m.group(1)), int(m.group(2))
fig = plt.figure(figsize=(7.2, 7.8)); gs = fig.add_gridspec(3, 2, hspace=1.0, wspace=0.45, height_ratios=[1, 1.7, 0.8])
a0 = fig.add_subplot(gs[0, 0]); a1 = fig.add_subplot(gs[0, 1]); a2 = fig.add_subplot(gs[1, :]); a3 = fig.add_subplot(gs[2, :])
groups = [('same lab', ['first MISSED | same lab | same site', 'first MISSED | same lab | other site'], ['first recorded | same lab | same site', 'first recorded | same lab | other site']),
          ('new lab,\nsame site', ['first MISSED | other lab | same site'], ['first recorded | other lab | same site']),
          ('new lab,\nnew site', ['first MISSED | other lab | other site'], ['first recorded | other lab | other site'])]
for i, (lab, fm, fr) in enumerate(groups):
    for j, keys in enumerate((fr, fm)):
        k = [frac(q) for q in keys]; m, n = sum(v[0] for v in k), sum(v[1] for v in k)
        a0.bar(i + (j - 0.5) * 0.38, m / n * 100 if n else 0, 0.36, color=MISS if j else REC); a0.text(i + (j - 0.5) * 0.38, (m / n * 100 if n else 0) + 1.5, f'{m}/{n}', ha='center', fontsize=6.3)
a0.set_xticks(range(3)); a0.set_xticklabels([g[0] for g in groups]); a0.set_ylabel('% of later uses missed'); a0.set_ylim(0, 90)
a0.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=REC), plt.Rectangle((0, 0), 1, 1, color=MISS)], labels=["code's first bonded use recorded", "code's first bonded use missed"], fontsize=6, frameon=False, loc='upper left', bbox_to_anchor=(0, -0.32), ncol=1)
a0.set_title('A  Later uses of the same code (in-window)')
bl = [('same code', 'lab after first miss, same code'), ('other code,\nsame drug', 'lab after first miss, other code, same drug'), ('other code,\nother drug', 'lab after first miss, other code, different drug'), ('labs with no\nearlier miss', 'baseline labs (earlier units, none missed)')]
for i, (lab, key) in enumerate(bl):
    m, n = frac(key); a1.bar(i, m / n * 100, 0.6, color=[MISS, REC, PDBJ, GREY][i]); a1.text(i, m / n * 100 + 1.2, f'{m}/{n}', ha='center', fontsize=6.5)
a1.set_xticks(range(4)); a1.set_xticklabels([b[0] for b in bl], fontsize=6.5); a1.set_ylabel('% of later units missed'); a1.set_ylim(0, 75)
a1.set_title("B  A lab's later units after its first miss (in-window)")
DR = pd.read_csv(f'{R}/ligand_code_effect/step6/window_results.tsv', sep='\t')
subs = ['strict', 'drop HUR/80I', 'drop Li J. lab', 'drop PDBj']
terms = [('orig', 'joint GLM c_first_miss', "code's first use missed: bond-like set (as in v1.0)", '#e8a39b'), ('win', 'joint GLM c_first_miss', "code's first use missed: in-window set", MISS),
         ('win', 'joint GLM lab_prior_miss_other', 'lab missed earlier on another code (in-window)', GREY), ('win', 'joint GLM pdbj', 'processed at PDBj (in-window)', DARK)]
for si, sb in enumerate(subs):
    for ti, (st, t, lab, col) in enumerate(terms):
        r = DR[(DR['set'] == st) & (DR.subset == sb) & (DR.design == t)]
        if len(r) == 0: continue
        o, lo, hi = ci(r.stat.iloc[0]); y = si * 5 + ti
        a2.plot([lo, hi], [y, y], color=col, lw=1.4); a2.plot(o, y, 'o', color=col, ms=3.6); a2.text(170, y, f'{o:.2f} ({lo:.2f}–{hi:.2f})', va='center', fontsize=5.4, color=col)
a2.set_xscale('log'); a2.axvline(1, color='k', lw=0.6, ls=':'); a2.set_xlim(0.05, 800)
a2.set_yticks([si * 5 + 1.5 for si in range(4)]); a2.set_yticklabels(['all strict units', 'without HUR and 80I', 'without the Li J. lab', 'without PDBj']); a2.invert_yaxis()
a2.set_xlabel('adjusted odds ratio for a later miss (95% CI clustered by code; same lab and deposition year also in model)', fontsize=6.5)
a2.legend(handles=[Line2D([], [], color=c, marker='o', label=l) for _, _, l, c in terms], fontsize=5.6, frameon=False, loc='upper left', ncol=2, bbox_to_anchor=(0, -0.25))
a2.set_title('C  Joint logistic model on later uses: bond-like set against in-window set')
AT = pd.read_csv(f'{R}/ligand_code_effect/step6/window_attribution.tsv', sep='\t')
lbls = {'all in-window units': 'all in-window units', 'drug-like in-window units': 'drug-like', 'in-window strict set': 'strict set'}
for i, r in AT.iterrows():
    tot = r.misses; parts = [(r.first_day_misses, '#9fb3c8', "on the code's first deposition day"), (r.later_misses_first_recorded, GREY, 'later use, first use recorded'),
                             (r.expected_if_first_recorded, '#e8a39b', 'later use, first use missed: expected anyway'), (r.excess, MISS, 'later use, first use missed: excess')]
    left = 0
    for v, col, lab in parts:
        a3.barh(i, v / tot * 100, left=left, color=col, label=lab if i == 0 else None); left += v / tot * 100
    a3.text(101, i, f'{int(r.later_misses_first_missed)} / {int(tot)}', va='center', fontsize=6)
a3.set_yticks(range(len(AT))); a3.set_yticklabels([lbls[x] for x in AT['set']], fontsize=6.5); a3.invert_yaxis(); a3.set_xlim(0, 100)
a3.set_xlabel('% of in-window (code, entry) misses'); a3.legend(fontsize=5.6, frameon=False, ncol=2, loc='upper left', bbox_to_anchor=(0, -0.42))
a3.set_title('D  Share of in-window misses that follow a missed first use')
fig.text(0, -0.02, 'In-window strict set: drug-like codes with ≤10 unbonded uses and ≥1 contact under 2.0 Å inside a MAXIT window (all ten pairs of Table 1 of the distance-table note v1.1);\n'
         '5,611 (code, entry) units, 299 missed; misses counted on in-window contacts only. Lab = last author. "Same drug" = same InChIKey first block or skeleton key.', fontsize=6.3, color='#555')
save(fig, 'fig4_carry_forward')

# ---------- Figure 5: example timelines ----------
T = pd.read_csv(f'{R}/archive/code_timelines.tsv', sep='\t', dtype={'entry': str}, low_memory=False)
SU = pd.read_csv(f'{R}/ligand_code_effect/step4/strict_units.tsv', sep='\t', dtype={'entry': str}, low_memory=False)
rows = [('MKZ', 'MKZ: 6T5V (Phillips, PDBe, 2019) → 7RP3, 7RP4 (Wang, RCSB, 2021)'), ('XVE', 'XVE: 2XYP (Mittl, PDBe, 2010) → 4Q24 (Gondry, RCSB, 2014)'),
        ('59H', '59H: 5DF7 (Owens, PDBe, 2015) → 8VEQ (Davies, RCSB, 2023)'), ('HUR', 'HUR: 7E18 (Senda) → Zhang, Li (all PDBj, 2021–22)'),
        ('ALD', 'ALD (escaped): 1BP4 (RCSB, 1998) → 4 labs, 2020 on'), ('FFQ', 'FFQ (escaped): 1UAE (RCSB, 1996) → 3 labs')]
fig, ax = plt.subplots(figsize=(7.2, 3.0))
for i, (c, lab) in enumerate(rows):
    g = SU[SU.lig_comp == c].copy(); g['d'] = pd.to_datetime(g.dep)
    u = T[(T.lig_comp == c) & (T.status == 'unbonded')].copy(); u['d'] = pd.to_datetime(u.dep)
    ax.plot([g.d.min(), g.d.max()], [i, i], color='#ddd', lw=1, zorder=0)
    ax.scatter(u.d, [i] * len(u), marker='o', s=14, facecolors='none', edgecolors=GREY, lw=0.7)
    for st, m, col in [('recorded', 'o', REC), ('missed', 'X', MISS), ('mixed', 'D', PDBJ)]:
        q = g[g.status == st]; ax.scatter(q.d, [i] * len(q), marker=m, s=26, color=col, zorder=3)
    f = g.sort_values('d').iloc[0]; ax.scatter([f.d], [i], marker='s', s=70, facecolors='none', edgecolors='k', lw=0.9, zorder=4)
    nb = len(g); nm = int((g.status != 'recorded').sum())
    ax.text(1.01, i, f'{nm}/{nb} bonded uses not recorded', transform=ax.get_yaxis_transform(), va='center', fontsize=6.5)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[1] for r in rows], fontsize=6.5); ax.invert_yaxis(); ax.set_xlabel('deposition date')
ax.legend(handles=[Line2D([], [], marker='o', ls='', color=REC, label='recorded'), Line2D([], [], marker='X', ls='', color=MISS, label='not recorded (close contact only)'),
                   Line2D([], [], marker='o', ls='', mfc='none', mec=GREY, label='use without a nucleophile bond'), Line2D([], [], marker='s', ls='', mfc='none', mec='k', label="code's first bonded use")],
          fontsize=6, frameon=False, ncol=4, loc='upper left', bbox_to_anchor=(0, 1.13))
save(fig, 'fig5_example_codes')

# ---------- Figure 6: hypotheses tested ----------
S = pd.read_csv(f'{R}/ligand_code_effect/step3/attribution_sensitivity.tsv', sep='\t')
fig = plt.figure(figsize=(7.2, 6.4)); gs = fig.add_gridspec(2, 2, hspace=0.6, wspace=0.6, height_ratios=[1, 1.25])
a0 = fig.add_subplot(gs[0, 0]); a1 = fig.add_subplot(gs[1, :]); a2 = fig.add_subplot(gs[0, 1])
lbl = {'drug-like': 'all drug-like', 'drug-like, specific codes (<=10 unbonded uses)': 'specific codes', 'drug-like, bond-like distance': 'bond-like distance', 'drug-like, specific + bond-like': 'specific + bond-like'}
for i, r in S.iterrows():
    for k, col, off in [('unbonded', GREY, -0.15), ('missed', MISS, 0.15)]:
        lo, hi = (float(v) for v in r[f'{k}_CI'].split('-')); a0.plot([lo, hi], [i + off] * 2, color=col, lw=1.3); a0.plot(r[f'{k}_OR'], i + off, 'o', color=col, ms=3.5)
a0.set_xscale('log'); a0.axvline(1, color='k', ls=':', lw=0.6); a0.set_yticks(range(len(S))); a0.set_yticklabels([lbl[s] for s in S.subset], fontsize=6.5); a0.invert_yaxis()
a0.set_xlabel('OR for a later miss vs code born recorded'); a0.set_title('A  "Born unbonded" vs "born on a missed bond"')
a0.legend(handles=[Line2D([], [], color=GREY, marker='o', label='code first defined on an unbonded model'), Line2D([], [], color=MISS, marker='o', label='code first defined on a missed bond')], fontsize=5.8, frameon=False, loc='upper left', bbox_to_anchor=(-0.5, -0.22))
t27 = open(f'{R}/ligand_code_effect/step5/27_stdout.txt').read().split('miss ~ c_first_miss + C(c_first_era)')[1].split('miss ~')[0]
nm = {'C(c_first_era)[T.2014-19]': "code's first use 2014–19", 'C(c_first_era)[T.2020+]': "code's first use 2020+", 'C(c_first_era)[T.<=2007]': "code's first use ≤2007",
      'C(ref, Treatment("PHENIX"))[T.BUSTER]': 'refinement: BUSTER', 'C(ref, Treatment("PHENIX"))[T.CNS]': 'refinement: CNS', 'C(ref, Treatment("PHENIX"))[T.REFMAC]': 'refinement: REFMAC',
      'C(ref, Treatment("PHENIX"))[T.other]': 'refinement: other', 'C(era)[T.2008-13]': 'deposited 2008–13', 'C(era)[T.2014-19]': 'deposited 2014–19', 'C(era)[T.2020+]': 'deposited 2020+',
      'c_first_miss': "code's first bonded use missed", 'pdbj': 'processed at PDBj'}
ys = []
for line in t27.strip().splitlines():
    m = re.match(r'\s*(.+?): OR ([\d.]+) \(([\d.]+)-([\d.]+)\) p=(\S+)', line)
    if m: ys.append((nm[m.group(1)], *map(float, m.group(2, 3, 4)), m.group(5)))
for i, (n, o, lo, hi, p) in enumerate(ys):
    col = MISS if 'missed' in n else DARK if 'PDBj' in n else GREY
    a1.plot([lo, hi], [i, i], color=col, lw=1.2); a1.plot(o, i, 'o', color=col, ms=3.5)
a1.set_yticks(range(len(ys))); a1.set_yticklabels([f'{n}  (p={p})' for n, *_, p in ys], fontsize=6.3); a1.invert_yaxis(); a1.set_xscale('log'); a1.axvline(1, color='k', ls=':', lw=0.6)
a1.set_xlabel('OR for a later miss (1,703 later uses, 115 missed)'); a1.set_title('C  Refinement software and era, adjusted for code history and site')
cl = [('carried\n(first missed,\nlater missed)', 4, 12), ('escaped\n(first missed,\nlater recorded)', 9, 9), ('always\nrecorded', 402, 406)]
for i, (n, k, m) in enumerate(cl):
    a2.bar(i, k / m * 100, 0.6, color=[MISS, REC, GREY][i]); a2.text(i, k / m * 100 + 3, f'{k}/{m}', ha='center', fontsize=6.5)
a2.set_xticks(range(3)); a2.set_xticklabels([c[0] for c in cl], fontsize=6.3); a2.set_ylim(0, 118); a2.set_ylabel('% of codes with PCM today')
a2.set_title('B  PCM annotation today, by code history')
fig.text(0, 0.0, 'A: drug-like units, archive-wide (step 3). B: 479 codes with later strict uses; 3,454 PCM modify dates fall on 2024-09-27. C: strict set; reference levels PHENIX,\n'
         'deposited ≤2007, first use 2008–13; CI clustered by code.', fontsize=6.5, color='#555')
save(fig, 'fig6_hypotheses')
print('figures written to', F)
