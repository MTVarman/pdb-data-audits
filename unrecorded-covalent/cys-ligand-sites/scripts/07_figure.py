import os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'); O = f'{D}/out'
a = pd.read_csv(f'{O}/contacts_all.tsv', sep='\t', low_memory=False, keep_default_na=False, na_values=[''])
for c in ['recorded', 'nearest_lig_for_cys']: a[c] = a[c].astype(str).eq('True')
a['process_site'] = a.process_site.fillna('(blank)')
base = a[a.nearest_lig_for_cys & a.dist.between(1.6, 2.4)]
dl = base[(base.el_class == 'C') & base.category.isin(['other_ligand', 'named_drug'])]
hc = base[(base.el_class == 'C') & (base.category == 'heme_porphyrin')]
ss = base[(base.el_class == 'S') & base.category.isin(['other_ligand', 'additive_or_thiol', 'named_drug'])]
w = 0.05; e = np.round(np.arange(1.6, 2.4001, w), 3); cen = e[:-1] + w / 2
def g(x):
    t = x.groupby(pd.cut(x.dist, e, right=False), observed=False).recorded.agg(['sum', 'count']); t['pct'] = 100 * t['sum'] / t['count'].replace(0, np.nan); return t
G, H, S = g(dl), g(hc), g(ss)
BLUE, ORANGE, INK, GREY, PURP = '#2b6fd6', '#e8642c', '#111', '#777', '#8a5cc2'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'axes.spines.top': False, 'axes.spines.right': False})
fig = plt.figure(figsize=(9, 12), facecolor='#fbfaf7')
lo = dl[dl.dist < 2.0]; hi = dl[dl.dist >= 2.0]; b1 = dl[(dl.dist >= 2.0) & (dl.dist < 2.1)]; b2 = dl[dl.dist >= 2.1]
fig.text(0.07, 0.979, 'Unrecorded covalent bonds to cysteine,\nby wwPDB processing site', fontsize=19, weight='bold', color=INK, linespacing=1.15, va='top')
fig.text(0.07, 0.901, f'Cys SG to ligand carbon, {len(dl):,} contacts in {dl.entry.nunique():,} entries\n(inhibitors, drugs, substrates; heme, bilins, flavins, retinal, PYP, additives, lipids, metals excluded)', fontsize=10.5, color='#444')
ax1 = fig.add_axes([0.10, 0.625, 0.84, 0.255], facecolor='#fbfaf7')
ax1.axvspan(1.6, 2.0, color=BLUE, alpha=0.06, lw=0); ax1.axvspan(2.0, 2.4, color=ORANGE, alpha=0.06, lw=0)
ok = H['count'] > 0
ax1.step(e[:-1][ok.values], H.pct[ok].values, where='post', color=GREY, lw=1.2, alpha=0.8, label=f'heme c / porphyrin, same sample (n={len(hc):,})')
ok = S['count'] > 4
ax1.plot(cen[ok.values], S.pct[ok].values, ls='--', marker='s', ms=4, color=PURP, lw=1.2, label=f'ligand S (disulfide-type), n={len(ss):,}; bins n≥5')
ok = G['count'] > 0
ax1.step(np.r_[e[:-1][ok.values], 2.4], np.r_[G.pct[ok].values, G.pct[ok].values[-1]], where='post', color=INK, lw=2.2, label='drug-like ligands (C)')
ax1.scatter(cen[ok.values], G.pct[ok], s=np.clip(G['count'][ok] / 5, 10, 220), c=[BLUE if x < 2.0 else ORANGE for x in cen[ok.values]], zorder=3, edgecolor='white', lw=0.6)
for x, (s_, n) in zip(cen, zip(G['sum'], G['count'])):
    if n: ax1.text(x, -11, f'{int(n)}', ha='center', fontsize=7, color=GREY)
ax1.axvline(2.0, ls='--', color=INK, lw=1.3); ax1.set_xlim(1.6, 2.4); ax1.set_ylim(-16, 112)
ax1.set_ylabel('% recorded as covale in _struct_conn', fontsize=10.5); ax1.set_xlabel('Cys SG to ligand atom distance (Å); 0.05 Å bins, contacts per bin in grey', fontsize=10)
ax1.text(1.79, 50, f'{100*lo.recorded.mean():.1f}% recorded\n{lo.recorded.sum():,} of {len(lo):,}', ha='center', fontsize=14, weight='bold', color=BLUE)
ax1.text(2.205, 93, f'2.00–2.10 Å: {b1.recorded.sum()} of {len(b1)} ({100*b1.recorded.mean():.0f}%)\n2.10–2.40 Å: {b2.recorded.sum()} of {len(b2)} ({100*b2.recorded.mean():.0f}%)', ha='center', fontsize=11.5, weight='bold', color=ORANGE)
ax1.legend(loc='lower left', fontsize=8.5, frameon=False, bbox_to_anchor=(0.0, 0.08))
ax1.set_title('A  Share of contacts recorded as a bond (dot size = contacts)', loc='left', fontsize=11, color=INK)
# panel B
ax2 = fig.add_axes([0.10, 0.315, 0.84, 0.205], facecolor='#fbfaf7')
sites = ['RCSB', 'PDBE', 'PDBC', 'PDBJ', 'BNL']
rows = []
for s_ in sites:
    x = lo[lo.process_site == s_]; y = hi[hi.process_site == s_]
    ent = x.groupby('entry').recorded.agg(lambda r: (~r).any())
    rows.append((s_, x.recorded.sum(), len(x), (~x.recorded).sum(), ent.sum(), len(ent), y.recorded.sum(), len(y)))
X = np.arange(len(sites)); bw = 0.38
p1 = [100 * r[1] / r[2] if r[2] else np.nan for r in rows]; p2 = [100 * r[6] / r[7] if r[7] else np.nan for r in rows]
ax2.bar(X - bw / 2, p1, bw, color=BLUE, label='contacts < 2.0 Å recorded')
ax2.bar(X + bw / 2, p2, bw, color=ORANGE, label='contacts 2.0–2.4 Å recorded')
for i, r in enumerate(rows):
    ax2.text(i - bw / 2, p1[i] + 1.5, f'{r[1]}/{r[2]}', ha='center', fontsize=8.5, color=BLUE)
    if r[7]: ax2.text(i + bw / 2, p2[i] + 1.5, f'{r[6]}/{r[7]}', ha='center', fontsize=8.5, color=ORANGE)
    ax2.text(i, -24, f'{r[3]} normal-length misses\nentries hit: {r[4]}/{r[5]} ({100*r[4]/max(r[5],1):.0f}%)', ha='center', fontsize=8, color=INK)
ax2.set_xticks(X); ax2.set_xticklabels(sites, fontsize=11, weight='bold'); ax2.tick_params(axis='x', pad=4)
ax2.set_ylim(0, 125); ax2.set_yticks([0,20,40,60,80,100]); ax2.set_ylabel('% recorded in _struct_conn', fontsize=10.5)
ax2.legend(loc='upper center', fontsize=9.5, frameon=False, ncol=2)
pj = [r for r in rows if r[0] == 'PDBJ'][0]; tot = sum(r[3] for r in rows)
ax2.set_title(f'B  By wwPDB processing site: PDBj has {100*pj[2]/len(lo):.0f}% of contacts under 2.0 Å\n    but {pj[3]} of the {tot} unrecorded ones', loc='left', fontsize=11, color=INK)

# panel C: the two matched cases, chemistry held fixed
ax3 = fig.add_axes([0.10, 0.100, 0.84, 0.095]); ax3.axis('off')
def mrow(df, comp, site):
    x = lo[(lo.lig_comp == comp) & (lo.process_site == site)]
    return int(x.recorded.sum()), len(x)
cases = [('Nirmatrelvir (4WI) on Mpro Cys145', [('RCSB',) + mrow(lo, '4WI', 'RCSB'), ('PDBC',) + mrow(lo, '4WI', 'PDBC'),
                                               ('PDBe',) + mrow(lo, '4WI', 'PDBE'), ('PDBj',) + mrow(lo, '4WI', 'PDBJ')]),
         ('Ubiquinone-0 (UQ0) on glyoxylate carboligase Cys588', [('PDBe 8BEO',) + mrow(lo, 'UQ0', 'PDBE'),
                                                                  ('PDBj 8I01/05/08',) + mrow(lo, 'UQ0', 'PDBJ')])]
ax3.text(0, 1.12, 'C  Same ligand, same residue, same protein — recorded or not by site', fontsize=11, color=INK, transform=ax3.transAxes)
for k, (lab, cells) in enumerate(cases):
    y = 0.80 - k * 0.50
    ax3.text(0, y, lab, fontsize=9.5, color=INK, weight='bold', transform=ax3.transAxes)
    for j, (nm, r_, n_) in enumerate(cells):
        col = ORANGE if n_ and r_ / n_ < 0.5 else BLUE
        ax3.text(0.02 + j * 0.15, y - 0.26, f'{nm}  {r_}/{n_}', fontsize=8.5, color=col,
                 weight='bold', transform=ax3.transAxes)
fig.text(0.07, 0.048, 'Sample: RCSB full-text "covalent cysteine" or "covalent inhibitor", plus entries with named covalent drugs:\n24,144 entries, mmCIF downloaded 14 Nov 2026. Nearest ligand per Cys; heavy atoms, occupancy > 0, matching alt confs,\nsymmetry-aware. Recorded = a covale row joining that Cys and that ligand residue.', fontsize=8.5, color=GREY)
fig.text(0.07, 0.013, 'Michael T. Varman · Morphogenesis Analytics · morphalytica.com', fontsize=10.5, color=INK, weight='bold')
fig.savefig(f'{O}/cys_covalent_two_rulers.png', dpi=200, facecolor=fig.get_facecolor())
fig.savefig(f'{O}/cys_covalent_two_rulers.pdf', facecolor=fig.get_facecolor())
fig.savefig(f'{O}/cys_covalent_two_rulers.svg', facecolor=fig.get_facecolor())
print(rows)
