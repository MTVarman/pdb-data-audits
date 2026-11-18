"""Recorded fraction by distance, processing site, year; drug subset. Writes out/summary_numbers.txt and bin tables."""
import os, numpy as np, pandas as pd
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'); O = f'{D}/out'
d = pd.read_csv(f'{O}/contacts_all.tsv', sep='\t', low_memory=False, keep_default_na=False, na_values=[''])
d['recorded'] = d.recorded.astype(str).eq('True'); d['nearest_lig_for_cys'] = d.nearest_lig_for_cys.astype(str).eq('True')
d['bonded_elsewhere'] = d.bonded_elsewhere.astype(str).eq('True')
d['process_site'] = d.process_site.fillna('(blank)')
L = []
def P(*a): s = ' '.join(str(x) for x in a); print(s); L.append(s)
def frac(x): return f"{int(x.sum())}/{len(x)} ({100*x.mean():.1f}%)" if len(x) else '0/0'
def band(x):
    return pd.cut(x.dist, [1.6, 2.0, 2.1, 2.4000001], right=False, labels=['1.60-<2.00', '2.00-<2.10', '2.10-2.40'])
base = d[(d.el_class == 'C') & d.nearest_lig_for_cys & d.dist.between(1.6, 2.4)]
pops = {'DRUGLIKE (other_ligand + named_drug)': base[base.category.isin(['other_ligand', 'named_drug'])],
        'named_drug only': base[base.category == 'named_drug'],
        'additive_or_thiol (C partner)': base[base.category == 'additive_or_thiol'],
        'lipid_acyl_coa': base[base.category == 'lipid_acyl_coa'],
        'flavin': base[base.category == 'flavin'],
        'heme_porphyrin (in-sample control)': base[base.category == 'heme_porphyrin'],
        'bilin (in-sample control)': base[base.category == 'bilin'],
        'pyp_chromophore': base[base.category == 'pyp_chromophore']}
P('# Cys SG - ligand C contacts, nearest ligand per Cys, 1.60-2.40 A')
for k, x in pops.items():
    b = band(x)
    P(f'{k}: contacts {len(x)}, entries {x.entry.nunique()} | ' + ' | '.join(f'{lab}: {frac(x.recorded[b == lab])}' for lab in b.cat.categories))
dl = pops['DRUGLIKE (other_ligand + named_drug)']
for lab, x in [('excl. bonded_elsewhere', dl[~dl.bonded_elsewhere]), ('X-ray only', dl[dl.method.fillna('').str.contains('X-RAY')]),
               ('full occupancy both', dl[(dl.lig_occ >= 1) & (dl.cys_occ >= 1)]), ('no altloc either', dl[(dl.cys_alt.fillna('') == '') & (dl.lig_alt.fillna('') == '')]),
               ('same asymmetric unit (sym_image 0)', dl[dl.sym_image == 0])]:
    b = band(x); P(f'  DRUGLIKE sensitivity [{lab}]: n={len(x)} | ' + ' | '.join(f'{c}: {frac(x.recorded[b == c])}' for c in b.cat.categories))
# 0.05 A bins
edges = np.round(np.arange(1.6, 2.4001, 0.05), 2)
def bins(x, w=0.05):
    e = np.round(np.arange(1.6, 2.4001, w), 3)
    c = pd.cut(x.dist, e, right=False, include_lowest=True)
    g = x.groupby(c, observed=False).recorded.agg(['sum', 'count']); g['pct'] = (100 * g['sum'] / g['count']).round(1)
    return g
g = bins(dl); g.to_csv(f'{O}/bins_druglike_0.05A.tsv', sep='\t')
P('\n# DRUGLIKE 0.05 A bins (recorded/count)'); [P(f'  {i}: {int(r["sum"])}/{int(r["count"])}  {r["pct"]}%') for i, r in g.iterrows()]
for k in ['heme_porphyrin (in-sample control)', 'bilin (in-sample control)']:
    bins(pops[k]).to_csv(f'{O}/bins_{k.split()[0]}_0.05A.tsv', sep='\t')
g2 = bins(dl[dl.dist.between(1.9, 2.1)], 0.01) if False else None
x = dl[dl.dist.between(1.9, 2.1)]; e = np.round(np.arange(1.90, 2.1001, 0.01), 2)
gg = x.groupby(pd.cut(x.dist, e, right=False), observed=False).recorded.agg(['sum', 'count'])
P('\n# DRUGLIKE 0.01 A bins around the line'); [P(f'  {i}: {int(r["sum"])}/{int(r["count"])}') for i, r in gg.iterrows()]
P('  exactly 2.000:', frac(dl.recorded[dl.dist == 2.0]), '| 1.990-1.999:', frac(dl.recorded[dl.dist.between(1.99, 1.9995)]), '| 2.001-2.010:', frac(dl.recorded[dl.dist.between(2.0005, 2.01)]))
# processing site
P('\n# DRUGLIKE by processing site')
nl = dl[dl.dist < 2.0]; lg = dl[dl.dist >= 2.0]
tot_miss = (~nl.recorded).sum()
for s, x in dl.groupby('process_site'):
    a = x[x.dist < 2.0]; b2 = x[x.dist >= 2.0]
    P(f'  {s}: <2.0 {frac(a.recorded)}; >=2.0 {frac(b2.recorded)}; share of all <2.0 contacts {100*len(a)/max(len(nl),1):.1f}%; share of <2.0 misses {(~a.recorded).sum()}/{tot_miss}')
site_tab = dl.assign(band=band(dl).astype(str)).groupby(['process_site', 'band']).recorded.agg(['sum', 'count']).reset_index()
site_tab.to_csv(f'{O}/site_druglike.tsv', sep='\t', index=False)
# controls by site too
for k in ['heme_porphyrin (in-sample control)', 'bilin (in-sample control)']:
    x = pops[k]; a = x[x.dist < 2.0]; m = (~a.recorded).sum()
    P(f'\n# {k} normal-length (<2.0) misses by site: total {m}')
    for s, y in a.groupby('process_site'): P(f'  {s}: <2.0 {frac(y.recorded)}; misses {(~y.recorded).sum()}')
# year
P('\n# DRUGLIKE <2.0 recorded by deposition year group')
dy = pd.to_numeric(nl.deposit_year, errors='coerce')
for lab, m in [('<2000', dy < 2000), ('2000-2009', dy.between(2000, 2009)), ('2010-2014', dy.between(2010, 2014)), ('2015-2019', dy.between(2015, 2019)), ('2020-2022', dy.between(2020, 2022)), ('2023-2026', dy >= 2023)]:
    P(f'  {lab}: {frac(nl.recorded[m])}')
P('\n# DRUGLIKE >=2.0 recorded by deposition year group')
dy2 = pd.to_numeric(lg.deposit_year, errors='coerce')
for lab, m in [('<2010', dy2 < 2010), ('2010-2019', dy2.between(2010, 2019)), ('2020-2026', dy2 >= 2020)]:
    P(f'  {lab}: {frac(lg.recorded[m])}')
# drug classes
P('\n# Named drugs (all named-drug contacts, C partner, nearest, 1.6-2.4)')
nd = pops['named_drug only']
for (cl, dr), x in nd.groupby(['drug_class', 'drug']):
    b = band(x); P(f'  {cl} | {dr} ({x.lig_comp.iloc[0]}): n={len(x)} entries={x.entry.nunique()} | ' + ' | '.join(f'{c}: {frac(x.recorded[b == c])}' for c in b.cat.categories))
# S partners
s = d[(d.el_class == 'S') & d.nearest_lig_for_cys & d.dist.between(1.6, 2.4)]
P('\n# Cys SG - ligand S (disulfide-type) contacts, nearest per Cys, by category')
for k, x in s.groupby('category'):
    b = band(x); P(f'  {k}: n={len(x)} | ' + ' | '.join(f'{c}: {frac(x.recorded[b == c])}' for c in b.cat.categories))
x = s; b = pd.cut(x.dist, [1.6, 1.95, 2.0, 2.05, 2.1, 2.2, 2.4000001], right=False)
P('  all S, finer:', ' | '.join(f'{i}: {frac(x.recorded[b == i])}' for i in b.cat.categories))
# below 1.6
lo = d[(d.el_class == 'C') & d.nearest_lig_for_cys & (d.dist < 1.6) & d.category.isin(['other_ligand', 'named_drug'])]
P(f'\n# DRUGLIKE contacts below 1.6 A (not in window): {frac(lo.recorded)}')
E2 = pd.read_csv(f'{O}/entries_scanned.tsv', sep='\t')
P(f'\nEntries scanned: {len(E2)}; with any SG contact <=2.6 A to a non-polymer non-metal atom: {(E2.n_pairs > 0).sum()}; symmetry fallback (ASU only): {(E2.symmetry_checked == False).sum()}')
P('Processing sites of scanned entries:', E2.process_site.fillna('(blank)').value_counts().to_dict())
open(f'{O}/summary_numbers.txt', 'w').write('\n'.join(L) + '\n')
dl.to_csv(f'{O}/contacts_druglike_C_1.6-2.4.tsv', sep='\t', index=False)
