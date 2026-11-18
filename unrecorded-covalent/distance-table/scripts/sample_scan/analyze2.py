"""Evidence-sample tables.

Usage: python3 -I analyze2.py out_prefix
Reads <pre>_contacts_all.csv and <pre>_window_evidence.csv.
Writes the binned tables, the classified unrecorded-site list, and prints every
number that goes into the report.
"""
import collections, csv, json, statistics, sys

BAND = 0.40
RESIDUES = ['SER', 'THR', 'TYR', 'LYS', 'HIS', 'ASP', 'GLU', 'SEC', 'CYS']
HB = {('O', 'O'), ('O', 'N'), ('N', 'O'), ('N', 'N')}


def num(v, d=None):
    try:
        return float(v)
    except Exception:
        return d


def tf(v):
    return str(v).strip().lower() in ('true', '1', 'yes')


def load(p):
    rows = list(csv.DictReader(open(p)))
    for r in rows:
        for k in ('dist_A', 'probe_line_A', 'offset_from_line_A', 'resolution'):
            r[k] = num(r[k])
        for k in ('prot_occ', 'lig_occ'):
            r[k] = num(r[k], 0)
        for k in ('recorded_covale', 'in_band', 'past_line', 'residue_already_bonded',
                  'hbond_capable_pair', 'same_comp_recorded_covalent_in_this_entry',
                  'clash_listed', 'in_report', 'clash_elsewhere_same_residue_pair'):
            if k in r:
                r[k] = tf(r[k])
    return rows


def verdict(r):
    """What the contact most plausibly is. Conservative: only the last class is
    presented as a candidate unrecorded covalent bond."""
    if r['hbond_capable_pair']:
        return 'hydrogen_bond_or_salt_bridge'
    if r['family'] == 'additive_cryoprotectant':
        return 'cryoprotectant_or_additive'
    if r['family'] == 'class1_cofactor_ptm':
        return 'cofactor_or_PTM_chemistry'
    if r['family'] == 'class3_nucleotide_phosphate':
        return 'phosphoryl_transfer_or_nucleotide'
    if r['residue_already_bonded']:
        return 'residue_already_bonded_through_another_atom'
    if r['ligand_atom_tier'] == 'C':
        return 'no_bond_as_modelled_unreactive_atom'
    if r['ligand_atom_tier'] in ('A', 'B') and r['lig_el'] in ('C', 'P', 'S', 'B'):
        if r['ligand_atom_kind'].startswith('C=') or r['ligand_atom_kind'] in ('C#N / C#C',):
            return 'candidate_bond_or_near_attack_complex'
        return 'candidate_unrecorded_covalent_bond'
    return 'other'


def rec_dist_table(allr, pre):
    print('\n' + '=' * 100)
    print('RECORDED (_struct_conn covale) BOND LENGTHS vs PROBE\'S BOND LINE, by element pair')
    print('  how much headroom there is between the way the bond is actually modelled')
    print('  and the distance at which the validator stops perceiving it')
    rows = []
    print(f'{"pair":6s} {"line":>6s} {"n_rec":>7s} {"min":>6s} {"med":>6s} {"p95":>6s} {"p99":>6s} '
          f'{"max":>6s} {"%>=line":>8s}')
    for pair in sorted({r['pair'] for r in allr if r['recorded_covale']}):
        b = [r['dist_A'] for r in allr if r['pair'] == pair and r['recorded_covale']]
        if len(b) < 5:
            continue
        line = next(r['probe_line_A'] for r in allr if r['pair'] == pair)
        b.sort()
        p = lambda f: b[min(len(b) - 1, int(f * len(b)))]
        over = sum(1 for x in b if x >= line)
        rows.append({'pair': pair, 'probe_line_A': line, 'n_recorded': len(b), 'min': b[0],
                     'median': round(statistics.median(b), 3), 'p95': p(0.95), 'p99': p(0.99),
                     'max': b[-1], 'pct_recorded_at_or_past_line': round(100 * over / len(b), 2),
                     'headroom_median_to_line': round(line - statistics.median(b), 3)})
        print(f'{pair:6s} {line:6.2f} {len(b):7d} {b[0]:6.3f} {statistics.median(b):6.3f} '
              f'{p(0.95):6.3f} {p(0.99):6.3f} {b[-1]:6.3f} {100*over/len(b):8.2f}')
    with open(f'{pre}_recorded_bond_lengths.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return rows


def main():
    pre = sys.argv[1]
    allr = load(f'{pre}_contacts_all.csv')
    win = load(f'{pre}_window_evidence.csv')
    for r in allr:
        r['verdict'] = verdict(r)
    for r in win:
        r['verdict'] = verdict(r)
    ev = {(r['pdbid'], r['ligand'], r['protein_atom']): r for r in win}
    usable = [r for r in allr if r['occupancy_class'] == 'full' and r['prot_occ'] > 0 and r['lig_occ'] > 0]

    print('=' * 100)
    print('CONTACT INVENTORY  (closest ligand atom per entry / ligand instance / nucleophile atom,')
    print('                    organic non-polymer or branched ligand, model 1, image 0, d < 2.6 A)')
    print(f'rows {len(allr)}  entries {len({r["pdbid"] for r in allr})}')
    print(f'full occupancy, no altloc crossing: {len(usable)} rows, {len({r["pdbid"] for r in usable})} entries')
    print(f'\n{"res":5s} {"rows":>8s} {"entries":>8s} {"recorded":>9s} {"rec%":>6s}')
    inv = []
    for res in RESIDUES:
        b = [r for r in usable if r['residue'] == res]
        if b:
            rec = sum(r['recorded_covale'] for r in b)
            print(f'{res:5s} {len(b):8d} {len({r["pdbid"] for r in b}):8d} {rec:9d} {100*rec/len(b):6.1f}')
            inv.append({'residue': res, 'rows': len(b), 'entries': len({r['pdbid'] for r in b}),
                        'recorded': rec, 'recorded_pct': round(100 * rec / len(b), 1)})
    with open(f'{pre}_inventory.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(inv[0].keys()))
        w.writeheader()
        w.writerows(inv)

    rec_dist_table(usable, pre)

    # ---- 1. in-band counts ----
    print('\n' + '=' * 100)
    print('1. IN-BAND CONTACTS (probe line <= d < line + 0.40 A), recorded vs not; full occupancy')
    band_tab = []
    print(f'{"residue":8s} {"pair":6s} {"line":>5s} {"band":>11s} {"n":>6s} {"rec":>5s} {"unrec":>6s} '
          f'{"rec%":>6s} {"entries":>8s}')
    for res in RESIDUES:
        for pair in sorted({r['pair'] for r in usable if r['residue'] == res and r['in_band']}):
            b = [r for r in usable if r['residue'] == res and r['pair'] == pair and r['in_band']]
            line = b[0]['probe_line_A']
            rec = sum(r['recorded_covale'] for r in b)
            band_tab.append({'residue': res, 'pair': pair, 'probe_line_A': line,
                             'band': f'{line:.2f}-{line+BAND:.2f}', 'n': len(b), 'recorded': rec,
                             'unrecorded': len(b) - rec, 'recorded_pct': round(100 * rec / len(b), 1),
                             'entries': len({r['pdbid'] for r in b}),
                             'hbond_capable': b[0]['hbond_capable_pair']})
            print(f'{res:8s} {pair:6s} {line:5.2f} {line:5.2f}-{line+BAND:<5.2f} {len(b):6d} {rec:5d} '
                  f'{len(b)-rec:6d} {100*rec/len(b):6.1f} {len({r["pdbid"] for r in b}):8d}'
                  + ('  [H-bond-capable]' if b[0]['hbond_capable_pair'] else ''))
    with open(f'{pre}_band_counts.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(band_tab[0].keys()))
        w.writeheader()
        w.writerows(band_tab)

    # ---- 2. 0.05 A bins ----
    print('\n' + '=' * 100)
    print('2. RECORDED FRACTION BY DISTANCE, 0.05 A bins. Non-H-bond-capable element pairs only')
    print('   ("record line" = highest 0.05 A bin with n>=3 at or above 50% recorded, with every')
    print('    n>=3 bin beyond it below 50%)')
    bins, lines_tab = [], []

    def do_bins(label, data, res, pair):
        line = data[0]['probe_line_A']
        seq = []
        lo = 1.10
        while lo < 2.60:
            g = [r for r in data if lo <= r['dist_A'] < lo + 0.05]
            if g:
                rec = sum(r['recorded_covale'] for r in g)
                row = {'set': label, 'residue': res, 'pair': pair, 'probe_line_A': line,
                       'bin_lo': round(lo, 2), 'bin_hi': round(lo + 0.05, 2), 'n': len(g),
                       'recorded': rec, 'recorded_pct': round(100 * rec / len(g), 1),
                       'contains_probe_line': lo <= line < lo + 0.05}
                seq.append(row)
                bins.append(row)
            lo = round(lo + 0.05, 2)
        cross = None
        big = [x for x in seq if x['n'] >= 3]
        for i, x in enumerate(big):
            if x['recorded_pct'] >= 50 and all(y['recorded_pct'] < 50 for y in big[i + 1:]):
                cross = x
        n = len(data)
        rec = sum(r['recorded_covale'] for r in data)
        past = [r for r in data if r['past_line']]
        print(f'\n  {label}  {res} {pair}   probe line {line:.2f} A   n={n}  recorded {rec} ({100*rec/n:.1f}%)')
        print(f'    {"bin":>12s} {"n":>7s} {"rec":>6s} {"rec%":>6s}')
        for x in seq:
            mark = '  <== PROBE LINE' if x['contains_probe_line'] else ''
            if cross and x is cross:
                mark += '  <== RECORD LINE'
            print(f'    {x["bin_lo"]:.2f}-{x["bin_hi"]:.2f} {x["n"]:7d} {x["recorded"]:6d} '
                  f'{x["recorded_pct"]:6.1f}{mark}')
        rl = f'{cross["bin_lo"]:.2f}-{cross["bin_hi"]:.2f}' if cross else 'none (no bin reaches 50%)'
        print(f'    record line {rl}   |  past the probe line: {len(past)} rows, '
              f'{sum(r["recorded_covale"] for r in past)} recorded')
        lines_tab.append({'set': label, 'residue': res, 'pair': pair, 'probe_line_A': line,
                          'n': n, 'recorded': rec, 'recorded_pct': round(100 * rec / n, 1),
                          'record_line_bin': rl,
                          'record_line_hi': cross['bin_hi'] if cross else None,
                          'record_minus_probe_line': round(cross['bin_hi'] - line, 2) if cross else None,
                          'n_past_probe_line': len(past),
                          'recorded_past_probe_line': sum(r['recorded_covale'] for r in past)})

    for res in RESIDUES:
        for pair in sorted({r['pair'] for r in usable if r['residue'] == res}):
            b = [r for r in usable if r['residue'] == res and r['pair'] == pair]
            if len(b) < 25 or (b[0]['lig_el'], b[0]['prot_el']) in HB:
                continue
            do_bins('per-residue', b, res, pair)

    print('\n  POOLED over all non-Cys nucleophiles:')
    for pair in sorted({r['pair'] for r in usable if r['residue'] != 'CYS'}):
        b = [r for r in usable if r['pair'] == pair and r['residue'] != 'CYS']
        if len(b) < 25 or (b[0]['lig_el'], b[0]['prot_el']) in HB:
            continue
        do_bins('pooled-nonCys', b, 'ALL_nonCys', pair)
    for pair in sorted({r['pair'] for r in usable if r['residue'] == 'CYS'}):
        b = [r for r in usable if r['pair'] == pair and r['residue'] == 'CYS']
        if len(b) < 25:
            continue
        do_bins('cys-control', b, 'CYS', pair)

    with open(f'{pre}_record_line_bins.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(bins[0].keys()))
        w.writeheader()
        w.writerows(bins)
    with open(f'{pre}_record_lines.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(lines_tab[0].keys()))
        w.writeheader()
        w.writerows(lines_tab)

    # ---- 3/4 ----
    print('\n' + '=' * 100)
    print('3/4. UNRECORDED IN-BAND CONTACTS (non-Cys, full occupancy): what they are,')
    print('     density support, and the clash flip')
    print('     density thresholds fixed before any report was read:')
    print('     both partners RSCC >= 0.80 and RSR <= 0.30 and occupancy > 0')
    unrec = [r for r in usable if r['in_band'] and not r['recorded_covale'] and r['residue'] != 'CYS']
    print(f'\n  total {len(unrec)} rows, {len({r["pdbid"] for r in unrec})} entries')
    print('\n  by what the contact most plausibly is:')
    vt = []
    print(f'  {"verdict":46s} {"rows":>6s} {"entries":>8s} {"inrep":>6s} {"densOK":>7s} {"clash":>6s}')
    for v, n in collections.Counter(r['verdict'] for r in unrec).most_common():
        b = [r for r in unrec if r['verdict'] == v]
        e = [ev.get((r['pdbid'], r['ligand'], r['protein_atom'])) for r in b]
        e = [x for x in e if x and x['in_report']]
        dens = [x for x in e if x['density_class'] == 'density_supported']
        cl = [x for x in e if x['clash_listed']]
        vt.append({'verdict': v, 'rows': n, 'entries': len({r['pdbid'] for r in b}),
                   'in_validation_report': len(e), 'density_supported': len(dens),
                   'density_supported_entries': len({x['pdbid'] for x in dens}),
                   'clash_outlier_at_pair': len(cl)})
        print(f'  {v:46s} {n:6d} {len({r["pdbid"] for r in b}):8d} {len(e):6d} {len(dens):7d} {len(cl):6d}')
    with open(f'{pre}_verdict_summary.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(vt[0].keys()))
        w.writeheader()
        w.writerows(vt)

    print('\n  by residue and element pair:')
    rows3 = []
    print(f'  {"res":5s} {"pair":6s} {"unrec":>6s} {"inrep":>6s} {"densOK":>7s} {"clash":>6s} {"clash%":>7s} {"entries":>8s}')
    for res in RESIDUES:
        for pair in sorted({r['pair'] for r in unrec if r['residue'] == res}):
            b = [r for r in unrec if r['residue'] == res and r['pair'] == pair]
            e = [ev.get((r['pdbid'], r['ligand'], r['protein_atom'])) for r in b]
            e = [x for x in e if x and x['in_report']]
            dens = [x for x in e if x['density_class'] == 'density_supported']
            cl = [x for x in e if x['clash_listed']]
            rows3.append({'residue': res, 'pair': pair, 'unrecorded_in_band': len(b),
                          'in_validation_report': len(e), 'density_supported': len(dens),
                          'clash_outlier_at_pair': len(cl),
                          'clash_pct_of_in_report': round(100 * len(cl) / len(e), 1) if e else None,
                          'entries': len({r['pdbid'] for r in b}),
                          'density_supported_entries': len({x['pdbid'] for x in dens})})
            print(f'  {res:5s} {pair:6s} {len(b):6d} {len(e):6d} {len(dens):7d} {len(cl):6d} '
                  f'{(100*len(cl)/len(e) if e else 0):7.1f} {len({r["pdbid"] for r in b}):8d}')
    with open(f'{pre}_unrecorded_evidence_summary.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows3[0].keys()))
        w.writeheader()
        w.writerows(rows3)

    tot = [ev.get((r['pdbid'], r['ligand'], r['protein_atom'])) for r in unrec]
    tot = [x for x in tot if x and x['in_report']]
    ds = [x for x in tot if x['density_class'] == 'density_supported']
    print(f'\n  in a wwPDB validation report with both partners present: {len(tot)}')
    print(f'  DENSITY-SUPPORTED: {len(ds)} rows, {len({x["pdbid"] for x in ds})} entries')
    print('  density class:', dict(collections.Counter(x['density_class'] for x in tot)))
    print(f'  clash outlier at the flagged atom pair: {sum(1 for x in tot if x["clash_listed"])}')
    print(f'  clash somewhere between the same residue pair: '
          f'{sum(1 for x in tot if x["clash_elsewhere_same_residue_pair"])}')
    print('\n  density-supported, excluding H-bond-capable pairs and additives:')
    dsc = [x for x in ds if x['verdict'] not in ('hydrogen_bond_or_salt_bridge', 'cryoprotectant_or_additive')]
    print(f'    {len(dsc)} rows, {len({x["pdbid"] for x in dsc})} entries;',
          dict(collections.Counter(x['verdict'] for x in dsc)))

    # clash flip vs offset
    print('\n  CLASH FLIP: unrecorded contacts by offset from the probe line (report present, full occ)')
    flip = []
    for label, data in (('non-Cys, non-H-bond pairs',
                         [x for x in win if x['residue'] != 'CYS' and not x['recorded_covale']
                          and x['occupancy_class'] == 'full' and x['in_report']
                          and not x['hbond_capable_pair']]),
                        ('non-Cys, H-bond-capable pairs',
                         [x for x in win if x['residue'] != 'CYS' and not x['recorded_covale']
                          and x['occupancy_class'] == 'full' and x['in_report'] and x['hbond_capable_pair']]),
                        ('Cys control',
                         [x for x in win if x['residue'] == 'CYS' and not x['recorded_covale']
                          and x['occupancy_class'] == 'full' and x['in_report']])):
        print(f'\n   {label}: n={len(data)}')
        lo = -0.45
        while lo < 0.45:
            g = [x for x in data if lo <= x['offset_from_line_A'] < lo + 0.05]
            if g:
                c = sum(1 for x in g if x['clash_listed'])
                print(f'    offset {lo:+.2f}..{lo+0.05:+.2f} {len(g):6d}  clash {c:5d}  {100*c/len(g):6.1f}%')
                flip.append({'set': label, 'offset_lo': round(lo, 2), 'offset_hi': round(lo + 0.05, 2),
                             'n': len(g), 'clash_listed': c, 'clash_pct': round(100 * c / len(g), 1)})
            lo = round(lo + 0.05, 2)
        bel = [x for x in data if x['offset_from_line_A'] < 0]
        ab = [x for x in data if x['offset_from_line_A'] >= 0]
        print(f'    below the line {sum(1 for x in bel if x["clash_listed"])}/{len(bel)}   '
              f'at or past the line {sum(1 for x in ab if x["clash_listed"])}/{len(ab)}')
    with open(f'{pre}_clash_flip.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(flip[0].keys()))
        w.writeheader()
        w.writerows(flip)

    # ---- 5. sites ----
    print('\n' + '=' * 100)
    print('5. SITES')
    ds.sort(key=lambda x: (x['residue'], x['dist_A']))
    with open(f'{pre}_density_supported_sites.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(ds[0].keys()))
        w.writeheader()
        w.writerows(ds)
    cand = [x for x in ds if x['verdict'] in ('candidate_unrecorded_covalent_bond',
                                              'candidate_bond_or_near_attack_complex',
                                              'cofactor_or_PTM_chemistry',
                                              'phosphoryl_transfer_or_nucleotide')]
    cand.sort(key=lambda x: (x['residue'], x['dist_A']))
    with open(f'{pre}_candidate_sites.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(cand[0].keys()) if cand else ['none'])
        w.writeheader()
        w.writerows(cand)
    print(f'  density-supported sites: {len(ds)} rows / {len({x["pdbid"] for x in ds})} entries')
    print(f'  of those, chemistry-plausible (cofactor/PTM, phosphoryl, candidate bond or NAC): '
          f'{len(cand)} rows / {len({x["pdbid"] for x in cand})} entries')
    print('\n  top components among the chemistry-plausible density-supported set:')
    for comp, n in collections.Counter(x['ligand_comp'] for x in cand).most_common(40):
        nm = next(x['ligand_name'] for x in cand if x['ligand_comp'] == comp)
        print(f'    {comp:8s} {n:4d}  {nm[:72]}')
    print('\n  every chemistry-plausible density-supported site:')
    for x in cand:
        print(f"    {x['pdbid']} {x['residue']:4s} {x['protein_atom']:20s} {x['ligand_comp']:7s} "
              f"{x['ligand_atom']:5s} {x['pair']:5s} d={x['dist_A']:.3f} line={x['probe_line_A']:.2f} "
              f"res={x['resolution']} clash={x['clash_listed']} v={x['verdict'][:28]} | {x['title'][:58]}")


if __name__ == '__main__':
    main()
