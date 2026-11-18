"""Tabulate the allowed-bond distance windows in each MAXIT source release.
Usage: python3 -I maxit_table.py versions_dir out.csv
versions_dir/<version>/... must hold BondUtil_global.h (v10+) or _allowed_bond.h (v8) unpacked from
https://sw-tools.rcsb.org/apps/MAXIT/maxit-v<version>-prod-src.tar.gz"""
import re, glob, sys, csv
REL = {'8.120': '2008-02-12', '8.125': '2013-03-01', '8.200': '2015-03-18', '10.000': '2017-06-07', '10.100': '2019-12-04',
       '10.200': '2020-07-29', '11.000': '2020-07-29', '11.100': '2020-08-25', '11.200': '2024-04-21', '11.300': '2024-10-11', '11.400': '2025-10-16'}
PAIRS = ['C_S', 'C_O', 'C_N', 'B_O', 'B_N', 'S_S', 'O_P', 'N_P', 'O_S', 'N_S', 'C_SE', 'O_O', 'N_O', 'N_N', 'C_C']
out = []
for v, date in REL.items():
    fs = glob.glob(f'{sys.argv[1]}/{v}/**/BondUtil_global.h', recursive=True) + glob.glob(f'{sys.argv[1]}/{v}/**/_allowed_bond.h', recursive=True)
    t = open(fs[0], errors='replace').read()
    d = {m.group(1): (m.group(2), m.group(3), m.group(4)) for m in re.finditer(r'\{\s*"([A-Z_]+)",\s*"\w+",\s*"\w+",\s*([0-9.]+),\s*([0-9.]+),\s*(\d)', t)}
    for p in PAIRS:
        lo, hi, ty = d.get(p, ('', '', ''))
        out.append({'maxit_release': v, 'release_date': date, 'source_file': fs[0].split('/')[-1], 'rows_in_table': len(d), 'pair': p, 'lower_limit_A': lo, 'upper_limit_A': hi, 'type': ty})
with open(sys.argv[2], 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
for p in PAIRS:
    vals = {(o['lower_limit_A'], o['upper_limit_A']) for o in out if o['pair'] == p}
    print(p, sorted(vals), 'identical in all releases' if len(vals) == 1 else 'DIFFERS')
