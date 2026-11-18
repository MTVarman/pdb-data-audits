"""Two ways of locating each pair's record line, from the 0.05 A bin table."""
import csv, collections

P = '/home/claude/out/lane3'
bins = list(csv.DictReader(open(f'{P}_record_line_bins.csv')))
for b in bins:
    b['lo'] = float(b['bin_lo']); b['hi'] = float(b['bin_hi'])
    b['n'] = int(b['n']); b['pct'] = float(b['recorded_pct']); b['line'] = float(b['probe_line_A'])
g = collections.defaultdict(list)
for b in bins:
    g[(b['set'], b['residue'], b['pair'])].append(b)

print(f'{"set":14s} {"res":11s} {"pair":5s} {"probe":>6s} {"collapse":>10s} {"last>=50%":>10s} '
      f'{"collapse-probe":>15s}  n   rec')
rows = []
for k, v in sorted(g.items()):
    v.sort(key=lambda x: x['lo'])
    big = [x for x in v if x['n'] >= 3]
    if not big:
        continue
    # collapse point: the first bin (n>=3) below 50% that has no later n>=3 bin at >=80%
    collapse = None
    for i, x in enumerate(big):
        if x['pct'] < 50 and not any(y['pct'] >= 80 for y in big[i + 1:]):
            collapse = x
            break
    last50 = None
    for i, x in enumerate(big):
        if x['pct'] >= 50 and all(y['pct'] < 50 for y in big[i + 1:]):
            last50 = x
    n = sum(x['n'] for x in v)
    rec = sum(x['n'] * x['pct'] / 100 for x in v)
    line = v[0]['line']
    cs = f'{collapse["lo"]:.2f}' if collapse else 'none'
    ls = f'{last50["lo"]:.2f}-{last50["hi"]:.2f}' if last50 else 'none'
    diff = f'{collapse["lo"]-line:+.2f}' if collapse else ''
    print(f'{k[0]:14s} {k[1]:11s} {k[2]:5s} {line:6.2f} {cs:>10s} {ls:>10s} {diff:>15s} {n:5d} {rec:6.0f}')
    rows.append({'set': k[0], 'residue': k[1], 'pair': k[2], 'probe_line_A': line,
                 'record_collapse_bin_lo': collapse['lo'] if collapse else None,
                 'record_collapse_minus_probe_line': round(collapse['lo'] - line, 2) if collapse else None,
                 'last_bin_at_or_above_50pct': ls, 'n': n, 'recorded': round(rec)})
with open(f'{P}_record_lines_two_ways.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
