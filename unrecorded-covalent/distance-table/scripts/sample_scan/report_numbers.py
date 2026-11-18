"""Consolidated counts of the evidence-sample scan, as first reported."""
import collections, csv, json, statistics

P = '/home/claude/out/lane3'
arms = json.load(open('/home/claude/work/lane3/arms.json'))
run = json.load(open(f'{P}_window_evidence_run.json'))


def tf(v):
    return str(v).strip().lower() == 'true'


import sys
sys.path.insert(0, '/home/claude/work/lane3')
from analyze2 import verdict as _verdict

allr = list(csv.DictReader(open(f'{P}_contacts_all.csv')))
win = list(csv.DictReader(open(f'{P}_window_evidence.csv')))
for r in allr + win:
    r['hbond_capable_pair'] = tf(r['hbond_capable_pair'])
    r['residue_already_bonded'] = tf(r['residue_already_bonded'])
    r['verdict'] = _verdict(r)
for r in allr + win:
    r['d'] = float(r['dist_A'])
    r['line'] = float(r['probe_line_A'])
    r['rec'] = tf(r['recorded_covale'])
    r['ib'] = tf(r['in_band'])
    r['hb'] = tf(r['hbond_capable_pair'])
for r in win:
    r['cl'] = tf(r['clash_listed'])
    r['ir'] = tf(r['in_report'])
usable = [r for r in allr if r['occupancy_class'] == 'full']
ev = {(r['pdbid'], r['ligand'], r['protein_atom']): r for r in win}

print('## SELECTION')
print('arms:', arms['counts'], 'union', arms['union'])
print('retrieved_utc', arms['retrieved_utc'], 'baseline pool', arms['baseline_pool_size'],
      'stride', arms['baseline_stride'])
print('warhead counts', arms['warhead_counts'])
print('excluded warheads', arms['excluded_warheads'])
print('report failures', run['report_failures'])

print('\n## SCALE')
print('entries parsed 48072 (all ok). contact rows', len(allr), 'entries', len({r['pdbid'] for r in allr}))
print('full-occupancy rows', len(usable), 'entries', len({r['pdbid'] for r in usable}))
print('window rows', len(win), 'in report', sum(r['ir'] for r in win))

print('\n## PER-RESIDUE, IN BAND, NON-H-BOND PAIRS ONLY')
tot = collections.Counter()
for res in ['SER', 'THR', 'TYR', 'LYS', 'HIS', 'ASP', 'GLU', 'SEC']:
    b = [r for r in usable if r['residue'] == res and r['ib'] and not r['hb']]
    rec = sum(r['rec'] for r in b)
    tot['n'] += len(b)
    tot['rec'] += rec
    print(f'{res:4s} n={len(b):4d} rec={rec:4d} unrec={len(b)-rec:4d} '
          f'rec%={100*rec/len(b) if b else 0:5.1f} entries={len({r["pdbid"] for r in b})}')
print(f'ALL  n={tot["n"]:4d} rec={tot["rec"]:4d} unrec={tot["n"]-tot["rec"]:4d} '
      f'rec%={100*tot["rec"]/tot["n"]:5.1f}')
hb = [r for r in usable if r['residue'] != 'CYS' and r['ib'] and r['hb']]
print(f'H-bond-capable in-band (O-O, O-N, N-O, N-N): n={len(hb)} recorded={sum(r["rec"] for r in hb)} '
      f'entries={len({r["pdbid"] for r in hb})}')
cy = [r for r in usable if r['residue'] == 'CYS' and r['ib']]
print(f'CYS control in-band: n={len(cy)} recorded={sum(r["rec"] for r in cy)} '
      f'({100*sum(r["rec"] for r in cy)/len(cy):.1f}%) entries={len({r["pdbid"] for r in cy})}')

print('\n## DENSITY / CLASH, UNRECORDED IN-BAND, NON-CYS')
unrec = [r for r in usable if r['ib'] and not r['rec'] and r['residue'] != 'CYS']
e = [ev.get((r['pdbid'], r['ligand'], r['protein_atom'])) for r in unrec]
e = [x for x in e if x and x['ir']]
ds = [x for x in e if x['density_class'] == 'density_supported']
print(f'unrecorded in-band {len(unrec)} rows / {len({r["pdbid"] for r in unrec})} entries; '
      f'in report {len(e)}; density-supported {len(ds)} rows / {len({x["pdbid"] for x in ds})} entries')
nh = [x for x in e if not x['hb']]
nhds = [x for x in nh if x['density_class'] == 'density_supported']
print(f'  excluding H-bond-capable pairs: {len(nh)} rows / {len({x["pdbid"] for x in nh})} entries; '
      f'density-supported {len(nhds)} rows / {len({x["pdbid"] for x in nhds})} entries')
print('  clash at pair:', sum(x['cl'] for x in e), 'of', len(e),
      '| non-H-bond pairs:', sum(x['cl'] for x in nh), 'of', len(nh))
keep = {'candidate_unrecorded_covalent_bond', 'candidate_bond_or_near_attack_complex',
        'cofactor_or_PTM_chemistry', 'phosphoryl_transfer_or_nucleotide'}
cand = [x for x in nhds if x['verdict'] in keep]
print(f'  chemistry-plausible subset: {len(cand)} rows / {len({x["pdbid"] for x in cand})} entries')
print('  verdicts in that subset:', dict(collections.Counter(x['verdict'] for x in cand)))
print('  residues:', dict(collections.Counter(x['residue'] for x in cand)))
print('  pairs:', dict(collections.Counter(x['pair'] for x in cand)))

print('\n## B-O SUB-RESULT')
bo = [r for r in usable if r['pair'] == 'B-O' and r['ib'] and not r['rec']]
boe = [ev.get((r['pdbid'], r['ligand'], r['protein_atom'])) for r in bo]
boe = [x for x in boe if x and x['ir']]
print(f'Ser/Thr B-O unrecorded in-band: {len(bo)} rows / {len({r["pdbid"] for r in bo})} entries; '
      f'density-supported {sum(1 for x in boe if x["density_class"]=="density_supported")}; '
      f'clash-listed {sum(x["cl"] for x in boe)}')

print('\n## RECORD LINE SUMMARY (pooled non-Cys + Cys control)')
for r in csv.DictReader(open(f'{P}_record_lines.csv')):
    if r['set'] in ('pooled-nonCys', 'cys-control'):
        print(f"  {r['set']:14s} {r['pair']:5s} probe {r['probe_line_A']:>5s}  record line {r['record_line_bin']:>11s}  "
              f"diff {r['record_minus_probe_line']:>6s}  n={r['n']:>5s} rec={r['recorded']:>5s} "
              f"past-line n={r['n_past_probe_line']:>4s} rec={r['recorded_past_probe_line']:>4s}")

print('\n## RECORDED BOND LENGTH vs LINE')
for r in csv.DictReader(open(f'{P}_recorded_bond_lengths.csv')):
    print(f"  {r['pair']:5s} line {r['probe_line_A']:>5s} n={r['n_recorded']:>5s} median {r['median']:>6s} "
          f"p99 {r['p99']:>6s} max {r['max']:>6s} headroom {r['headroom_median_to_line']:>6s} "
          f"%at/past line {r['pct_recorded_at_or_past_line']:>5s}")

print('\n## WITHIN-ENTRY SPLIT')
cl = list(csv.DictReader(open(f'{P}_within_entry_split.csv')))
nc = [x for x in cl if x['residue'] != 'CYS']
print('mixed groups', len(cl), 'non-Cys', len(nc), 'clean split', sum(1 for x in nc if tf(x['split_clean'])))

print('\n## CLOSE-CONTACT TABLE')
cc = list(csv.DictReader(open(f'{P}_close_contact_check.csv')))
print('in-band unrecorded rows', len(cc), 'in close-contact table',
      sum(1 for x in cc if tf(x['in_close_contact_table'])))

print('\n## DRUGS')
an = list(csv.DictReader(open(f'{P}_density_supported_sites_annotated.csv')))
ac = [r for r in an if r['verdict'] in keep and not tf(r['hbond_capable_pair'])]
apr = {r['ligand_comp']: (r['drugbank_name'], r['drug_groups']) for r in ac
       if 'approved' in (r['drug_groups'] or '')}
print('components mapping to an approved DrugBank drug:', apr)
print('targets:', dict(collections.Counter((r['uniprot_names'].split('|')[0].strip() or
                                            r['entity_descriptions'][:40]) for r in ac
                                           if r['verdict'].startswith('candidate'))))
