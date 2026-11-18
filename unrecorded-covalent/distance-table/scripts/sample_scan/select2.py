"""Evidence-sample candidate set.

Arms (none selected on _struct_conn content; all X-RAY DIFFRACTION, resolution <= 2.5 A):
  A_text     covalent-chemistry words in struct.title or struct_keywords.text
  B_warhead  entry contains a non-polymer component whose CCD structure contains a
             specific electrophilic-warhead substructure (RCSB chemical service,
             sub-struct-graph-strict). Over-broad patterns that matched most of the
             archive (O=CC "aldehyde" 163,006 entries; COP(=O)O 32,219) are excluded.
  C_baseline every 4th entry, by entry ID, of all X-ray <= 2.5 A entries that have a
             non-polymer entity (pool 139,519) -- unenriched reference arm
ids.txt is written in a seed-fixed shuffled order so any prefix of it is a uniform
random subsample of the union.
"""
import json, random, sys, time
sys.path.insert(0, '/home/claude/work/lane3')
from rcsb3 import XRAY, res_le, text, andq, orq, ids

STRIDE = 4
TERMS = 'covalent covalently adduct acyl-enzyme irreversible warhead thioester'
WARHEADS = {
    'boronic_acid': 'B(O)O',
    'beta_lactam': 'O=C1CCN1',
    'epoxide': 'C1CO1',
    'sulfonyl_fluoride': 'FS(=O)=O',
    'chloroacetamide': 'ClCC(=O)N',
    'vinyl_sulfone': 'C=CS(=O)=O',
    'aziridine': 'C1CN1',
    'acrylamide': 'NC(=O)C=C',
    'nitrile': 'C#N',
    'carbamate': 'O=C(O)N',
}
EXCLUDED_WARHEADS = {'aldehyde O=CC': 163006, 'phosphonate_ester COP(=O)O': 32219}


def chem(smi):
    return {'type': 'terminal', 'service': 'chemical',
            'parameters': {'value': smi, 'type': 'descriptor', 'descriptor_type': 'SMILES',
                           'match_type': 'sub-struct-graph-strict'}}


arms, qs, wh_counts = {}, {}, {}
qs['A_text'] = andq(XRAY, res_le(2.5),
                    orq(text('struct.title', 'contains_words', TERMS),
                        text('struct_keywords.text', 'contains_words', TERMS)))
got, total = ids(qs['A_text'])
arms['A_text'] = got
print(f'A_text      {len(got)}', flush=True)

wh = set()
for name, smi in WARHEADS.items():
    try:
        g, t = ids(andq(XRAY, res_le(2.5), chem(smi)))
        wh |= set(g)
        wh_counts[name] = len(g)
        print(f'  {name:20s} {smi:14s} {len(g):6d}', flush=True)
    except Exception as e:
        wh_counts[name] = f'FAILED {str(e)[:60]}'
        print(f'  {name:20s} {smi:14s} FAILED', flush=True)
arms['B_warhead'] = sorted(wh)
print(f'B_warhead   {len(wh)}', flush=True)

base = json.load(open('/home/claude/work/lane3/baseline_pool.json'))
qs['C_baseline'] = base['query']
arms['C_baseline'] = base['ids'][::STRIDE]
print(f'C_baseline  {len(arms["C_baseline"])} of pool {base["total"]} (every {STRIDE}th by entry ID)', flush=True)

union = sorted(set().union(*arms.values()))
print('UNION', len(union), flush=True)
json.dump({'retrieved_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
           'baseline_pool_retrieved_utc': base['retrieved_utc'],
           'terms': TERMS, 'warheads': WARHEADS, 'excluded_warheads': EXCLUDED_WARHEADS,
           'warhead_counts': wh_counts, 'queries': qs,
           'counts': {k: len(v) for k, v in arms.items()},
           'baseline_pool_size': base['total'], 'baseline_stride': STRIDE,
           'union': len(union), 'shuffle_seed': 20261115, 'arms': arms},
          open('/home/claude/work/lane3/arms.json', 'w'))
sh = list(union)
random.Random(20261115).shuffle(sh)
open('/home/claude/work/lane3/ids.txt', 'w').write('\n'.join(sh) + '\n')
print('ids.txt written, seed-shuffled', len(sh))
