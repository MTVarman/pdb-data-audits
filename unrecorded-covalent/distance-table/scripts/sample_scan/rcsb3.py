"""RCSB search helpers for the evidence sample."""
import json, time, urllib.request

URL = 'https://search.rcsb.org/rcsbsearch/v2/query'


def post(qy, tries=4, timeout=300):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(URL, json.dumps(qy).encode(), {'Content-Type': 'application/json'})
            return json.loads(urllib.request.urlopen(req, timeout=timeout).read())
        except Exception as e:  # noqa
            last = e
            time.sleep(3 * (i + 1))
    raise RuntimeError(f'search failed: {last}')


XRAY = {'type': 'terminal', 'service': 'text',
        'parameters': {'attribute': 'exptl.method', 'operator': 'exact_match', 'value': 'X-RAY DIFFRACTION'}}


def res_le(v):
    return {'type': 'terminal', 'service': 'text',
            'parameters': {'attribute': 'rcsb_entry_info.resolution_combined', 'operator': 'less_or_equal', 'value': v}}


def text(attr, op, val):
    return {'type': 'terminal', 'service': 'text', 'parameters': {'attribute': attr, 'operator': op, 'value': val}}


def chem(smiles):
    return {'type': 'terminal', 'service': 'chemical',
            'parameters': {'value': smiles, 'type': 'descriptor', 'descriptor_type': 'SMILES',
                           'match_type': 'graph-relaxed'}}


def andq(*nodes):
    return {'type': 'group', 'logical_operator': 'and', 'nodes': list(nodes)}


def orq(*nodes):
    return {'type': 'group', 'logical_operator': 'or', 'nodes': list(nodes)}


def count(node):
    return post({'query': node, 'return_type': 'entry', 'request_options': {'return_counts': True}})['total_count']


def ids(node, cap=400000):
    out, start, rows, total = [], 0, 10000, 0
    while True:
        r = post({'query': node, 'return_type': 'entry',
                  'request_options': {'paginate': {'start': start, 'rows': rows},
                                      'results_content_type': ['experimental'],
                                      'sort': [{'sort_by': 'rcsb_entry_container_identifiers.entry_id',
                                                'direction': 'asc'}]}})
        got = [x['identifier'] for x in r.get('result_set', [])]
        out += got
        total = r.get('total_count', 0)
        if not got or len(out) >= min(total, cap):
            break
        start += rows
    return out, total
