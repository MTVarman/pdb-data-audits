"""Fetch CCD name/type/formula for every component seen in the scan (RCSB GraphQL)."""
import json, os, urllib.request, sys
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
comps = set()
for fn in ('scan.jsonl', 'scan_big.jsonl'):
    for l in open(f'{D}/out/{fn}'):
        r = json.loads(l)
        for p in r.get('pairs', []): comps.add(p['lig_comp'])
comps = sorted(comps); out = {}
for i in range(0, len(comps), 100):
    ids = comps[i:i+100]
    q = '{chem_comps(comp_ids:%s){chem_comp{id name type formula}}}' % json.dumps(ids)
    req = urllib.request.Request('https://data.rcsb.org/graphql', data=json.dumps({'query': q}).encode(), headers={'Content-Type': 'application/json'})
    for c in json.load(urllib.request.urlopen(req, timeout=120))['data']['chem_comps']:
        if c: cc = c['chem_comp']; out[cc['id']] = {'name': cc['name'], 'type': cc['type'], 'formula': cc['formula']}
json.dump(out, open(f'{D}/out/ccd_info.json', 'w'), indent=0)
print(len(comps), len(out))
