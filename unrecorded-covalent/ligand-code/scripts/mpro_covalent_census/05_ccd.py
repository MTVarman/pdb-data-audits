"""CCD name/formula/synonyms/DrugBank name for every non-polymer component in the Mpro entries and every partner comp."""
import json, urllib.request
D='/workspace/mpro_covalent_census'
comps=set()
for l in open(f'{D}/out/scan.jsonl'):
    r=json.loads(l); comps|=set(r.get('nonpoly_comps',{})); comps|={p['lig_comp'] for p in r['pairs']}
comps=sorted(comps); out={}
for i in range(0,len(comps),100):
    q='{chem_comps(comp_ids:%s){chem_comp{id name type formula formula_weight} rcsb_chem_comp_synonyms{name} drugbank{drugbank_info{name}} rcsb_chem_comp_descriptor{SMILES}}}'%json.dumps(comps[i:i+100])
    req=urllib.request.Request('https://data.rcsb.org/graphql',data=json.dumps({'query':q}).encode(),headers={'Content-Type':'application/json'})
    for c in json.load(urllib.request.urlopen(req,timeout=120))['data']['chem_comps']:
        if not c: continue
        cc=c['chem_comp']
        out[cc['id']]={'name':cc['name'],'type':cc['type'],'formula':cc['formula'],'fw':cc['formula_weight'],
            'synonyms':[s['name'] for s in (c.get("rcsb_chem_comp_synonyms") or [])][:15],
            'drugbank':((c.get('drugbank') or {}).get('drugbank_info') or {}).get('name'),
            'smiles':(c.get('rcsb_chem_comp_descriptor') or {}).get('SMILES')}
json.dump(out,open(f'{D}/out/ccd_info.json','w'),indent=0); print(len(comps),len(out))
