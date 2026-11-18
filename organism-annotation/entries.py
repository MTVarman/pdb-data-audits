import json, urllib.request, urllib.parse, random, sys
def search(taxid, rows=None):
    q={"query":{"type":"group","logical_operator":"and","nodes":[
        {"type":"terminal","service":"text","parameters":{"attribute":"rcsb_entity_source_organism.taxonomy_lineage.id","operator":"exact_match","value":str(taxid)}},
        {"type":"terminal","service":"text","parameters":{"attribute":"entity_poly.rcsb_entity_polymer_type","operator":"exact_match","value":"Protein"}},
        {"type":"terminal","service":"text","parameters":{"attribute":"rcsb_entry_info.structure_determination_methodology","operator":"exact_match","value":"experimental"}}]},
       "return_type":"entry","request_options":{"paginate":{"start":0,"rows":10000},"results_content_type":["experimental"]}}
    ids=[]; start=0
    while True:
        q["request_options"]["paginate"]["start"]=start
        url="https://search.rcsb.org/rcsbsearch/v2/query?json="+urllib.parse.quote(json.dumps(q))
        d=json.load(urllib.request.urlopen(url,timeout=120))
        got=[r["identifier"] for r in d.get("result_set",[])]
        ids+=got; start+=10000
        if start>=d["total_count"] or not got: break
    return ids
groups={"dog":9615,"cat":9685,"horse":9796,"human":9606}
out={}
for g,t in groups.items():
    ids=search(t); print(g,len(ids),file=sys.stderr); out[g]=ids
random.seed(42)
out["human_sample"]=random.sample(out["human"],2000)
json.dump(out,open("entries.json","w"))
