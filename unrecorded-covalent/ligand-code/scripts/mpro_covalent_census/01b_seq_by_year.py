"""The RCSB sequence service caps at 2000 hits; partition the sequence search by initial release date to get all hits."""
import json, urllib.request, time
D='/workspace/mpro_covalent_census'
n=open(f'{D}/data/nsp5_sars2.fasta').read().split()[1]
def q(lo,hi):
    body={'query':{'type':'group','logical_operator':'and','nodes':[
      {'type':'terminal','service':'sequence','parameters':{'evalue_cutoff':1e-10,'identity_cutoff':0.85,'sequence_type':'protein','value':n}},
      {'type':'terminal','service':'text','parameters':{'attribute':'rcsb_accession_info.initial_release_date','operator':'range','value':{'from':lo,'to':hi,'include_lower':True,'include_upper':False}}}]},
      'return_type':'polymer_entity','request_options':{'return_all_hits':True,'results_verbosity':'compact'}}
    r=urllib.request.Request('https://search.rcsb.org/rcsbsearch/v2/query',data=json.dumps(body).encode(),headers={'Content-Type':'application/json'})
    s=urllib.request.urlopen(r,timeout=300).read(); return json.loads(s)['result_set'] if s else []
bounds=['1990-01-01','2020-01-01','2020-07-01','2021-01-01','2021-07-01','2022-01-01','2022-07-01','2023-01-01','2023-07-01','2024-01-01','2024-07-01','2025-01-01','2025-07-01','2026-01-01','2026-07-01','2027-01-01']
allh=set()
for lo,hi in zip(bounds,bounds[1:]):
    h=q(lo,hi); print(lo,hi,len(h)); assert len(h)<2000; allh|=set(h); time.sleep(0.5)
print('total entities',len(allh),'entries',len({x.split('_')[0] for x in allh}))
d=json.load(open(f'{D}/data/search_raw.json')); d['A_partitioned']=sorted(allh); json.dump(d,open(f'{D}/data/search_raw.json','w'))
