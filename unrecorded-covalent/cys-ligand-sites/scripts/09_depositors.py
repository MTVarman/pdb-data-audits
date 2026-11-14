"""Depositing-group and target metadata for the entries behind the drug-like contacts under 2.0 A.

Reads out/contacts_druglike_C_1.6-2.4.tsv, asks the RCSB Data API (GraphQL) for each entry's
audit authors, primary citation, title and source organisms, and writes

  out/depositor_groups.tsv      one row per entry: two grouping keys and its contact counts
  out/nirmatrelvir_targets.tsv  one row per entry with a nirmatrelvir (4WI) contact under 2.0 A

No author names are written. citation_group is the primary citation's DOI or PubMed ID; an entry
without either gets a hash of its author list, so entries from one unpublished set still group.
last_author_group is a hash of the last audit author's name.

usage: python3 scripts/09_depositors.py        (standard library only; ~25 requests)
"""
import csv, hashlib, json, os, time, urllib.request
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'); O = f'{D}/out'
API = 'https://data.rcsb.org/graphql'

def h(s): return hashlib.sha1(s.encode('utf-8')).hexdigest()[:12]

def fetch(ids):
    q = ('{entries(entry_ids:%s){rcsb_id struct{title} audit_author{name} '
         'rcsb_primary_citation{pdbx_database_id_DOI pdbx_database_id_PubMed} '
         'polymer_entities{rcsb_entity_source_organism{scientific_name}}}}' % json.dumps(ids))
    req = urllib.request.Request(API, data=json.dumps({'query': q}).encode(), headers={'Content-Type': 'application/json'})
    for attempt in range(4):
        try:
            return json.load(urllib.request.urlopen(req, timeout=90))['data']['entries']
        except Exception:
            if attempt == 3: raise
            time.sleep(3)

rows = [r for r in csv.DictReader(open(f'{O}/contacts_druglike_C_1.6-2.4.tsv', newline=''), delimiter='\t') if float(r['dist']) < 2.0]
ents = sorted({r['entry'] for r in rows})
info = {}
for i in range(0, len(ents), 100):
    for e in fetch(ents[i:i + 100]): info[e['rcsb_id']] = e
missing = [e for e in ents if e not in info]
if missing: raise SystemExit(f'no metadata returned for {missing[:10]}')

site, n, miss = {}, {}, {}
for r in rows:
    e = r['entry']; site[e] = r['process_site']; n[e] = n.get(e, 0) + 1
    miss[e] = miss.get(e, 0) + (r['recorded'] != 'True')

with open(f'{O}/depositor_groups.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['entry', 'process_site', 'citation_group', 'last_author_group', 'n_contacts_lt2', 'n_unrecorded_lt2'])
    for e in ents:
        au = [a['name'] for a in (info[e].get('audit_author') or [])]
        c = info[e].get('rcsb_primary_citation') or {}
        if c.get('pdbx_database_id_DOI'): cg = 'doi:' + c['pdbx_database_id_DOI'].lower()
        elif c.get('pdbx_database_id_PubMed'): cg = 'pmid:%s' % c['pdbx_database_id_PubMed']
        else: cg = 'authors:' + h('|'.join(au))
        w.writerow([e, site[e], cg, 'la:' + h(au[-1] if au else '?'), n[e], miss[e]])

with open(f'{O}/nirmatrelvir_targets.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['entry', 'process_site', 'source_organism', 'title'])
    for e in sorted({r['entry'] for r in rows if r['lig_comp'] == '4WI'}):
        org = sorted({o['scientific_name'] for pe in (info[e].get('polymer_entities') or [])
                      for o in (pe.get('rcsb_entity_source_organism') or []) if o.get('scientific_name')})
        w.writerow([e, site[e], '; '.join(org), (info[e].get('struct') or {}).get('title') or ''])
print(f'{len(ents)} entries -> out/depositor_groups.tsv, out/nirmatrelvir_targets.tsv')
