"""Drug and target annotation for the evidence-sample site list (RCSB Data API GraphQL).

Usage: python3 -I annotate.py sites.csv out.csv
Adds, per row: DrugBank name / drug groups for the ligand component, the polymer
entity descriptions for the entry, and the UniProt protein names.
"""
import csv, json, sys, time, urllib.request

GQL = 'https://data.rcsb.org/graphql'


def q(query, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(GQL, json.dumps({'query': query}).encode(),
                                         {'Content-Type': 'application/json'})
            return json.loads(urllib.request.urlopen(req, timeout=180).read())
        except Exception as e:  # noqa
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(str(last))


def comps(ids):
    out = {}
    for i in range(0, len(ids), 50):
        b = ids[i:i + 50]
        lst = ','.join(f'"{x}"' for x in b)
        r = q(f'{{chem_comps(comp_ids:[{lst}]){{rcsb_id chem_comp{{name formula type}} '
              f'drugbank{{drugbank_info{{name drug_groups drug_categories}}}} '
              f'rcsb_chem_comp_synonyms{{name provenance_source}}}}}}')
        if r.get('errors'):
            print('GraphQL errors:', str(r['errors'])[:300], flush=True)
        for c in (r.get('data', {}).get('chem_comps') or []):
            if not c:
                continue
            db = ((c.get('drugbank') or {}).get('drugbank_info') or {})
            syn = [s['name'] for s in (c.get('rcsb_chem_comp_synonyms') or [])
                   if s.get('provenance_source') in ('DrugBank', 'CHEBI', 'PDB-CARE', 'ChEBI')]
            out[c['rcsb_id']] = {
                'comp_name': (c.get('chem_comp') or {}).get('name') or '',
                'comp_formula': (c.get('chem_comp') or {}).get('formula') or '',
                'drugbank_name': db.get('name') or '',
                'drug_groups': '|'.join(db.get('drug_groups') or []),
                'drug_categories': '|'.join((db.get('drug_categories') or [])[:4]),
                'synonyms': '|'.join(syn[:3]),
            }
    return out


def entries(ids):
    out = {}
    for i in range(0, len(ids), 40):
        b = ids[i:i + 40]
        lst = ','.join(f'"{x}"' for x in b)
        r = q(f'{{entries(entry_ids:[{lst}]){{rcsb_id struct{{title}} '
              f'polymer_entities{{rcsb_polymer_entity{{pdbx_description}} '
              f'rcsb_polymer_entity_container_identifiers{{uniprot_ids}} '
              f'uniprots{{rcsb_uniprot_protein{{name{{value}}}}}}}}}}}}')
        for e in (r.get('data', {}).get('entries') or []):
            if not e:
                continue
            desc, up, upn = [], [], []
            for pe in (e.get('polymer_entities') or []):
                d = ((pe.get('rcsb_polymer_entity') or {}).get('pdbx_description') or '')
                if d:
                    desc.append(d)
                up += ((pe.get('rcsb_polymer_entity_container_identifiers') or {}).get('uniprot_ids') or [])
                for u in (pe.get('uniprots') or []):
                    nm = (((u or {}).get('rcsb_uniprot_protein') or {}).get('name') or {}).get('value')
                    if nm:
                        upn.append(nm)
            out[e['rcsb_id']] = {'entity_descriptions': ' | '.join(dict.fromkeys(desc))[:200],
                                 'uniprot_ids': ' '.join(dict.fromkeys(up)),
                                 'uniprot_names': ' | '.join(dict.fromkeys(upn))[:200]}
    return out


def main():
    src, dst = sys.argv[1], sys.argv[2]
    rows = list(csv.DictReader(open(src)))
    if not rows:
        print('no rows')
        return
    cs = sorted({r['ligand_comp'] for r in rows})
    es = sorted({r['pdbid'] for r in rows})
    print('components', len(cs), 'entries', len(es), flush=True)
    cann, eann = comps(cs), entries(es)
    print('annotated components', len(cann), 'entries', len(eann), flush=True)
    keys = list(rows[0].keys()) + ['comp_name', 'comp_formula', 'drugbank_name', 'drug_groups',
                                   'drug_categories', 'synonyms', 'entity_descriptions',
                                   'uniprot_ids', 'uniprot_names']
    blank_c = {k: '' for k in ('comp_name', 'comp_formula', 'drugbank_name', 'drug_groups',
                               'drug_categories', 'synonyms')}
    blank_e = {k: '' for k in ('entity_descriptions', 'uniprot_ids', 'uniprot_names')}
    with open(dst, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({**r, **blank_c, **blank_e,
                        **cann.get(r['ligand_comp'], {}), **eann.get(r['pdbid'], {})})
    napp = sum(1 for r in rows if 'approved' in (cann.get(r['ligand_comp'], {}).get('drug_groups') or ''))
    print('rows whose ligand component maps to an approved drug in DrugBank:', napp)
    print('rows whose ligand component maps to any DrugBank entry:',
          sum(1 for r in rows if cann.get(r['ligand_comp'], {}).get('drugbank_name')))


if __name__ == '__main__':
    main()
