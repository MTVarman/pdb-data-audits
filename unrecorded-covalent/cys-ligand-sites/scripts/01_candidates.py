"""Build the candidate entry list (reproducible: RCSB Search API, queried 2026-11-14).
Arms:
  A  full_text "covalent cysteine"
  B  full_text "covalent inhibitor"
  C  every entry containing a named covalent drug component (comp IDs verified against CCD names)
Writes candidates.tsv (entry, arms) and drug_comps.tsv.
"""
import sys, json, os
sys.path.insert(0, os.path.dirname(__file__)); from rcsb_q import search, t
DRUGS = {  # comp_id: (drug, class) -- checked via data.rcsb.org chemcomp names/synonyms
 'MOV':('sotorasib (AMG 510, bound form)','KRAS G12C'), 'A1B7W':('adagrasib','KRAS G12C'), 'M1X':('MRTX849/adagrasib (6UT0 ligand)','KRAS G12C'),
 '91S':('ARS-1620','KRAS G12C'), '5UT':('ARS-853','KRAS G12C'), 'A1AWR':('divarasib (bound form)','KRAS G12C'),
 '1E8':('ibrutinib','BTK'), '8E8':('ibrutinib (bound form)','BTK'), 'A1H14':('ibrutinib (unbound form)','BTK'),
 'XQQ':('acalabrutinib','BTK'), 'BA0':('zanubrutinib','BTK'),
 '0WM':('afatinib','EGFR/HER'), '0WN':('afatinib (bound form)','EGFR/HER'), 'YY3':('osimertinib','EGFR/HER'),
 'HKI':('neratinib (bound form)','EGFR/HER'), '1C9':('dacomitinib','EGFR/HER'),
 '4WI':('nirmatrelvir (bound form)','SARS-CoV-2 Mpro'), 'ZGW':('nirmatrelvir','SARS-CoV-2 Mpro'), 'OW1':('nirmatrelvir (reacted form)','SARS-CoV-2 Mpro'),
}
arms = {}
for name, q in [('A', {"type":"terminal","service":"full_text","parameters":{"value":"covalent cysteine"}}),
                ('B', {"type":"terminal","service":"full_text","parameters":{"value":"covalent inhibitor"}})]:
    r = search(q, 'entry', True)
    for e in r['result_set']: arms.setdefault(e, set()).add(name)
    print(name, len(r['result_set']))
for c in DRUGS:
    r = search(t("rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id","exact_match",c), 'entry', True)
    ids = r.get('result_set', [])
    print(c, DRUGS[c][0], len(ids))
    for e in ids: arms.setdefault(e, set()).add('C')
out = os.path.join(os.path.dirname(__file__), '..', 'out')
with open(f'{out}/candidates.tsv','w') as f:
    f.write('entry\tarms\n')
    for e in sorted(arms): f.write(f"{e}\t{''.join(sorted(arms[e]))}\n")
with open(f'{out}/drug_comps.tsv','w') as f:
    f.write('comp_id\tdrug\tclass\n')
    for c,(d,k) in DRUGS.items(): f.write(f'{c}\t{d}\t{k}\n')
print('total candidates', len(arms))
