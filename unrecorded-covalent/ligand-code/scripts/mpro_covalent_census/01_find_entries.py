"""Find every SARS-CoV-2 Mpro (nsp5) entry in the PDB via the RCSB Search API.
Three routes, unioned, then each polymer entity is checked locally against reference nsp5 sequences.
 A) sequence search (mmseqs2) with SARS-CoV-2 nsp5 (P0DTD1 3264-3569), identity >= 0.85, e<=1e-10
 B) UniProt accession P0DTD1 or P0DTC1 on any polymer entity (pp1ab / pp1a); then filtered by alignment to nsp5
 C) full text 'main protease' / '3CL' / 'Mpro' / 'nsp5' AND organism SARS-CoV-2 (taxid 2697049)
"""
import json, urllib.request, os, time
D = os.path.dirname(os.path.abspath(__file__)) + '/..'
def fetch(u, data=None, hdr=None):
    for i in range(5):
        try:
            r = urllib.request.Request(u, data=data, headers=hdr or {})
            with urllib.request.urlopen(r, timeout=300) as f: return f.read()
        except Exception as e: err = e; time.sleep(3*(i+1))
    raise err
pp = json.loads(fetch('https://rest.uniprot.org/uniprotkb/P0DTD1.json'))['sequence']['value']
NSP5 = pp[3263:3569]; assert NSP5.startswith('SGFRKMAF') and NSP5.endswith('SGVTFQ'), NSP5
open(f'{D}/data/nsp5_sars2.fasta','w').write('>P0DTD1_3264-3569\n'+NSP5+'\n')
def search(q, rt='polymer_entity'):
    body = {'query': q, 'return_type': rt, 'request_options': {'return_all_hits': True, 'results_verbosity': 'compact'}}
    s = fetch('https://search.rcsb.org/rcsbsearch/v2/query', json.dumps(body).encode(), {'Content-Type': 'application/json'})
    return json.loads(s)['result_set'] if s else []
def t(attr, op, val): return {'type': 'terminal', 'service': 'text', 'parameters': {'attribute': attr, 'operator': op, 'value': val}}
A = search({'type': 'terminal', 'service': 'sequence', 'parameters': {'evalue_cutoff': 1e-10, 'identity_cutoff': 0.85, 'sequence_type': 'protein', 'value': NSP5}})
B = search({'type': 'group', 'logical_operator': 'and', 'nodes': [
    t('rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession', 'in', ['P0DTD1', 'P0DTC1']),
    t('rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name', 'exact_match', 'UniProt')]})
C = search({'type': 'group', 'logical_operator': 'and', 'nodes': [
    {'type': 'group', 'logical_operator': 'or', 'nodes': [{'type': 'terminal', 'service': 'full_text', 'parameters': {'value': v}} for v in ['"main protease"', '3CLpro', 'Mpro', 'nsp5', '"3C-like"']]},
    t('rcsb_entity_source_organism.taxonomy_lineage.id', 'exact_match', '2697049')]}, rt='entry')
print('A seq', len(A), 'B uniprot', len(B), 'C text(entries)', len(C))
json.dump({'A': A, 'B': B, 'C': C, 'nsp5': NSP5, 'queried_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}, open(f'{D}/data/search_raw.json', 'w'))
