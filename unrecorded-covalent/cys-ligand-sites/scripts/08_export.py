"""Write the per-contact export tables (unrecorded <2.0 A, 2.0-2.4 A, named drugs, site entry-level).
Optional: flag rows already listed in Varman's strict-cut-491 file if its path is given as argv[1]."""
import os, sys, pandas as pd
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'); O = f'{D}/out'
d = pd.read_csv(f'{O}/contacts_druglike_C_1.6-2.4.tsv', sep='\t', keep_default_na=False, na_values=[''])
cols = ['entry','process_site','deposit_year','method','resolution','lig_comp','lig_name','drug','drug_class','cys_chain','cys_seq','cys_alt','lig_chain','lig_seq','lig_atom','lig_alt','lig_occ','dist','recorded','rec_atom_exact','rec_dist','bonded_elsewhere']
m = d[(d.dist < 2.0) & (~d.recorded)][cols].sort_values(['process_site', 'entry']).copy()
if len(sys.argv) > 1:
    s = pd.read_csv(sys.argv[1]); s = s[s.residue_name == 'CYS']
    keys = set(f'{a}|{b}|{c}' for a, b, c in zip(s.entry, s.chain, s.residue_number))
    m['in_strict_cut_491'] = [f'{a}|{b}|{c}' in keys for a, b, c in zip(m.entry, m.cys_chain, m.cys_seq)]
m.to_csv(f'{O}/unrecorded_under_2.0A_druglike.tsv', sep='\t', index=False)
d[d.category == 'named_drug'][cols].sort_values(['drug_class', 'lig_comp', 'entry']).to_csv(f'{O}/named_drug_contacts.tsv', sep='\t', index=False)
d[d.dist >= 2.0][cols].sort_values('dist').to_csv(f'{O}/druglike_2.0-2.4A_contacts.tsv', sep='\t', index=False)
a = d[d.dist < 2.0]
e = a.groupby(['entry', 'process_site']).recorded.agg(lambda x: (~x).any()).reset_index()
e.groupby('process_site').recorded.agg(['sum', 'count']).rename(columns={'sum': 'entries_with_miss', 'count': 'entries'}).to_csv(f'{O}/site_entry_level_druglike.tsv', sep='\t')
print(len(m), m.get('in_strict_cut_491', pd.Series(dtype=bool)).sum())
