"""One row per unrecorded main-set Cys145 contact with metadata, dataset membership, close-contact flag, version history."""
import pandas as pd
D='/workspace/mpro_covalent_census/out'
X=pd.read_csv(f'{D}/contacts_main_datasets.tsv',sep='\t',keep_default_na=False,low_memory=False)
X=X[X.recorded.astype(str)=='False']
A=pd.read_csv(f'{D}/annotation_fields.tsv',sep='\t',keep_default_na=False)[['entry','nuc_chain','lig_comp','dist','flagged_close_contact']]
V=pd.read_csv(f'{D}/unrecorded_version_history.tsv',sep='\t',keep_default_na=False)[['entry','nuc_chain','lig_comp','dist','covale_by_version']]
E=pd.read_csv(f'{D}/entries.tsv',sep='\t',keep_default_na=False)[['entry','deposit_date','authors','citation_journal','citation_year','citation_doi','title']]
P=pd.read_csv(f'{D}/protein_modification_flag.tsv',sep='\t',keep_default_na=False)[['entry','has_protein_modification']]
for t in (A,V): t['dist']=t.dist.astype(float)
X['dist']=X.dist.astype(float)
X=X.merge(A,on=['entry','nuc_chain','lig_comp','dist'],how='left').merge(V,on=['entry','nuc_chain','lig_comp','dist'],how='left').merge(E,on='entry',how='left').merge(P,on='entry',how='left')
X['last_author']=X.authors.str.split(';').str[-1]
cols=['entry','process_site','deposit_date','release_year','lig_comp','drug','lig_heavy_atoms','nuc_chain','nuc_seq','lig_atom','dist','band','flagged_close_contact','has_protein_modification','covale_by_version',
      'plinder_is_covalent','legacy_split','rnp_ligand','biolip_ligand_site','last_author','citation_journal','citation_year','citation_doi','title']
X[cols].sort_values(['band','process_site','entry','nuc_chain']).to_csv(f'{D}/unrecorded_contacts.tsv',sep='\t',index=False)
print(len(X))
