"""Read _pdbx_entry_details.has_protein_modification and _pdbx_modification_feature rows touching Cys145
for every covalent-complex entry; cross-tab against recorded/unrecorded."""
import gemmi, pandas as pd
D='/workspace/mpro_covalent_census'
EN=pd.read_csv(f'{D}/out/entries_covalent.tsv',sep='\t')
rows=[]
for e in EN.entry:
    b=gemmi.cif.read(f'{D}/data/cif/{e}.cif.gz').sole_block()
    v=b.find_value('_pdbx_entry_details.has_protein_modification')
    mf=b.get_mmcif_category('_pdbx_modification_feature.')
    n145=0
    if mf:
        for i in range(len(mf.get('ordinal',[]))):
            if (mf.get('modified_residue_auth_seq_id') or [None]*999)[i]=='145' or (mf.get('label_seq_id') and False): n145+=1
    rows.append({'entry':e,'has_protein_modification':gemmi.cif.as_string(v) if v else None,'n_modification_feature':len(mf.get('ordinal',[])) if mf else 0,'n_modfeat_145':n145})
X=EN[['entry','process_site','has_unrec_lt2','all_unrec']].merge(pd.DataFrame(rows),on='entry')
X.to_csv(f'{D}/out/protein_modification_flag.tsv',sep='\t',index=False)
print(pd.crosstab([X.all_unrec],X.has_protein_modification.fillna('absent')))
print(pd.crosstab([X.all_unrec],X.n_modfeat_145>0))
