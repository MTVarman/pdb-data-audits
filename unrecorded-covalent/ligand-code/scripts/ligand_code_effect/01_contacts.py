"""Merge covalent-contact rows (ligand carbon 1.6-2.4 A from Cys SG) from the Mpro census (Cys145) and the
Cys covalent scan (drug-like, all Cys). One row per (entry, cys chain, cys seq, ligand comp/chain/seq).
Mpro rows take precedence where both exist (Mpro census is the more complete for Mpro)."""
import pandas as pd, json
M='/workspace/mpro_covalent_census/out'; C='/workspace/cys_covalent_scan/out'; O='/workspace/ligand_code_effect/out'
m=pd.read_csv(f'{M}/contacts_main_cys145_C_1.6-2.4.tsv',sep='\t',dtype=str,keep_default_na=False)
m=m.rename(columns={'nuc_chain':'cys_chain','nuc_seq':'cys_seq','nuc_ins':'cys_ins'}); m['source']='mpro_census'
c=pd.read_csv(f'{C}/contacts_druglike_C_1.6-2.4.tsv',sep='\t',dtype=str,keep_default_na=False); c['source']='cys_scan'
cols=['entry','cys_chain','cys_seq','cys_ins','lig_comp','lig_chain','lig_seq','lig_ins','lig_atom','dist','recorded','rec_atom_exact','process_site','deposit_year','method','resolution','category','drug','lig_name','source']
d=pd.concat([m[cols],c[cols]]); key=['entry','cys_chain','cys_seq','cys_ins','lig_comp','lig_chain','lig_seq','lig_ins']
d=d.drop_duplicates(key,keep='first')
# deposit dates
dd={}
for l in open(f'{C}/entries_scanned.tsv'):
    f=l.rstrip('\n').split('\t')
    if f[0]!='entry': dd[f[0]]=(f[3],f[6],None)
for fn in ('/workspace/mpro_covalent_census/out/scan.jsonl',):
    for l in open(fn):
        r=json.loads(l); au=(r.get('authors') or '').split(';')
        dd[r['entry']]=(r.get('deposit_date'),r.get('initial_release'),au[-1] if au and au[-1] else None)
d['deposit_date']=d.entry.map(lambda e:dd.get(e,(None,))[0])
d['release_date']=d.entry.map(lambda e:dd.get(e,(None,None))[1])
d['last_author']=d.entry.map(lambda e:dd.get(e,(None,None,None))[2])
d['dist']=d.dist.astype(float); d['recorded']=d.recorded.eq('True')
d.to_csv(f'{O}/contacts_merged.tsv',sep='\t',index=False)
print(len(d), d.source.value_counts().to_dict(), d.entry.nunique(), d.lig_comp.nunique(), (d.dist<2).sum())
