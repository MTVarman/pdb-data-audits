"""Flatten scan.jsonl into contacts and entries tables with ligand categories and drug names."""
import json, sys, re
import pandas as pd
sys.path.insert(0,'/workspace/cys_covalent_scan/scripts')
D='/workspace/mpro_covalent_census'
ccd=json.load(open(f'{D}/out/ccd_info.json'))
THIOL_ADD={'BME','DTT','DTU','DTV','GSH','CYS','2ME','MRY','CME','TCE','DTD','MPO','CSD','DCY','GTT','GDS','GSF','GSN','MTL','ACM','MPD','EDO','GOL','PEG','PGE','1PE','PG4','ACT','FMT','DMS','EPE','MES','TRS','BU3','CIT','SO4','PO4','IMD','NHE','MLI','IPA','EOH','MOH','BCN','B3P','ACY','AKG','TLA','TAR','PGO','P6G','SIN','12P','15P','PE8','PE4','TMO','DOD'}
METALS=re.compile(r'\b(Fe|Ni|Co|Cu|Zn|Mn|Mo|W|V|Mg|Ca|Cd|Hg|Pt|Au|Ag|Pd|Ru|Rh|Ir|Os|Re|Ga|Al|Sn|Pb|Cr|Na|K|Li|Tl|Sb|Bi|Te|Ti|Zr|Y|La|Ce|Eu|Gd|Tb|Yb|Lu|Sm|U|As|Se)\d*\b')
DRUGS={'4WI':'nirmatrelvir (bound, thioimidate)','ZGW':'nirmatrelvir (nitrile)','OW1':'nirmatrelvir (reacted form)',
 '7YY':'ensitrelvir (non-covalent control)','K36':'GC376 (bisulfite adduct)','UED':'GC376/GC373 (bound form)',
 'U5G':'boceprevir (bound)','SV6':'telaprevir (bound)','V2M':'PF-00835231 (bound)','80I':'lufotrelvir (PF-07304814)',
 'YDL':'ibuzatrelvir (bound)','7ON':'leritrelvir (bound)','A1AA0':'zevotrelvir/EDP-235 (bound)','NNA':'narlaprevir (bound)',
 'FHR':'bofutrelvir/11a','ZQB':'pomotrelvir (bound)','9JT':'ebselen (bound)','VR4':'GC376 alpha-ketoamide analog','A1EN0':'leritrelvir (unreacted-form CCD)','PQL':'simnotrelvir (nitrile CCD)','HUR':'YH-53'}
def cat(c):
    if c in DRUGS: return 'named_drug'
    i=ccd.get(c,{}); f=i.get('formula') or ''
    if METALS.search(f): return 'metal_or_Se_containing'
    if c in THIOL_ADD: return 'additive_or_thiol'
    return 'other_ligand'
rows=[]; ents=[]
for l in open(f'{D}/out/scan.jsonl'):
    r=json.loads(l)
    e={k:r.get(k) for k in ['entry','process_site','deposit_site','deposit_date','initial_release','last_revision','n_revisions','method','resolution','title','deposit_group','deposit_group_title','deposit_group_type','citation_title','citation_journal','citation_year','citation_doi','citation_pubmed','authors','software','n_mpro_chains','pos145_residues','n_covale_rows_total','symmetry_checked']}
    e['n_covale_rows_pos145']=len(r['covale_rows_pos145'])
    e['covale_pos145_partners']=';'.join(sorted({(c['ptnr2_auth_comp_id'] if c['ptnr1_auth_comp_id'] in ('CYS','SER') else c['ptnr1_auth_comp_id']) or '' for c in r['covale_rows_pos145']}))
    e['nonpoly_comps']=';'.join(sorted(r['nonpoly_comps']))
    e['drugs_present']=';'.join(sorted(DRUGS[c] for c in r['nonpoly_comps'] if c in DRUGS))
    ents.append(e)
    for p in r['pairs']: rows.append({'entry':r['entry']}|p)
E=pd.DataFrame(ents); C=pd.DataFrame(rows)
E['deposit_year']=E.deposit_date.str[:4]; E['release_year']=E.initial_release.str[:4]
E['xchem_like_id']=E.entry.str.match(r'^(5R|5S|7G|7H|7I|5Q)')
C=C.merge(E[['entry','process_site','deposit_year','release_year','method','resolution','deposit_group','xchem_like_id']],on='entry')
C['category']=[('peptide_polymer' if t=='polymer_nonMpro' else cat(c)) for c,t in zip(C.lig_comp,C.partner_type)]
C['drug']=C.lig_comp.map(lambda c:DRUGS.get(c,''))
C['lig_name']=C.lig_comp.map(lambda c:(ccd.get(c,{}).get('name') or '')[:90])
C['size_class']=C.lig_heavy_atoms.map(lambda n:'fragment(<=20 heavy)' if n is not None and n==n and n<=20 else 'drug-like(>20 heavy)')
C['dist_band']=pd.cut(C.dist,[0,1.6,2.0,2.4,2.6],right=False,labels=['<1.6','1.6-<2.0','2.0-<2.4','2.4-2.6'])
# nearest partner per nucleophile atom & element class (C vs other)
C['el_class']=C.lig_el.map(lambda x:'C' if x=='C' else x)
C=C.sort_values('dist'); C['nearest_for_nuc']=~C.duplicated(['entry','nuc_chain','nuc_seq','nuc_ins','nuc_atom','el_class'])
C.to_csv(f'{D}/out/contacts_all.tsv',sep='\t',index=False); E.to_csv(f'{D}/out/entries.tsv',sep='\t',index=False)
print(len(E),len(C)); print(C.groupby(['nuc_res','nsp5_pos','nuc_atom']).size().sort_values().tail(20)); print(C.category.value_counts())
