"""Flatten scan output into contacts TSV with ligand category labels."""
import json, os, re
import pandas as pd
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
ccd = json.load(open(f'{D}/out/ccd_info.json'))
drugs = pd.read_csv(f'{D}/out/drug_comps.tsv', sep='\t').set_index('comp_id')
HEME = {'HEM','HEC','HEA','HEB','HEO','HDD','HDM','HAS','HEV','HE5','DHE','SRM','COH','CCH','HCO','HNI','ZNH','1CP','FEC','HIF','VEA','VER','HMG','HE6','MHM','PNI','HKL','6CO','6CQ','CLN','ISW','MD9','HEG','HES','HP5','UFE','T2Z','HCB','HNN','HME','HRM','HFM','HEW'}
BILIN = {'BLA','LBV','LBW','CYC','PEB','PUB','3Q8','EL5','BLR','RBF','PCB','URO','DBV','15P','BV','PVN','BPB'}
FLAVIN = {'FAD','FMN','RBF','FDA','FAE','FMA','6FA','8FM','FNR','FAS','FSH','FAO','C1L','LFN'}
THIOL_ADD = {'BME','DTT','DTU','DTV','GSH','CYS','2ME','MRY','CME','TCE','DTD','MPO','CSD','DCY','GTT','GDS','GSF','GSN','MTL','ACM','MPD','EDO','GOL','PEG','PGE','1PE','PG4','ACT','FMT','DMS','EPE','MES','TRS','BU3','CIT','SO4','PO4','IMD','NHE','MLI','IPA','EOH','MOH','BCN','B3P','ACY','AKG','TLA','TAR','PGO','P6G','SIN','12P','15P'}
LIPID_ACYL = {'PLM','MYR','STE','OLA','PAM','HXA','DAO','DKA','UND','LAU','PLC','COA','ACO','SCA','MLC','CAA','HXC','DCC','MCA'}
METALS = re.compile(r'\b(Fe|Ni|Co|Cu|Zn|Mn|Mo|W|V|Mg|Ca|Cd|Hg|Pt|Au|Ag|Pd|Ru|Rh|Ir|Os|Re|Ga|Al|Sn|Pb|Cr|Na|K|Li|Tl|Sb|Bi|Te|Ti|Zr|Y|La|Ce|Eu|Gd|Tb|Yb|Lu|Sm|U|As)\d*\b')
def cat(c):
    if c in drugs.index: return 'named_drug'
    info = ccd.get(c, {}); n = (info.get('name') or '').upper(); f = info.get('formula') or ''
    if c in HEME or 'HEME' in n or 'PORPHYRIN' in n or 'CHLOROPHYLL' in n: return 'heme_porphyrin'
    if c in BILIN or 'BILIN' in n or 'BILIVERDIN' in n or 'PHYCO' in n: return 'bilin'
    if c == 'RET' or 'RETINAL' in n: return 'retinal'
    if c == 'HC4': return 'pyp_chromophore'
    if c in FLAVIN or 'FLAVIN' in n: return 'flavin'
    if METALS.search(f): return 'metal_containing'
    if c in THIOL_ADD: return 'additive_or_thiol'
    if c in LIPID_ACYL or 'COENZYME A' in n: return 'lipid_acyl_coa'
    return 'other_ligand'
rows = []; ent = []
for fn in ('scan.jsonl', 'scan_big.jsonl'):
    for l in open(f'{D}/out/{fn}'):
        r = json.loads(l)
        if fn == 'scan.jsonl' and r['status'] == 'skipped_too_large': continue
        ent.append({k: r.get(k) for k in ('entry','status','process_site','deposit_date','method','resolution','initial_release','last_revision','n_atoms','symmetry_checked')} | {'n_pairs': len(r.get('pairs', []))})
        for p in r.get('pairs', []):
            rows.append({'entry': r['entry'], 'process_site': r.get('process_site'), 'deposit_year': (r.get('deposit_date') or '')[:4],
                         'method': r.get('method'), 'resolution': r.get('resolution')} | p)
df = pd.DataFrame(rows); E = pd.DataFrame(ent)
df['category'] = df.lig_comp.map(cat)
df['lig_name'] = df.lig_comp.map(lambda c: (ccd.get(c, {}).get('name') or '')[:90])
df['drug'] = df.lig_comp.map(lambda c: drugs.drug.get(c, ''))
df['drug_class'] = df.lig_comp.map(lambda c: drugs['class'].get(c, ''))
df['bonded_elsewhere'] = (~df.recorded) & ((df.cys_cov_other.fillna('') != '') | (df.lig_cov_other_protein_res.fillna('') != ''))
df.to_csv(f'{D}/out/contacts_all.tsv', sep='\t', index=False)
E.to_csv(f'{D}/out/entries_scanned.tsv', sep='\t', index=False)
print(len(E), E.status.value_counts().to_dict(), len(df))
print(df[df.el_class=='C'].category.value_counts())
