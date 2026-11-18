"""Step 6a. Classify nucleophile-ligand contacts against the MAXIT 11.400 _allowed_bonds link windows in
Varman, distance-table note v1.1 (doi:10.5281/zenodo.24295951), Table 1, bounds inclusive ("at or between"):
C-S 1.55-2.00, C-O 1.10-1.55, C-N 1.20-1.60, B-O 1.20-1.65, B-N 1.20-1.65, S-S 1.85-2.20, O-P 1.40-1.75,
N-P 1.54-1.77, O-S 1.35-1.80, N-S 1.45-1.90, C-Se 1.75-2.15 (C-Se cannot occur with an S/O/N nucleophile).
Pairs are unordered: O-S covers ligand S-Ser OG and ligand O-Cys SG. Other element pairs are 'no window'.
Universe: residue pairs (nucleophile residue, ligand residue) as in script 14 -- recorded = covale row joins them
(distance = pdbx_dist_value, closest row); unrecorded = no covale row and _pdbx_validate_close_contact lists a
nucleophile atom within <=2.2 A of a ligand atom of element C, B, P, S, Si, Se or As (closest listed pair).
Script 14 stopped unrecorded pairs at <2.0 A; here they are extended to the parser's 2.2 A limit so that contacts at
exactly 2.00 A (inside the C-S window) and 2.00-2.2 A are counted. Script 14 also ignored unrecorded contacts with a
ligand O or N; for Cys SG these are O-S and N-S pairs, which have windows, so unrecorded SG-ligand O/N close contacts
(<=2.2 A) are added for residue pairs not already in the table (src 'cc_SG_O/N'). Existing rows are not changed.
Writes out/step6/contacts_windowed.tsv and out/step6/29_stdout.txt."""
import json, pandas as pd, numpy as np
O='/workspace/ligand_code_effect'; A=f'{O}/archive'; OUT=f'{O}/out/step6'
C=json.load(open(f'{O}/out/ccd_parsed_archive.json'))
MAIN={('CYS','SG'),('SER','OG'),('LYS','NZ'),('HIS','NE2'),('HIS','ND1')}
EL={'C','B','P','S','Si','Se','As'}
NUCEL={'SG':'S','OG':'O','NZ':'N','NE2':'N','ND1':'N'}
T1={('C','S'):(1.55,2.00),('C','O'):(1.10,1.55),('C','N'):(1.20,1.60),('B','O'):(1.20,1.65),('B','N'):(1.20,1.65),('S','S'):(1.85,2.20),
    ('O','P'):(1.40,1.75),('N','P'):(1.54,1.77),('O','S'):(1.35,1.80),('N','S'):(1.45,1.90),('C','Se'):(1.75,2.15)}   # Table 1 as published
WIN={}
for (a,b),w in T1.items(): WIN[(a,b)]=w; WIN[(b,a)]=w      # keyed (nucleophile element, ligand element), both orders
def elem(c,a):
    r=C.get(c)
    return r['el'].get(a) if r and a in r['el'] else None
B=pd.read_csv(f'{A}/bonds_all.tsv',sep='\t',dtype={'entry':str},low_memory=False)
B['src']='bonds_all'
# extension: unrecorded pairs with closest listed contact 2.0-2.2 A (not in bonds_all, which stops at <2.0)
ext=[]
for l in open(f'{A}/parsed.jsonl'):
    r=json.loads(l)
    if 'err' in r: continue
    rec={(x[0],x[2],x[3],x[4] or '',x[5],x[7],x[8],x[9] or '') for x in r['cov'] if (x[0],x[1]) in MAIN}
    best={}
    for x in r['cc']:
        if (x[0],x[1]) not in MAIN: continue
        el=elem(x[5],x[6])
        if el not in EL: continue
        k=(x[0],x[2],x[3],x[4] or '',x[5],x[7],x[8],x[9] or '')
        if k in rec: continue
        if k not in best or x[10]<best[k][0]: best[k]=(x[10],x[1],x[6],el)
    for k,v in best.items():
        if v[0]<2.0: continue   # already in bonds_all
        ext.append(dict(entry=r['e'],nuc_res=k[0],nuc_chain=k[1],nuc_seq=k[2],nuc_ins=k[3],lig_comp=k[4],lig_chain=k[5],lig_seq=k[6],lig_ins=k[7],
                        dist=v[0],nuc_atom=v[1],lig_atom=v[2],lig_el=v[3],recorded=False,ps=r.get('ps'),dep=r.get('dep'),au_last=r.get('au_last'),src='cc_2.0-2.2'))
    have=set(best)
    son={}
    for x in r['cc']:
        if (x[0],x[1])!=('CYS','SG'): continue
        el=elem(x[5],x[6])
        if el not in ('O','N'): continue
        k=(x[0],x[2],x[3],x[4] or '',x[5],x[7],x[8],x[9] or '')
        if k in rec or k in have: continue
        if k not in son or x[10]<son[k][0]: son[k]=(x[10],x[1],x[6],el)
    for k,v in son.items():
        ext.append(dict(entry=r['e'],nuc_res=k[0],nuc_chain=k[1],nuc_seq=k[2],nuc_ins=k[3],lig_comp=k[4],lig_chain=k[5],lig_seq=k[6],lig_ins=k[7],
                        dist=v[0],nuc_atom=v[1],lig_atom=v[2],lig_el=v[3],recorded=False,ps=r.get('ps'),dep=r.get('dep'),au_last=r.get('au_last'),src='cc_SG_O/N'))
X=pd.concat([B,pd.DataFrame(ext)],ignore_index=True)
X['nuc_el']=X.nuc_atom.map(NUCEL)
X['pair']=X.lig_el.fillna('?')+'-'+X.nuc_el      # ligand element first, e.g. C-S
def wstat(r):
    w=WIN.get((r.nuc_el,r.lig_el))
    if w is None: return 'no window'
    if pd.isna(r.dist): return 'no distance'
    return 'below' if r.dist<w[0] else ('above' if r.dist>w[1] else 'inside')
X['window']=X.apply(wstat,axis=1)
X['inside']=X.window=='inside'
CD=pd.read_csv(f'{O}/out/archive_codes.tsv',sep='\t',index_col=0)
X['cat']=X.lig_comp.map(CD['cat']).fillna('unclassified')
X.to_csv(f'{OUT}/contacts_windowed.tsv',sep='\t',index=False)
L=[]; P=lambda *a: L.append(' '.join(map(str,a)))
P('residue-pair contacts: bonds_all',len(B),'+ unrecorded 2.0-2.2 A extension',int((X.src=='cc_2.0-2.2').sum()),'+ unrecorded Cys SG-ligand O/N',int((X.src=='cc_SG_O/N').sum()),'=',len(X))
P('unrecorded contacts at exactly 2.00 A (extension):',int(((X.src=='cc_2.0-2.2')&(X.dist==2.0)).sum()))
t=X.groupby(['pair','window','recorded']).size().unstack(fill_value=0); t.columns=['unrecorded' if not c else 'recorded' for c in t.columns]
P('\n## contacts by element pair (ligand-nucleophile) and window'); P(t.to_string())
for lab,s in [('ALL inside-window',X[X.inside]),('inside, < 2.0 A',X[X.inside&(X.dist<2.0)])]:
    n=len(s); u=int((~s.recorded).sum()); P(f'\n{lab}: {n} contacts, unrecorded {u} ({100*u/n:.2f}%), recorded {n-u} ({100*(n-u)/n:.2f}%)')
P('\n## inside-window by pair: total, unrecorded, %')
for p,s in X[X.inside].groupby('pair'): P(f'   {p}: {len(s)} {int((~s.recorded).sum())} {100*(~s.recorded).mean():.2f}%')
lt=X[(X.dist<2.0)&(X.src=='bonds_all')]
P('\n## bonds < 2.0 A (the 43,488 of the addendum) by window status: total, unrecorded')
for w,s in lt.groupby('window'): P(f'   {w}: {len(s)} {int((~s.recorded).sum())}')
P('unrecorded < 2.0 A total',int((~lt.recorded).sum()))
P('\n## no-window pairs < 2.0 A: total, unrecorded')
for p,s in lt[lt.window=='no window'].groupby('pair'): P(f'   {p}: {len(s)} {int((~s.recorded).sum())}')
I=X[X.inside]
P('\n## inside-window unrecorded by ligand category: total, unrecorded, %')
for c,s in I.groupby('cat'): P(f'   {c}: {len(s)} {int((~s.recorded).sum())} {100*(~s.recorded).mean():.1f}%')
P('\n## inside-window unrecorded by processing site: total, unrecorded, %')
for c,s in I.groupby(I.ps.fillna('NA')): P(f'   {c}: {len(s)} {int((~s.recorded).sum())} {100*(~s.recorded).mean():.1f}%')
nb=I[~I.cat.isin(['bilin','heme'])]; P(f'\ninside-window, excluding bilins and hemes: {len(nb)} contacts, unrecorded {int((~nb.recorded).sum())} ({100*(~nb.recorded).mean():.2f}%)')
dl=I[I.cat=='drug-like']; P(f'inside-window, drug-like: {len(dl)} contacts, unrecorded {int((~dl.recorded).sum())} ({100*(~dl.recorded).mean():.2f}%)')
peb=I[I.lig_comp=='PEB']; P(f'inside-window PEB: {len(peb)} contacts, unrecorded {int((~peb.recorded).sum())}')
P(f'inside-window unrecorded where the same nucleophile residue has a covale row to another ligand residue: {int((~I.recorded & (I.nuc_cov_other==True)).sum())}')
P('inside-window unrecorded entries:',I[~I.recorded].entry.nunique(),' codes:',I[~I.recorded].lig_comp.nunique())

# boundary sensitivity: close-contact distances are printed to 2 decimals, pdbx_dist_value to 3
def onb(r):
    w=WIN.get((r.nuc_el,r.lig_el)); return w is not None and (abs(r.dist-w[0])<1e-9 or abs(r.dist-w[1])<1e-9)
bd=I[I.apply(onb,axis=1)]
P(f'\ninside-window contacts lying exactly on a window bound: {len(bd)} (unrecorded {int((~bd.recorded).sum())}); at 2.00 A C-S: {int(((bd.pair=="C-S")&(bd.dist==2.0)).sum())}')
ex=I.drop(bd.index); P(f'inside-window with bounds exclusive: {len(ex)} contacts, unrecorded {int((~ex.recorded).sum())} ({100*(~ex.recorded).mean():.2f}%)')
# rule cross-check: windowed pairs outside the window
wo=X[X.window.isin(['below','above'])]; P(f'windowed-pair contacts outside the window: {len(wo)}, recorded {int(wo.recorded.sum())} ({100*wo.recorded.mean():.1f}%)')
P('\n## inside-window unrecorded by site, bilins vs other: site, bilin unrecorded, other unrecorded, other total')
for c,s in I.groupby(I.ps.fillna('NA')):
    b=s[s.cat=='bilin']; o=s[s.cat!='bilin']; P(f'   {c}: {int((~b.recorded).sum())} {int((~o.recorded).sum())} {len(o)}')
nbl=I[I.cat!='bilin']; u=nbl[~nbl.recorded]; P(f'inside-window non-bilin: {len(nbl)} contacts, unrecorded {len(u)} ({100*len(u)/len(nbl):.2f}%), at PDBj {int((u.ps=="PDBJ").sum())}')

P('\n## in-window contacts by source: total, unrecorded')
for c,s in I.groupby('src'): P(f'   {c}: {len(s)} {int((~s.recorded).sum())}')
P('in-window contacts at >= 2.0 A (outside the unit analysis):',int((I.dist>=2.0).sum()),' unrecorded:',int(((I.dist>=2.0)&~I.recorded).sum()))
P('\n## in-window by Table 1 pair (unordered): total, unrecorded, %')
I2=I.copy(); I2['t1']=[('-'.join(k) if k in T1 else '-'.join(k[::-1])) for k in zip(I2.lig_el,I2.nuc_el)]
for c,s in I2.groupby('t1'): P(f'   {c}: {len(s)} {int((~s.recorded).sum())} {100*(~s.recorded).mean():.2f}%')
# leritrelvir
lr=X[X.lig_comp.isin(['A1EN0','7ON'])&X.entry.isin(['23LW','23LX','23LY','23LZ','23MA','23MB','23MC','23MD','23ME','23MF','23MG','23MH','23MI','23MJ','23MK','23ML'])]
for c,s in lr.groupby('lig_comp'): P(f'leritrelvir 2026 batch {c}: contacts {len(s)}, pair {sorted(set(s.pair))}, header distance {s.dist.min():.3f}-{s.dist.max():.3f}, window {s.window.value_counts().to_dict()}, recorded {int(s.recorded.sum())}')
open(f'{OUT}/29_stdout.txt','w').write('\n'.join(L)); print('\n'.join(L))
