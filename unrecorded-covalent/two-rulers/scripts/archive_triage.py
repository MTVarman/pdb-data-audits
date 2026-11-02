"""Triage the whole-PDB sweep.

Usage: python3 -I archive_triage.py archive.jsonl ccd_dir out_prefix

Rows: one per (entry, ligand instance, protein nucleophile atom), closest ligand atom,
pairs with no _struct_conn record for that atom pair. Columns added:
  band            <2.0 (comparable to the published scan floor) or 2.0-2.5 (new)
  occupancy_class full / partial / altloc_crossing (same rule as analyze.py)
  ligand_atom_tier, ligand_atom_kind   from the CCD definition of the ligand atom
      A  reactive centre: non-aromatic C=C / C=O / C=N / C#N carbon, open valence,
         B, P, Se, or thiol S
      B  sp3 carbon carrying a heteroatom (O/N/S/halogen), not a methyl: covers
         reacted forms drawn as their product (hemithioacetal as alcohol etc.)
      C  methyl, aromatic, or sp3 carbon with only C/H neighbours: no bond as modelled
  family          class1_cofactor_ptm / class3_nucleotide_phosphate / '' (approximation
                  of the original triage families by CCD id; listed in COF below)
  probe_line_A    covalent radius sum + 0.2 A (probe's bond-perception distance)
"""
import json, os, sys, urllib.request, csv, collections
from concurrent.futures import ThreadPoolExecutor
import gemmi

COV = {'C': 0.77, 'N': 0.70, 'O': 0.66, 'S': 1.04, 'B': 0.86, 'P': 1.10, 'Se': 1.17}
VAL = {'C': 4, 'B': 3, 'P': 5, 'S': 2, 'Se': 2}
ORD = {'SING': 1, 'DOUB': 2, 'TRIP': 3, 'AROM': 1.5}
# class 1: covalent cofactors and PTM-like attachments (approximating the 11 families)
COF = set('''PEB PUB PCB CYC PVB BLA BLV LBV BLR RBV
HEC HEM HEB HAS HEO
FAD FMN FDA FAE
PLP PMP P1T PLR PDP
BTN BTI LPA LPB
RET
PLM MYR STE OLA LAU DAO PAM MYA
FAR GER GRG FPP GPP
PNS PPS 4PS
GSH GTT GSF GDS
ACE FOR'''.split())
NUCLEOTIDE_LIKE = set('''ATP ADP AMP ANP ACP APC AGS ADX GTP GDP GMP GNP GCP GSP GAV
CTP CDP CMP UTP UDP UMP TTP TDP TMP DTP DCP DGT DUT ITP IDP
NAD NAP NDP NAI FAD FMN COA ACO SAM SAH'''.split())
CLASS3_RES = {'ASP', 'GLU', 'LYS', 'SER', 'THR'}

def coexist(pa, la, po, lo):
    if pa and la and pa != la:
        return 'altloc_crossing'
    return 'full' if po >= 0.99 and lo >= 0.99 else 'partial'

def ccd_fetch(comp, d):
    p = os.path.join(d, f'{comp}.cif')
    if not os.path.exists(p):
        try:
            data = urllib.request.urlopen(f'https://files.rcsb.org/ligands/download/{comp}.cif', timeout=60).read()
            open(p, 'wb').write(data)
        except Exception:
            return None
    return p

def ccd_atoms(path):
    b = gemmi.cif.read_file(path).sole_block()
    atoms = {r[0].strip('"'): (r[1].capitalize() if len(r[1]) > 1 else r[1], r[2])
             for r in b.find('_chem_comp_atom.', ['atom_id', 'type_symbol', 'pdbx_aromatic_flag'])}
    nb = collections.defaultdict(list)
    for r in b.find('_chem_comp_bond.', ['atom_id_1', 'atom_id_2', 'value_order', 'pdbx_aromatic_flag']):
        a1, a2 = r[0].strip('"'), r[1].strip('"')
        nb[a1].append((a2, r[2].upper(), r[3])); nb[a2].append((a1, r[2].upper(), r[3]))
    name = (b.find_value('_chem_comp.name') or '').strip('"; \n')
    return atoms, nb, name

def classify(atoms, nb, atom):
    if atom not in atoms:
        return '?', 'atom not in CCD'
    el, arom = atoms[atom]
    bonds = nb.get(atom, [])
    val = sum(ORD.get(o, 1) for _, o, _ in bonds)
    els = [atoms.get(o, ('?', ''))[0] for o, _, _ in bonds]
    nH = sum(1 for e in els if e in ('H', 'D'))
    heavy = [e for e in els if e not in ('H', 'D')]
    orders = [o for _, o, _ in bonds]
    if el in ('B', 'P', 'Se'):
        return 'A', f'{el}'
    if el == 'S':
        if nH >= 1 or val < VAL['S']:
            return 'A', 'thiol / open S'
        return 'C' if val <= 2 else 'B', ('thioether / ring S' if val <= 2 else 'oxidised S (sulfonyl etc.)')
    if el != 'C':
        return 'C', f'{el}'
    if val < 4:
        return 'A', 'open valence C'
    if arom == 'Y' or any(f == 'Y' for _, _, f in bonds):
        return 'C', 'aromatic C'
    if 'TRIP' in orders:
        return 'A', 'C#N / C#C'
    if 'DOUB' in orders:
        dbl = [atoms.get(o, ('?', ''))[0] for o, ord_, _ in bonds if ord_ == 'DOUB']
        return 'A', 'C=' + ('O' if 'O' in dbl else 'N' if 'N' in dbl else 'S' if 'S' in dbl else 'C')
    if any(e in ('O', 'N', 'S', 'F', 'Cl', 'CL', 'Br', 'BR', 'I', 'Se', 'SE', 'B', 'P') for e in heavy):
        return 'B', 'sp3 C with heteroatom'
    # reacted warheads are often drawn in product form with H in place of the protein
    # bond (6T5V MKZ: acrylamide drawn as propan-1-one, bonded carbon = 'methyl').
    # Keep sp3 carbons within two bonds of an atom carrying a double bond to O/N/S.
    def has_hetero_double(a):
        return any(o == 'DOUB' and atoms.get(x, ('?', ''))[0] in ('O', 'N', 'S') for x, o, _ in nb.get(a, []))
    n1 = [x for x, _, _ in bonds if atoms.get(x, ('?', ''))[0] not in ('H', 'D')]
    n2 = {y for x in n1 for y, _, _ in nb.get(x, []) if y != atom and atoms.get(y, ('?', ''))[0] not in ('H', 'D')}
    if any(has_hetero_double(x) for x in n1):
        return 'B', 'sp3 C alpha to C=O/S=O/C=N (reacted-warhead pattern)'
    if any(has_hetero_double(y) for y in n2):
        return 'B', 'sp3 C beta to C=O/S=O/C=N (reacted-warhead pattern)'
    if nH >= 3:
        return 'C', 'methyl'
    return 'C', 'sp3 C, C/H only'

def main():
    src, ccd_dir, pre = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(ccd_dir, exist_ok=True)
    st = collections.Counter(); rows = {}; meta = {}
    for line in open(src):
        r = json.loads(line)
        st[r.get('status')] += 1
        meta[r['pdbid']] = (r.get('method', ''), r.get('resolution'))
        for c in r.get('contacts', []):
            if c['conn_pair']:
                continue
            if c.get('lig_entity') == 'branched':
                st['branched_contact_rows_set_aside'] += 1
                continue  # glycans: reported separately, not part of the organic-ligand scan
            k = (r['pdbid'], c['lig'], c['prot'])
            if k not in rows or c['dist'] < rows[k]['dist']:
                rows[k] = {**c, 'pdbid': r['pdbid']}
    comps = sorted({k[1].split('_')[0] for k in rows})
    with ThreadPoolExecutor(16) as ex:
        paths = dict(zip(comps, ex.map(lambda c: ccd_fetch(c, ccd_dir), comps)))
    cache = {}
    out = []
    for (pdbid, lig, prot), c in rows.items():
        comp = lig.split('_')[0]
        if comp not in cache:
            cache[comp] = ccd_atoms(paths[comp]) if paths.get(comp) else None
        if cache[comp]:
            tier, kind = classify(cache[comp][0], cache[comp][1], c['lig_atom'])
            cname = cache[comp][2]
        else:
            tier, kind, cname = '?', 'CCD unavailable', ''
        pres = prot.split(':')[1][:3]
        pel = 'Se' if prot.endswith(':SE') else prot.split(':')[2][0]
        fam = ''
        if comp in COF:
            fam = 'class1_cofactor_ptm'
        elif (comp in NUCLEOTIDE_LIKE or c['lig_el'] == 'P') and pres in CLASS3_RES:
            fam = 'class3_nucleotide_phosphate'
        method, res = meta[pdbid]
        d = c['dist']
        out.append({'pdbid': pdbid, 'ligand': lig, 'ligand_name': cname[:80], 'protein_atom': prot,
                    'residue': pres, 'ligand_atom': c['lig_atom'], 'pair': f"{c['lig_el']}-{pel}",
                    'dist_A': d, 'band': '<2.0' if d < 2.0 else '2.0-2.5',
                    'probe_line_A': round(COV.get(c['lig_el'], 0) + COV.get(pel, 0) + 0.2, 2),
                    'occupancy_class': coexist(c['prot_altloc'], c['lig_altloc'], c['prot_occ'], c['lig_occ']),
                    'prot_altloc': c['prot_altloc'], 'lig_altloc': c['lig_altloc'],
                    'prot_occ': c['prot_occ'], 'lig_occ': c['lig_occ'],
                    'ligand_atom_tier': tier, 'ligand_atom_kind': kind, 'family': fam,
                    'residue_level_conn': ' '.join(c['conn_residues']),
                    'method': method, 'resolution': res})
    out.sort(key=lambda x: (x['dist_A'], x['pdbid']))
    with open(f'{pre}_contacts.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
    print('entries', sum(st.values()), dict(st))
    for band in ('<2.0', '2.0-2.5'):
        b = [x for x in out if x['band'] == band]
        print(f'\n== {band}: rows {len(b)}, entries {len({x["pdbid"] for x in b})}')
        print(' occupancy', collections.Counter(x['occupancy_class'] for x in b))
        f = [x for x in b if x['occupancy_class'] == 'full']
        print(' full: tier', collections.Counter(x['ligand_atom_tier'] for x in f), 'family', collections.Counter(x['family'] for x in f))
        cand = [x for x in f if x['ligand_atom_tier'] in ('A', 'B') and not x['family']]
        print(' full, tier A/B, not class1/3:', len(cand), 'entries', len({x['pdbid'] for x in cand}))
        print('   by residue', collections.Counter(x['residue'] for x in cand).most_common(12))
        xr = [x for x in cand if 'X-RAY' in x['method'] and x['resolution'] and x['resolution'] <= 2.5]
        print('   ...X-ray <=2.5 A:', len(xr), 'entries', len({x['pdbid'] for x in xr}),
              '| Cys:', sum(1 for x in xr if x['residue'] == 'CYS'), 'entries', len({x['pdbid'] for x in xr if x['residue'] == 'CYS'}))

if __name__ == '__main__':
    main()
