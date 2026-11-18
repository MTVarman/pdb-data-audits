"""For every full-occupancy serine/threonine-to-boron contact beyond the annotation limit (1.65 A) with no record:
does the validation report list a clash on the boron atom, and does it list the same hydroxyl oxygen (or its
hydrogen) against another atom of the same boronate?  Is that atom bonded to the boron in the Chemical
Component Dictionary?
Usage: python3 -I boron_neighbour_clashes.py data/sample_contacts.csv.gz val_dir cif_dir ccd_dir out.csv
val_dir holds <id>_validation.xml.gz, cif_dir <ID>.cif, ccd_dir <COMP>.cif (fetched from files.rcsb.org when absent)."""
import csv, sys, gzip, os, collections, re, urllib.request, gemmi
import xml.etree.ElementTree as ET
src, vdir, cdir, ccd, outp = sys.argv[1:6]
def op(path): return gzip.open(path, 'rt', newline='') if path.endswith('.gz') else open(path, newline='')
def fetch(url, path):
    if not os.path.exists(path) or os.path.getsize(path) < 500:
        os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
        open(path, 'wb').write(urllib.request.urlopen(url, timeout=120).read())
    return path
rows = [r for r in csv.DictReader(op(src)) if r['pair'] == 'B-O' and r['residue'] in ('SER', 'THR')
        and r['occupancy_class'] == 'full' and r['recorded_covale'] != 'True' and float(r['dist_A']) > 1.65]
cache, bondcache = {}, {}
def load(pid):
    if pid not in cache:
        p = pid.lower()
        root = ET.parse(gzip.open(fetch(f'https://files.rcsb.org/pub/pdb/validation_reports/{p[1:3]}/{p}/{p}_validation.xml.gz', f'{vdir}/{p}_validation.xml.gz'))).getroot()
        cl = collections.defaultdict(list)
        for sg in root.iter('ModelledSubgroup'):
            for c in sg.findall('clash'):
                cl[c.get('cid')].append((sg.get('chain'), sg.get('resname'), sg.get('resnum'), c.get('atom'), c.get('dist')))
        st = gemmi.read_structure(fetch(f'https://files.rcsb.org/download/{pid}.cif', f'{cdir}/{pid}.cif')); st.setup_entities()
        el = {(ch.name, res.name, str(res.seqid.num), a.name): a.element.name for ch in st[0] for res in ch for a in res}
        cache[pid] = (cl, el)
    return cache[pid]
def bonded(comp):
    if comp not in bondcache:
        g = collections.defaultdict(set)
        for r in gemmi.cif.read(fetch(f'https://files.rcsb.org/ligands/download/{comp}.cif', f'{ccd}/{comp}.cif')).sole_block().find('_chem_comp_bond.', ['atom_id_1', 'atom_id_2']):
            a, c = gemmi.cif.as_string(r[0]), gemmi.cif.as_string(r[1]); g[a].add(c); g[c].add(a)
        bondcache[comp] = g
    return bondcache[comp]
out = []
for r in rows:
    cl, el = load(r['pdbid'])
    pch, presx, patom = r['protein_atom'].split(':'); pres = r['residue']; pnum = presx[len(pres):]
    comp, lch, lnum = r['ligand'].rsplit('_', 2)
    on_boron = 0; partners = []; nuc = []
    for cid, v in cl.items():
        a = [x for x in v if x[0] == pch and x[1] == pres and x[2] == pnum]
        b = [x for x in v if x[0] == lch and x[1] == comp and x[2] == lnum]
        if a and b:
            e = el.get((lch, comp, lnum, b[0][3]), 'H' if b[0][3].startswith('H') else '?')
            if b[0][3] == r['ligand_atom']: on_boron += 1
            partners.append(f"{a[0][3]}~{b[0][3]}({e}) {b[0][4]}A")
            if re.match(r'^(OG1?|HG1?)$', a[0][3]): nuc.append(b[0][3])
    d = float(r['dist_A'])
    out.append({'pdbid': r['pdbid'], 'protein_atom': r['protein_atom'], 'ligand': r['ligand'], 'boron_atom': r['ligand_atom'], 'dist_A': r['dist_A'],
                'zone': 'between_limit_and_probe_line' if d <= 1.72 else ('past_probe_line_to_2.12' if d <= 2.12 else 'beyond_2.12'),
                'clash_records_naming_the_boron': on_boron, 'n_clash_records_residue_to_same_boronate': len(partners),
                'clash_partners': '; '.join(partners), 'resolution': r['resolution'], 'deposited': r['deposited'],
                'n_records_nucleophile_O_or_H_to_same_boronate': len(nuc),
                'n_of_those_naming_an_atom_bonded_to_the_boron': sum(1 for a in nuc if a in bonded(comp)[r['ligand_atom']])})
with open(outp, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
for z in ('between_limit_and_probe_line', 'past_probe_line_to_2.12', 'beyond_2.12'):
    s = [o for o in out if o['zone'] == z]
    print(z, 'contacts', len(s), 'entries', len({o['pdbid'] for o in s}), '| clash naming boron', sum(o['clash_records_naming_the_boron'] for o in s),
          '| nucleophile O or H listed against another boronate atom:', sum(1 for o in s if o['n_records_nucleophile_O_or_H_to_same_boronate']),
          '| records', sum(o['n_records_nucleophile_O_or_H_to_same_boronate'] for o in s), '| naming an atom bonded to the boron', sum(o['n_of_those_naming_an_atom_bonded_to_the_boron'] for o in s))
