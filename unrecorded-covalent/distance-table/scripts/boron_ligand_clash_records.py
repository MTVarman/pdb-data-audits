"""Every boron-containing ligand residue in the entries of the boron table: clash records on the residue,
and how many of them name a boron atom.  Elements are taken from the entry's own atom_site records.
Usage: python3 -I boron_ligand_clash_records.py boron_clash_check.csv cif_dir val_dir out.csv"""
import csv, sys, os, gzip, collections, urllib.request, gemmi
import xml.etree.ElementTree as ET
ids = sorted({r['pdbid'] for r in csv.DictReader(open(sys.argv[1]))})
cdir, vdir = sys.argv[2], sys.argv[3]
os.makedirs(cdir, exist_ok=True); os.makedirs(vdir, exist_ok=True)
def get(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 1000: return True
    for i in range(3):
        try:
            with urllib.request.urlopen(url, timeout=90) as r: open(path, 'wb').write(r.read()); return True
        except Exception as e: err = e
    print('  FAIL', url, err); return False
tot = collections.Counter(); partner = collections.Counter(); closest = []; table = []
for pid in ids:
    p = pid.lower()
    if not get(f'https://files.rcsb.org/download/{pid}.cif', f'{cdir}/{pid}.cif'): tot['cif missing'] += 1; continue
    if not get(f'https://files.rcsb.org/pub/pdb/validation_reports/{p[1:3]}/{p}/{p}_validation.xml.gz', f'{vdir}/{p}_validation.xml.gz'): tot['validation report missing'] += 1; continue
    st = gemmi.read_structure(f'{cdir}/{pid}.cif'); st.setup_entities()
    bres = {}
    for ch in st[0]:
        for res in ch:
            if res.het_flag == 'H' and any(a.element.name == 'B' for a in res):
                bres[(ch.name, str(res.seqid.num), res.name)] = {a.name: a.element.name for a in res}
    if not bres: tot['entries with no boron residue found'] += 1; continue
    tot['entries'] += 1; tot['boron-containing ligand residues'] += len(bres)
    root = ET.parse(gzip.open(f'{vdir}/{p}_validation.xml.gz')).getroot()
    per = {k: collections.Counter() for k in bres}
    for sg in root.iter('ModelledSubgroup'):
        k = (sg.get('chain'), sg.get('resnum'), sg.get('resname'))
        if k in bres:
            for cl in sg.findall('clash'):
                tot['clash records on those residues'] += 1
                an = cl.get('atom')
                el = bres[k].get(an, 'H' if an.startswith('H') else '?')
                partner[el] += 1; per[k][el] += 1
                if el == 'B': tot['clash records naming a boron atom'] += 1; closest.append((pid, k, cl.get('atom'), cl.get('dist')))
    for k in sorted(per):
        table.append({'pdbid': pid, 'chain': k[0], 'resnum': k[1], 'component': k[2], 'boron_atoms': ' '.join(a for a, e in bres[k].items() if e == 'B'),
                      'clash_records': sum(per[k].values()), 'clash_records_naming_boron': per[k]['B'],
                      'clash_records_by_element': ' '.join(f'{e}:{n}' for e, n in sorted(per[k].items()))})
with open(sys.argv[4], 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(table[0])); w.writeheader(); w.writerows(table)
print(dict(tot)); print('clash records by element of the ligand atom named:', dict(partner)); print(closest[:10])
