import csv, sys
src = sys.argv[1]
rows = [x for x in csv.DictReader(open(src)) if x['organism'] == 'Severe acute respiratory syndrome coronavirus 2']
def assign(x):
    r = (x['residue_name'], x['residue_number'])
    if r == ('CYS', '145'):
        return 'Mpro (nsp5)', 'Catalytic Cys145 thiol bonded to the inhibitor warhead carbon'
    if r == ('CYS', '111'):
        return 'PLpro (nsp3)', 'Catalytic Cys111 thiol bonded to the inhibitor warhead carbon (irreversible PLpro inhibitor series)'
    if r == ('GLU', '203'):
        return 'PLpro (nsp3), C111S', 'PanDDA fragment screen; fragment has no electrophile for an ester to Glu203; same 1.375-1.379 A contact for three unrelated fragments; weak ligand density. Likely modelling artifact'
    if r == ('ASP', '164'):
        return 'PLpro (nsp3), C111S', 'PanDDA fragment screen; phenol fragment has no electrophile for an ester to Asp164; density supported. Chemically unexpected, check the map'
order = {'Mpro (nsp5)': 0, 'PLpro (nsp3)': 1, 'PLpro (nsp3), C111S': 2}
out = []
for x in rows:
    prot, note = assign(x)
    out.append({
        'pdb_id': x['entry'], 'protein': prot, 'chain': x['chain'],
        'residue': f"{x['residue_name']} {x['residue_number']}", 'protein_atom': x['protein_atom'],
        'ligand': x['ligand_comp'], 'ligand_residue_number': x['ligand_number'], 'ligand_atom': x['ligand_atom'],
        'distance_A': x['distance_A'], 'resolution_A': f"{float(x['resolution_A']):.2f}", 'release_date': x['release_date'],
        'protein_rscc': x['protein_rscc'], 'protein_rsr': x['protein_rsr'],
        'ligand_rscc': x['ligand_rscc'], 'ligand_rsr': x['ligand_rsr'], 'ligand_occupancy': x['ligand_occupancy'],
        'density_class': x['density_class'], 'clash_class': x['clash_class'], 'interpretation': note,
    })
out.sort(key=lambda d: (order[d['protein']], d['pdb_id'], d['chain']))
w = csv.DictWriter(open('sars_cov2_unrecorded_covalent_appendix.csv', 'w', newline=''), fieldnames=list(out[0]))
w.writeheader(); w.writerows(out)
import collections
print(len(out), len({d['pdb_id'] for d in out}))
for p in order:
    s = [d for d in out if d['protein'] == p]
    print(p, len(s), 'sites', len({d['pdb_id'] for d in s}), 'entries', collections.Counter(d['density_class'] for d in s))
print('density_supported sites', sum(d['density_class']=='density_supported' for d in out), 'entries', len({d['pdb_id'] for d in out if d['density_class']=='density_supported'}))
mp=[d for d in out if d['protein']=='Mpro (nsp5)']
print('Mpro dist range', min(float(d['distance_A']) for d in mp), max(float(d['distance_A']) for d in mp))
print('Mpro years', collections.Counter(d['release_date'][:4] for d in mp))
print('Mpro ligands', collections.Counter(d['ligand'] for d in mp).most_common())
