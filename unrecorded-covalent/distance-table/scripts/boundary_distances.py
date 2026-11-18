"""Exact distances for the pairs whose rounded distance (0.001 A) lies exactly on a boundary used in the note:
a link-table limit, the validator's bond line, the line + 0.40 A, or the close-contact cut-off 2.195 A.
The scan tables store distances rounded to 0.001 A; these pairs are re-measured from the deposited coordinates.
Archive pairs: model 1, symmetry copies included, matching conformers (as in the archive scan).
Sample pairs: asymmetric unit only, the two conformers named in the row.
Usage: python3 -I boundary_distances.py archive_residue_pairs.csv.gz sample_contacts.csv.gz cif_dir out.csv"""
import csv, gzip, os, re, sys, urllib.request
import gemmi
arch, samp, cifdir, dst = sys.argv[1:5]
MAXIT = {'C-S': (1.55, 2.00), 'C-O': (1.10, 1.55), 'C-N': (1.20, 1.60), 'B-O': (1.20, 1.65), 'B-N': (1.20, 1.65), 'S-S': (1.85, 2.20),
         'O-P': (1.40, 1.75), 'O-S': (1.35, 1.80), 'N-S': (1.45, 1.90), 'N-P': (1.54, 1.77), 'C-Se': (1.75, 2.15)}
COV = {'C': 0.77, 'N': 0.70, 'O': 0.66, 'S': 1.04, 'P': 1.10, 'B': 0.86, 'Se': 1.17}
def canon(a, b): return '-'.join(sorted([a.capitalize(), b.capitalize()]))
def bounds(p):
    a, b = p.split('-'); pl = round(COV[a] + COV[b] + 0.2, 2)
    return {MAXIT[p][0], MAXIT[p][1], pl, round(pl + 0.40, 2), 2.195}
todo = []
for r in csv.DictReader(gzip.open(arch, 'rt', newline='')):
    p = canon(r['prot_el'], r['lig_el'])
    if p in MAXIT and float(r['dist']) in bounds(p):
        todo.append(('archive', r['entry'], (r['prot_chain'], r['prot_res'], r['prot_seq'], r['prot_atom']), (r['lig_chain'], r['lig_comp'], r['lig_seq'], r['lig_atom']), r['dist'], True, None))
for r in csv.DictReader(gzip.open(samp, 'rt', newline='')):
    p = canon(r['lig_el'], r['prot_el'])
    if p in MAXIT and r['hbond_capable_pair'] != 'True' and float(r['dist_A']) in bounds(p):
        ch, resx, atom = r['protein_atom'].split(':')
        comp, lch, lseq = r['ligand'].rsplit('_', 2)
        todo.append(('sample', r['pdbid'], (ch, r['residue'], resx[len(r['residue']):], atom), (lch, comp, re.match(r'-?\d+', lseq).group(0), r['ligand_atom']), r['dist_A'], False,
                     (r['prot_altloc'], r['lig_altloc'])))
out = []
for e in sorted({t[1] for t in todo}):
    path = os.path.join(cifdir, f'{e}.cif.gz')
    if not os.path.exists(path):
        req = urllib.request.Request(f'https://files.rcsb.org/download/{e}.cif.gz', headers={'User-Agent': 'pdb-struct-audit/1.0'})
        with urllib.request.urlopen(req, timeout=300) as u, open(path, 'wb') as f: f.write(u.read())
    st = gemmi.read_structure(path); st.setup_entities(); m = st[0]
    ns = gemmi.NeighborSearch(m, st.cell, 4).populate()
    for src, _, pa, la, dr, sym, alts in [t for t in todo if t[1] == e]:
        best = None
        for ch in m:
            if ch.name != pa[0]: continue
            for res in ch:
                if res.name != pa[1] or str(res.seqid.num) != pa[2]: continue
                for at in res:
                    if at.name != pa[3] or at.occ <= 0: continue
                    for mk in ns.find_atoms(at.pos, '\0', radius=2.7):
                        if not sym and mk.image_idx != 0: continue
                        c = mk.to_cra(m)
                        if (c.chain.name, c.residue.name, str(c.residue.seqid.num), c.atom.name) != la or c.atom.occ <= 0: continue
                        if alts is None:      # archive scan: matching conformers, or one of the two has none
                            if at.altloc != '\0' and c.atom.altloc != '\0' and at.altloc != c.atom.altloc: continue
                        elif (at.altloc.strip('\0'), c.atom.altloc.strip('\0')) != alts: continue   # sample: the conformers named in the row
                        d = st.cell.find_nearest_pbc_image(at.pos, c.atom.pos, mk.image_idx).dist() if sym else at.pos.dist(c.atom.pos)
                        if best is None or d < best: best = d
        out.append({'source': src, 'entry': e, 'protein_chain': pa[0], 'protein_residue': pa[1], 'protein_seq': pa[2], 'protein_atom': pa[3],
                    'ligand_chain': la[0], 'ligand_comp': la[1], 'ligand_seq': la[2], 'ligand_atom': la[3],
                    'protein_altloc': alts[0] if alts else '', 'ligand_altloc': alts[1] if alts else '',
                    'dist_rounded': dr, 'dist_exact': f'{best:.5f}' if best is not None else ''})
with open(dst, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
print(len(out), 'pairs in', len({o['entry'] for o in out}), 'entries;', sum(1 for o in out if not o['dist_exact']), 'not found')
