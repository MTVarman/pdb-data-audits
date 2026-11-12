import gemmi, sys, gzip
st = gemmi.read_structure(sys.argv[1]); st.setup_entities()
m = st[0]
# find the Lys NZ - ATP O1B pair
best = None
ns = gemmi.NeighborSearch(m, st.cell, 4).populate()
for ch in m:
    for r in ch:
        if r.name == 'ATP':
            a = r.find_atom('O1B', '*')
            if not a: continue
            for mk in ns.find_atoms(a.pos, '\0', radius=1.6):
                cra = mk.to_cra(m)
                if cra.residue.name == 'LYS' and cra.atom.name == 'NZ':
                    d = cra.atom.pos.dist(a.pos)
                    if best is None or d < best[0]: best = (d, cra.chain.name, cra.residue.seqid.num, ch.name, r.seqid.num)
print('pair', best)
d, lc, ln, ac, an = best
import numpy as np
def write(shift, name):
    out = gemmi.Structure(); out.cell = st.cell; out.spacegroup_hm = st.spacegroup_hm
    mm = gemmi.Model('1')
    c1 = gemmi.Chain('A'); 
    for r in m[lc]:
        if abs(r.seqid.num - ln) <= 0 and r.entity_type == gemmi.EntityType.Polymer: c1.add_residue(r)
    c2 = gemmi.Chain('B'); atp = [r for r in m[ac] if r.name == 'ATP' and r.seqid.num == an][0].clone() if hasattr(gemmi.Residue,'clone') else None
    if atp is None:
        import copy; atp = copy.deepcopy([r for r in m[ac] if r.name == 'ATP' and r.seqid.num == an][0])
    nz = [r for r in m[lc] if r.seqid.num == ln][0].find_atom('NZ', '*').pos
    o = atp.find_atom('O1B', '*').pos
    v = gemmi.Position(o.x - nz.x, o.y - nz.y, o.z - nz.z); L = v.length()
    for a in atp:
        a.pos = gemmi.Position(a.pos.x + v.x / L * shift, a.pos.y + v.y / L * shift, a.pos.z + v.z / L * shift)
    c2.add_residue(atp); mm.add_chain(c1); mm.add_chain(c2); out.add_model(mm); out.setup_entities()
    out.write_pdb(name)
    return L + shift
for tgt in (None, 2.0, 2.1, 2.15, 2.2, 2.25, 2.29, 2.4, 2.8, 3.1):
    sh = 0 if tgt is None else tgt - d
    nm = 'frag_dep.pdb' if tgt is None else f'frag_{round(tgt*100)}.pdb'
    print(nm, round(write(sh, nm), 3))
