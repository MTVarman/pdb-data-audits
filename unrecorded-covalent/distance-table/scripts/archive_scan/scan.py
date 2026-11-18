"""Contact scanner: ligand atoms within 2.6 A of a protein nucleophile, and whether
`_struct_conn` records a covalent bond between that residue pair.

Method follows unrecorded-covalent/cys-ligand-sites/README.md and strict-cut-491/README.md:
  model 1 only, hydrogens removed, occupancy > 0 on both atoms, altlocs must match or be
  absent, symmetry-aware nearest periodic image via gemmi NeighborSearch, ligand = non-polymer
  non-water residue with > 1 heavy atom containing carbon, contacted ligand atom is never a
  metal, protein partner must be a polymer residue.
Recorded = a _struct_conn row with conn_type_id starting 'covale' joining that protein residue
to that ligand residue in either order, under any symmetry operator.
"""
import gemmi

CUT = 2.6

SIDECHAIN_NUC = {"SG", "OG", "OG1", "OG2", "OH", "NZ", "ND1", "NE2", "OD1", "OD2",
                 "OE1", "OE2", "ND2", "NE", "NH1", "NH2", "SD", "SE", "SE1"}
BACKBONE_NUC = {"N", "O", "OXT"}
NUC = SIDECHAIN_NUC | BACKBONE_NUC

METALS = {"LI", "BE", "NA", "MG", "AL", "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO",
          "NI", "CU", "ZN", "GA", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD",
          "AG", "CD", "IN", "SN", "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD",
          "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT",
          "AU", "HG", "TL", "PB", "BI", "PO", "FR", "RA", "AC", "TH", "PA", "U", "NP", "PU",
          "SB", "AS", "GE", "SI", "B"}

WATER = {"HOH", "DOD", "WAT", "H2O", "D2O"}


def _covale_keys(block):
    """Exact residue-pair keys for every `covale*` _struct_conn row.

    Returns (auth_pairs, label_pairs): sets of frozenset({(asym, comp, seq), (asym, comp, seq)}).
    Keys are exact residue identities; no chain- or sequence-wildcard matching, which would
    over-match in entries holding several copies of the same component in one chain.
    """
    auth, label = set(), set()
    cat = block.get_mmcif_category("_struct_conn.")
    if not cat or "conn_type_id" not in cat:
        return auth, label
    n = len(cat["conn_type_id"])

    def col(name):
        c = cat.get(name)
        return c if c else [None] * n

    ct = col("conn_type_id")
    cols = {
        "auth": (col("ptnr1_auth_asym_id"), col("ptnr1_auth_comp_id"), col("ptnr1_auth_seq_id"),
                 col("ptnr2_auth_asym_id"), col("ptnr2_auth_comp_id"), col("ptnr2_auth_seq_id")),
        "label": (col("ptnr1_label_asym_id"), col("ptnr1_label_comp_id"), col("ptnr1_label_seq_id"),
                  col("ptnr2_label_asym_id"), col("ptnr2_label_comp_id"), col("ptnr2_label_seq_id")),
    }
    for i in range(n):
        if not str(ct[i] or "").lower().startswith("covale"):
            continue
        for space, (xa, xc, xs, ya, yc, ys) in cols.items():
            if xa[i] in (None, "?", ".") or ya[i] in (None, "?", "."):
                continue
            def nz(v):
                v = "" if v is None else str(v)
                return "." if v in ("?", "") else v
            p = (nz(xa[i]), nz(xc[i]), nz(xs[i]))
            q = (nz(ya[i]), nz(yc[i]), nz(ys[i]))
            (auth if space == "auth" else label).add(frozenset([p, q]))
    return auth, label


def scan_path(path, entry_id):
    doc = gemmi.cif.read(path)
    block = doc.sole_block()
    cov_auth, cov_label = _covale_keys(block)
    st = gemmi.make_structure_from_block(block)
    st.setup_entities()
    st.remove_hydrogens()
    while len(st) > 1:
        del st[1]
    if len(st) == 0:
        return []
    model = st[0]

    poly_ids = set()
    for ch in model:
        for res in ch.get_polymer():
            poly_ids.add((ch.name, res.seqid.num, res.name))

    lig = []
    for ch in model:
        for res in ch:
            if res.name in WATER:
                continue
            if (ch.name, res.seqid.num, res.name) in poly_ids:
                continue
            hv = [a for a in res if a.element.name not in ("H", "D")]
            if len(hv) <= 1:
                continue
            if not any(a.element.name == "C" for a in hv):
                continue
            lig.append((ch, res))
    if not lig:
        return []

    ns = gemmi.NeighborSearch(model, st.cell, CUT + 0.5).populate()
    seen = {}
    for ch, res in lig:
        for atom in res:
            if atom.occ <= 0 or atom.element.name.upper() in METALS:
                continue
            for m in ns.find_atoms(atom.pos, "\0", radius=CUT):
                cra = m.to_cra(model)
                if cra.atom.name not in NUC or cra.atom.occ <= 0:
                    continue
                if (cra.chain.name, cra.residue.seqid.num, cra.residue.name) not in poly_ids:
                    continue
                if cra.chain.name == ch.name and cra.residue.seqid.num == res.seqid.num:
                    continue
                if atom.altloc and cra.atom.altloc and atom.altloc != cra.atom.altloc:
                    continue
                d = st.cell.find_nearest_pbc_image(cra.atom.pos, atom.pos, m.image_idx).dist()
                if d > CUT or d < 0.5:
                    continue
                key = (cra.chain.name, cra.residue.name, cra.residue.seqid.num, cra.atom.name,
                       cra.atom.element.name, ch.name, res.name, res.seqid.num, atom.name,
                       atom.element.name)
                if key in seen and seen[key][0] <= d:
                    continue
                p = (cra.chain.name, cra.residue.name, str(cra.residue.seqid.num))
                q = (ch.name, res.name, str(res.seqid.num))
                rec = frozenset([p, q]) in cov_auth
                if not rec and cov_label:
                    lp = (cra.residue.subchain or ".", cra.residue.name,
                          str(cra.residue.label_seq) if cra.residue.label_seq is not None else ".")
                    lq = (res.subchain or ".", res.name,
                          str(res.label_seq) if res.label_seq is not None else ".")
                    rec = frozenset([lp, lq]) in cov_label
                seen[key] = (d, rec)

    return [{"entry": entry_id, "prot_chain": k[0], "prot_res": k[1], "prot_seq": k[2],
             "prot_atom": k[3], "prot_el": k[4], "lig_chain": k[5], "lig_comp": k[6],
             "lig_seq": k[7], "lig_atom": k[8], "lig_el": k[9],
             "dist": round(v[0], 3), "recorded": int(v[1])}
            for k, v in seen.items()]
