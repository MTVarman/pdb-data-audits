#!/usr/bin/env python3
"""
Weekly PDB release auditor.

Usage:  release_audit.py <YYYY-MM-DD> [workdir]

Pulls every entry released on that date, then runs three checks that need
no agreement about where a bond ends:

  PARITY   heme c (HEC) attaches through exactly two cysteine thioethers
           (CXXCH). A heme with two cysteines in contact and fewer than two
           covale records is under-recorded. No distance threshold decides
           this -- the count does.

  CONTACT  nucleophilic protein atom within 2.3 A of a ligand heavy atom
           with no covale/disulf record. Banded, not judged.

  GEOMETRY contacts under 1.4 A: closer than any real single bond, i.e.
           atoms modelled on top of each other.

Emits <date>_findings.csv and <date>_brief.md.
"""
import sys, os, json, csv, gzip, urllib.request, collections
import gemmi
import concurrent.futures as cf

NUC = {'CYS':['SG'],'SER':['OG'],'THR':['OG1'],'LYS':['NZ'],'TYR':['OH'],
       'HIS':['NE2','ND1'],'ASP':['OD1','OD2'],'GLU':['OE1','OE2'],
       'MET':['SD'],'ARG':['NH1','NH2','NE'],'ASN':['ND2'],'GLN':['NE2'],'TRP':['NE1']}
METALS = set('LI NA K RB CS MG CA SR BA MN FE CO NI CU ZN CD HG PT AU AG PB V CR MO W '
             'AL GA IN TL SN SB BI RU RH PD OS IR LA CE PR ND SM EU GD TB DY HO ER TM '
             'YB LU U ZR HF Y SC TI'.split())
SKIP = {'HOH','DOD','WAT','UNX','UNL'}
SEARCH = 'https://search.rcsb.org/rcsbsearch/v2/query'


def release_ids(date):
    q = {"query":{"type":"terminal","service":"text","parameters":{
            "attribute":"rcsb_accession_info.initial_release_date",
            "operator":"equals","value":date}},
         "return_type":"entry",
         "request_options":{"paginate":{"start":0,"rows":10000},
                            "results_content_type":["experimental"]}}
    req = urllib.request.Request(SEARCH,
            data=json.dumps(q).encode(), headers={'Content-Type':'application/json'})
    d = json.loads(urllib.request.urlopen(req, timeout=120).read())
    return [r['identifier'] for r in d.get('result_set', [])]


def fetch(args):
    pdbid, cifdir = args
    p = os.path.join(cifdir, pdbid + '.cif.gz')
    if os.path.exists(p) and os.path.getsize(p) > 0:
        return pdbid, 'have'
    for _ in range(3):
        try:
            d = urllib.request.urlopen(
                f'https://files.rcsb.org/download/{pdbid}.cif.gz', timeout=120).read()
            with open(p, 'wb') as fh:
                fh.write(d)
            return pdbid, 'ok'
        except Exception as e:
            err = str(e)
    return pdbid, 'fail ' + err


BOND_MAX = 2.30      # generous upper edge for a C-S thioether
SITE_MAX = 12.0      # a heme belonging to this motif is at least this close
HIS_FE_MAX = 2.8     # axial His NE2-Fe; beyond this the His is not coordinating


def _records(block):
    rec = set()
    for r in block.find('_struct_conn.', ['conn_type_id',
                                          'ptnr1_auth_asym_id','ptnr1_auth_comp_id','ptnr1_auth_seq_id',
                                          'ptnr2_auth_asym_id','ptnr2_auth_comp_id','ptnr2_auth_seq_id']):
        if r[0] not in ('covale', 'disulf'):
            continue
        a = (r[1], r[2], r[3]); c = (r[4], r[5], r[6])
        rec.add((a, c)); rec.add((c, a))
    return rec


def heme_c_sites(block, st):
    rec = _records(block)
    m = st[0]

    hemes = [(ch.name, res) for ch in m for res in ch if res.name == 'HEC']
    out = []

    for ch in m:
        resl = [r for r in ch if r.name not in ('HOH', 'DOD')]
        by_num = {r.seqid.num: r for r in resl}
        for r in resl:
            if r.name != 'CYS':
                continue
            n = r.seqid.num
            c2, h = by_num.get(n + 3), by_num.get(n + 4)
            if not (c2 is not None and c2.name == 'CYS'
                    and h is not None and h.name == 'HIS'):
                continue                                   # not C-x-x-C-H

            sgs = [(cy, next((a for a in cy if a.name == 'SG' and a.occ > 0), None))
                   for cy in (r, c2)]
            anchor = next((sg for _, sg in sgs if sg is not None), None)
            if anchor is None:
                continue

            # the heme this motif belongs to: nearest HEC iron to the motif's
            # histidine NE2 (falls back to a cysteine SG if the His has no NE2)
            ne2 = next((a for a in h if a.name == 'NE2' and a.occ > 0), None)
            probe = ne2 if ne2 is not None else anchor
            best = None
            for hch, hres in hemes:
                fe = next((a for a in hres if a.element.name == 'Fe' and a.occ > 0), None)
                if fe is None:
                    continue
                d = probe.pos.dist(fe.pos)
                if best is None or d < best[0]:
                    best = (d, hch, hres)
            if best is None or best[0] > SITE_MAX:
                out.append((ch.name, n, None, 'NO_HEME', [], None))
                continue
            his_fe = round(best[0], 2) if ne2 is not None else None

            _, hch, hres = best
            carb = [a for a in hres if a.element.name == 'C' and a.occ > 0]
            det = []
            for cy, sg in sgs:
                if sg is None:
                    det.append((f'{ch.name}:CYS{cy.seqid.num}', None, False)); continue
                d = min(sg.pos.dist(a.pos) for a in carb)
                key = ((ch.name, 'CYS', str(cy.seqid.num)),
                       (hch, 'HEC', str(hres.seqid.num)))
                det.append((f'{ch.name}:CYS{cy.seqid.num}', round(d, 3), key in rec))

            bonded = [x for x in det if x[1] is not None and x[1] <= BOND_MAX]
            if len(bonded) == 2 and all(x[2] for x in bonded):
                state = 'OK'
            elif len(bonded) == 2:
                state = 'UNRECORDED'
            elif his_fe is not None and his_fe <= HIS_FE_MAX:
                state = 'UNBONDED_HIS_ON_FE'      # His holds the iron; thioethers not built
            else:
                state = 'UNBONDED_DISPLACED'      # whole motif sits away from the heme
            out.append((ch.name, n, f'{hch}:HEC{hres.seqid.num}', state, det, his_fe))
    return out


def ok_alt(a, b):
    return (not a.has_altloc()) or (not b.has_altloc()) or a.altloc == b.altloc


def scan(path):
    """Returns (entry, meta, contacts, hemes). Imported lazily so the parent
    process never holds gemmi state."""
    import gemmi
    entry = os.path.basename(path)[:4]
    doc = gemmi.cif.read(path)
    b = doc.sole_block()
    st = gemmi.make_structure_from_block(b)
    st.setup_entities()
    st.remove_hydrogens()

    meta = {
        'entry': entry,
        'method': b.find_value('_exptl.method') or '',
        'res': b.find_value('_refine.ls_d_res_high')
               or b.find_value('_em_3d_reconstruction.resolution') or '',
        'title': (b.find_value('_struct.title') or '').strip("'\"")[:160],
    }

    rec = set()
    for r in b.find('_struct_conn.', ['conn_type_id',
                                      'ptnr1_auth_asym_id','ptnr1_auth_comp_id','ptnr1_auth_seq_id',
                                      'ptnr2_auth_asym_id','ptnr2_auth_comp_id','ptnr2_auth_seq_id']):
        if r[0] not in ('covale', 'disulf'):
            continue
        a = (r[1], r[2], r[3]); c = (r[4], r[5], r[6])
        rec.add((a, c)); rec.add((c, a))

    model = st[0]
    ns = gemmi.NeighborSearch(model, st.cell, 3).populate()

    contacts = []

    for ch in model:
        for res in ch:
            if res.name in NUC and res.het_flag != 'H':
                for a in res:
                    if a.name not in NUC[res.name] or a.occ <= 0:
                        continue
                    best = None
                    for m in ns.find_atoms(a.pos, '\0', radius=2.3):
                        c = m.to_cra(model)
                        if c.residue.het_flag != 'H' or c.residue.name in SKIP:
                            continue
                        if (c.atom.element.name.upper() in METALS
                                or c.atom.occ <= 0 or not ok_alt(a, c.atom)):
                            continue
                        d = a.pos.dist(c.atom.pos)
                        if best is None or d < best[0]:
                            best = (d, c)
                    if best:
                        d, c = best
                        key = ((ch.name, res.name, str(res.seqid.num)),
                               (c.chain.name, c.residue.name, str(c.residue.seqid.num)))
                        is_rec = key in rec
                        contacts.append([entry, ch.name, res.name, res.seqid.num, a.name,
                                         c.residue.name, c.chain.name, c.residue.seqid.num,
                                         c.atom.name, c.atom.element.name, round(d, 3),
                                         'recorded' if is_rec else 'unrecorded'])
    return entry, meta, contacts, heme_c_sites(b, st)


def main():
    date = sys.argv[1]
    work = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.getcwd(), 'release_' + date)
    cifdir = os.path.join(work, 'cif')
    os.makedirs(cifdir, exist_ok=True)

    ids = release_ids(date)
    print(f'{date}: {len(ids)} entries released', flush=True)
    if not ids:
        return

    with cf.ThreadPoolExecutor(12) as ex:
        res = list(ex.map(fetch, [(i, cifdir) for i in ids]))
    bad = [r for r in res if r[1].startswith('fail')]
    print(f'downloaded {len(res)-len(bad)}/{len(res)}', flush=True)

    all_contacts, all_hemes, metas, failed = [], {}, {}, []
    for n, pdbid in enumerate(ids, 1):
        p = os.path.join(cifdir, pdbid + '.cif.gz')
        if not os.path.exists(p):
            continue
        try:
            e, meta, contacts, sites = scan(p)
            metas[e] = meta
            all_contacts += contacts
            all_hemes[e] = sites
        except Exception as ex:
            failed.append((pdbid, str(ex)[:120]))
            if len(failed) <= 3:
                print(f'  !! {pdbid}: {ex}', flush=True)
        if n % 50 == 0:
            print(f'  scanned {n}/{len(ids)}', flush=True)

    findings = []

    # --- PARITY: every CXXCH heme c site carries exactly two thioethers ---
    NOTE = {
        'UNRECORDED':          'both thioethers in bonding range; covale record missing',
        'UNBONDED_HIS_ON_FE':  'motif His coordinates the iron; cysteines modelled beyond bonding range',
        'UNBONDED_DISPLACED':  'motif present with heme c nearby; neither His nor Cys in contact',
    }
    for entry, sites in sorted(all_hemes.items()):
        for chain, seq, heme, state, det, his_fe in sites:
            if state in ('OK', 'NO_HEME'):
                continue
            desc = '; '.join(
                f'{c} {("%.3f" % d) if d is not None else "noSG"} A'
                f'{" (recorded)" if r else " (unrecorded)"}' for c, d, r in det)
            hf = f'His-Fe {his_fe} A' if his_fe is not None else 'His-Fe n/a'
            findings.append(['PARITY', entry,
                             f'CXXCH {chain}{seq} -> {heme} [{state}]',
                             f'{NOTE[state]} :: {hf}; {desc}',
                             metas.get(entry, {}).get('title', '')])

    # --- GEOMETRY: closer than any real bond ---
    for c in all_contacts:
        if c[10] < 1.40:
            findings.append(['GEOMETRY', c[0],
                             f'{c[1]}:{c[2]}{c[3]} {c[4]} - {c[5]} {c[8]} {c[10]} A',
                             'closer than any single bond; atoms overlap',
                             metas.get(c[0], {}).get('title', '')])

    # --- CONTACT: unrecorded, in or above the ruler band ---
    for c in all_contacts:
        if c[11] == 'unrecorded' and c[10] >= 1.40:
            band = '<2.0' if c[10] < 2.0 else '2.0-2.3'
            findings.append(['CONTACT', c[0],
                             f'{c[1]}:{c[2]}{c[3]} {c[4]} - {c[5]} {c[6]}{c[7]} {c[8]} {c[10]} A',
                             f'no covale record; {band} A band',
                             metas.get(c[0], {}).get('title', '')])

    out_csv = os.path.join(work, f'{date}_findings.csv')
    with open(out_csv, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['check', 'entry', 'detail', 'note', 'title'])
        w.writerows(findings)

    n = collections.Counter(f[0] for f in findings)
    ent = len({f[1] for f in findings})
    with open(os.path.join(work, f'{date}_brief.md'), 'w') as fh:
        fh.write(f'# PDB release audit — {date}\n\n')
        fh.write(f'{len(ids)} entries released. {len(findings)} findings across {ent} entries.\n\n')
        fh.write(f'- **Parity violations (heme c):** {n["PARITY"]}\n')
        fh.write(f'- **Impossible geometry (<1.40 A):** {n["GEOMETRY"]}\n')
        fh.write(f'- **Unrecorded contacts:** {n["CONTACT"]}\n\n')
        if failed:
            fh.write(f'{len(failed)} entries failed to parse.\n\n')
        for cat in ('PARITY', 'GEOMETRY', 'CONTACT'):
            rows = [f for f in findings if f[0] == cat][:15]
            if not rows:
                continue
            fh.write(f'## {cat}\n\n')
            for r in rows:
                fh.write(f'- `{r[1]}` {r[2]} — {r[3]}\n')
            fh.write('\n')

    if failed:
        print(f'\n!! {len(failed)}/{len(ids)} entries failed to parse'
              f' -- first: {failed[0][0]} {failed[0][1]}')
    print(f'\n{len(findings)} findings / {ent} entries')
    for k, v in n.most_common():
        print(f'  {k}: {v}')
    print('->', out_csv)


if __name__ == '__main__':
    main()
