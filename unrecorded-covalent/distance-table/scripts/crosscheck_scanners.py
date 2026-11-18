"""Cross-check the two independently written scanners on the entries and atom pairs they share.
Usage: python3 -I crosscheck_scanners.py sample_contacts.csv data/archive_atom_pairs.csv.gz   (first argument uncompressed)"""
import csv, gzip, sys, collections, re
l3p, l2p = sys.argv[1], sys.argv[2]
l3 = {}
ents3 = set()
for r in csv.DictReader(open(l3p)):
    ch, resx, atom = r['protein_atom'].split(':')
    m = re.match(r'([A-Z0-9]{1,3}?)(-?\d+)([A-Za-z]?)$', resx)
    # residue name is known from the 'residue' column
    resn = r['residue']; seq = resx[len(resn):]
    comp, lch, lseq = r['ligand'].rsplit('_', 2)
    key = (r['pdbid'], ch, resn, seq, atom, lch, comp, lseq, r['ligand_atom'])
    l3[key] = (float(r['dist_A']), r['recorded_covale'] == 'True', 'Covale' in (r['residue_level_conn'] or ''), r['occupancy_class'])
    ents3.add(r['pdbid'])
print('sample scan rows', len(l3), 'entries', len(ents3))
n2 = 0; both = 0; dd = collections.Counter(); rec = collections.Counter(); bad = []
seen2 = set()
for r in csv.DictReader(gzip.open(l2p, 'rt')):
    if r['entry'] not in ents3: continue
    n2 += 1
    key = (r['entry'], r['prot_chain'], r['prot_res'], r['prot_seq'], r['prot_atom'], r['lig_chain'], r['lig_comp'], r['lig_seq'], r['lig_atom'])
    seen2.add(key)
    if key not in l3: continue
    both += 1
    d3, rec_atom, rec_res, occ = l3[key]
    diff = abs(d3 - float(r['dist']))
    dd['<=0.001' if diff <= 0.0011 else ('<=0.01' if diff <= 0.01 else '>0.01')] += 1
    r2 = r['recorded'] == '1'
    r3 = rec_atom or rec_res
    rec[(r2, r3)] += 1
    if r2 != r3 and len(bad) < 15: bad.append((key, d3, r['dist'], r2, rec_atom, rec_res))
print('archive scan rows in shared entries', n2, '| atom pairs present in both', both)
print('distance agreement', dict(dd))
print('recorded (archive scan residue-level, sample scan atom-or-residue-level):', {str(k): v for k, v in rec.items()})
for b in bad: print('  disagree', b)
# sample-scan rows missing from the archive scan (same entries), by reason
miss = collections.Counter()
ents2 = {k[0] for k in seen2}
for k, v in l3.items():
    if k in seen2: continue
    if k[0] not in ents2: miss['entry not in archive table']+=1; continue
    miss['pair absent']+=1
print('sample scan rows not matched in the archive scan:', dict(miss))
