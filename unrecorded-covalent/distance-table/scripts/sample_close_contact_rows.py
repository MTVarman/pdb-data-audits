"""The close-contact tables of the sample entries, as the scan read them from each file.
Writes one row per row of _pdbx_validate_close_contact, for every sample entry that has at least one contact.
Usage: python3 -I sample_close_contact_rows.py scan_dir sample_contacts.csv.gz out_rows.csv.gz
scan_dir holds the scan output (shard_*.jsonl): one line per entry, with the file's close-contact rows under 'close'."""
import csv, glob, gzip, json, sys
scan_dir, contacts, out_rows = sys.argv[1:4]
need = {r['pdbid'] for r in csv.DictReader(gzip.open(contacts, 'rt', newline=''))}
close = {}
for p in sorted(glob.glob(f'{scan_dir}/shard_*.jsonl')):
    for line in open(p):
        r = json.loads(line)
        if r.get('pdbid') in need: close[r['pdbid']] = r.get('close') or []
print('entries with a contact', len(need), 'found in the scan output', len(close))
with gzip.open(out_rows, 'wt', newline='') as f:
    w = csv.writer(f)
    w.writerow(['pdbid', 'model', 'chain_1', 'comp_1', 'seq_1', 'atom_1', 'chain_2', 'comp_2', 'seq_2', 'atom_2', 'dist'])
    n = 0
    for e in sorted(close):
        for c in close[e]:
            w.writerow([e] + list(c)); n += 1
print('close-contact rows written', n, 'in', sum(1 for e in close if close[e]), 'entries')
