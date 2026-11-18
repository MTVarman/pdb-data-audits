"""Download + scan one shard of the at-risk population. Resumable.

usage: python3 -I 06_scan_shard.py SHARD NSHARDS
writes data/scan/contacts_<SHARD>.tsv   (one row per contact)
       data/scan/done_<SHARD>.tsv       (entry, n_contacts, status)  = the scanned-entry log
"""
import json, os, queue, sys, threading, time, traceback, urllib.request

sys.path.insert(0, "/home/claude/work")
import scan

SHARD = int(sys.argv[1])
NSHARDS = int(sys.argv[2])
D = "/home/claude/work/data"
OUT = f"{D}/scan"
TMP = f"/home/claude/work/cif_{SHARD}"
os.makedirs(OUT, exist_ok=True)
os.makedirs(TMP, exist_ok=True)

atrisk = [l.strip().split("\t")[0] for l in open(f"{D}/atrisk_ids.txt") if l.strip()]
mine = [e for i, e in enumerate(atrisk) if i % NSHARDS == SHARD]

done_path = f"{OUT}/done_{SHARD}.tsv"
cont_path = f"{OUT}/contacts_{SHARD}.tsv"
already = set()
if os.path.exists(done_path):
    for l in open(done_path):
        already.add(l.split("\t")[0])
todo = [e for e in mine if e not in already]
print(f"shard {SHARD}: {len(mine)} entries, {len(already)} done, {len(todo)} to do", flush=True)

COLS = ["entry", "prot_chain", "prot_res", "prot_seq", "prot_atom", "prot_el",
        "lig_chain", "lig_comp", "lig_seq", "lig_atom", "lig_el", "dist", "recorded"]
if not os.path.exists(cont_path):
    with open(cont_path, "w") as fh:
        fh.write("\t".join(COLS) + "\n")

q = queue.Queue(maxsize=200)
STOP = object()


FIRST_ERR = [None]


def fetch(ids):
    for e in ids:
        p = f"{TMP}/{e}.cif.gz"
        req = urllib.request.Request(f"https://files.rcsb.org/download/{e}.cif.gz",
                                     headers={"User-Agent": "pdb-struct-audit/1.0"})
        ok = False
        for a in range(4):
            try:
                with urllib.request.urlopen(req, timeout=180) as r:
                    data = r.read()
                with open(p, "wb") as fh:
                    fh.write(data)
                ok = True
                break
            except Exception as ex:
                if a == 0 and FIRST_ERR[0] is None:
                    FIRST_ERR[0] = repr(ex)[:200]
                    print("  first download error:", FIRST_ERR[0], flush=True)
                time.sleep(2 + 3 * a)
        q.put((e, p if ok else None))


NTHREAD = 12
chunks = [todo[i::NTHREAD] for i in range(NTHREAD)]
threads = [threading.Thread(target=fetch, args=(c,), daemon=True) for c in chunks]
for t in threads:
    t.start()


def closer():
    for t in threads:
        t.join()
    q.put(STOP)


threading.Thread(target=closer, daemon=True).start()

n = 0
t0 = time.time()
with open(cont_path, "a") as cf, open(done_path, "a") as df:
    while True:
        item = q.get()
        if item is STOP:
            break
        e, p = item
        if p is None:
            df.write(f"{e}\t0\tdownload_failed\n")
            n += 1
            continue
        try:
            if os.path.getsize(p) > 25_000_000:
                df.write(f"{e}\t0\ttoo_large_skipped\n")
                n += 1
                os.remove(p)
                continue
            rows = scan.scan_path(p, e)
            for r in rows:
                cf.write("\t".join(str(r[c]) for c in COLS) + "\n")
            df.write(f"{e}\t{len(rows)}\tok\n")
        except Exception as ex:
            df.write(f"{e}\t0\tscan_error:{type(ex).__name__}\n")
        finally:
            try:
                os.remove(p)
            except OSError:
                pass
        n += 1
        if n % 2000 == 0:
            cf.flush(); df.flush()
            r = n / (time.time() - t0)
            print(f"  shard {SHARD}: {n}/{len(todo)} {r:.1f}/s eta "
                  f"{(len(todo)-n)/r/60:.0f} min", flush=True)
print(f"shard {SHARD} finished {n} entries in {(time.time()-t0)/60:.1f} min", flush=True)
