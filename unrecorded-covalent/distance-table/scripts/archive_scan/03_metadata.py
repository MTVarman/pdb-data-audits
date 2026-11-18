"""Per-entry metadata for every experimental PDB entry, from the RCSB Data API (GraphQL).

Writes shards to data/meta/NNNN.jsonl (resumable), one JSON object per entry:
  entry, site, deposit_site, dep_date, rel_date, methods, resolution,
  n_np (nonpolymer entity count), n_prot (protein polymer entity count),
  cit (doi:/pmid:/authors:hash), la (hashed last audit author), comps [[id, formula], ...]
"""
import hashlib, json, os, sys, threading, time
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rcsb

OUT = "/home/claude/work/data/meta"
os.makedirs(OUT, exist_ok=True)
IDS = [l.strip() for l in open("/home/claude/work/data/experimental_ids.txt") if l.strip()]
BATCH = 100
Q = """{entries(entry_ids:%s){
 rcsb_id
 pdbx_database_status{process_site deposit_site recvd_initial_deposition_date}
 rcsb_accession_info{initial_release_date}
 exptl{method}
 rcsb_entry_info{resolution_combined nonpolymer_entity_count polymer_entity_count_protein}
 rcsb_primary_citation{pdbx_database_id_DOI pdbx_database_id_PubMed}
 audit_author{name}
 nonpolymer_entities{nonpolymer_comp{chem_comp{id formula}}}
}}"""

lock = threading.Lock()
done = 0


def h(s):
    return hashlib.sha1(s.encode()).hexdigest()[:12]


def row(e):
    st = e.get("pdbx_database_status") or {}
    ai = e.get("rcsb_accession_info") or {}
    ei = e.get("rcsb_entry_info") or {}
    c = e.get("rcsb_primary_citation") or {}
    au = [a["name"] for a in (e.get("audit_author") or [])]
    if c.get("pdbx_database_id_DOI"):
        cit = "doi:" + c["pdbx_database_id_DOI"].lower()
    elif c.get("pdbx_database_id_PubMed"):
        cit = "pmid:%s" % c["pdbx_database_id_PubMed"]
    else:
        cit = "authors:" + h("|".join(au))
    res = ei.get("resolution_combined") or []
    comps = [[x["nonpolymer_comp"]["chem_comp"]["id"],
              x["nonpolymer_comp"]["chem_comp"].get("formula") or ""]
             for x in (e.get("nonpolymer_entities") or [])
             if x.get("nonpolymer_comp") and x["nonpolymer_comp"].get("chem_comp")]
    return {"entry": e["rcsb_id"], "site": st.get("process_site"),
            "deposit_site": st.get("deposit_site"),
            "dep_date": (st.get("recvd_initial_deposition_date") or "")[:10],
            "rel_date": (ai.get("initial_release_date") or "")[:10],
            "methods": sorted({m["method"] for m in (e.get("exptl") or []) if m.get("method")}),
            "resolution": res[0] if res else None,
            "n_np": ei.get("nonpolymer_entity_count"),
            "n_prot": ei.get("polymer_entity_count_protein"),
            "cit": cit, "la": "la:" + h(au[-1] if au else "?"),
            "comps": comps}


def work(i):
    global done
    path = f"{OUT}/{i:05d}.jsonl"
    if os.path.exists(path) and os.path.getsize(path) > 0:
        with lock:
            done += 1
        return
    chunk = IDS[i * BATCH:(i + 1) * BATCH]
    d = rcsb.gql(Q % json.dumps(chunk))
    ents = [e for e in (d.get("entries") or []) if e]
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        for e in ents:
            fh.write(json.dumps(row(e)) + "\n")
    os.replace(tmp, path)
    with lock:
        done += 1
        if done % 100 == 0:
            print(f"  {done} batches", flush=True)


nb = (len(IDS) + BATCH - 1) // BATCH
print("batches", nb)
t = time.time()
with ThreadPoolExecutor(max_workers=10) as ex:
    list(ex.map(work, range(nb)))
print("elapsed %.1f min" % ((time.time() - t) / 60))
