"""Shared RCSB helpers (stdlib only)."""
import json, time, urllib.request, urllib.error

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"


def _post(url, payload, timeout=240, tries=5):
    data = json.dumps(payload).encode()
    last = None
    for a in range(tries):
        try:
            req = urllib.request.Request(url, data=data,
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            body = e.read()[:400].decode("utf8", "replace")
            last = f"HTTP {e.code} {body}"
            if e.code in (400, 404):
                raise RuntimeError(last)
        except Exception as e:  # noqa
            last = repr(e)
        time.sleep(2 + 3 * a)
    raise RuntimeError(f"failed after {tries}: {last}")


def search(query, return_type="entry", rows=0, start=0, extra=None):
    if rows == 0:
        opts = {"return_counts": True}
    else:
        opts = {"paginate": {"start": start, "rows": rows},
                "results_verbosity": "compact"}
    if extra:
        opts.update(extra)
    q = {"query": query, "return_type": return_type, "request_options": opts}
    return _post(SEARCH, q)


def count(query, return_type="entry"):
    d = search(query, return_type=return_type, rows=0)
    return d.get("total_count", 0)


def ids(query, return_type="entry", page=10000):
    out, start = [], 0
    while True:
        d = search(query, return_type=return_type, rows=page, start=start)
        got = d.get("result_set") or []
        got = [g["identifier"] if isinstance(g, dict) else g for g in got]
        out.extend(got)
        total = d.get("total_count", len(out))
        start += page
        if not got or len(out) >= total:
            return out, total


def term(attribute, value, operator="exact_match", service="text"):
    return {"type": "terminal", "service": service,
            "parameters": {"attribute": attribute, "operator": operator, "value": value}}


def rng(attribute, lo, hi, service="text"):
    return {"type": "terminal", "service": service,
            "parameters": {"attribute": attribute, "operator": "range",
                           "value": {"from": lo, "to": hi,
                                     "include_lower": True, "include_upper": False}}}


def AND(*nodes):
    return {"type": "group", "logical_operator": "and", "nodes": list(nodes)}


def OR(*nodes):
    return {"type": "group", "logical_operator": "or", "nodes": list(nodes)}


def gql(query):
    d = _post(GRAPHQL, {"query": query})
    if "errors" in d and not d.get("data"):
        raise RuntimeError(str(d["errors"])[:400])
    return d["data"]
