import json,urllib.request,sys
def search(q, rt="entry", allrows=False):
    body={"query":q,"return_type":rt,"request_options":{"return_all_hits":allrows,"results_verbosity":"compact"} if allrows else {"paginate":{"start":0,"rows":5}}}
    r=urllib.request.Request("https://search.rcsb.org/rcsbsearch/v2/query",data=json.dumps(body).encode(),headers={"Content-Type":"application/json"})
    try:
        with urllib.request.urlopen(r,timeout=120) as f:
            s=f.read()
            return json.loads(s) if s else {"total_count":0}
    except urllib.error.HTTPError as e: return {"err":e.read().decode()[:500]}
def t(attr,op,val): return {"type":"terminal","service":"text","parameters":{"attribute":attr,"operator":op,"value":val}}
def AND(*n): return {"type":"group","logical_operator":"and","nodes":list(n)}
