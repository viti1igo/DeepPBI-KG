"""Where do saved gene-to-GO annotations sit in the GO DAG? Run from the repo root.

Reads the gene/GO pairs printed by cell 3 of code/富集分析data_process.ipynb and
walks the is_a and part_of edges in code/GO/go-basic.obo.
"""
import collections
import json
import re
import statistics as st
from functools import lru_cache

terms, alt, cur = {}, {}, None
for line in open("code/GO/go-basic.obo"):
    line = line.rstrip("\n")
    if line == "[Term]":
        cur = {"parents": set(), "ns": None, "obsolete": False}
        continue
    if line.startswith("["):
        cur = None
        continue
    if cur is None or ": " not in line:
        continue
    key, value = line.split(": ", 1)
    if key == "id":
        terms[value] = cur
        cur["id"] = value
    elif key == "namespace":
        cur["ns"] = value
    elif key == "is_obsolete":
        cur["obsolete"] = True
    elif key == "alt_id":
        alt[value] = cur["id"]
    elif key == "is_a":
        cur["parents"].add(value.split(" ! ")[0])
    elif key == "relationship" and value.split()[0] == "part_of":
        cur["parents"].add(value.split()[1])

children = collections.defaultdict(set)
for term_id, term in terms.items():
    for parent in term["parents"]:
        children[parent].add(term_id)


@lru_cache(None)
def ancestors(term_id):
    found = set()
    for parent in terms[term_id]["parents"]:
        if parent in terms:
            found |= {parent} | ancestors(parent)
    return frozenset(found)


@lru_cache(None)
def depth(term_id):
    parents = [p for p in terms[term_id]["parents"] if p in terms]
    return 1 + max(map(depth, parents)) if parents else 0


nb = json.load(open("code/富集分析data_process.ipynb"))
text = "".join("".join(o.get("text", [])) for o in nb["cells"][3].get("outputs", []))
pairs = re.findall(r"^(\S.*?)\t(GO:\d{7})\s*$", text, flags=re.M)
genes = collections.defaultdict(set)
for gene, go in pairs:
    go = alt.get(go, go)
    if go in terms:
        genes[gene].add(go)

annotations = [go for s in genes.values() for go in s]
live = [t for t, v in terms.items() if not v["obsolete"]]
with_ancestor = [sum(1 for go in s if s & ancestors(go)) for s in genes.values()]
print("pairs", len(pairs), "alt-id pairs", sum(g in alt for _, g in pairs))
print("genes", len(genes), "annotations", len(annotations))
print("terms per gene, median", st.median(map(len, genes.values())))
print("annotations on a leaf", round(sum(not children[g] for g in annotations) / len(annotations), 3))
print("leaf share of all live terms", round(sum(not children[t] for t in live) / len(live), 3))
print("annotation depth median", st.median(map(depth, annotations)), "max", max(map(depth, annotations)))
print("live term depth median", st.median(map(depth, live)))
print("annotations with an ancestor on the same gene", sum(with_ancestor))
print("genes with such a term", sum(n > 0 for n in with_ancestor))
print("live terms with more than one parent", round(sum(len(terms[t]["parents"]) > 1 for t in live) / len(live), 3))
