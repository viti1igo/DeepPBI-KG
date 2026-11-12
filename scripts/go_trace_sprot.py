"""Trace each saved gene-to-GO label back to the Swiss-Prot entries with exactly its GO IDs.

The DIAMOND output behind the notebook isn't published, so this is evidence, not the original hit:
a label with many matching entries doesn't say which one was hit. Standard library only.

  python scripts/go_trace_sprot.py code/富集分析data_process.ipynb uniprot_sprot.dat.gz 2022_01 trace.csv

Also compares the few DIAMOND hits the notebook still shows in its cell 0 output.
"""
import collections
import csv
import gzip
import json
import re
import sys

from go_rollup import read_pairs
from go_uniprot_per_protein import parse_sprot

HIT_CELL = 0


def visible_hits(notebook):
    """{label: accession} for the DIAMOND rows printed in the notebook."""
    with open(notebook, encoding="utf-8") as handle:
        outputs = json.load(handle)["cells"][HIT_CELL].get("outputs", [])
    text = "".join("".join(o.get("data", {}).get("text/plain", [])) for o in outputs)
    return dict(re.findall(r"^\d+\s+(\S+=\S+)\s+([A-Z0-9]+)\.\d+\s", text, flags=re.M))


def trace(pairs, entries, release):
    """One row per label: the entries whose GO IDs equal the label's."""
    by_terms = collections.defaultdict(list)
    for e in entries:
        by_terms[frozenset(e["go"])].append(e)
    rows = []
    for label, listed in sorted(pairs.items()):
        found = sorted(by_terms.get(frozenset(listed), []), key=lambda e: e["ac"])
        rows.append({"label": label, "terms": len(listed), "release": release, "exact_entries": len(found),
                     "examples": "; ".join(f'{e["ac"]} ({e["name"].rstrip(".")})' for e in found[:3])})
    return rows


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    notebook, sprot, release, out = sys.argv[1:]
    pairs = read_pairs(notebook)
    with gzip.open(sprot, "rt", encoding="utf-8", errors="replace") as handle:
        entries = list(parse_sprot(handle))
    rows = trace(pairs, entries, release)
    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    by_accession = {e["ac"]: set(e["go"]) for e in entries}
    hits = {k: v for k, v in visible_hits(notebook).items() if k in pairs}
    same = sum(by_accession.get(ac) == pairs[label] for label, ac in hits.items())
    counts = [r["exact_entries"] for r in rows]
    print(f"{release}: {len(entries)} entries, {len(rows)} labels")
    print(f"labels equal to at least one entry: {sum(c > 0 for c in counts)}, to exactly one: {counts.count(1)}, "
          f"to none: {counts.count(0)}")
    print(f"DIAMOND hits shown in the notebook whose entry has the same GO IDs: {same}/{len(hits)}")
    print("written to", out)


if __name__ == "__main__":
    main()
