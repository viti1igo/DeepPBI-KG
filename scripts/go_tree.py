"""Draw the GO tree for example genes: one HTML page per gene, with no extra packages.

  python scripts/go_tree.py go-basic.obo code/富集分析data_process.ipynb out_folder --level 3 \
      --gene gene=dnaA --gene gene=aarA

Each row is one depth from the domain root. Boxes show the terms the gene lists, which are its
deepest, the ancestors worked out from the GO graph, and the terms it rolls up to at --level.
"""
import argparse
import csv
import hashlib
import html
import os
import re
import textwrap

from go_dag import Dag, RELATION_SETS, obo_version, parse_obo
from go_rollup import deepest, read_pairs, resolve_genes

BOX_W, BOX_H, GAP_X, ROW = 200, 62, 16, 104
LEFT, TOP = 86, 12
DOMAIN_NAMES = {"F": "Molecular function", "P": "Biological process", "C": "Cellular component"}
EDGE_DASH = {"is_a": "", "part_of": "6 4"}  # Other relations are dotted.

STYLE = """
.viz { color-scheme: light; --surface: #fcfcfb; --ink: #0b0b0b; --ink2: #52514e; --line: #a3a29c; --band: #f0efec;
  --deepest: #2a78d6; --listed: #1baf7a; --target: #eb6834;
  background: var(--surface); color: var(--ink); font: 14px/1.4 system-ui, sans-serif; padding: 20px; min-height: 100vh; box-sizing: border-box; }
@media (prefers-color-scheme: dark) { .viz { color-scheme: dark; --surface: #1a1a19; --ink: #f0efec; --ink2: #c3c2b7;
  --line: #6d6c66; --band: #262624; --deepest: #3987e5; --listed: #199e70; --target: #d95926; } }
body { margin: 0; } h1 { font-size: 20px; margin: 0 0 4px; } h2 { font-size: 16px; margin: 28px 0 4px; }
p { margin: 0 0 8px; color: var(--ink2); } .scroll { overflow-x: auto; }
.source { font-size: 12px; font-family: ui-monospace, monospace; }
text { fill: var(--ink); font-size: 12px; } .small { fill: var(--ink2); font-size: 11px; } .chosen { font-weight: 700; }
.box { fill: var(--surface); stroke: var(--line); stroke-width: 1; stroke-dasharray: 4 3; }
.deepest { fill: color-mix(in srgb, var(--deepest) 18%, var(--surface)); stroke: var(--deepest); stroke-width: 2; stroke-dasharray: none; }
.listed { fill: color-mix(in srgb, var(--listed) 18%, var(--surface)); stroke: var(--listed); stroke-width: 2; stroke-dasharray: none; }
.ring { fill: none; stroke: var(--target); stroke-width: 3; }
.edge { fill: none; stroke: var(--line); stroke-width: 1.25; } .band { fill: var(--band); }
.legend { display: flex; flex-wrap: wrap; gap: 6px 18px; margin: 8px 0 0; color: var(--ink2); font-size: 13px; }
.key { display: inline-block; width: 18px; height: 12px; border-radius: 3px; margin-right: 6px; vertical-align: -1px; box-sizing: border-box; }
table { border-collapse: collapse; margin-top: 10px; font-size: 13px; } th, td { text-align: left; padding: 3px 12px 3px 0; border-bottom: 1px solid var(--band); }
th { color: var(--ink2); font-weight: 600; } details { margin-top: 8px; }
"""


def provenance(obo_path, pairs_path, relations):
    """Where a drawing's data came from, so a page can be traced back to its inputs."""
    with open(pairs_path, "rb") as handle:
        digest = hashlib.sha256(handle.read()).hexdigest()[:12]
    return (f"Terms from {os.path.basename(pairs_path)} (sha256 {digest}); "
            f"graph from GO {obo_version(obo_path)}, following {relations}.")


def trace_note(row):
    """One line from go_trace_sprot.py's table: which Swiss-Prot entries have this gene's GO IDs."""
    count = int(row["exact_entries"])
    if not count:
        return f"No Swiss-Prot {row['release']} entry has exactly these GO IDs."
    prefix = "e.g. " if count > 3 else ""
    return f"These GO IDs equal {count} Swiss-Prot {row['release']} entries: {prefix}{row['examples']}."


def tree_nodes(listed, dag, level, depth="longest"):
    """Return {term: {listed, deepest, target}} for the listed terms and all their ancestors."""
    tips = deepest(listed, dag)
    targets = set().union(*(dag.at_level(t, level, depth) for t in tips))
    every = listed | set().union(*(dag.ancestors(t) for t in listed))
    return {t: {"listed": t in listed, "deepest": t in tips, "target": t in targets} for t in every}


def layout(nodes, dag, depth="longest"):
    """Return {term: (x, y)} with one row per depth, ordered to keep edges short."""
    measure = dag.longest_depth if depth == "longest" else dag.shortest_depth
    rows = {}
    for term in sorted(nodes):
        rows.setdefault(measure(term), []).append(term)
    widest = max(map(len, rows.values()))
    pitch = BOX_W + GAP_X
    pos = {}

    def place():
        for number, row in rows.items():
            start = LEFT + (widest - len(row)) * pitch / 2
            for i, term in enumerate(row):
                pos[term] = (start + i * pitch, TOP + number * ROW)

    place()
    children = {t: [c for c in nodes if t in dag.parents[c]] for t in nodes}
    # ponytail: fixed barycenter sweeps can leave some crossings; use Graphviz if big trees need better.
    for sweep in range(8):
        linked = dag.parents if sweep % 2 == 0 else children
        for number in sorted(rows, reverse=sweep % 2 == 1):
            def centre(term):
                near = [pos[x][0] for x in linked[term] if x in pos]
                return sum(near) / len(near) if near else pos[term][0]
            rows[number].sort(key=centre)
        place()
    return pos


def edge_relation(terms, child, parent):
    found = [r for r, p in terms[child]["edges"] if p == parent]
    return "is_a" if "is_a" in found else found[0]


def role_words(info, is_leaf):
    words = ["deepest listed" if info["deepest"] else "listed" if info["listed"] else "from GO graph"]
    if is_leaf:
        words.append("GO leaf")
    return words


def render_svg(nodes, terms, dag, level, depth="longest"):
    measure = dag.longest_depth if depth == "longest" else dag.shortest_depth
    pos = layout(nodes, dag, depth)
    width = max(x for x, _ in pos.values()) + BOX_W + 12
    height = max(y for _, y in pos.values()) + BOX_H + 12
    out = [f'<svg width="{width:.0f}" height="{height}" viewBox="0 0 {width:.0f} {height}" role="img">']
    rows = sorted({measure(t) for t in nodes})
    for row in rows:
        y = TOP + row * ROW
        if row == level:
            out.append(f'<rect class="band" x="0" y="{y - 8}" width="{width:.0f}" height="{BOX_H + 16}" rx="4"/>')
        chosen = ' class="chosen"' if row == level else ' class="small"'
        out.append(f'<text x="6" y="{y + BOX_H / 2 + 4}"{chosen}>Level {row}</text>')
    for child in sorted(nodes):
        for parent in sorted(dag.parents[child]):
            x1, y1 = pos[parent][0] + BOX_W / 2, pos[parent][1] + BOX_H
            x2, y2 = pos[child][0] + BOX_W / 2, pos[child][1]
            relation = edge_relation(terms, child, parent)
            dash = EDGE_DASH.get(relation, "2 3")
            mid = (y1 + y2) / 2
            attr = f' stroke-dasharray="{dash}"' if dash else ""
            out.append(f'<path class="edge"{attr} d="M{x1:.0f} {y1} C{x1:.0f} {mid} {x2:.0f} {mid} {x2:.0f} {y2}">'
                       f'<title>{html.escape(terms[child]["name"])} {relation} {html.escape(terms[parent]["name"])}</title></path>')
    for term in sorted(nodes):
        info, (x, y) = nodes[term], pos[term]
        name = terms[term]["name"]
        kind = "deepest" if info["deepest"] else "listed" if info["listed"] else ""
        words = role_words(info, dag.is_leaf(term)) + (["roll-up target"] if info["target"] else [])
        tip = html.escape(f'{term} {name}: {", ".join(words)}; longest depth {dag.longest_depth(term)}, '
                          f'shortest depth {dag.shortest_depth(term)}')
        out.append(f'<g><title>{tip}</title><rect class="box {kind}" x="{x:.0f}" y="{y}" width="{BOX_W}" height="{BOX_H}" rx="6"/>')
        if info["target"]:
            out.append(f'<rect class="ring" x="{x - 3:.0f}" y="{y - 3}" width="{BOX_W + 6}" height="{BOX_H + 6}" rx="8"/>')
        lines = textwrap.wrap(name, 30, max_lines=2, placeholder="…")
        for i, line in enumerate(lines):
            out.append(f'<text x="{x + 8:.0f}" y="{y + 18 + i * 15}">{html.escape(line)}</text>')
        tag = "leaf" if dag.is_leaf(term) else ""
        out.append(f'<text class="small" x="{x + 8:.0f}" y="{y + 53}">{term}{" · " + tag if tag else ""}</text></g>')
    out.append("</svg>")
    return "".join(out)


def render_table(nodes, terms, dag, level, depth="longest"):
    rows = ["<tr><th>GO ID</th><th>Name</th><th>Role</th><th>Longest depth</th><th>Shortest depth</th>"
            f"<th>Level {level} target</th></tr>"]
    for term in sorted(nodes, key=lambda t: (dag.longest_depth(t), t)):
        role = ", ".join(role_words(nodes[term], dag.is_leaf(term)))
        rows.append(f'<tr><td>{term}</td><td>{html.escape(terms[term]["name"])}</td><td>{role}</td>'
                    f'<td>{dag.longest_depth(term)}</td><td>{dag.shortest_depth(term)}</td>'
                    f'<td>{"yes" if nodes[term]["target"] else ""}</td></tr>')
    return f"<details><summary>Table view</summary><table>{''.join(rows)}</table></details>"


LEGEND = """<div class="legend">
<span><span class="key" style="background:color-mix(in srgb,var(--deepest) 18%,var(--surface));border:2px solid var(--deepest)"></span>Deepest listed term</span>
<span><span class="key" style="background:color-mix(in srgb,var(--listed) 18%,var(--surface));border:2px solid var(--listed)"></span>Listed, with a listed descendant</span>
<span><span class="key" style="border:1px dashed var(--line)"></span>Not listed; from the GO graph</span>
<span><span class="key" style="border:3px solid var(--target)"></span>Roll-up target at the chosen level</span>
<span>Lines: solid is_a, dashed part_of, dotted regulates</span>
</div>"""


def render_page(gene, by_domain, terms, dag, level, depth, dropped, source=""):
    parts = [f"<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\"><title>{html.escape(gene)}</title><style>{STYLE}</style></head>"
             f"<body><div class=\"viz\"><h1>GO tree for {html.escape(gene)}</h1>"
             f"<p>Rolled up to level {level}, using {depth}-path depth. Hover a box for details."
             + (f" {dropped} listed IDs are obsolete or unknown, so they aren't drawn." if dropped else "") + "</p>"
             + (f"<p class=\"source\">{html.escape(source)}</p>" if source else "") + LEGEND]
    for domain in "FPC":
        listed = by_domain.get(domain, {}).get(gene)
        if not listed:
            continue
        nodes = tree_nodes(listed, dag, level, depth)
        tips = sum(n["deepest"] for n in nodes.values())
        parts.append(f"<h2>{DOMAIN_NAMES[domain]}</h2><p>{len(listed)} listed, {tips} deepest, {len(nodes) - len(listed)} added from the GO graph.</p>"
                     f'<div class="scroll">{render_svg(nodes, terms, dag, level, depth)}</div>'
                     + render_table(nodes, terms, dag, level, depth))
    parts.append("</div></body></html>")
    return "".join(parts)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("obo")
    parser.add_argument("pairs")
    parser.add_argument("out")
    parser.add_argument("--gene", action="append", required=True)
    parser.add_argument("--level", type=int, default=3)
    parser.add_argument("--trace", help="table from go_trace_sprot.py")
    parser.add_argument("--depth", choices=["longest", "shortest"], default="longest")
    parser.add_argument("--relations", choices=list(RELATION_SETS), default="is_a+part_of")
    args = parser.parse_args()

    terms, alt = parse_obo(args.obo)
    dag = Dag(terms, alt, RELATION_SETS[args.relations])
    genes = read_pairs(args.pairs)
    missing = [g for g in args.gene if g not in genes]
    if missing:
        parser.error(f"no saved GO terms for: {', '.join(missing)}")
    os.makedirs(args.out, exist_ok=True)
    source = provenance(args.obo, args.pairs, args.relations)
    traced = {}
    if args.trace:
        with open(args.trace, encoding="utf-8") as handle:
            traced = {row["label"]: trace_note(row) for row in csv.DictReader(handle)}
    for gene in args.gene:
        by_domain, dropped = resolve_genes({gene: genes[gene]}, terms, dag)
        path = os.path.join(args.out, re.sub(r"\W+", "_", gene) + ".html")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(render_page(gene, by_domain, terms, dag, args.level, args.depth, dropped,
                                     f"{source} {traced.get(gene, '')}".strip()))
        print(path)


if __name__ == "__main__":
    main()
