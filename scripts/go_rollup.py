"""Move each gene's GO terms up the GO tree to a GO slim, and to a chosen level.

Input is the gene-to-GO pairs saved in a notebook's cell 3, or a table of gene<TAB>GO IDs
(separated by "," or ";"). Standard library only. Run from the repo root:

  python scripts/go_rollup.py go-basic.obo code/富集分析data_process.ipynb out_folder \
      --genes host_key=model/host_key_gene_0.00004.csv --genes phage_key=model/phage_key_gene_0.0001.csv

Writes to out_folder:
  gene_terms.tsv     every listed term, with its depths and whether it is a deepest term
  rollup_SLIM.tsv    each deepest term moved up to the slim's most specific terms, or "other"
  slim_summary.csv   for every gene group: slim coverage, distinct terms, terms per gene
  rollup_level_N.tsv each deepest term moved up to level N
  levels.csv         for every level: distinct terms, terms per gene, and terms already above it
  measurements.csv   branch or independent-node measures, for every edge set and gene group

Levels count edges from the domain root (0). A term already above the level stays as it is
and is marked. A term can reach several terms at one level, or several slim terms.

The default is goslim_metagenomics with level 3 as a secondary view; the reasons are in
docs/GO_LITERATURE_REVIEW.md.
"""
import argparse
import collections
import csv
import json
import os
import re
import statistics as st

from go_dag import Dag, RELATION_SETS, parse_obo
from go_uniprot_per_protein import background, measure

CELL = 3  # The notebook cell that printed the gene-to-GO pairs.


def read_pairs(path):
    """Return {gene: set of GO IDs as written in the file}."""
    genes = collections.defaultdict(set)
    if path.endswith(".ipynb"):
        with open(path, encoding="utf-8") as handle:
            outputs = json.load(handle)["cells"][CELL].get("outputs", [])
        text = "".join("".join(o.get("text", [])) for o in outputs)
        for gene, go in re.findall(r"^(\S.*?)\t(GO:\d{7})\s*$", text, flags=re.M):
            genes[gene].add(go)
    else:
        with open(path, encoding="utf-8") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if len(row) > 1:
                    found = set(re.findall(r"GO:\d{7}|GO:[A-Z]\w*", row[1]))
                    if found:
                        genes[row[0]] |= found
    return dict(genes)


def read_gene_list(path):
    with open(path, encoding="utf-8") as handle:
        return {row["gene"] for row in csv.DictReader(handle)}


def resolve_genes(genes, terms, dag):
    """Return ({domain: {gene: live terms}}, count of obsolete or unknown IDs dropped)."""
    by_domain = collections.defaultdict(dict)
    dropped = 0
    for gene, listed in genes.items():
        for raw in listed:
            term = dag.resolve(raw)
            if term is None:
                dropped += 1
            else:
                by_domain[terms[term]["domain"]].setdefault(gene, set()).add(term)
    return by_domain, dropped


def deepest(listed, dag):
    """The listed terms that are not an ancestor of another listed term."""
    return listed - set().union(*(dag.ancestors(t) for t in listed))


def level_rows(genes, dag, depth, levels):
    """Summarize {gene: terms} at each level."""
    measure_depth = dag.longest_depth if depth == "longest" else dag.shortest_depth
    tips = {gene: deepest(listed, dag) for gene, listed in genes.items()}
    rows = []
    for level in levels:
        rolled = [set().union(*(dag.at_level(t, level, depth) for t in ts)) for ts in tips.values()]
        all_tips = [t for ts in tips.values() for t in ts]
        rows.append({
            "level": level,
            "genes": len(tips),
            "distinct_terms": len(set().union(*rolled)),
            "median_terms_per_gene": st.median(map(len, rolled)),
            "pct_deepest_terms_above_level": 100 * sum(measure_depth(t) < level for t in all_tips) / len(all_tips),
        })
    return rows


def slim_row(genes, dag, slim):
    """Summarize {gene: terms} on a GO slim."""
    tips = [deepest(listed, dag) for listed in genes.values()]
    all_tips = [t for ts in tips for t in ts]
    mapped = [set().union(*(dag.to_slim(t, slim) for t in ts)) for ts in tips]
    return {
        "genes": len(tips),
        "deepest_terms": len(all_tips),
        "pct_deepest_terms_covered": 100 * sum(bool(dag.to_slim(t, slim)) for t in all_tips) / len(all_tips),
        "distinct_slim_terms": len(set().union(*mapped)),
        "median_slim_terms_per_gene": st.median(map(len, mapped)),
        "genes_only_other": sum(not m for m in mapped),
    }


def write_table(path, rows, delimiter):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("obo")
    parser.add_argument("pairs")
    parser.add_argument("out")
    parser.add_argument("--slim", default="goslim_metagenomics", help="a subset named in the OBO file")
    parser.add_argument("--level", type=int, default=3)
    parser.add_argument("--depth", choices=["longest", "shortest"], default="longest")
    parser.add_argument("--relations", choices=list(RELATION_SETS), default="is_a+part_of")
    parser.add_argument("--genes", action="append", default=[], metavar="NAME=CSV",
                        help="a gene group, from a CSV with a 'gene' column")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    terms, alt = parse_obo(args.obo)
    genes = read_pairs(args.pairs)
    groups = {"all": set(genes)}
    for spec in args.genes:
        name, path = spec.split("=", 1)
        members = read_gene_list(path)
        groups[name] = members & set(genes)
        print(f"{name}: {len(groups[name])} of {len(members)} genes have saved GO terms")

    dags = {name: Dag(terms, alt, rel) for name, rel in RELATION_SETS.items()}
    dag = dags[args.relations]
    try:
        slim = dag.slim_terms(terms, args.slim)
    except ValueError as error:
        parser.error(str(error))
    by_domain, dropped = resolve_genes(genes, terms, dag)
    total = sum(map(len, genes.values()))
    print(f"{len(genes)} genes, {total} gene-term pairs, {dropped} dropped as obsolete or unknown")

    listed_rows, rollup_rows, slim_rows = [], [], []
    for domain, members in sorted(by_domain.items()):
        for gene, listed in sorted(members.items()):
            tips = deepest(listed, dag)
            for term in sorted(listed):
                listed_rows.append({"gene": gene, "domain": domain, "go_id": term, "name": terms[term]["name"],
                                    "deepest": "yes" if term in tips else "no",
                                    "longest_depth": dag.longest_depth(term),
                                    "shortest_depth": dag.shortest_depth(term),
                                    "go_leaf": "yes" if dag.is_leaf(term) else "no"})
            measure_depth = dag.longest_depth if args.depth == "longest" else dag.shortest_depth
            for term in sorted(tips):
                for target in sorted(dag.at_level(term, args.level, args.depth)):
                    rollup_rows.append({"gene": gene, "domain": domain, "deepest_go_id": term,
                                        "deepest_name": terms[term]["name"], "deepest_depth": measure_depth(term),
                                        "level": args.level, "level_go_id": target, "level_name": terms[target]["name"],
                                        "above_level": "yes" if measure_depth(term) < args.level else "no"})
                for target in sorted(dag.to_slim(term, slim)) or ["other"]:
                    slim_rows.append({"gene": gene, "domain": domain, "deepest_go_id": term,
                                      "deepest_name": terms[term]["name"], "slim_go_id": target,
                                      "slim_name": terms[target]["name"] if target in terms else "other"})
    write_table(os.path.join(args.out, "gene_terms.tsv"), listed_rows, "\t")
    write_table(os.path.join(args.out, f"rollup_level_{args.level}.tsv"), rollup_rows, "\t")

    write_table(os.path.join(args.out, f"rollup_{args.slim}.tsv"), slim_rows, "\t")

    slim_table, level_table = [], []
    for group, members in groups.items():
        for domain, domain_genes in sorted(by_domain.items()):
            chosen = {g: t for g, t in domain_genes.items() if g in members}
            if chosen:
                slim_table.append({"group": group, "domain": domain, "slim": args.slim, "relations": args.relations,
                                   **slim_row(chosen, dag, slim)})
                max_level = max(dag.longest_depth(t) for ts in chosen.values() for t in ts)
                for row in level_rows(chosen, dag, args.depth, range(max_level + 1)):
                    level_table.append({"group": group, "domain": domain, "depth": args.depth,
                                        "relations": args.relations, **row})
    write_table(os.path.join(args.out, "slim_summary.csv"), slim_table, ",")
    write_table(os.path.join(args.out, "levels.csv"), level_table, ",")

    measures = []
    for relation, relation_dag in dags.items():
        relation_domains, _ = resolve_genes(genes, terms, relation_dag)
        for group, members in groups.items():
            for domain, domain_genes in sorted(relation_domains.items()):
                chosen = [t for g, t in domain_genes.items() if g in members]
                if chosen:
                    measures.append({"relations": relation, "group": group, "domain": domain,
                                     **measure(chosen, relation_dag), **background(relation_dag, terms, domain)})
    write_table(os.path.join(args.out, "measurements.csv"), measures, ",")
    print(f"wrote {len(listed_rows)} listed terms, {len(rollup_rows)} level rows and {len(slim_rows)} slim rows "
          f"to {args.out}")


if __name__ == "__main__":
    main()
