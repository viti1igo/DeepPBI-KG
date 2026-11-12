"""Are the GO IDs UniProt records for a protein a branch, or independent deepest nodes?

Measured per protein and per GO domain from a Swiss-Prot flat file and go-basic.obo.
Standard library only.

  python scripts/go_uniprot_per_protein.py go-basic.obo uniprot_sprot.txt.gz out.csv
  python scripts/go_uniprot_per_protein.py --selftest

Depth is the longest path from the domain root.
Branch measures leave out the domain root, which every term reaches. Ancestors are plain
graph reachability over the chosen edges, without GO's relation-chain rules.
"""
import collections
import csv
import gzip
import re
import statistics as st
import sys
import tempfile
import os

from go_dag import Dag, RELATION_SETS, parse_obo

EXPERIMENTAL = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"}


def parse_sprot(lines):
    """Yield one dict per entry: accession, kingdom, name, taxid, hosts, go {id: evidence set}."""
    entry = None
    for line in lines:
        tag, rest = line[:2], line[5:].rstrip("\n")
        if tag == "ID":
            entry = {"ac": None, "kingdom": None, "name": "", "taxid": None, "hosts": set(), "go": {}}
        elif tag == "AC" and entry["ac"] is None:
            entry["ac"] = rest.split(";")[0]
        elif tag == "OS":
            entry["name"] += rest
        elif tag == "OC" and entry["kingdom"] is None:
            entry["kingdom"] = rest.split(";")[0].rstrip(".")
        elif tag == "OX":
            entry["taxid"] = re.search(r"NCBI_TaxID=(\d+)", rest).group(1)
        elif tag == "OH":
            entry["hosts"].update(re.findall(r"NCBI_TaxID=(\d+)", rest))
        elif tag == "DR" and rest.startswith("GO; "):
            fields = rest.split("; ")
            entry["go"].setdefault(fields[1], set()).add(fields[3].split(":")[0])
        elif tag == "//":
            yield entry


def assign_groups(entries):
    """Return the group names per entry, and how many phages each rule found."""
    bacterial_taxa = {e["taxid"] for e in entries if e["kingdom"] == "Bacteria"}
    routes = collections.Counter()
    groups = []
    for e in entries:
        names = ["all"]
        if e["kingdom"] == "Bacteria":
            names.append("bacteria")
        elif e["kingdom"] == "Viruses":
            by_name = "phage" in e["name"].lower()
            by_host = bool(e["hosts"] & bacterial_taxa)
            if by_name or by_host:
                names.append("phages")
                routes["name and host" if by_name and by_host else "name only" if by_name else "host only"] += 1
        groups.append(names)
    return groups, routes


def measure(protein_terms, dag):
    """protein_terms holds one set of GO IDs per protein, all in one domain and all live."""
    m = collections.Counter()
    terms_per, depths, tip_fraction, closure_share = [], [], [], []
    for listed in protein_terms:
        ancestors = {t: dag.ancestors(t) for t in listed}
        above = set().union(*ancestors.values())
        with_ancestor = [t for t in listed if ancestors[t] & listed]
        tips = listed - above

        m["proteins"] += 1
        m["annotations"] += len(listed)
        m["single_term_proteins"] += len(listed) == 1
        m["proteins_with_pair"] += bool(with_ancestor)
        m["annotations_with_listed_ancestor"] += len(with_ancestor)
        m["annotations_that_are_listed_ancestors"] += len(listed & above)
        m["annotations_on_leaf"] += sum(map(dag.is_leaf, listed))
        m["proteins_all_leaves"] += all(map(dag.is_leaf, listed))
        m["tips"] += len(tips)
        m["tips_on_leaf"] += sum(map(dag.is_leaf, tips))
        terms_per.append(len(listed))
        depths.extend(map(dag.longest_depth, listed))

        for tip in tips:
            branch = ancestors[tip] - dag.roots
            if branch:
                shown = len(branch & listed) / len(branch)
                tip_fraction.append(shown)
                m["tips_with_branch"] += 1
                m["tips_full_branch"] += shown == 1
                m["tips_empty_branch"] += shown == 0
        closure = (listed | above) - dag.roots
        if closure:
            closure_share.append(len(closure & listed) / len(closure))

    if not m["proteins"]:
        return None
    n, a, t = m["proteins"], m["annotations"], m["tips"]
    multi = n - m["single_term_proteins"]
    tb = m["tips_with_branch"] or 1
    return {
        "proteins": n,
        "annotations": a,
        "median_terms": st.median(terms_per),
        "pct_single_term": 100 * m["single_term_proteins"] / n,
        "pct_with_term_and_ancestor": 100 * m["proteins_with_pair"] / n,
        "pct_with_term_and_ancestor_of_multi": 100 * m["proteins_with_pair"] / multi if multi else 0,
        "pct_only_independent": 100 * (n - m["proteins_with_pair"]) / n,
        "pct_only_independent_multi": 100 * (multi - m["proteins_with_pair"]) / n,
        "pct_annotations_with_listed_ancestor": 100 * m["annotations_with_listed_ancestor"] / a,
        "pct_annotations_that_are_listed_ancestors": 100 * m["annotations_that_are_listed_ancestors"] / a,
        "mean_pct_branch_listed": 100 * st.mean(tip_fraction) if tip_fraction else 0,
        "median_pct_branch_listed": 100 * st.median(tip_fraction) if tip_fraction else 0,
        "pct_tips_full_branch": 100 * m["tips_full_branch"] / tb,
        "pct_tips_empty_branch": 100 * m["tips_empty_branch"] / tb,
        "median_pct_closure_listed": 100 * st.median(closure_share) if closure_share else 0,
        "pct_annotations_on_leaf": 100 * m["annotations_on_leaf"] / a,
        "pct_tips_on_leaf": 100 * m["tips_on_leaf"] / t,
        "pct_proteins_all_leaves": 100 * m["proteins_all_leaves"] / n,
        "median_depth": st.median(depths),
    }


def background(dag, terms, domain):
    live = [t for t in dag.live if terms[t]["domain"] == domain]
    return {"go_live_terms": len(live),
            "go_pct_leaf": 100 * sum(map(dag.is_leaf, live)) / len(live),
            "go_median_depth": st.median(map(dag.longest_depth, live)),
            "go_pct_multi_parent": 100 * sum(len(dag.parents[t]) > 1 for t in live) / len(live)}


def run(obo_path, sprot_lines, out_path):
    terms, alt = parse_obo(obo_path)
    entries = list(parse_sprot(sprot_lines))
    groups, routes = assign_groups(entries)

    stats = collections.Counter(entries=len(entries))
    prepared = []  # (group names, {domain: {evidence subset: set of terms}})
    for e, names in zip(entries, groups):
        by_domain = {}
        stats["entries_with_go"] += bool(e["go"])
        for raw, evidence in e["go"].items():
            stats["raw_annotations"] += 1
            term = alt.get(raw, raw)
            stats["alt_id_annotations"] += term != raw
            if term not in terms:
                stats["unknown_term_annotations"] += 1
            elif terms[term]["obsolete"]:
                stats["obsolete_annotations"] += 1
            else:
                subsets = by_domain.setdefault(terms[term]["domain"], collections.defaultdict(set))
                subsets["all"].add(term)
                if evidence & EXPERIMENTAL:
                    subsets["experimental"].add(term)
                if "IEA" in evidence:
                    subsets["electronic"].add(term)
        prepared.append((names, by_domain))

    rows = []
    for relation_name, relations in RELATION_SETS.items():
        dag = Dag(terms, alt, relations)
        pools = collections.defaultdict(list)
        for names, by_domain in prepared:
            for domain, subsets in by_domain.items():
                for evidence, listed in subsets.items():
                    for name in names:
                        pools[name, evidence, domain].append(listed)
        for (group, evidence, domain), proteins in sorted(pools.items()):
            result = measure(proteins, dag)
            rows.append({"relations": relation_name, "group": group, "evidence": evidence,
                         "domain": domain, **result, **background(dag, terms, domain)})

    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return stats, routes, rows


def selftest():
    # R <- A <- B <- D (is_a), B <-part_of- C, E regulates D, F is_a R; X is obsolete; old -> B.
    obo = """format-version: 1.2

[Term]
id: GO:R
namespace: biological_process

[Term]
id: GO:A
namespace: biological_process
is_a: GO:R ! root

[Term]
id: GO:B
alt_id: GO:old
namespace: biological_process
is_a: GO:A ! a

[Term]
id: GO:C
namespace: biological_process
relationship: part_of GO:B ! b

[Term]
id: GO:D
namespace: biological_process
is_a: GO:B ! b

[Term]
id: GO:E
namespace: biological_process
is_a: GO:R ! root
relationship: regulates GO:D ! d

[Term]
id: GO:F
namespace: biological_process
is_a: GO:R ! root

[Term]
id: GO:X
namespace: biological_process
is_obsolete: true

[Typedef]
id: part_of
is_a: GO:shouldBeIgnored
"""
    flat = """ID   P1 Reviewed; 1 AA.
AC   P1; Q1;
OS   Escherichia coli.
OC   Bacteria; Pseudomonadota;
OC   Enterobacterales.
OX   NCBI_TaxID=562;
DR   GO; GO:D; P:d; IDA:X.
DR   GO; GO:old; P:b; IEA:X.
DR   GO; GO:X; P:gone; IEA:X.
DR   GO; GO:missing; P:gone; IEA:X.
//
ID   P2 Reviewed; 1 AA.
AC   P2;
OS   Enterobacteria virus T.
OC   Viruses; Duplodnaviria.
OX   NCBI_TaxID=9;
OH   NCBI_TaxID=562; Escherichia coli.
DR   GO; GO:D; P:d; IEA:X.
DR   GO; GO:F; P:f; IEA:X.
//
ID   P3 Reviewed; 1 AA.
AC   P3;
OS   Influenza A.
OC   Viruses.
OX   NCBI_TaxID=11;
OH   NCBI_TaxID=9606; Homo sapiens.
DR   GO; GO:C; P:c; IEA:X.
//
"""
    with tempfile.TemporaryDirectory() as folder:
        obo_path, out = os.path.join(folder, "t.obo"), os.path.join(folder, "t.csv")
        open(obo_path, "w").write(obo)
        stats, routes, rows = run(obo_path, flat.splitlines(True), out)

    assert stats["raw_annotations"] == 7 and stats["alt_id_annotations"] == 1
    assert stats["obsolete_annotations"] == 1 and stats["unknown_term_annotations"] == 1
    assert routes == {"host only": 1}

    assert [e["kingdom"] for e in parse_sprot(flat.splitlines(True))] == ["Bacteria", "Viruses", "Viruses"]
    row = {(r["relations"], r["group"], r["evidence"]): r for r in rows}
    p1 = row["is_a", "bacteria", "all"]  # P1 lists D and its ancestor B.
    assert p1["proteins"] == 1 and p1["pct_with_term_and_ancestor"] == 100
    assert p1["pct_annotations_with_listed_ancestor"] == 50
    assert p1["pct_annotations_that_are_listed_ancestors"] == 50
    assert p1["mean_pct_branch_listed"] == 50  # Tip D: branch {B, A}, B listed.
    assert p1["pct_annotations_on_leaf"] == 50 and p1["median_depth"] == 2.5

    p2 = row["is_a", "phages", "all"]  # P2 lists D and F, unrelated: F's branch is only the root.
    assert p2["pct_only_independent_multi"] == 100 and p2["pct_tips_empty_branch"] == 100

    # Under is_a, C has no parent, so it counts as a root; under part_of it hangs below B.
    assert row["is_a+part_of", "all", "all"]["proteins"] == 3
    assert row["is_a", "bacteria", "experimental"]["annotations"] == 1
    assert row["is_a", "bacteria", "electronic"]["annotations"] == 1

    print("selftest passed")


if __name__ == "__main__":
    if sys.argv[1:] == ["--selftest"]:
        selftest()
    elif len(sys.argv) != 4:
        sys.exit(__doc__)
    else:
        obo_file, sprot_file, out_file = sys.argv[1:]
        stats, routes, rows = run(obo_file, gzip.open(sprot_file, "rt", encoding="utf-8"), out_file)
        print(dict(stats))
        print("phage rule:", dict(routes))
        print(len(rows), "rows written to", out_file)
