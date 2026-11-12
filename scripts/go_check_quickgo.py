"""Check go_dag.py's ancestors against QuickGO on a reproducible random sample.

  python scripts/go_check_quickgo.py go-basic.obo quickgo_check.csv
  python scripts/go_check_quickgo.py go-basic.obo out.csv --sample 300 --seed 0

The same OBO file and seed always pick the same terms, so two runs are comparable.
The sample is split evenly across F, P and C. Standard library only.

QuickGO serves its own GO release, so a difference can be a release difference rather
than a bug: the script prints both releases and lists every differing term.
`go-basic.obo` also leaves out edges that cross GO domains, which QuickGO keeps, so
some `part_of` differences are expected by design.
"""
import argparse
import csv
import json
import random
import urllib.error
import urllib.request

from go_dag import Dag, RELATION_SETS, obo_version, parse_obo

API = "https://www.ebi.ac.uk/QuickGO/services/ontology/go"
QUICKGO_RELATIONS = {"is_a": "is_a", "is_a+part_of": "is_a,part_of"}
BATCH = 20


def sample_terms(dag, terms, count, seed):
    """Count terms, split evenly across the three domains, chosen by seed."""
    rng = random.Random(seed)
    chosen = []
    for domain in "FPC":
        pool = sorted(t for t in dag.live - dag.roots if terms[t]["domain"] == domain)
        chosen += rng.sample(pool, min(count // 3, len(pool)))
    return sorted(chosen)


def fetch(path):
    with urllib.request.urlopen(API + path, timeout=60) as response:
        return json.load(response)


def quickgo_ancestors(ids, relations):
    """{term: set of ancestors} from QuickGO, without the term itself."""
    found = {}
    for start in range(0, len(ids), BATCH):
        batch = ids[start:start + BATCH]
        data = fetch(f"/terms/{','.join(batch)}/ancestors?relations={relations}")
        for result in data["results"]:
            found[result["id"]] = set(result.get("ancestors", [])) - {result["id"]}
    missing = set(ids) - set(found)
    if missing:
        raise SystemExit(f"QuickGO returned nothing for {len(missing)} terms, e.g. {sorted(missing)[:3]}")
    return found


def rows(chosen, terms, dag, theirs):
    """One row per term: our ancestors against QuickGO's."""
    out = []
    for term in chosen:
        ours, their = dag.ancestors(term), theirs[term]
        out.append({"go_id": term, "domain": terms[term]["domain"], "name": terms[term]["name"],
                    "ours": len(ours), "quickgo": len(their), "match": "yes" if ours == their else "no",
                    "only_ours": " ".join(sorted(ours - their)),
                    "only_quickgo": " ".join(sorted(their - ours))})
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("obo")
    parser.add_argument("out")
    parser.add_argument("--sample", type=int, default=300)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--relations", choices=list(QUICKGO_RELATIONS), default="is_a+part_of")
    args = parser.parse_args()

    terms, alt = parse_obo(args.obo)
    dag = Dag(terms, alt, RELATION_SETS[args.relations])
    chosen = sample_terms(dag, terms, args.sample, args.seed)
    try:
        theirs = quickgo_ancestors(chosen, QUICKGO_RELATIONS[args.relations])
        about = fetch("/about")
    except urllib.error.URLError as error:
        raise SystemExit(f"QuickGO request failed, so nothing was checked: {error}")

    table = rows(chosen, terms, dag, theirs)
    with open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)

    differ = [r for r in table if r["match"] == "no"]
    print(f"our GO release: {obo_version(args.obo)}")
    print(f"QuickGO release: {about.get('go', {}).get('version', 'unknown')}")
    print(f"{len(table)} terms checked with {args.relations}, seed {args.seed}: "
          f"{len(table) - len(differ)} match, {len(differ)} differ")
    for row in differ:
        print(f"  {row['go_id']} {row['name']} ({row['domain']}): "
              f"only ours [{row['only_ours']}] only QuickGO [{row['only_quickgo']}]")
    print("written to", args.out)


if __name__ == "__main__":
    main()
