"""Validate go_extract.py on DeepPBI-KG's annotated example genomes. Specific to this repo.

  python scripts/go_extract_validate.py inputs go_data/extract_test
  python scripts/go_extract_validate.py compare go_data/extract_test code/富集分析data_process.ipynb \
      go_data/extract_test/out_2022_01 go_data/extract_test/out_2026_03

`inputs` gathers every Prokka protein from example_host_annotation and example_phage_annotation
into examples.faa, and writes locus_gene.tsv (locus tag, gene name, genome, host or phage).
`compare` prints Markdown tables for docs/GO_EXTRACT_VALIDATION.md.

The authors' query sequences aren't published, so the comparison is by Prokka gene name: a name
counts as reproduced if any of our proteins with that name has exactly the GO IDs of any of the
authors' labels with that name. That is lenient, and the example genomes may not be in their table.
"""
import collections
import glob
import json
import os
import re
import statistics as st
import sys

from go_extract import best_hits, read_hits
from go_rollup import read_pairs

GROUPS = {"host": "example_host_annotation", "phage": "example_phage_annotation"}
THRESHOLDS = [(0, 0), (30, 50), (40, 50), (50, 70)]


def inputs(out):
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "examples.faa"), "w") as faa, open(os.path.join(out, "locus_gene.tsv"), "w") as tsv:
        for group, folder in GROUPS.items():
            for genome in sorted(glob.glob(os.path.join(folder, "prokka_*"))):
                if os.path.basename(genome).startswith("._"):
                    continue
                for path in glob.glob(os.path.join(genome, "*.faa")):
                    if not os.path.basename(path).startswith("._"):
                        with open(path) as handle:
                            faa.write(handle.read())
                for path in glob.glob(os.path.join(genome, "*.tsv")):
                    if os.path.basename(path).startswith("._"):
                        continue
                    with open(path) as handle:
                        next(handle)
                        for line in handle:
                            fields = line.rstrip("\n").split("\t")
                            tsv.write(f"{fields[0]}\t{fields[2]}\t{os.path.basename(genome)}\t{group}\n")


def base(name):
    return re.sub(r"_\d+$", "", name)


def compare(folder, notebook, runs):
    group, gene = {}, {}
    with open(os.path.join(folder, "locus_gene.tsv")) as handle:
        for line in handle:
            locus, name, _, kind = line.rstrip("\n").split("\t")
            group[locus], gene[locus] = kind, name
    with open(os.path.join(folder, "examples.faa")) as handle:
        proteins = [line[1:].split()[0] for line in handle if line.startswith(">")]
    total = collections.Counter(group[p] for p in proteins)
    authors = collections.defaultdict(list)
    for label, terms in read_pairs(notebook).items():
        if label.startswith("gene="):
            authors[base(label[5:])].append(terms)
    print(f"Proteins: {total['host']} host, {total['phage']} phage, "
          f"from {len(glob.glob(os.path.join(GROUPS['host'], 'prokka_*')))} host and "
          f"{len(glob.glob(os.path.join(GROUPS['phage'], 'prokka_*')))} phage genomes.\n")

    print("| Run | Identity / coverage | Host with a hit | Phage with a hit |\n|---|---|---|---|")
    for run in runs:
        raw = list(read_hits(os.path.join(run, "diamond_hits.tsv")))
        for identity, coverage in THRESHOLDS:
            found = collections.Counter(group[q] for q in best_hits(raw, identity, coverage, 1e-5))
            print(f"| {os.path.basename(run)} | {identity}% / {coverage}% | "
                  f"{100 * found['host'] / total['host']:.1f}% | {100 * found['phage'] / total['phage']:.1f}% |")

    print("\n| Run | Database | Proteins with GO (host / phage) | Pairs dropped | "
          "Gene names in both | Reproduced exactly | Median best Jaccard |\n|---|---|---|---|---|---|---|")
    for run in runs:
        with open(os.path.join(run, "provenance.json")) as handle:
            prov = json.load(handle)
        ours = read_pairs(os.path.join(run, "pairs.tsv"))
        with_go = collections.Counter(group[q] for q in ours)
        named = collections.defaultdict(list)
        for q, terms in ours.items():
            if gene.get(q):
                named[base(gene[q])].append(terms)
        shared = [n for n in named if n in authors]
        exact = sum(any(a == o for a in authors[n] for o in named[n]) for n in shared)
        jaccard = [max(len(a & o) / len(a | o) for a in authors[n] for o in named[n]) for n in shared]
        print(f"| {os.path.basename(run)} | {prov['database']['sprot_file']} | "
              f"{with_go['host']} / {with_go['phage']} | {prov['pairs_dropped_retired_or_unknown']} | "
              f"{len(shared)} | {exact} ({100 * exact / len(shared):.1f}%) | {st.median(jaccard):.2f} |")


if __name__ == "__main__":
    if sys.argv[1:2] not in (["inputs"], ["compare"]) or len(sys.argv) < (3 if sys.argv[1] == "inputs" else 5):
        sys.exit(__doc__)
    if sys.argv[1] == "inputs":
        inputs(sys.argv[2])
    else:
        compare(sys.argv[2], sys.argv[3], sys.argv[4:])
