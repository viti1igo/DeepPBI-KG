"""Extract GO terms for any set of proteins, or take them from another tool, and summarise them.

Build the search database once, from the Swiss-Prot file that is also used to look up GO:

  python scripts/go_extract.py build-db uniprot_sprot.dat.gz sprot_db

Then run on proteins, or on another tool's output:

  python scripts/go_extract.py run proteins.faa out --db sprot_db --obo go-basic.obo
  python scripts/go_extract.py run sample.emapper.annotations out --format emapper --obo go-basic.obo
  python scripts/go_extract.py run interproscan.tsv out --format interproscan --obo go-basic.obo
  python scripts/go_extract.py run pairs.tsv out --format table --obo go-basic.obo

Writes to out: pairs.tsv (gene, GO ID), hits.tsv (FASTA input only), provenance.json, and the
summaries from go_rollup.py. Needs DIAMOND for FASTA input; otherwise standard library only.

FASTA input keeps each protein's best Swiss-Prot hit that passes the thresholds and copies that
entry's GO terms, so the result is a homology prediction, not an experimental annotation.
"""
import argparse
import datetime
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

from go_dag import Dag, RELATION_SETS, obo_version, parse_obo
from go_rollup import read_pairs, resolve_genes
from go_uniprot_per_protein import parse_sprot

HERE = os.path.dirname(os.path.abspath(__file__))
FIELDS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore", "qcovhsp"]
GO_ID = re.compile(r"GO:\d{7}")


def sprot_fasta(lines):
    """Yield (accession, sequence) from a Swiss-Prot flat file."""
    accession, seq, in_seq = None, [], False
    for line in lines:
        tag = line[:2]
        if tag == "ID":
            accession, seq, in_seq = None, [], False
        elif tag == "AC" and accession is None:
            accession = line[5:].split(";")[0].strip()
        elif tag == "SQ":
            in_seq = True
        elif tag == "//":
            yield accession, "".join(seq)
        elif in_seq:
            seq.append(line.replace(" ", "").strip())


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def diamond():
    path = shutil.which("diamond") or os.path.join(os.path.dirname(sys.executable), "diamond")
    if not os.path.exists(path):
        raise SystemExit("DIAMOND not found; install it, e.g. conda install -c bioconda diamond")
    return path


def build_db(sprot, out, threads):
    os.makedirs(out, exist_ok=True)
    fasta = os.path.join(out, "swissprot.faa")
    with gzip.open(sprot, "rt", encoding="utf-8") as handle, open(fasta, "w") as fa:
        count = 0
        for accession, seq in sprot_fasta(handle):
            fa.write(f">{accession}\n{seq}\n")
            count += 1
    subprocess.run([diamond(), "makedb", "--in", fasta, "-d", os.path.join(out, "swissprot"),
                    "--threads", str(threads), "--quiet"], check=True)
    os.remove(fasta)
    info = {"sprot_file": os.path.basename(sprot), "sprot_sha256": sha256(sprot), "entries": count}
    with open(os.path.join(out, "db.json"), "w") as handle:
        json.dump(info, handle, indent=2)
    print(f"{count} Swiss-Prot entries indexed in {out}")


def best_hits(rows, min_identity, min_coverage, max_evalue):
    """{query: best row by bit score} among rows that pass every threshold."""
    best = {}
    for row in rows:
        if (row["pident"] >= min_identity and row["qcovhsp"] >= min_coverage and row["evalue"] <= max_evalue
                and (row["qseqid"] not in best or row["bitscore"] > best[row["qseqid"]]["bitscore"])):
            best[row["qseqid"]] = row
    return best


def load_json(path):
    with open(path) as handle:
        return json.load(handle)


def lines(path):
    with open(path, encoding="utf-8") as handle:
        yield from handle


def read_hits(path):
    with open(path) as handle:
        for line in handle:
            values = line.rstrip("\n").split("\t")
            row = dict(zip(FIELDS, values))
            for key in FIELDS[2:]:
                row[key] = float(row[key])
            yield row


def go_by_accession(sprot):
    with gzip.open(sprot, "rt", encoding="utf-8") as handle:
        return {e["ac"]: set(e["go"]) for e in parse_sprot(handle)}


def from_fasta(args, out):
    db = load_json(os.path.join(args.db, "db.json"))
    raw = os.path.join(out, "diamond_hits.tsv")
    subprocess.run([diamond(), "blastp", "-q", args.input, "-d", os.path.join(args.db, "swissprot"),
                    "-o", raw, "--outfmt", "6", *FIELDS, "-k", "5", "--evalue", str(args.max_evalue),
                    "--threads", str(args.threads), "--quiet"], check=True)
    hits = best_hits(read_hits(raw), args.min_identity, args.min_coverage, args.max_evalue)
    go = go_by_accession(args.sprot)
    pairs = {q: go.get(h["sseqid"], set()) for q, h in hits.items()}
    with open(os.path.join(out, "hits.tsv"), "w") as handle:
        handle.write("\t".join(FIELDS + ["go_terms"]) + "\n")
        for q, h in sorted(hits.items()):
            handle.write("\t".join(str(h[f]) for f in FIELDS) + f"\t{len(pairs[q])}\n")
    queries = sum(1 for line in lines(args.input) if line.startswith(">"))
    return pairs, {"database": db, "proteins": queries, "proteins_with_hit": len(hits),
                   "thresholds": {"min_identity": args.min_identity, "min_query_coverage": args.min_coverage,
                                  "max_evalue": args.max_evalue},
                   "diamond": subprocess.run([diamond(), "version"], capture_output=True, text=True).stdout.strip()}


def read_emapper(path):
    """eggNOG-mapper annotations: the query and GOs columns."""
    pairs, header = {}, None
    for line in lines(path):
        if line.startswith("#query"):
            header = line[1:].rstrip("\n").split("\t")
        elif header and not line.startswith("#"):
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            pairs.setdefault(row["query"], set()).update(GO_ID.findall(row.get("GOs", "")))
    if header is None:
        raise SystemExit(f"{path}: no '#query' header, so not an eggNOG-mapper annotations file")
    return pairs


def read_interproscan(path):
    """InterProScan TSV: GO IDs from column 14, for every match of each protein."""
    pairs = {}
    for line in lines(path):
        fields = line.rstrip("\n").split("\t")
        pairs.setdefault(fields[0], set()).update(GO_ID.findall(fields[13]) if len(fields) > 13 else [])
    return pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build-db", help="index a Swiss-Prot flat file for DIAMOND")
    b.add_argument("sprot")
    b.add_argument("out")
    b.add_argument("--threads", type=int, default=4)
    r = sub.add_parser("run", help="extract or read GO terms, then summarise them")
    r.add_argument("input")
    r.add_argument("out")
    r.add_argument("--format", choices=["fasta", "emapper", "interproscan", "table"], default="fasta")
    r.add_argument("--obo", required=True)
    r.add_argument("--db", help="folder from build-db (FASTA input)")
    r.add_argument("--sprot", help="the Swiss-Prot flat file the database was built from (FASTA input)")
    # ponytail: one fixed set of thresholds, chosen on 2026-11-12 after the DeepPBI-KG validation;
    # per-dataset tuning would need a labelled set for that dataset.
    r.add_argument("--min-identity", type=float, default=30)
    r.add_argument("--min-coverage", type=float, default=50)
    r.add_argument("--max-evalue", type=float, default=1e-5)
    r.add_argument("--threads", type=int, default=4)
    r.add_argument("--slim", default="goslim_metagenomics")
    r.add_argument("--level", type=int, default=3)
    args = parser.parse_args()

    if args.command == "build-db":
        return build_db(args.sprot, args.out, args.threads)

    os.makedirs(args.out, exist_ok=True)
    details = {}
    if args.format == "fasta":
        if not (args.db and args.sprot):
            parser.error("FASTA input needs --db and --sprot")
        if load_json(os.path.join(args.db, "db.json"))["sprot_sha256"] != sha256(args.sprot):
            parser.error("--sprot is not the file the database was built from")
        pairs, details = from_fasta(args, os.path.abspath(args.out))
    else:
        pairs = {"emapper": read_emapper, "interproscan": read_interproscan, "table": read_pairs}[args.format](args.input)
    pairs = {g: t for g, t in pairs.items() if t}
    if not pairs:
        raise SystemExit("no GO terms found in the input")

    table = os.path.join(args.out, "pairs.tsv")
    with open(table, "w") as handle:
        for gene, terms in sorted(pairs.items()):
            handle.writelines(f"{gene}\t{t}\n" for t in sorted(terms))
    terms, alt = parse_obo(args.obo)
    by_domain, dropped = resolve_genes(pairs, terms, Dag(terms, alt, RELATION_SETS["is_a+part_of"]))
    kept = set().union(*(genes.keys() for genes in by_domain.values()))
    provenance = {"command": sys.argv, "date": datetime.date.today().isoformat(), "input": args.input,
                  "input_sha256": sha256(args.input), "format": args.format, "go_release": obo_version(args.obo),
                  "genes_with_go": len(pairs), "pairs": sum(map(len, pairs.values())),
                  "pairs_dropped_retired_or_unknown": dropped, "genes_left_without_go": len(pairs) - len(kept),
                  **details}
    with open(os.path.join(args.out, "provenance.json"), "w") as handle:
        json.dump(provenance, handle, indent=2)
    print(f"{provenance['genes_with_go']} genes with GO terms, {provenance['pairs']} pairs, "
          f"{dropped} pairs dropped as retired or unknown in GO {provenance['go_release']}")
    subprocess.run([sys.executable, os.path.join(HERE, "go_rollup.py"), args.obo, table, args.out,
                    "--slim", args.slim, "--level", str(args.level)], check=True)


if __name__ == "__main__":
    main()
