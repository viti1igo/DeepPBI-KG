# GO scripts

Scripts for working with Gene Ontology (GO) annotations: extract GO terms for any set of proteins, clean the IDs, move terms up the GO graph to a GO slim or a level, draw per-gene trees, and record where every result came from. They were built and tested on the DeepPBI-KG dataset. Results and the reasons for the defaults are in [docs/GO_TASK_STATUS.md](../docs/GO_TASK_STATUS.md), [docs/GO_LITERATURE_REVIEW.md](../docs/GO_LITERATURE_REVIEW.md) and [docs/GO_EXTRACT_VALIDATION.md](../docs/GO_EXTRACT_VALIDATION.md).

## Requirements

- Python 3.9 or later. The scripts use only the standard library.
- [DIAMOND](https://github.com/bbuchfink/diamond), only for protein FASTA input. Tested with 2.2.8: `conda install -c bioconda diamond`.

Run every command from the repo root.

## Data

| Data | Used for | Download | Version used here |
|---|---|---|---|
| GO graph, `go-basic.obo` | Every script | https://purl.obolibrary.org/obo/go/go-basic.obo (about 32 MB) | releases/2026-07-26. `code/GO/go-basic.obo` holds releases/2022-03-22 |
| Swiss-Prot flat file, `uniprot_sprot.dat.gz` | GO lookup, DIAMOND database, measurements | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz (about 0.7 GB) | 2026_03 |
| Swiss-Prot, an older release | Matching the DeepPBI-KG authors' annotations | `https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2022_01/knowledgebase/uniprot_sprot-only2022_01.tar.gz` (about 1.5 GB; extract `uniprot_sprot.dat.gz`) | 2022_01 |
| DeepPBI-KG gene-to-GO pairs | Test input | Cell 3 of `code/富集分析data_process.ipynb` in this repo | 11,680 labels, 54,931 pairs |

Both links above marked `current_release` and `go-basic.obo` always serve the newest release, so record the release you used. `go-basic.obo` states its release in its `data-version` line; Swiss-Prot states it in `reldate.txt` next to the download.

```bash
mkdir -p go_data
curl -L -o go_data/go-basic.obo https://purl.obolibrary.org/obo/go/go-basic.obo
curl -L -o go_data/uniprot_sprot_2026_03.dat.gz https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
```

GO and UniProt data are released under CC BY 4.0. Credit the [Gene Ontology Consortium](https://geneontology.org/docs/go-citation-policy/) and [UniProt](https://www.uniprot.org/help/publications) when you publish results.

`go_data/` is ignored by git, so downloads and outputs stay local.

## Main command: `go_extract.py`

Turns proteins, or another tool's GO output, into cleaned GO terms and summaries.

### 1. Build the search database once

```bash
python scripts/go_extract.py build-db go_data/uniprot_sprot_2026_03.dat.gz go_data/sprot_2026_03_db
```

This builds the DIAMOND database from the same Swiss-Prot file that is later used to look up GO, so the two can't come from different releases. It takes about 40 seconds and writes `swissprot.dmnd` and `db.json`, which records the file name, its SHA-256 and the number of entries.

### 2. Run on your data

| Input | Command |
|---|---|
| Protein FASTA | `python scripts/go_extract.py run proteins.faa out --db go_data/sprot_2026_03_db --sprot go_data/uniprot_sprot_2026_03.dat.gz --obo go_data/go-basic.obo` |
| eggNOG-mapper | `python scripts/go_extract.py run sample.emapper.annotations out --format emapper --obo go_data/go-basic.obo` |
| InterProScan TSV | `python scripts/go_extract.py run sample.tsv out --format interproscan --obo go_data/go-basic.obo` |
| Plain table | `python scripts/go_extract.py run pairs.tsv out --format table --obo go_data/go-basic.obo` |

A plain table has one gene per line: the gene ID, a tab, then GO IDs separated by commas or semicolons. A gene can also appear on several lines.

**How FASTA input is annotated:** DIAMOND finds up to 5 Swiss-Prot hits per protein. The script keeps the hit with the highest bit score among those that pass every threshold, and copies that entry's GO terms. The result is a prediction by similarity, not an experimental annotation.

| Option | Default | Meaning |
|---|---|---|
| `--min-identity` | 30 | Minimum % identical amino acids in the aligned region |
| `--min-coverage` | 50 | Minimum % of the query protein covered by the alignment |
| `--max-evalue` | 1e-5 | Maximum E-value |
| `--threads` | 4 | DIAMOND threads |
| `--slim` | `goslim_metagenomics` | GO slim to summarise with; any subset named in the OBO file works |
| `--level` | 3 | Level for the secondary summary, counted from the domain root |

The defaults were chosen after the DeepPBI-KG validation. Swiss-Prot has few phage proteins, so only about 5% of phage proteins get GO this way; for phage data, eggNOG-mapper or InterProScan output is a better input.

### 3. Outputs

| File | Contents |
|---|---|
| `pairs.tsv` | One line per gene and GO ID, as found, before any clean-up |
| `hits.tsv` | FASTA input only: the hit kept for each protein, with identity, coverage, E-value and bit score |
| `diamond_hits.tsv` | FASTA input only: all DIAMOND hits, before the thresholds |
| `provenance.json` | Command, date, input file and its SHA-256, GO release, database file and its SHA-256, thresholds, DIAMOND version, protein counts, and the number of pairs dropped because their GO ID is retired or unknown in this GO release |
| `gene_terms.tsv` | Every live term per gene, its depths, whether it is a deepest term and whether it is a GO leaf |
| `rollup_goslim_metagenomics.tsv` | Each deepest term with the most specific GO slim terms above it, or `other` |
| `slim_summary.csv` | Per GO domain: genes, deepest terms, % covered by the slim, distinct slim terms, slim terms per gene |
| `rollup_level_3.tsv` | Each deepest term with its ancestors at the chosen level |
| `levels.csv` | For every level: distinct terms, terms per gene, and how many deepest terms are already above it |
| `measurements.csv` | Whether genes list whole branches or independent terms, for three sets of GO edges |

**How the terms are cleaned:** alternate IDs become their current ID, and retired or unknown IDs are dropped and counted. Each gene's terms are split by GO domain, and only its deepest terms are kept: those that aren't an ancestor of another listed term. Moving up follows `is_a` and `part_of` edges.

## Other scripts

| Script | What it does | Example |
|---|---|---|
| `go_dag.py` | Library: reads `go-basic.obo` and gives parents, ancestors, longest and shortest depth, terms at a level, GO slim terms and alternate IDs | Imported by the others |
| `go_rollup.py` | Runs the clean-up and summaries on a gene-to-GO table or the DeepPBI-KG notebook. `go_extract.py run` calls it | `python scripts/go_rollup.py go_data/go-basic.obo code/富集分析data_process.ipynb go_rollup_out --genes host_key=model/host_key_gene_0.00004.csv --genes phage_key=model/phage_key_gene_0.0001.csv` |
| `go_tree.py` | One HTML page per gene showing its GO terms in the GO graph, with the data sources stated on the page | `python scripts/go_tree.py go_data/go-basic.obo code/富集分析data_process.ipynb go_data/trees --gene gene=dnaA --gene gene=aarA --trace go_data/sprot_2022_01/trace.csv` |
| `go_trace_sprot.py` | For each gene, lists the Swiss-Prot entries that have exactly its GO IDs, to trace where annotations came from | `python scripts/go_trace_sprot.py code/富集分析data_process.ipynb go_data/sprot_2022_01/uniprot_sprot.dat.gz 2022_01 go_data/sprot_2022_01/trace.csv` |
| `go_check_quickgo.py` | Checks `go_dag.py`'s ancestors against EBI's QuickGO on a random sample, the same for a given seed | `python scripts/go_check_quickgo.py go_data/go-basic.obo go_data/quickgo_check.csv --seed 0` |
| `go_uniprot_per_protein.py` | Measures, per Swiss-Prot protein, whether its GO terms form branches or independent nodes | `python scripts/go_uniprot_per_protein.py go_data/go-basic.obo go_data/uniprot_sprot_2026_03.dat.gz go_data/go_per_protein.csv` |
| `go_extract_validate.py` | Reproduces the DeepPBI-KG validation tables. Specific to this repo | See [docs/GO_EXTRACT_VALIDATION.md](../docs/GO_EXTRACT_VALIDATION.md) |
| `download_data.sh` | Restores DeepPBI-KG's genomes and BLAST databases; `--with-swissprot` adds the authors' SwissProt DIAMOND database | `scripts/download_data.sh` |
| `check_feature_similarity.py` | Rechecks the WGS vs key gene feature comparison | `python scripts/check_feature_similarity.py` |

`go_tree.py` options: `--gene` (repeat for each gene), `--level` (default 3), `--trace` (a table from `go_trace_sprot.py`), `--depth longest|shortest` and `--relations`. `go_check_quickgo.py` needs internet access to `www.ebi.ac.uk`.

## Tests

```bash
python tests/test_go_dag.py
python tests/test_go_extract.py
```

`test_go_extract.py` uses only small inline data and `tests/fixtures/mini.obo`, and runs one real DIAMOND search when DIAMOND is installed. `test_go_dag.py` uses `code/GO/go-basic.obo` and the DeepPBI-KG notebook from this repo. `python scripts/go_uniprot_per_protein.py --selftest` runs that script's own check.

## Limits

- GO from FASTA is copied from the most similar Swiss-Prot protein, so it can be wrong where proteins are similar but differ in function.
- Results depend on the GO and Swiss-Prot releases, so always keep `provenance.json` with the results.
- GO slims don't fit every organism. `goslim_metagenomics` covers only about half of phage cellular component terms.
