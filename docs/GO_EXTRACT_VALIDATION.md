# Validation of `go_extract.py`

Run on 2026-11-12. The tool turns protein sequences, or another tool's GO output, into cleaned GO terms and GO slim and level summaries. This report checks both routes.

## Summary

- **Protein route:** on DeepPBI-KG's 17 annotated example genomes, 78.9% of gene names found in both our output and the authors' saved table get exactly the authors' GO terms, when the database is Swiss-Prot 2022_01. With Swiss-Prot 2026_03 it's 17.3%, because UniProt has changed its GO annotations since then.
- **Phage proteins rarely get GO from Swiss-Prot:** 5.3% with the default thresholds, against 55.1% of host proteins.
- **eggNOG-mapper reader:** on PBIP's real file, every protein's GO terms match the file and PBIP's own table.
- **Retired IDs:** 23% of eggNOG's GO pairs and 3.8% of Swiss-Prot 2022_01's are retired or unknown in GO 2026-07-26, so they're dropped. `provenance.json` now records this count.

## Data and versions

| Item | Version |
|---|---|
| GO graph | `go-basic.obo` releases/2026-07-26 |
| Swiss-Prot, older | 2022_01 (`uniprot_sprot.dat.gz` from the release-2022_01 archive, sha256 starts `0405ebee8d0d`) |
| Swiss-Prot, current | 2026_03 (sha256 starts `c87208822a56`) |
| DIAMOND | 2.2.8 |
| Thresholds | identity ≥ 30%, query coverage ≥ 50%, E-value ≤ 1e-5, best hit by bit score |
| Proteins | 20,481 from Prokka: 19,844 from 8 host genomes, 637 from 9 phage genomes |

Swiss-Prot 2022_01 is used because the authors' saved GO sets match it: 11,663 of their 11,680 labels have exactly the GO IDs of an entry in that release (`scripts/go_trace_sprot.py`).

## Protein route: comparison with the authors' table

The authors' query sequences aren't published, so the comparison is by Prokka gene name. A name counts as reproduced if any of our proteins with that name has exactly the GO IDs of any of the authors' labels with that name.

| Database | Proteins with GO (host / phage) | Pairs dropped | Gene names in both | Reproduced exactly | Median best Jaccard |
|---|---|---|---|---|---|
| Swiss-Prot 2022_01 | 10,404 / 28 | 1,707 | 2,639 | 2,081 (78.9%) | 1.00 |
| Swiss-Prot 2026_03 | 10,408 / 30 | 60 | 2,629 | 456 (17.3%) | 0.60 |

Jaccard is the overlap of two GO sets: shared terms divided by all terms. For each gene name we took the best pair.

## Protein route: effect of the thresholds

Share of proteins with a hit that passes the thresholds (E-value ≤ 1e-5 in every row):

| Database | Identity / coverage | Host | Phage |
|---|---|---|---|
| 2026_03 | none | 61.3% | 8.6% |
| 2026_03 | **30% / 50% (default)** | **55.1%** | **5.3%** |
| 2026_03 | 40% / 50% | 43.6% | 1.3% |
| 2026_03 | 50% / 70% | 32.6% | 0.2% |

The 2022_01 database gives figures within 0.3 points of these.

## GO slim coverage

Share of deepest terms that `goslim_metagenomics` covers (F / P / C):

| Run | Coverage |
|---|---|
| DeepPBI-KG examples, Swiss-Prot 2026_03 | 95.9 / 94.1 / 95.3% |
| DeepPBI-KG examples, Swiss-Prot 2022_01 | 97.0 / 91.9 / 96.3% |
| PBIP phage proteins, eggNOG-mapper | 87.7 / 80.4 / 50.6% |

## eggNOG-mapper reader on a real file

PBIP's `phage_eggnog.emapper.annotations` (eggNOG-mapper 2.1.15) was used only to test the reader. Its results aren't part of the DeepPBI-KG evaluation.

| Check | Result |
|---|---|
| Proteins in the file | 8,404 |
| Proteins with GO, read directly from the GOs column | 2,877 |
| Proteins with GO, from the reader | 2,877, each with the same GO set |
| PBIP's own table (`unified_go_kegg.tsv`): eggNOG IDs with GO | 1,230, all the same GO set |
| Pairs dropped as retired or unknown | 31,195 of 133,981 (23%) |

## Limits

- **Lenient comparison:** matching by gene name lets one good protein stand for a name. The authors' table may not even include these 17 genomes.
- **Not a test of biological truth:** reproducing the authors' table shows the method matches theirs, not that the GO terms are correct. Both copy GO from the nearest Swiss-Prot protein.
- **Few phage proteins:** 637 phage proteins from 9 genomes. The phage result is clear, but small.
- **InterProScan reader:** tested only on a small made-up sample, not on a real file.
- **Thresholds:** chosen once from this run, not tuned against labelled data.

## Reproduce

From the repo root:

```bash
python scripts/go_extract.py build-db go_data/sprot_2022_01/uniprot_sprot.dat.gz go_data/sprot_2022_01/db
python scripts/go_extract.py build-db go_data/uniprot_sprot_2026_03.dat.gz go_data/sprot_2026_03_db
python scripts/go_extract_validate.py inputs go_data/extract_test
python scripts/go_extract.py run go_data/extract_test/examples.faa go_data/extract_test/out_2022_01 --db go_data/sprot_2022_01/db --sprot go_data/sprot_2022_01/uniprot_sprot.dat.gz --obo go_data/go-basic.obo
python scripts/go_extract.py run go_data/extract_test/examples.faa go_data/extract_test/out_2026_03 --db go_data/sprot_2026_03_db --sprot go_data/uniprot_sprot_2026_03.dat.gz --obo go_data/go-basic.obo
python scripts/go_extract_validate.py compare go_data/extract_test code/富集分析data_process.ipynb go_data/extract_test/out_2022_01 go_data/extract_test/out_2026_03
python scripts/go_extract.py run ../PBIP/results/eggnog/phage_eggnog.emapper.annotations go_data/extract_test/out_pbip_eggnog --format emapper --obo go_data/go-basic.obo
```

Downloads: Swiss-Prot 2022_01 from `https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2022_01/knowledgebase/uniprot_sprot-only2022_01.tar.gz`, and the current release and `go-basic.obo` as listed in `docs/GO_TASK_STATUS.md`. Outputs in `go_data/` aren't committed.
