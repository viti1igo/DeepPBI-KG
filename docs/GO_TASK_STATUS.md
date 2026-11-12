# GO task status

Last checked: 2026-11-12. How to run every script is in [scripts/README.md](../scripts/README.md).

## Original task

From the user's notes:

> The Gene Ontology is structured as a directed acyclic graph (DAG) divided into three independent domains: Molecular Function F, Biological Process P, and Cellular Component C.
> [GO - Ontology Lookup Service]
>
> Next tasks
> - **Are GO IDs recorded by UniProt a branch or independent node of highest depth (furthest from root)**
> - Visualize GO mappings extracted from UniProt
>   - AmiGO 2 Manual: Visualize - GO Wiki
> - How far should we traverse up the tree?
>   - A script to auto-traverse to a certain level
>   - Literature review of where past studies traverse to
>     - Detection of Gene Ontology Clusters Using Biclustering Algorithms (SN Computer Science, Springer Nature)
> - Extra task: look into workflow orchestration tools like Snakemake or Nextflow

KEGG was set aside, and the work covers GO only.

## Summary

| Question | Answer | Where |
|---|---|---|
| Branch or deepest node? | Mostly independent nodes, not branches, for molecular function (F) and biological process (P). Cellular component (C) is weaker: about a third of its terms sit below another listed term in genes with two or more C terms | [Question 1](#question-1-branch-or-independent-nodes) |
| Visualize GO mappings | Trees for 4 example genes, with their data sources on each page. Checked against QuickGO (300 terms) and, by hand, AmiGO (1 gene) | [Question 2](#question-2-visualize-go-mappings) |
| How far up to traverse? | No standard level exists. Default: `goslim_metagenomics`, with level 3 as a secondary view | [Question 3](#question-3-how-far-up-to-traverse), [GO_LITERATURE_REVIEW.md](GO_LITERATURE_REVIEW.md) |
| Snakemake or Nextflow | Not started | |
| GO for any dataset (added later) | `scripts/go_extract.py`, validated on DeepPBI-KG | [GO_EXTRACT_VALIDATION.md](GO_EXTRACT_VALIDATION.md) |

## Data used

| Data | Version | Source |
|---|---|---|
| GO graph | `go-basic.obo` releases/2026-07-26 | https://purl.obolibrary.org/obo/go/go-basic.obo |
| GO graph, older | `go-basic.obo` releases/2022-03-22, in `code/GO/` | Shipped with the upstream repo |
| Swiss-Prot, current | 2026_03 (575,748 entries) | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz |
| Swiss-Prot, older | 2022_01 (566,996 entries) | https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2022_01/knowledgebase/uniprot_sprot-only2022_01.tar.gz |
| DeepPBI-KG gene-to-GO pairs | 11,680 labels, 54,931 pairs, printed in cell 3 of `code/富集分析data_process.ipynb` | Upstream notebook |
| DeepPBI-KG key genes | `model/host_key_gene_0.00004.csv`, `model/phage_key_gene_0.0001.csv` | Upstream repo |

The `current_release` link always points to the newest release, so a later download gives different numbers. Downloads go in `go_data/`, which isn't committed.

**Where DeepPBI-KG's gene-to-GO pairs come from:** the notebook built them by matching proteins to Swiss-Prot with DIAMOND, but that output isn't published. 11,663 of the 11,680 labels have exactly the GO IDs of one Swiss-Prot 2022_01 entry, and so do all 10 DIAMOND hits still printed in the notebook (`scripts/go_trace_sprot.py`). The 17 that don't match are labels such as `protein=putative`, which seem to merge several proteins. So each label carries the GO terms of a single Swiss-Prot entry from about that release. 2022_01 fits, but it isn't confirmed as the release the authors used.

## Question 1: branch or independent nodes?

### On Swiss-Prot

Measured per protein and per GO domain on Swiss-Prot 2026_03 (551,801 proteins with GO terms) and GO 2026-07-26 (`scripts/go_uniprot_per_protein.py`).

- **Not branches:** a protein's deepest listed term shows only 0.1–2.7% of its branch on average, not counting the root.
- **Independent:** about 87–94% of proteins never list a term with its ancestor, and about 99% in bacteria.
- **Not always the bottom of GO:** the listed terms are where the annotation stops, not necessarily GO leaves.

So upward traversal has to be computed from the GO graph; the annotations don't already contain it.

Using `is_a` and `part_of` edges, each cell F / P / C:

| Measure | All Swiss-Prot | Bacteria | Phages |
|---|---|---|---|
| Proteins with GO terms in the domain | 482,781 / 449,845 / 467,132 | 301,770 / 269,694 / 260,506 | 813 / 1,012 / 781 |
| Lists a term and its ancestor | 5.8 / 8.8 / 13.0% | 1.0 / 1.2 / 0.8% | 6.3 / 2.5 / 2.4% |
| Only one term | 35.3 / 53.9 / 65.2% | 34.4 / 63.6 / 84.3% | 50.1 / 44.2 / 79.0% |
| Several terms, all independent | 58.9 / 37.3 / 21.8% | 64.5 / 35.3 / 14.9% | 43.7 / 53.4 / 18.6% |
| Average share of branch listed | 1.0 / 0.9 / 2.7% | 0.1 / 0.1 / 0.2% | 0.6 / 0.3 / 0.9% |
| Deepest listed terms with no ancestor listed | 95.7 / 91.0 / 83.8% | 99.4 / 99.1 / 99.3% | 96.2 / 97.7 / 97.8% |
| Annotations on a GO leaf | 60.3 / 29.1 / 15.1% | 68.0 / 32.8 / 15.2% | 30.7 / 35.0 / 27.9% |
| Median annotation depth | 5 / 7 / 4 | 5 / 7 / 3 | 4 / 7 / 3 |

For comparison, 80.3 / 49.5 / 68.6% of all live GO terms are leaves, and their median depth is 5 / 7 / 6.

The edges used change the result, mainly for component and process terms (all Swiss-Prot):

| Edges | Lists a term and its ancestor | Annotations on a GO leaf | Live terms with more than one parent |
|---|---|---|---|
| `is_a` | 5.7 / 7.6 / 4.1% | 60.3 / 35.6 / 34.1% | 18.9 / 50.9 / 13.9% |
| Adding `part_of` | 5.8 / 8.8 / 13.0% | 60.3 / 29.1 / 15.1% | 19.0 / 64.8 / 50.0% |
| Adding `regulates` | 5.8 / 10.6 / 13.0% | 60.3 / 21.6 / 15.1% | 19.0 / 72.5 / 50.0% |

Ancestor overlap comes mostly from non-electronic annotations in eukaryotes: 0.0 / 0.2 / 0.0% of proteins with electronic annotations only, and 16.0 / 22.8 / 31.5% with experimental annotations only. Eukaryotes hold 96.6% of all term-and-ancestor pairs.

### On DeepPBI-KG's saved genes

Same conclusion. Measured with `scripts/go_rollup.py` on the 11,680 labels, using GO 2026-07-26 and `is_a` and `part_of` edges. Each cell F / P / C.

| Measure | All saved genes | Host key genes | Phage key genes |
|---|---|---|---|
| Genes | 10,399 / 8,820 / 7,854 | 1,609 / 1,349 / 1,159 | 152 / 130 / 92 |
| Only one term | 33.0 / 51.9 / 57.0% | 32.3 / 49.4 / 55.2% | 21.1 / 40.0 / 67.4% |
| Several terms, all independent | 56.0 / 37.9 / 11.7% | 57.3 / 39.4 / 12.8% | 67.8 / 49.2 / 20.7% |
| Lists a term and its ancestor | 11.0 / 10.2 / 31.3% | 10.4 / 11.1 / 32.0% | 11.2 / 10.8 / 12.0% |
| Average share of branch listed | 1.5 / 0.7 / 8.2% | 1.3 / 0.9 / 8.4% | 0.9 / 1.2 / 2.8% |
| Annotations on a GO leaf | 48.9 / 28.7 / 8.2% | 50.6 / 29.7 / 6.9% | 46.5 / 30.4 / 13.8% |
| Median annotation depth | 5 / 7 / 3 | 5 / 7 / 3 | 5 / 7 / 3 |

**Why overlap is higher than in Swiss-Prot 2026_03 bacteria:**
- **The labels carry 2022-era annotations.** Each label matches a single Swiss-Prot 2022_01 entry, so the overlap isn't from merging several hits per gene, as an earlier version of this document suggested. It comes from what those entries listed at the time.
- **One retired term inflates the component number.** `GO:0016021` (integral component of membrane) is now an alternate ID for `membrane`, so it can look like an ancestor of `plasma membrane`. 1,691 genes list both IDs. Leaving out `GO:0016021` lowers component overlap from 31.3% to 9.2%.

### Checks on the measure

| Check | Result |
|---|---|
| Synthetic: add every ancestor to the P terms of 3,000 saved genes that have two or more | Share of terms that are ancestors of another listed term rises to 87.0% |
| Synthetic: keep only each gene's deepest terms | Falls to 0.0% |
| Known branches: PBIP's eggNOG-mapper lists | 75.0 / 80.5 / 64.1% |
| DeepPBI-KG's saved genes, all genes | 3.8 / 3.7 / 14.9% |
| DeepPBI-KG's saved genes, only genes with two or more terms | 5.6 / 7.6 / 34.7% |
| Second measure: share of all ancestors that are listed, DeepPBI-KG vs PBIP | 1.9 / 1.4 / 15.8% vs 86.7 / 92.6 / 92.7% |

The measure separates branches from independent nodes, and a second measure agrees. Genes with a single term can't list an ancestor, so they lower the average: among genes with two or more terms, component reaches 34.7%. **The independent-node answer holds for F and P, but only partly for C.** These checks were run once and aren't saved as a script.

**Limits:**
- Swiss-Prot only, so it doesn't describe unreviewed proteins.
- 3,821 Swiss-Prot annotations (0.12%), and 2,273 of DeepPBI-KG's pairs (4.1%), use IDs that are retired or unknown in GO 2026-07-26 and were dropped.
- Phages in Swiss-Prot are a small group: 2,097 proteins, matched by a name containing "phage" or by a bacterial host.
- Only 1,758 of 3,102 host key genes and 160 of 274 phage key genes have saved GO terms.
- Ancestors are plain reachability over the chosen edges, without GO's relation-chain rules. This mainly affects the `regulates` row.

## Question 2: visualize GO mappings

`scripts/go_tree.py` draws one HTML page per gene. Each GO domain is drawn from its root, one row per depth. The page marks the terms the gene lists, the ancestors added from the GO graph, and the terms it rolls up to at level 3. Each page also states the input file and its hash, the GO release, the edges followed, and which Swiss-Prot entries have exactly that gene's GO IDs.

Trees exist for `gene=dnaA`, `gene=aarA`, `gene=acm` and `gene=marA`, chosen as examples:
- **Main set:** `go_data/trees/`, on GO 2026-07-26, chosen by the user on 2026-11-06.
- **For comparison:** `go_data/trees_2022/`, on GO 2022-03-22, the same era as the annotations.

The annotations are older than the main graph, so retired IDs are updated or dropped. For example, `aarA`'s `GO:0016021` becomes `membrane`.

| Check | Result |
|---|---|
| Ancestors computed by `scripts/go_dag.py` against QuickGO, 300 random terms (100 per domain, seed 0), `is_a` only | 300 of 300 match |
| The same with `is_a` and `part_of` | 278 match. All 22 differences are molecular function terms where QuickGO adds terms from another domain, because `go-basic.obo` leaves out edges between domains. QuickGO served GO 2026-09-27 (`go-plus`) |
| The 4 trees against the notebook, the OBO file and QuickGO | No term missing or extra; every edge exists in the OBO file; deepest terms recomputed match. Run once, not saved as a script |
| `gene=aarA` by hand in AmiGO 2.5.17 (data loaded 2026-08-06) | All 14 terms and 16 edges appear. AmiGO also draws three UBERON anatomy terms, which aren't GO terms. Saved as `go_data/trees/amigo_aarA.png` |

**Limits:** only 4 example genes, not chosen by any biological criterion. QuickGO and AmiGO check the graph, not the depths or the layout. The layout can cross edges where AmiGO's doesn't.

## Question 3: how far up to traverse?

The full review is in [GO_LITERATURE_REVIEW.md](GO_LITERATURE_REVIEW.md), with 27 references.

- **No standard level exists,** and a level is a poor proxy for how specific a term is. On the saved genes, Spearman correlation between level and information content is 0.11 / 0.11 / 0.33.
- **Default, approved by the user:** start from each gene's deepest terms, move up with `is_a` and `part_of`, and summarise with `goslim_metagenomics`. Level 3 stays as a secondary view.
- **Coverage on the saved genes:** `goslim_metagenomics` covers 97.0 / 92.5 / 94.4% of deepest terms with 35 / 39 / 31 terms. Level 3 needs 233 / 197 / 83 distinct terms.

`scripts/go_rollup.py` applies the default and can also report any level, by longest or shortest path.

How far each level reaches on the saved genes (longest path):

| Level | Distinct terms | Deepest terms already above this level |
|---|---|---|
| 1 | 19 / 14 / 3 | 0 / 0 / 0% |
| 2 | 69 / 63 / 40 | 1.1 / 0.0 / 0.5% |
| 3 | 233 / 197 / 83 | 10.8 / 2.0 / 7.6% |
| 4 | 540 / 447 / 139 | 23.3 / 7.6 / 65.2% |
| 5 | 1,483 / 703 / 192 | 35.0 / 21.7 / 88.6% |
| 6 | 1,857 / 924 / 235 | 62.7 / 33.2 / 93.2% |

## Extension: GO for any dataset

Added on 2026-11-12, with DeepPBI-KG as the test dataset. `scripts/go_extract.py` takes protein FASTA through DIAMOND and Swiss-Prot, or reads eggNOG-mapper, InterProScan or a plain table, then applies the steps above. Results and limits are in [GO_EXTRACT_VALIDATION.md](GO_EXTRACT_VALIDATION.md). In short: with Swiss-Prot 2022_01, 78.9% of gene names found in both reproduce the authors' GO terms exactly, and phage proteins rarely get GO from Swiss-Prot (5.3%).

## Decisions

- Scope of the first measurement: all reviewed proteins, with bacteria and phages also measured separately.
- KEGG set aside; GO only.
- Default summary: `goslim_metagenomics`, with level 3 as a secondary view (approved 2026-11-06).
- Main trees on GO 2026-07-26 (chosen 2026-11-06).
- Extraction thresholds: identity 30%, query coverage 50%, E-value 1e-5 (chosen 2026-11-12).
- PBIP is used only to test the eggNOG-mapper reader on a real file.

## Open questions

- Which dataset this task is for: confirm with the supervisor.
- Whether the goal is interpreting results or feeding GO into the model: assumed to be interpreting, so confirm with the supervisor. A model would need its own validation.
- Snakemake or Nextflow: not started.
