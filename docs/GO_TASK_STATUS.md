# GO task status

Last checked: 2026-11-06.

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

## Status: the GO hierarchy work has not been done

The repo has flat GO/KEGG enrichment from the upstream paper, but nothing that walks the GO hierarchy. `code/parse_go_obofile.py` keeps only ID, name and namespace, and drops the `is_a` and `relationship` edges. A search of the repo found no traversal, AmiGO, biclustering, Snakemake or Nextflow code. Searches of `PBIP`, `sphae`, `tools`, `datasets` and `archive` under `/data/axnguyen` found no GO work either; the rest of that directory wasn't searched.

| Task | Status | Evidence |
|---|---|---|
| Branch or deepest node? | Preliminary answer below | `scripts/go_annotation_depth.py`, added 2026-11-06 |
| Visualize GO mappings | Not started | `code/GO_KEGG.R` makes enrichment bar and dot plots only, with no DAG view |
| Auto-traverse script | Not started | No ancestor or level code |
| Literature review of traversal depth | Not started in this repo | The paper is linked in the notes only |
| Snakemake or Nextflow | Not started in this repo | `/data/axnguyen/sphae` has `.snakemake` folders, so it may be a useful example |

## Preliminary answer: mostly neither

UniProt-derived annotations here are usually several separate nodes at different depths. They are often not the deepest possible terms, and they don't list a whole branch.

Measured on 54,912 annotations for 11,680 genes, using `is_a` and `part_of` edges:

| Measure | Result | What it suggests |
|---|---|---|
| Annotations on a leaf term | 31.7% (58.2% of all live GO terms are leaves) | Usually not the deepest node |
| Median depth | 5 (maximum 14; all live terms: 7) | Shallower than the average term |
| Terms per gene, median | 4 | Several separate nodes per gene |
| Annotations with an ancestor on the same gene | 8.2% (4,500), on 19.5% of genes (2,283) | Rarely a listed branch |
| Live terms with more than one parent | 54.5% | "Up the tree" is often several paths, so a level isn't one fixed node |

So upward traversal would have to be computed from the DAG; the annotations don't already contain it.

**Limits of this result:**
- It uses GO IDs saved in notebook cell 3, not a fresh UniProt query.
- Those IDs come from DIAMOND hits to SwissProt, combined per gene. Merging several hits may create some of the ancestor overlap.
- Genes are keyed by name, such as `gene=mdtE_1`, so the notebook's naming may affect the per-gene counts.
- `go-basic.obo` is from 2022-03-22, and 24 annotations point to terms now marked obsolete.

To confirm it, repeat the measurement per protein on current UniProt or GOA data.

## Suggested order

1. Confirm the answer above per protein using a current GO release and UniProt data.
2. Fix `code/UniProt2GO_annotate.py:15` before generating new mappings (see `docs/CODE_REVIEW.md`).
3. Build a small GO DAG module with parents, ancestors and depth, keeping alternate IDs. Add tests using a few hand-checked terms.
4. Add an upward-traversal script. Decide first how to define a "level", since a term can have several paths of different lengths: shortest path, longest path, or a GO slim.
5. Visualize mappings for a few example genes, using AmiGO or the local DAG module.
6. Do the literature review and use it to choose a default level.
7. Compare Snakemake and Nextflow only once the steps above are worth automating.

## Questions for the user

- Should the analysis use all genes, or only key genes?
- Is the goal to feed GO levels into the model, or to interpret results?
- Which organisms matter most?
