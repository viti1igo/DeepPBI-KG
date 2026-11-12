# How far up the Gene Ontology should we traverse?

Literature review for question 3 of the GO task ([GO_TASK_STATUS.md](GO_TASK_STATUS.md)). Written 2026-11-06.

## Summary

**The literature gives no standard GO level to stop at.** None of the papers reviewed here justifies a fixed number, and several warn against using level as a measure of specificity [3]–[6]. What past studies do depends on the goal:

| Goal | Common approach | Examples |
|---|---|---|
| Summarize or interpret genes | Map terms to a GO slim, or choose terms for the question | [6], [12], [13], [15] |
| Reduce an enrichment result | Keep all terms for testing, then remove redundant ones | [18], [19] |
| Predict function or interactions | Propagate terms upward, then filter or select features | [20]–[22] |
| Browse a fixed tree level | Offered by plotting tools, with known drawbacks | [17] |

**Recommendation for DeepPBI-KG, if the goal is interpreting results:**
1. Start from each gene's deepest listed terms.
2. Move up with `is_a` and `part_of` only.
3. Summarize with `goslim_metagenomics`, and report any uncovered terms as "other".
4. Keep level 3 as a secondary view, not as the biological unit.

The user approved this as the default on 2026-11-06, and `scripts/go_rollup.py` now applies it. It is still a recommendation, not a tested result. If the goal is model features, don't choose a level by hand; select features using training data only [20], [21].

## Scope and method

I searched Consensus, PubMed and Google Scholar, read papers through UTS Library access where needed, and followed references from key papers. This is a focused review, not a systematic one. There was no screening protocol, and studies may have been missed.

Citation details were checked against Europe PMC, PubMed or the publisher's page. Papers I could only read in abstract are marked in [Papers read](#papers-read).

## The paper linked in the task notes

Ali *et al.* [1] **don't traverse the GO tree**, so their paper doesn't directly answer this question. They:
- encode human proteins from STRING as binary vectors of about 17,900 GO terms and 271 KEGG pathways (the paper gives both 17,931 and 17,921 terms);
- keep associations with a mutual information score above 0.3, leaving 172;
- bicluster the result with QUBIC.

So they select terms by a data-driven score, not by level. Their highlighted terms are general ones. In our current GO file, `GO:1901363` (heterocyclic compound binding) and `GO:0043170` (macromolecule metabolic process) are both at level 3 by longest path, while `GO:0044428` and `GO:0044260` are now obsolete. These levels are my calculation, not the paper's claim.

Two limits matter if we reuse the method. Their negative examples are random protein pairs, and the paper doesn't say whether GO terms were propagated to ancestors.

DeepPBI-KG's own paper [2] uses GO enrichment on key genes to support interpretation. In the parts I read, it doesn't set a GO level.

## What the literature says

### A level isn't a measure of specificity

- Alterovitz *et al.* [3] report that DAG levels are poor indicators of specificity, and they partition GO by information content instead.
- Clark and Radivojac [4] note that terms at the same depth can have very different specificity.
- Pesquita *et al.* [5] explain why edge counts can mislead: terms and edges aren't evenly distributed, so edges at the same level needn't mean the same distance.
- Klopfenstein *et al.* [8] show how skewed the graph is. Among direct children of `biological_process`, descendant counts ran from over 18,000 to about 100.
- Gaudet and Dessimoz [9] describe the related "shallow annotation problem", and suggest information-theoretic measures as a partial remedy.

### Defining a level

Rhee *et al.* [6] define depth as the longest path from the root, and Manjang *et al.* [7] also treat the largest level as the correct one when a term has several. This supports the longest-path default already used in `scripts/go_rollup.py`.

Manjang *et al.* [7] also point out that a term's level alone is hard to interpret, because branches end at different depths.

### Which edges to follow

- Gaudet and Dessimoz [9] explain that `is_a` and `part_of` are transitive, but `regulates` isn't. Their example is a peptidase inhibitor reaching "proteolysis" through a regulation link, even though the inhibitor works against proteolysis.
- Hinderer *et al.* [10] treat `is_a` and `part_of` as relations of scope.
- GO's documentation [11] says annotations can be propagated in `go-basic`, which includes the regulation relations, so reachability alone doesn't settle this. The meaning of the annotation still changes when it crosses a regulation edge.

**So `is_a` + `part_of` is the safer default for grouping genes.** The `regulates` option should stay as a comparison only.

### GO slims and custom cuts

- The GO Consortium [12] describes subsets for summarizing gene sets and whole genomes.
- Davis *et al.* [13] generate dataset-specific slims from an information threshold, and report more statistical power than with full annotations. Jin and Lu [14] also select informative subsets automatically.
- MGnify [15] summarizes metagenomic GO annotations with `goslim_metagenomics`, which makes it the closest precedent here for microbial and viral sequence data.
- Rhee *et al.* [6] favour collapsing terms around the biological question before testing, and criticize tools limited to a fixed depth or to all nodes.
- Ye *et al.* [17] note that most WEGO users pick terms by tree level, and that this can bring in terms with no exact meaning.
- GOTrapper [16] lets users move up or down by a threshold on how many genes share a term.

### Enrichment results

These methods keep the hierarchy during analysis, then remove redundancy, rather than cutting at one level. GO Trimming [18] removes redundant terms along a path, and REVIGO [19] summarizes long, redundant term lists.

### GO terms as model features

- Wan and Freitas [20] compared hierarchical feature selection methods on GO features. HIP performed best overall on their aging-related gene datasets. Those aren't phage data, so the ranking may not carry over.
- DeepGO [21] propagates annotations, then keeps the most-annotated terms (932 BP, 589 MF, 436 CC). That is a frequency cut, not a level cut.
- Maetschke *et al.* [22] combine a GO-based inducer with a random forest for protein interactions, with AUC up to 0.88 across species.
- Bandyopadhyay and Mallick [23] build features from annotated terms and their ancestors, weighted by information content.
- TransformerGO [24] represents each protein as a set of GO terms.

None of these depends on one chosen level.

### Phage proteins

- PHROG [25] uses nine manually defined functional categories. GO was only one of the sources it compared against.
- GOPhage [26] reports that phage GO labels are very sparse: 912 terms across 1,173 UniProtKB phage proteins, or fewer than two proteins per term on average.

Sparse labels argue against a deep cut for phages, but I found no study that tests a GO level specifically on phage proteins. **This is a gap in the literature, or at least in this search.**

### Stability over time

Tomczak *et al.* [27] studied GO enrichment results from 2004 to 2015 and found they change as GO and its annotations evolve. Gaudet and Dessimoz [9] advise reporting releases and checking conclusions across several. Our own component overlap supports this: 13.9% on the 2022 GO file and 31.3% on the 2026 file.

## Checking this on our data

These are new one-off checks using `go-basic.obo` 2026-07-26, Swiss-Prot 2026_03 and DeepPBI-KG's saved genes, with `is_a` and `part_of`. Scratch code isn't in the repo. Each cell is F / P / C.

### Does a level predict specificity?

Barely. I computed information content from bacterial Swiss-Prot proteins, with annotations propagated upward.

| Measure | F | P | C |
|---|---|---|---|
| Spearman correlation, level vs information content | 0.11 | 0.11 | 0.33 |
| Level 3 information content, interquartile range (bits) | 8.20–15.62 | 7.90–15.23 | 8.69–14.99 |
| Level 3 information content, full range (bits) | 1.37–18.20 | 0.44–18.04 | 0.46–17.99 |

So two terms at level 3 can differ enormously in how specific they are, which matches [3]–[5]. Very rare terms also push values towards the maximum, so the top of each range says little.

### Which summary fits our genes?

Coverage is the share of each gene's deepest terms that reach at least one term in the subset, not counting roots.

| Subset | Terms available | All saved genes | Host key genes | Phage key genes |
|---|---|---|---|---|
| `goslim_metagenomics` | 35 / 40 / 33 | 97.0 / 92.5 / 94.4% | 96.8 / 91.6 / 93.6% | 98.9 / 93.4 / 90.8% |
| `goslim_prokaryote` | 34 / 49 / 14 | 65.2 / 82.9 / 87.0% | 65.6 / 81.9 / 86.6% | 64.1 / 73.0 / 80.7% |
| `goslim_generic` | 40 / 72 / 28 | 65.6 / 58.5 / 64.5% | 66.0 / 59.9 / 63.5% | 64.9 / 60.2 / 39.5% |
| `goslim_virus` | 0 / 13 / 7 | Not run | 0.0 / 1.6 / 0.1% | 0.0 / 8.0 / 0.0% |

`goslim_metagenomics` covers the most terms with a small vocabulary. `goslim_chembl` also covered 94.5 / 86.1 / 97.6% of deepest terms, but it has more terms and is aimed at drug discovery, so I wouldn't use it here.

For comparison, level 3 gives 233 / 197 / 83 distinct terms on all saved genes, while `goslim_metagenomics` uses 35 / 39 / 31.

## Recommendation

| Decision | Recommendation | Basis |
|---|---|---|
| Starting terms | Each gene's deepest listed terms | Our measurements in [GO_TASK_STATUS.md](GO_TASK_STATUS.md) |
| Edges | `is_a` + `part_of` | [9], [10] |
| Interpretation | `goslim_metagenomics`, with an "other" group | [12], [15] and the coverage above |
| Level 3 | Secondary view only | [3]–[6] and the correlations above |
| Level definition | Longest path | [6], [7] |
| Model features | Propagate, then select using training data only | [20], [21] |
| Reporting | Record GO and UniProt releases | [9], [27] |

**What would change this:**
- **Feeding GO into the model:** the slim recommendation wouldn't apply, and the choice would need its own validation.
- **A specific biological question:** a custom set of terms would be better than a general slim [6].
- **Rebuilt mappings:** the coverage above uses saved pairs that may include merged hits and the retired membrane ID, so it should be rechecked.

## Limits

- **Focused search:** it may have missed studies, especially ones that mention a level only in their methods.
- **Partial reading:** four papers were read in abstract only, and some others were checked for relevant passages rather than read in full.
- **Untested recommendation:** coverage shows the slim fits our terms, not that it gives better biological conclusions. A slim also isn't one term per gene: the median is 2–3 / 2 / 1 (F / P / C) for the key genes.
- **Chosen after seeing coverage:** I compared the subsets on the same genes the recommendation applies to.

## Papers read

**Full text, or the relevant sections:** [1], [2], [3]–[9], [10], [13], [15]–[18], [21], [24]–[27].

**Abstract or search-result summary only:** [14], [19], [22], [23]. For [20], I read the abstract and the publisher's summary of the results.

**PDFs in `papers/`** (ignored by git because of publisher copyright): Davis2010, Hassan2019, Hinderer2019, Hinderer2020, Jantzen2011, Klopfenstein2018, Manjang2020, Pesquita2009, Radivojac2013 (not cited here), Supek2011, Tomczak2018, Wan2018.

**Not downloaded:** Ali 2023 and Rhee 2008 didn't arrive through Chrome. Oxford Academic blocked automated downloads, and other papers failed too. All can be opened from the links below.

## References

[1] A. Ali, A. Ajil, A. Meenakshi Sundaram, and N. Joseph, "Detection of Gene Ontology clusters using biclustering algorithms," *SN Comput. Sci.*, vol. 4, Art. no. 217, 2023, doi: [10.1007/s42979-022-01624-w](https://doi.org/10.1007/s42979-022-01624-w).

[2] T. Wei *et al.*, "DeepPBI-KG: A deep learning method for the prediction of phage-bacteria interactions based on key genes," *Brief. Bioinform.*, vol. 25, no. 6, Art. no. bbae484, 2024, doi: [10.1093/bib/bbae484](https://doi.org/10.1093/bib/bbae484).

[3] G. Alterovitz, M. Xiang, M. Mohan, and M. F. Ramoni, "GO PaD: The Gene Ontology Partition Database," *Nucleic Acids Res.*, vol. 35, pp. D322–D327, 2007, doi: [10.1093/nar/gkl799](https://doi.org/10.1093/nar/gkl799).

[4] W. T. Clark and P. Radivojac, "Information-theoretic evaluation of predicted ontological annotations," *Bioinformatics*, vol. 29, no. 13, pp. i53–i61, 2013, doi: [10.1093/bioinformatics/btt228](https://doi.org/10.1093/bioinformatics/btt228).

[5] C. Pesquita, D. Faria, A. O. Falcão, P. Lord, and F. M. Couto, "Semantic similarity in biomedical ontologies," *PLoS Comput. Biol.*, vol. 5, no. 7, Art. no. e1000443, 2009, doi: [10.1371/journal.pcbi.1000443](https://doi.org/10.1371/journal.pcbi.1000443).

[6] S. Y. Rhee, V. Wood, K. Dolinski, and S. Draghici, "Use and misuse of the gene ontology annotations," *Nature Rev. Genet.*, vol. 9, no. 7, pp. 509–515, 2008, doi: [10.1038/nrg2363](https://doi.org/10.1038/nrg2363).

[7] K. Manjang, S. Tripathi, O. Yli-Harja, M. Dehmer, and F. Emmert-Streib, "Graph-based exploitation of gene ontology using GOxploreR for scrutinizing biological significance," *Sci. Rep.*, vol. 10, no. 1, Art. no. 16672, 2020, doi: [10.1038/s41598-020-73326-3](https://doi.org/10.1038/s41598-020-73326-3).

[8] D. V. Klopfenstein *et al.*, "GOATOOLS: A Python library for Gene Ontology analyses," *Sci. Rep.*, vol. 8, no. 1, Art. no. 10872, 2018, doi: [10.1038/s41598-018-28948-z](https://doi.org/10.1038/s41598-018-28948-z).

[9] P. Gaudet and C. Dessimoz, "Gene Ontology: Pitfalls, biases, and remedies," in *The Gene Ontology Handbook* (Methods in Molecular Biology, vol. 1446), C. Dessimoz and N. Škunca, Eds. New York, NY, USA: Humana Press, 2017, pp. 189–205, doi: [10.1007/978-1-4939-3743-1_14](https://doi.org/10.1007/978-1-4939-3743-1_14).

[10] E. W. Hinderer, R. M. Flight, R. Dubey, J. N. MacLeod, and H. N. B. Moseley, "Advances in gene ontology utilization improve statistical power of annotation enrichment," *PLoS One*, vol. 14, no. 8, Art. no. e0220728, 2019, doi: [10.1371/journal.pone.0220728](https://doi.org/10.1371/journal.pone.0220728).

[11] Gene Ontology Consortium, "Download ontology." Accessed: Nov. 6, 2026. [Online]. Available: <https://geneontology.org/docs/download-ontology/>

[12] Gene Ontology Consortium, "Guide to GO subsets." Accessed: Nov. 6, 2026. [Online]. Available: <https://geneontology.org/docs/go-subset-guide/>

[13] M. J. Davis, M. S. Sehgal, and M. A. Ragan, "Automatic, context-specific generation of Gene Ontology slims," *BMC Bioinform.*, vol. 11, Art. no. 498, 2010, doi: [10.1186/1471-2105-11-498](https://doi.org/10.1186/1471-2105-11-498).

[14] B. Jin and X. Lu, "Identifying informative subsets of the Gene Ontology with information bottleneck methods," *Bioinformatics*, vol. 26, no. 19, pp. 2445–2451, 2010, doi: [10.1093/bioinformatics/btq449](https://doi.org/10.1093/bioinformatics/btq449).

[15] A. L. Mitchell *et al.*, "MGnify: The microbiome analysis resource in 2020," *Nucleic Acids Res.*, vol. 48, no. D1, pp. D570–D578, 2020, doi: [10.1093/nar/gkz1035](https://doi.org/10.1093/nar/gkz1035).

[16] H. Hassan and S. Shanak, "GOTrapper: A tool to navigate through branches of gene ontology hierarchy," *BMC Bioinform.*, vol. 20, no. 1, Art. no. 20, 2019, doi: [10.1186/s12859-018-2581-8](https://doi.org/10.1186/s12859-018-2581-8).

[17] J. Ye *et al.*, "WEGO: A web tool for plotting GO annotations," *Nucleic Acids Res.*, vol. 34, pp. W293–W297, 2006, doi: [10.1093/nar/gkl031](https://doi.org/10.1093/nar/gkl031).

[18] S. G. Jantzen, B. J. Sutherland, D. R. Minkley, and B. F. Koop, "GO Trimming: Systematically reducing redundancy in large Gene Ontology datasets," *BMC Res. Notes*, vol. 4, Art. no. 267, 2011, doi: [10.1186/1756-0500-4-267](https://doi.org/10.1186/1756-0500-4-267).

[19] F. Supek, M. Bošnjak, N. Škunca, and T. Šmuc, "REVIGO summarizes and visualizes long lists of gene ontology terms," *PLoS One*, vol. 6, no. 7, Art. no. e21800, 2011, doi: [10.1371/journal.pone.0021800](https://doi.org/10.1371/journal.pone.0021800).

[20] C. Wan and A. A. Freitas, "An empirical evaluation of hierarchical feature selection methods for classification in bioinformatics datasets with gene ontology-based features," *Artif. Intell. Rev.*, vol. 50, pp. 201–240, 2018, doi: [10.1007/s10462-017-9541-y](https://doi.org/10.1007/s10462-017-9541-y).

[21] M. Kulmanov, M. A. Khan, and R. Hoehndorf, "DeepGO: Predicting protein functions from sequence and interactions using a deep ontology-aware classifier," *Bioinformatics*, vol. 34, no. 4, pp. 660–668, 2018, doi: [10.1093/bioinformatics/btx624](https://doi.org/10.1093/bioinformatics/btx624).

[22] S. R. Maetschke, M. Simonsen, M. J. Davis, and M. A. Ragan, "Gene Ontology-driven inference of protein-protein interactions using inducers," *Bioinformatics*, vol. 28, no. 1, pp. 69–75, 2012, doi: [10.1093/bioinformatics/btr610](https://doi.org/10.1093/bioinformatics/btr610).

[23] S. Bandyopadhyay and K. Mallick, "A new feature vector based on Gene Ontology terms for protein-protein interaction prediction," *IEEE/ACM Trans. Comput. Biol. Bioinform.*, vol. 14, no. 4, pp. 762–770, Jul./Aug. 2017, doi: [10.1109/TCBB.2016.2555304](https://doi.org/10.1109/TCBB.2016.2555304).

[24] I. Ieremie, R. M. Ewing, and M. Niranjan, "TransformerGO: Predicting protein-protein interactions by modelling the attention between sets of gene ontology terms," *Bioinformatics*, vol. 38, no. 8, pp. 2269–2277, 2022, doi: [10.1093/bioinformatics/btac104](https://doi.org/10.1093/bioinformatics/btac104).

[25] P. Terzian *et al.*, "PHROG: Families of prokaryotic virus proteins clustered using remote homology," *NAR Genom. Bioinform.*, vol. 3, no. 3, Art. no. lqab067, 2021, doi: [10.1093/nargab/lqab067](https://doi.org/10.1093/nargab/lqab067).

[26] J. Guan *et al.*, "GOPhage: Protein function annotation for bacteriophages by integrating the genomic context," *Brief. Bioinform.*, vol. 26, no. 1, Art. no. bbaf014, 2025, doi: [10.1093/bib/bbaf014](https://doi.org/10.1093/bib/bbaf014).

[27] A. Tomczak *et al.*, "Interpretation of biological experiments changes with evolution of the Gene Ontology and its annotations," *Sci. Rep.*, vol. 8, no. 1, Art. no. 5115, 2018, doi: [10.1038/s41598-018-23395-2](https://doi.org/10.1038/s41598-018-23395-2).
