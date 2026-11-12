# Code review

Reviewed 2026-11-06. Scope: the uncommitted WGS vs KG work, the GO/KEGG code, and the main pipeline. `code/key_gene_feature_generation.py`, `code/wgs_feature_generation.py`, `code/protein_features.py`, the key-gene helper scripts and `code/PHI_train.ipynb` were not read line by line.

**Verdict: request changes.** The prediction demo reproduces its reference results, but the WGS vs KG report must not be cited until it is corrected.

## Verified

Ran on 2026-11-06, using Python 3.9, torch 1.9.0+cpu, NumPy 1.26.4, pandas 2.3.3, scikit-learn 1.6.1 and Biopython 1.78:

- **Demo:** finished in about 2 minutes. All 72 pairs matched `example_output/result.csv` (largest score difference 2.5e-07), and key gene counts were identical.
- **Download script:** tested with `all_phage_db` only. The `.nsq` checksum matched the local copy. The full download wasn't run.

Not run: `code/train.py` and `code/GO_KEGG.R` (inputs missing), and `prokka_blast.sh`.

## Critical

### 1. The report's main conclusions don't match its own figures
`feature_similarity_report.md:39,67-68,76`

The report was written at 19:24 on 2026-06-24 and the figures at 19:42. Rerunning the script reproduces the figures, not the report.

| Claim | Report | Figure and rerun |
|---|---|---|
| Phage PC1 r | 0.9634, p < 1e-4 | **0.4203, p = 0.26** |
| Bacterium PC1 r | 0.9782 | 0.9854 |
| "Always matches" its own genome | Yes, strong diagonal | **No, see the table in item 2** |

The heatmap doesn't show the claimed "mathematical proof". In the phage heatmap, MG711466's own similarity is only 0.049 and the highest non-matching value is 0.449.

**Fix:** correct the script (item 2), regenerate the figures and rewrite sections 2-4.

### 2. Key gene features are scaled with the wrong scaler
`code/visualize_features.py:44-59,187,196`

The script scales both feature sets with `wgs_scaler.pkl`, but `DeepPBI-KG.py:1063` scales key gene features with `key_gene_scaler.pkl`. So the report's claim to match "the exact input seen by the classifier" is false. It also creates much of the "horizontal shift" that the report explains biologically: the phage WGS-KG gap on PC1 is 85.5 with the shared scaler and 0.9 with the correct ones.

| Measure | Phage, current | Phage, correct | Bacterium, current | Bacterium, correct |
|---|---|---|---|---|
| Matching-pair cosine, mean | 0.281 | 0.447 | 0.588 | 0.765 |
| Non-matching cosine, mean | 0.173 | 0.179 | 0.151 | 0.162 |
| PC1 r | 0.420 | 0.931 | 0.985 | 0.998 |
| Own genome ranks first, by WGS row | 8/9 | 7/9 | 6/8 | 7/8 |
| Own genome ranks first, by KG column | 2/9 | 8/9 | 4/8 | 8/8 |

The correct scalers support the idea better, but not "always": a few organisms still fail. With only 9 phages and 8 bacteria, this is a sanity check, not validation. To reproduce the table, run `scripts/check_feature_similarity.py`.

**Fix:** pass each feature set its own scaler, and report the sample size next to the conclusions.

## Important

### 3. UniProt filter checked the wrong column (fixed 2026-11-06)
`code/UniProt2GO_annotate.py:14-16`

It tested `lsplit[7]` but stored `lsplit[6]`. In UniProt's `idmapping_selected.tab` these are the UniRef100 and GO columns, so proteins without GO terms produced empty values, which notebook cells 4-5 then removed. Proteins with GO terms but no UniRef100 entry would have lost those terms, and rows with nothing after the GO column raised `IndexError`.

**Fixed:** it now tests `lsplit[6]` and strips only the newline. On UniProt's 1,000-row example file, output went from 1,000 genes (298 with blank GO values) to 702 genes with none. No row in that file has GO terms without UniRef100, so that loss wasn't seen there. Saved results weren't regenerated, because the inputs are missing.

### 4. GO parser drops alternate IDs and the hierarchy
`code/parse_go_obofile.py:38,46`

Line 46 replaces the table holding the alternate IDs. It also ignores `is_a` edges, which the GO task needs, and `DataFrame.append` was removed in pandas 2. The effect on existing results is small: 51 of 54,931 saved pairs use alternate IDs. **Fix:** `scripts/go_dag.py` now has a parser that keeps parents and alternate IDs. `code/parse_go_obofile.py` itself is unchanged, so `code/GO/go.tb` should move to the new parser the next time it is regenerated.

### 5. Enrichment plots don't do what they say
`code/GO_KEGG.R:26,33-38,43`

- The top-10 selection at lines 33-35 is overwritten at lines 36-38, so every significant term gets plotted.
- GO plots are labeled `-log10pvalue` but use `p.adjust`, while KEGG plots use raw `pvalue`.
- `arrange()` runs before `dplyr` is loaded.
- `ridgeplot(ego)` passes an `enricher` result to a GSEA plot.

The last two are from reading the code; the script couldn't be run.

### 6. Pipeline can't be rerun into existing folders
`DeepPBI-KG.py:1054,829-832`

`os.mkdir` raises `FileExistsError`, which also means the ReadMe command fails with the committed `template/` and `example_output/`. **Fix:** use `os.makedirs(..., exist_ok=True)`, and make sure stale files can't be mixed in.

### 7. Results can't be reproduced from this repo
`requirements.txt`, `code/GO_KEGG.R`, `code/train.py`

Training data and enrichment inputs are missing. `requirements.txt` omits `matplotlib`, `seaborn` and `scipy`, and the current server env has NumPy 2.0.2 with no Biopython or dominate, so the demo fails there. **Fix:** ask upstream for the training data, and replace `requirements.txt` with the versions in `GETTING_STARTED.md`.

### 8. No automated tests
GO DAG tests now exist: `python tests/test_go_dag.py`. **Fix:** add the demo comparison as a test.

## Suggestions

- **Duplicate pipeline:** `code/KGDeepPBI.py` closely mirrors `DeepPBI-KG.py`, and feature functions repeat across scripts. The knowledge graph's top four nodes are all copies of `protein_process()`. Keep one shared module.
- **Model selection:** `code/train.py` has no validation split and saves the last epoch only. Check how the paper's hyperparameters were chosen before relying on its test metrics.
- **Speed:** `DeepPBI-KG.py:1076-1081` searches the whole table for each pair, which grows quadratically with the number of pairs. A merge would be simpler and faster.
- **Typo:** `DeepPBI-KG.py:1057` has `0.0001,csv`. WGS runs appear not to use that path, so it looks harmless, but it should be fixed.
- **Shell quoting:** `prokka_blast.sh` loops over `ls` output and leaves `$val` unquoted, so paths with spaces break.
- **Unused code:** `StandardScaler` is imported but unused in `code/visualize_features.py:8`, and the report embeds a copy of the script that has already drifted.

## FYI

- **Security:** no secrets turned up in the scanned code and text files. `pickle.load` and `torch.load` can run code from a malicious file, so only load trusted models.
- **Earlier statement corrected:** dropped GO alternate IDs were first described as something that could change enrichment results. The measured effect is about 0.1% of pairs.
