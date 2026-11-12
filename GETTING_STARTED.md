# Getting started

DeepPBI-KG predicts phage-bacterium interactions with a deep neural network, using features from key genes (KG) and whole genomes (WGS). This copy adds two lines of work to the upstream project: a WGS vs KG feature comparison, and an open task on the Gene Ontology (GO) hierarchy.

## Read these first

| Read | For |
|---|---|
| [docs/GO_TASK_STATUS.md](docs/GO_TASK_STATUS.md) | The original task, what is done, and what to do next |
| [docs/CODE_REVIEW.md](docs/CODE_REVIEW.md) | Known bugs, and why `feature_similarity_report.md` can't be trusted yet |
| [scripts/README.md](scripts/README.md) | How to run the GO scripts, and which data and GO releases they use |
| [ReadMe.md](ReadMe.md) | Upstream install and usage guide |

Then check the real state with `git status` and `git log --oneline -5`.

## 1. Set up Python

Use Python 3.9. These versions ran the demo correctly on Linux:

```bash
conda create -n DeepPBI-KG python=3.9 -y
```

```bash
conda activate DeepPBI-KG
```

```bash
pip install torch==1.9.0+cpu -f https://download.pytorch.org/whl/torch_stable.html
```

```bash
pip install "numpy==1.26.4" "pandas==2.3.3" "scikit-learn==1.6.1" "scipy==1.13.1" "biopython==1.78" "dominate==2.9.1" "matplotlib==3.9.4" "seaborn==0.13.2"
```

Don't install from `requirements.txt`: it pins old versions and leaves out `matplotlib`, `seaborn` and `scipy`. NumPy must stay below 2, because PyTorch 1.9 fails with `RuntimeError: Numpy is not available` on NumPy 2. Biopython must stay below 1.80, because the code imports `GC` and `CodonUsage`.

The existing `DeepPBI-KG` conda env on the server needs a fix first: it has NumPy 2.0.2 and no Biopython or dominate.

`prokka` and `blast` are only needed to annotate new genomes (step 4). R and `clusterProfiler` are only needed for `code/GO_KEGG.R`.

## 2. Download the dataset

The original data is not stored in this repo. `git-lfs` must be installed. The script downloads a fixed upstream commit, skips folders you already have, and checks file counts.

```bash
scripts/download_data.sh
```

On the server, copy from the existing checkout instead:

```bash
scripts/download_data.sh --from /data/axnguyen/DeepPBI-KG
```

| Folder | Contents | Files | Size | Needed for |
|---|---|---|---|---|
| `phage_raw_data/` | Phage genomes (FASTA) | 3,513 | 375 MB | Prediction reports, rebuilding the BLAST DB |
| `host_raw_data/` | Bacterial genomes (FASTA) | 395 | 685 MB | Prediction reports, rebuilding the BLAST DB |
| `all_phage_db/` | Phage BLAST database | 8 | 65 MB | Aligning new phages |
| `all_host_db/` | Bacterial BLAST database | 8 | 446 MB | Aligning new bacteria |
| `code/swiss/` | SwissProt FASTA and DIAMOND DB | 2 | 0.5 GB | GO annotation only |

`code/swiss/` is optional:

```bash
scripts/download_data.sh --with-swissprot
```

The BLAST databases and SwissProt use upstream's Git LFS, so that download can fail if upstream runs out of LFS quota. You can rebuild the BLAST databases from the FASTA files if needed:

```bash
python integrate_seq.py
```

```bash
makeblastdb -dbtype nucl -in all_phage_seq.fasta -parse_seqids -out ./all_phage_db/all_phage_seq.db
```

```bash
makeblastdb -dbtype nucl -in all_host_seq.fasta -parse_seqids -out ./all_host_db/all_host_seq.db
```

**Not available anywhere we checked:**
- Training sets, such as `key_gene_train_random_set.xlsx`, so `code/train.py` can't be rerun.
- Enrichment inputs used by `code/GO_KEGG.R`: `my_degs*.csv`, `*_go_annotation.txt`, `*_kegg_annotation.txt`, `ko00001.json`, `*.gmt` and `*key_gene_importance.csv`.
- UniProt's `idmapping` file, used by `code/UniProt2GO_annotate.py`. It is a public UniProt download.

## 3. Check that it works

The demo takes about 2 minutes on CPU. Use new `--output` and `--template` folders every time, because the script stops with `FileExistsError` if its subfolders already exist. Don't point it at the committed `template/` or `example_output/`.

```bash
mkdir run_out run_tmp
```

```bash
python DeepPBI-KG.py --phage_annotation ./example_phage_annotation --bacterium_annotation ./example_host_annotation --output ./run_out --model ./model --template ./run_tmp --phage_raw_data ./phage_raw_data --bacterium_raw_data ./host_raw_data --phage_align_res ./phage_align_result --bacterium_align_res ./host_align_result
```

It worked if it prints `predict complete!` and writes `run_out/result.html`. All 72 scores in `run_out/result.csv` should match `example_output/result.csv` to within 0.000001.

## 4. Common commands

Annotate and align new genomes (needs `prokka`, `blast` and the BLAST databases). The output folders must already exist:

```bash
bash prokka_blast.sh -i ./my_phage -o ./my_phage_annotation -p "$(which prokka)" -d ./all_phage_db/all_phage_seq.db -b "$(which blastn)" -O ./my_phage_align
```

Recheck the WGS vs KG numbers:

```bash
python scripts/check_feature_similarity.py
```

Move saved GO annotations up the hierarchy to a GO slim and a chosen level, and measure where they sit (standard library only; `code/GO/go-basic.obo` is the 2022 release, so use a current one for new results):

```bash
python scripts/go_rollup.py code/GO/go-basic.obo code/富集分析data_process.ipynb go_rollup_out --genes host_key=model/host_key_gene_0.00004.csv --genes phage_key=model/phage_key_gene_0.0001.csv
```

Extract GO terms for any proteins, or from eggNOG-mapper or InterProScan output, and summarise them. See [scripts/README.md](scripts/README.md) for downloads and options:

```bash
python scripts/go_extract.py run proteins.faa out --db go_data/sprot_2026_03_db --sprot go_data/uniprot_sprot_2026_03.dat.gz --obo go_data/go-basic.obo
```

## 5. Where things are

| Path | Contents |
|---|---|
| `DeepPBI-KG.py` | Main prediction pipeline |
| `code/KGDeepPBI.py` | Older near-copy of the pipeline; edit `DeepPBI-KG.py` instead |
| `code/train.py`, `code/predict.py`, `code/RF.py` | Training, standalone prediction, random forest baseline |
| `code/GO_KEGG.R`, `code/富集分析data_process.ipynb` | GO/KEGG enrichment and its data preparation |
| `code/parse_go_obofile.py`, `code/UniProt2GO_annotate.py`, `code/kegg_json.py` | GO and KEGG file conversion |
| `code/GO/` | GO release 2022-03-22 (`go-basic.obo`, `go.tb`) |
| `code/visualize_features.py`, `feature_similarity_report.md` | WGS vs KG comparison |
| `model/` | Pretrained models, scalers, key gene lists |
| `template/`, `example_*`, `*_align_result/` | Saved intermediate files, demo inputs and demo results |
| `scripts/` | GO scripts (see `scripts/README.md`), data download and recheck scripts |

## 6. Running on a Mac

Copying the project is fine, but prediction always uses the CPU and `code/train.py` only looks for CUDA. Using the Mac GPU needs code changes and PyTorch 1.12 or later, and the versions above haven't been tested on a Mac.

## 7. Things to know

- **`feature_similarity_report.md` is out of date.** It disagrees with its own figures, so use the numbers in `docs/CODE_REVIEW.md`.
- **Upstream has no license.** Keep this repo private unless upstream gives permission.
- **Paths starting `.../`** in the R script and notebooks are placeholders for files we don't have.
- **Tests cover the GO code only:** run `python tests/test_go_dag.py` and `python tests/test_go_extract.py`. Step 3 is the only regression check for prediction.
