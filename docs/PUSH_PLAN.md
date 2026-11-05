# Plan: push to your GitHub

**Status: run on 2026-11-06.** This records how the repo was created.

## Result

- A new **private** repo, `viti1igo/DeepPBI-KG`. That name is currently free.
- Keeps all code, notes, results and upstream's 25 commits, including their authors.
- Excludes the original data from both the files and the history.
- Upload of about 110 MB, compared with 580 MB for the current history.

Your current folder isn't changed: the work happens in a new clone, so this folder is the fallback if anything goes wrong.

## Why not just `git rm` and push

- Raw data is already in the history, so a normal removal would still upload it.
- `origin` currently points to `Tongqing-Wei/DeepPBI-KG`, which is not yours. Step 6 sets up separate remotes.

## What goes where

| Excluded (original data) | Size | How to restore |
|---|---|---|
| `phage_raw_data/`, `host_raw_data/` | 1.06 GB | `scripts/download_data.sh` |
| `all_phage_db/`, `all_host_db/` | 511 MB | `scripts/download_data.sh` |
| `code/swiss/` | 0.5 GB | `scripts/download_data.sh --with-swissprot` |

| Kept | Why |
|---|---|
| `code/`, `scripts/`, `DeepPBI-KG.py`, `.agent/`, all docs | Code and notes |
| `figure/`, `feature_similarity_report.md`, `example_output/`, `*_align_result/`, `template/` | Results |
| `model/` | Needed for prediction |
| `example_phage/`, `example_host/`, `example_*_annotation/` | Needed for the demo. Treated as demo files, not original data |
| `code/GO/` (31 MB) | Needed for the GO task and fixes the GO version |
| `graphify-out/` reports, excluding its cache | Context for new sessions |

## Steps

### 1. Install the history tool

```bash
python3 -m pip install --user git-filter-repo
```

### 2. Make a clean clone and strip the data

```bash
GIT_LFS_SKIP_SMUDGE=1 git clone --no-local /data/axnguyen/DeepPBI-KG /data/axnguyen/DeepPBI-KG-clean
```

```bash
cd /data/axnguyen/DeepPBI-KG-clean
```

```bash
python3 -m git_filter_repo --invert-paths --path phage_raw_data/ --path host_raw_data/ --path all_phage_db/ --path all_host_db/ --path code/swiss/ --force
```

### 3. Copy in the uncommitted work

```bash
cd /data/axnguyen/DeepPBI-KG && git ls-files --others --exclude-standard -z | rsync -a --from0 --files-from=- ./ /data/axnguyen/DeepPBI-KG-clean/
```

```bash
cd /data/axnguyen/DeepPBI-KG-clean
```

### 4. Remove the Git LFS rules and commit

They only cover the excluded data.

```bash
git rm -q .gitattributes
```

```bash
git add -A
```

```bash
git status --short
```

```bash
git commit -m "Add GO task notes, feature comparison and session docs" -m "Exclude original data from the repo. scripts/download_data.sh restores it."
```

### 5. Check before uploading

The first command should print `0`. The second should show no file over 100 MB; the largest should be `go-basic.obo` at 31 MB.

```bash
git rev-list --objects --all | grep -cE ' (phage_raw_data|host_raw_data|all_phage_db|all_host_db|code/swiss)/'
```

```bash
git rev-list --objects --all | git cat-file --batch-check='%(objectsize) %(rest)' | sort -n | tail -3
```

### 6. Create the repo and push

```bash
gh repo create viti1igo/DeepPBI-KG --private --source=. --remote=origin
```

```bash
git remote add upstream https://github.com/Tongqing-Wei/DeepPBI-KG.git
```

```bash
git remote set-url --push upstream DISABLED
```

```bash
git push -u origin master
```

### 7. Test it like a new session

```bash
git clone https://github.com/viti1igo/DeepPBI-KG.git /data/axnguyen/DeepPBI-KG-verify
```

```bash
cd /data/axnguyen/DeepPBI-KG-verify && scripts/download_data.sh --from /data/axnguyen/DeepPBI-KG
```

Then run the demo in `GETTING_STARTED.md` step 3 and compare the scores. Delete `DeepPBI-KG-verify` afterwards. Keep the original folder until you're happy with the new repo, then use `DeepPBI-KG-clean` as your working copy.

## Before you run it

- **Visibility:** upstream has no license, so by default others have no right to reuse its code or data. Stay private, or get permission from the author before going public.
- **Commit IDs change:** rewriting the history gives every commit a new ID, so it won't line up with upstream's commits or merge from them cleanly.
- **Your email:** commits will show `popcapgames1234abcd@gmail.com` to anyone who can read the repo. Change `git config user.email` first if you don't want that, especially if the repo may become public.
- **Report:** `feature_similarity_report.md` will be pushed with its errors. `docs/CODE_REVIEW.md` records them, and correcting the report is a good next step.
