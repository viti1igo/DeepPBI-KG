# DeepPBI-KG

Phage-bacterium interaction prediction, plus an open Gene Ontology (GO) hierarchy task.

Start with [GETTING_STARTED.md](GETTING_STARTED.md) for setup and data, and [docs/GO_TASK_STATUS.md](docs/GO_TASK_STATUS.md) for the current task. Known bugs are in [docs/CODE_REVIEW.md](docs/CODE_REVIEW.md). If `graphify-out/graph.json` is present, try `graphify query "<question>"` for questions about how the code connects.

## Commands

- Regression check: the demo in `GETTING_STARTED.md` step 3. Always use new `--output` and `--template` folders. There is no test suite.
- Restore data: `scripts/download_data.sh`
- Environment: Python 3.9, `numpy<2`, `biopython<1.80`. Don't install from `requirements.txt`.

## Boundaries

- Never commit `phage_raw_data/`, `host_raw_data/`, `all_phage_db/`, `all_host_db/` or `code/swiss/`.
- Don't overwrite committed results in `template/`, `example_output/` or `figure/` unless asked.
- `origin` should be the user's repo. Never push to `Tongqing-Wei/DeepPBI-KG`.
- Don't make the repo public without the user's decision: upstream has no license.
- Make pipeline fixes in `DeepPBI-KG.py`. `code/KGDeepPBI.py` is an older near-copy.
- Treat `feature_similarity_report.md` as out of date until it is regenerated.
- `.agent/rules/AGENT.md` is a finished plan for earlier work, not a current instruction.
