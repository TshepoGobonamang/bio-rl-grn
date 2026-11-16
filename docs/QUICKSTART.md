# Bio-RL-GRN E1 quick start

This is an **exploratory expression-only baseline**, not a trained Bio-RL-GRN reinforcement learning policy and not a journal result.

## Run tests and a toy software smoke check

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
biorlgrn demo --output data/demo
biorlgrn run --data-dir data/demo --output-dir results/demo --max-edges 4
biorlgrn evaluate --data-dir data/demo --results-dir results/demo --k 3
```

## Real BEELINE dyn-LL-500 replicate 1

Official v4 Zenodo archive: https://zenodo.org/records/19009603 .
Download BEELINE-data.zip (~262 MB) into your local checkout or Downloads.

```bash
python scripts/extract_beeline_v4.py --archive BEELINE-data.zip --run dyn-LL-500-1 --output data/raw/beeline/dyn-LL-500-1
python scripts/run_e1.py --raw-dir data/raw/beeline/dyn-LL-500-1 --processed-dir data/processed/dyn_LL_500_1 --results-dir results/dyn_LL_500_1 --max-edges 19
```

The extraction checks the official v4 archive MD5. The runner converts genes-by-cells into cells-by-genes, fits expression-only Pearson, Spearman and graph-dependent greedy-ridge baselines, writes complete ranked edge lists, **then** reads the gold-reference GRN for final scoring.

To provide an independent TF list use the optional `--tf-list` flag. Without one, **all measured genes** are potential regulators (self loops excluded); this is a research limitation, not a data-driven TF discovery mechanism.

For GitHub-hosted execution, use **Actions → E1 BEELINE dyn-LL-500 exploratory benchmark → Run workflow** after the feature PR is merged. The workflow downloads the official archive on a temporary runner, does not commit original data, and uploads ranked predictions and result manifests. It is manual to avoid automatic 262 MB downloads on every change.

## Interpretation

- AUPRC/AUROC use a closed-world evaluation of the BEELINE simulated reference.
- The ridge optimizer repeatedly queries an internal scoring split: it may overfit and is not yet independently validated.
- Changes in expression reconstruction are predictive associations, **not** causal proof.
- The reliability module, RL agent, repeated-seed controls and independent biological validation are future milestones.

See [ROADMAP.md](ROADMAP.md) and [EXPERIMENT_PROTOCOL.md](EXPERIMENT_PROTOCOL.md).
