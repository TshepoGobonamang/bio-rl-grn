# GRN Reliability Benchmark

[![Research stage](https://img.shields.io/badge/stage-reproducible%20prototype-orange)](#scientific-status)

Clean-room research software for **graph-conditioned, cross-fitted reliability signals in directed gene-regulatory network inference**. This repository is a **new implementation**, not a claim that independent biological validation or reinforcement learning has succeeded.

## Scientific status

- **Implemented and locally tested:** Pearson/Spearman baselines, parent-conditional predictive gain, utility-based greedy graph editing, fixed and nonlearned source-weight controls, group-cross-fitted adaptive evidence scoring, fit-only feature selection, split and artifact hashes, separate reference evaluation, source-dropout/corruption, deterministic synthetic smoke demo.
- **Not completed:** independent unseen-cohort benchmarking; verified perturbation-derived direct causal edges; true independent biological replicates/controls for GSE216595; budget-matched trained RL policy; submission-ready results/manuscript.
- **Prior knowledge:** Previous research in `TshepoGobonamang/bio-rl-grn` demonstrated mixed/negative reliability results and showed source/reference overlap risk. It is not silently imported as confirmatory evidence here.

## Environment

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\\Scripts\\activate
python -m pip install -e '.[dev]'
python -m pytest -q
```

## Generate a deterministic **toy** smoke run

```bash
grnb smoke --output results/smoke
```

This creates a generated synthetic DAG and **fake randomly assigned prior scores**. Results are only a pipeline test, **not evidence for our paper**. Outputs include `fit/fit_manifest.json`, method rankings, `fit/greedy_trajectory.csv`, `fit/expression_test_results.csv`, and `evaluation/metrics.csv`.

## Real source contract

Prepare a cell × gene CSV with a unique `cell_id` index and numeric gene columns. Supply `tf.csv` with column `tf` for regulatory candidate selection. To use real biological evidence, add a source-backed `prior.csv` with columns `source,target,score` where `score ∈ [0,1]`. Missing prior rows remain **unobserved evidence**, not false regulatory edges.

```bash
grnb fit --expression data/real/expression.csv \
  --tf-list data/real/tfs.csv \
  --prior data/real/prior.csv \
  --output results/real/fit --max-genes 12 --top-k 10 --seed 42
```

**Reference is not passed to fit.** Only after the rankings and hashes are saved:

```bash
grnb evaluate --fitted results/real/fit \
  --reference data/real/reference.csv \
  --output results/real/evaluation --reference-complete --top-k 10
```

`--reference-complete` must mean the eligible reference candidate universe is genuinely exhaustively annotated (e.g. a fully specified synthetic DAG). **Do not use this option with an incomplete ChIP/curated/perturbation-positive-only reference**, because unknown edges cannot legitimately be scored as negatives. Use descriptive positive-only overlap instead, pending a suitable evaluation design.

To audit missing evidence or corrupted source strengths, run **separately frozen scenarios** with `--missing 0.25` and `--corrupt 0.10` (fractions apply to provided prior records; retain scenario in manifest). Prespecify the scenario matrix before inspecting evaluation labels.

## Structure

```text
src/grnbench/       Model, graph search, evidence masks, evaluation and CLI
configs/            Smoke and confirmatory-protocol declarations
tests/              Unit and end-to-end no-reference-leak checks
docs/               Methodology, source requirements and completion gates
paper/              Manuscript scaffold, no fabricated results
.github/workflows/  CI and synthetic smoke artifact workflow
scripts/            Convenience run helper
```

## Why no RL yet?

A trained add/remove/stop policy only belongs in the publication after a frozen benchmark shows what needs to improve. Any RL comparison must match the greedy policy's candidate set, graph utility, utility-query budget and wall time. Until that implementation and experiment exist, **RL is an extension, not a reported result**.

## Important limitations

- Pearson/Spearman are **symmetric** baseline scores; they cannot independently establish regulatory direction.
- Parent-conditional ridge gain is an **association and predictive improvement**, not direct TF binding or a causal intervention.
- Group-cross-fitted adaptive weights are **not calibrated biological truth probabilities**.
- `test_nmse` is untouched expression error after graph selection; reference AP is evaluated only after the fit manifest seals rankings.
- Gene filtering uses fitting cells only; candidate regulators still require a documented TF whitelist for non-toy experiments.
- Tests using toy synthetic priors are engineering verification only.
- Benchmarking must use genuinely new, independent network realizations, not previously inspected evaluation topologies presented as pristine holdouts.

See [protocol](docs/EXPERIMENT_PROTOCOL.md), [data contract](docs/DATA_CONTRACT.md), and [publication status](docs/STATUS_AND_DECISIONS.md).
