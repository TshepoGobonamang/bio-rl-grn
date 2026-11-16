# E1 prototype protocol (pre-registration draft)

**Goal:** validate data import, graph-conditioned ridge loss and static/greedy controls on BEELINE dyn-LL-500 replicate 1. Not yet a confirmatory experiment.

## Inputs

- BEELINE v4, dyn-LL-500-1 ExpressionData.csv, GroundTruthNetwork.csv, optional pseudotime.
- Unless an independent list is supplied, all measured genes are eligible regulators; disclose the all-gene assumption explicitly.
- Gold standard is never used in fit, model selection, prior construction, TF candidate filtering or reward.

## Split and modeling

- Deterministic shuffle: 75% fit, 25% internal graph scoring.
- Each target fits ridge regression from its selected parent set, with alpha=lambda_r*n_fit and an unpenalized intercept.
- Mean target score is held-out squared loss divided by fit-partition target variance.
- Marginal utility gain is reduction in graph loss minus normalized sparsity cost; greedy accepts maximum positive-gain addition.
- Pearson/Spearman score fitting-only absolute correlation.

## Evaluation and limits

- Rank every allowed directed candidate; only a separate evaluation operation reads gold.
- Report AUROC, AUPRC, prevalence baseline, Precision@K and Recall@K.
- This small, single-replicate study is exploratory, not evidence for RL superiority or biological validation.
- Repeatedly querying the internal graph-scoring partition may overfit; add independent untouched expression context for confirmatory analysis.

## Run official dataset

Official data: https://zenodo.org/records/19009603 (BEELINE-data.zip; approximately 262 MB).
The downloader/extractor checks publisher-provided MD5.

```bash
python scripts/extract_beeline_v4.py --archive BEELINE-data.zip --run dyn-LL-500-1 --output data/raw/beeline/dyn-LL-500-1
python scripts/run_e1.py --raw-dir data/raw/beeline/dyn-LL-500-1 --processed-dir data/processed/dyn_LL_500_1 --results-dir results/dyn_LL_500_1
```

Alternatively use GitHub Actions > E1 BEELINE dyn-LL-500 exploratory benchmark > Run workflow after merging the feature branch. The manual runner downloads data temporarily, runs models and uploads only result artifacts, not redistributed expression data.
