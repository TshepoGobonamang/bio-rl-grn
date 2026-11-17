# M2B — Learn and test a reliability gate without regulatory gold labels

**Milestone:** initial adaptive evidence reliability, not a reinforcement-learning agent. This is exploratory training on three **already-inspected** BEELINE topologies. No claim of biological causal discovery, calibration to edge truth, or journal-confirmatory superiority.

## 1. Pre-existing input: verified M2A cross-fitting dataset

Pinned [GitHub Actions M2A run #44089516809](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/44089516809) contains 1,170 fold-edge rows across dyn-LL, dyn-LI and dyn-BF with *separate* fitting-only predictors (`m2_features.csv`) and validation-only pseudo-targets (`m2_targets.csv`). The M2B loader checks allowlisted features, sample-isolation manifests, exact row-key alignment, missingness and the absence of reference labels.

**Target:** `y_temporal_reproducibility_TARGET` (0–1 soft target based on out-of-fold sign agreement). It is **not** the probability of a real regulatory edge. It is symmetric for reverse edges and cannot, on its own, establish orientation. The held-out expression gain target is not used to train the gate.

## 2. Leave-one-network-out design

For each network **N**, train on the other two networks and evaluate the entire withheld **N** at once. Folds/reversed edge pairs are correlated, so random-row splitting is prohibited. Fixed feature list is the 12-fitting-only columns in the M2A manifest, with explicitly masked missing evidence.

The adaptive model is intentionally minimal: constant-value imputation, z-score standardization and ridge regression with **alpha=30**, clipped to `[0,1]`. Every transform is fitted only on the two training network structures. Each training topology receives equal total sample weight to avoid letting the 306-edge long-linear network overwhelm the two 42-edge networks. Hyperparameter `alpha=30` is fixed in advance; no held-out topology used for selection.

Four comparable reliability weights:

- **fixed**: reliability 1 for observed evidence, 0 if unavailable.
- **global**: single source-wide reliability fitted as weighted mean of the **other networks'** soft targets.
- **stability**: direct, nonlearned within-fitting-fold bootstrap sign agreement.
- **adaptive**: cross-network regularized model using the allowed fitting-only features.

All weights are in `[0,1]` and exactly 0 if evidence is missing. The M2A feature `z_signed` is not divided by the weight in the graph utility; single-source reliability must genuinely attenuate temporal evidence.

## 3. Stage-one surrogate validation

Export `m2b_lono_metrics.csv`: pseudo-target MSE, MAE and mean gate value for each held-out network and method. Also export `m2b_meta_evaluation_only.csv` for diagnosis; it **contains** the held-out target and must never be used as inference input.

Export `m2b_inference_weights.csv`: **only** `dataset_id,source,target,method,rho`, aggregated over inner folds and free of any pseudo-target column. This is the only weight file used by GRN inference.

Do not interpret MSE improvement over fixed reliability=1 as proof of a useful learned gate. The primary sanity comparison is **adaptive versus learned global and nonlearned stability**. With three networks, any reported standard error or p-value treating fold-edges as independent is invalid.

## 4. Stage-two network reconstruction

Using the official BEELINE v4 expression and pseudotime, compare six methods on each held-out topology:

1. Expression-only graph-conditioned greedy.
2. Fixed temporal `rho=1` (the M1.2 baseline).
3. Globally weighted temporal greedy.
4. Nonlearned bootstrap-stability weighted temporal greedy.
5. **Adaptive** temporally weighted greedy using frozen leave-one-network-out gate.
6. Shuffled fit-cell pseudotime negative control.

All use the same independent 300/100/100 fit/graph-scoring/untouched-test split, candidate TF universe, ridge coefficient, **20-edit cap** and temporal contribution multiplier **0.015**. The frozen gate `rho_e` is computed from M2A fitting-fold-only features and multiplied with the signed temporal margin recomputed on the 300 outer fitting cells. Selected graphs are evaluated on previously untouched expression test cells before loading reference edges. Gold-standard edges are accessed **only** for final AUROC/AP/Precision@20/direction-concordance assessment.

**Limitations:** the pseudotime evidence from BoolODE may encode near-simulation time, gene direction does not follow causally from binned lagged association, and the same network topologies were studied earlier. These results are **exploratory**, not external biological validation.

## 5. Execution

After merge, `.github/workflows/m2b-reliability.yml` automatically runs once (workflow added to `main`) and can be dispatched manually thereafter. It retrieves the *pinned* completed M2A artifact, audits it, trains frozen leave-one-network-out gates, downloads official BEELINE v4 data on the runner and evaluates the six matched methods. No raw expression data are committed.

Local meta-learning from downloaded artifact:

```bash
python -m pip install -e ".[dev]"
python scripts/audit_m2_crossfit.py --root data/m2a
python scripts/run_m2b.py --m2a-root data/m2a --output-dir results/m2b
```

The output `m2b_manifest.json` lists the training feature schema and caveats.

## 6. Acceptance criteria for the next stage

- Does the adaptive model improve *held-out network* temporal reproducibility prediction compared with global and nonlearned stability weighting?
- Does its graph edit policy show higher gold-reference AP/direction concordance than equally informed fixed-temporal greedy, **without a material loss** in untouched expression prediction?
- Does shuffling pseudotime remove the apparent gain?
- Does performance stay stable across independent, **unexamined** topology realizations? (Not answered by these pilot datasets.)
- Before learning graph-conditioned reliability, verify that editing the graph changes the **full frozen-parameter utility**, not just newly added edge reward.

If improvements vanish under stronger comparisons, report the negative result and do not build an RL success claim upon them.
