# M1.1: First three BEELINE replicate results (2026-11-16)

**Grounded source:** [Successful GitHub Actions run #43898506683](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/43898506683), output artifact **m11-beeline-direction-holdout**, artifact ID `14239184402`. Official BEELINE v4 (`dyn-LL-500-1`, `-2`, `-3`).

**Status:** Exploratory, three expression replicates of the **same** regulatory topology. The source code contains **no trained RL agent** and **no adaptive biological prior-reliability learner**. Avoid attributing any observed gain to either mechanism.

## 1. Fixed experimental conditions

Per replicate: 500 simulated cells, 18 genes, 18 eligible gene regulators (all measured genes are eligible), 306 directed candidates after exclusion of self-regulation. Split into **300 fitting**, **100 internal graph-edit scoring**, and **100 untouched expression-test** cells. Ridge alpha=1.0, greedy edit budget=20, fixed reciprocal-edge penalty=0.01 for the experimental ablation, split seed 42.

**Reference audit resolved:** Each raw gold reference contains 19 interactions, exactly 18 eligible edges and **one self-loop `g18 → g18`** excluded by the candidate-universe rule. The source reference audit is archived with the run results.

## 2. Mean metrics across 3 expression replicates

| Method | AUROC | Average precision (AP) | AP SD | Precision@20 | Untouched-test loss ↓ | Reciprocal pairs among top 20 | Direction concordance |
|---|---:|---:|---:|---:|---:|---:|---:|
| Pearson | 0.8498 | 0.3852 | 0.0420 | **0.5000** | 0.6107 | 10.00 | 0.500 |
| Spearman | 0.7585 | 0.2647 | 0.0431 | 0.3333 | 0.7270 | 10.00 | 0.500 |
| **Graph-conditioned greedy** | **0.8906** | **0.5196** | 0.0523 | 0.4833 | **0.5806** | 5.67 | 0.519 |
| Reverse-edge penalty ablation | 0.6047 | 0.4171 | 0.0646 | 0.4500 | 0.6015 | **0.00** | 0.537 |

Average precision is `sklearn.metrics.average_precision_score` (noninterpolated AP). Trapezoidal PR-AUC is also in `m11_metrics.csv` and is a distinct measure. Random-ranking prevalence is `18/306 ≈ 0.0588`.

## 3. Per-replicate AP

| Replicate | Pearson | Spearman | Greedy | Reverse penalty |
|---|---:|---:|---:|---:|
| 1 | 0.3371 | 0.2374 | 0.4864 | 0.4917 |
| 2 | 0.4041 | 0.2424 | 0.4925 | 0.3822 |
| 3 | 0.4145 | 0.3144 | 0.5799 | 0.3775 |

Greedy has **34.9% larger mean AP** than Pearson ((0.519575 / 0.385210) - 1). This is descriptive, not a significance claim. **Top-20 true positives**: greedy (9, 10, 10) versus Pearson (10, 10, 10): the ranking advantage does **not** translate into a better top-20 recovery count.

## 4. Untouched expression-test loss

| Replicate | Null graph | Pearson | Spearman | Greedy | Reverse penalty |
|---|---:|---:|---:|---:|---:|
| 1 | 1.2618 | 0.8857 | 1.0119 | **0.8767** | 0.9252 |
| 2 | 0.9190 | 0.5080 | 0.6504 | **0.4736** | 0.4793 |
| 3 | 0.9034 | 0.4384 | 0.5186 | **0.3916** | 0.4000 |
| **Mean** | — | 0.6107 | 0.7270 | **0.5806** | 0.6015 |

The fitted ridge models underlying selected graphs are assessed only *after* graph selection on the third untouched cell subset. Greedy shows 4.9% lower **mean** untouched-test loss than Pearson, with improvement in all three replicates. This does not yet establish cross-topology or biological-context generalization.

## 5. Directionality is the key unresolved limitation

Pearson and Spearman are exactly symmetric under swapping source and target, producing **18 direction ties** each replicate. Greedy breaks these ties, but orientation accuracy is near chance: (0.500, 0.556, 0.500), mean **0.519**. Top-20 reciprocal pair counts are (7, 5, 5) for greedy versus 10 in every Pearson/Spearman replicate.

Applying a fixed reverse-edge penalty reduces reciprocal pairs to **zero** but simultaneously reduces mean AP from 0.5196 to 0.4171 and mean AUROC from 0.8906 to 0.6047. **Conclusion:** suppressing reciprocal edges by itself does not establish correct direction or benefit overall recovery. Experimental pseudotime or perturbation support should be evaluated independently; feedback edges can be genuine.

## 6. Interpretation and next acceptance gates

**Promising:** repeated improvement in AP and untouched expression prediction supports graph-conditioned regression as a useful engineering baseline.

**Unproven:** causal TF direction, superiority to well-matched graph-search alternatives, adaptive prior reliability, and any RL benefit.

Before moving to a journal claim:
- Examine additional **independent network topologies**, not just three replicates of one dyn-LL graph.
- Add temporally informed and/or experimentally supported regulatory orientation, with a strict provenance audit to ensure no evaluation-label leakage.
- Pre-specify whether expression test, early precision, or AP is the primary outcome; report confidence intervals at the network/replicate level.
- Evaluate sensitivity to edit budget, ridge coefficient, graph sparsity and independent seeds.
- Compare adaptive biological-prior weighting with **fixed and no-prior** baselines, then learned RL with same-utility greedy/local-search controls under compute-matched budgets.

The downloaded official artifact is the reproducibility record for full annotated rankings, 20-step trajectories, exact per-replicate metrics and gold-row exclusion audit. See [Issue #4](https://github.com/TshepoGobonamang/bio-rl-grn/issues/4) for follow-up tasks.
