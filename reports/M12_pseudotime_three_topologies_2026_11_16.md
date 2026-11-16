# M1.2 pseudotime-informed GRN orientation — empirical pilot (2026-11-16)

**Evidence:** [Official successful GitHub Actions run #43901680886](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/43901680886), artifact `m12-pseudotime-direction-pilot`, ID `14242545360`. Dataset: BEELINE v4 BoolODE, `dyn-LL-500-1` (long linear), `dyn-LI-500-1` (short linear), `dyn-BF-500-1` (bifurcating).

## Method / safeguards

Each topology contributes **one** 500-cell simulated dataset: 300 fit, 100 iterative graph-scoring and 100 untouched expression-test cells. Same candidate universe for all methods, self-loops excluded, no gold-label access until inference and independent expression testing are complete. Ridge alpha 1.0, temporal lag 1, 16 pseudotime bins per branch, pre-specified temporal weight **0.015**, 20 maximum edits.

Comparators: original graph-conditioned greedy; greedy + fitting-only lagged pseudotime evidence; greedy + **shuffled** fitting-cell pseudotime; and pseudotime-only ranking. These are **not** trained RL or adaptive-prior models, nor causal regulatory inference.

## Observed accuracy

| Topology | Method | AUROC | Average precision | TP@20 | Direction concordance | Reciprocal pairs@20 | Untouched expression loss ↓ |
|---|---|---:|---:|---:|---:|---:|---:|
| dyn-LL (18 eligible gold; 306 candidates) | Greedy | **.7809** | .4864 | 9 | .500 | 7 | **.8767** |
|  | Temporal greedy | .7741 | **.6566** | **12** | **.722** | **2** | .9365 |
|  | Shuffled time | .7730 | .5210 | 10 | .556 | 7 | .8843 |
|  | Temporal only | .7384 | .2886 | 7 | .722 | 0 | 1.0487 |
| dyn-LI (7 gold; 42 candidates) | Greedy | .7959 | .5801 | 5 | .429 | 8 | **.6861** |
|  | Temporal greedy | .8000 | **.6199** | 5 | .571 | 5 | .6892 |
|  | Shuffled time | **.8286** | .5587 | **7** | **.714** | 6 | .6879 |
|  | Temporal only | .7102 | .4561 | 5 | .714 | 0 | .7828 |
| dyn-BF (10 gold; 42 candidates; 2 branches) | Greedy | .5375 | .2918 | 6 | .375 | 9 | **.5696** |
|  | Temporal greedy | .6375 | .3385 | 7 | .750 | 8 | .5776 |
|  | Shuffled time | .5750 | .3080 | 6 | .625 | 9 | .5802 |
|  | Temporal only | **.8250** | **.7427** | **8** | **.875** | 0 | .7647 |

All 12 AUROC / noninterpolated AP values were **independently recomputed from evaluation-labeled edge-list CSVs** and match the workflow metrics.

## Results and interpretation

1. Temporal greedy improves AP over expression-only greedy on **all three** topologies: +35.0% relative (dyn-LL), +6.9% (dyn-LI), +16.0% (dyn-BF). Direction concordance improves from .500→.722, .429→.571 and .375→.750, respectively. Only one run per topology; **no statistical superiority claims**.
2. The mean expression-reconstruction trade-off fails the strict dual-success criterion: untouched test loss **increases** by 6.82%, 0.46%, and 1.41% on dyn-LL, dyn-LI and dyn-BF, respectively.
3. Shuffled time is not uniformly inferior in all metrics: e.g. dyn-LI directional concordance .714 exceeds actual temporal greedy .571. A few positive orientation comparisons are insufficient.
4. Temporal-only ranking is highly effective for the bifurcating simulation (AP .7427), but has poor untouched expression loss (.7647), and is worse in long-linear simulation (AP .2886). This is strong motivation to **test**, not assert, graph/context-dependent temporal evidence reliability.
5. Zero reciprocal edges in the top 20 for temporal-only follows in part mechanically from antisymmetric direction scores and near-half candidate selection; it does **not** validate causal orientation.

## Provenance and limitations

- dyn-LL: 19 raw gold edges, 18 eligible, 1 excluded self-loop; 18 genes and 1 pseudotime branch.
- dyn-LI: 8 raw, 7 eligible, 1 excluded self-loop; 7 genes, 1 branch.
- dyn-BF: 12 raw, 10 eligible, 2 excluded self-loops; 7 genes, **2 branches**.
- Pseudotime in BoolODE simulations may be simulation-time proxy, an **optimistic** directional-information scenario. Ordered scRNA-seq snapshots lack true within-cell temporal coupling (Scribe, Qiu et al., Cell Systems 2020).
- Cross-topology evidence consists of one expression realization per network family. Top-20 metrics cannot be compared across 42 and 306 candidate-edge universes without noting differing prevalence.
- Lagged predictive-gain targets use 16 aggregated bins and a fitting-only heuristic that needs further cross-fitting and sensitivity analyses. Exact temporal integration is **not** yet part of the frozen journal method.

## Next gate

Retain the expression-only greedy and pseudotime-greedy implementations. Do **not** declare pseudotime a solved causal direction mechanism. Next test (i) fitting-fold cross-fitted temporal evidence, (ii) independently estimated rather than simulator-time ordering, (iii) multiple seeds/network realizations, (iv) learned temporal-evidence reliability with no/fixed/shuffled controls and independent evaluation, and (v) matched compute budgets.

Research references: [BEELINE v1.1](https://murali-group.github.io/Beeline/reproducing-results.html), [BoolODE](https://github.com/Murali-group/BoolODE), [Scribe](https://doi.org/10.1016/j.cels.2020.02.003).
