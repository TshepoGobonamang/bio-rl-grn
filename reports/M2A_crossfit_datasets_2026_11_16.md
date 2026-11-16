# M2A official BEELINE cross-fit data audit — 2026-11-16

**Execution:** [GitHub Actions run #43906547502](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/43906547502), artifact `m2a-crossfit-reliability-dataset` (ID **14245383093**). Implemented through merged [PR #10](https://github.com/TshepoGobonamang/bio-rl-grn/pull/10), passing tests on Python 3.10, 3.11 and 3.12.

**Status:** Three actual **gold-blind meta-training datasets**, NOT trained reliability weights, GRN accuracy results, or an RL policy. Only BEELINE expression and pseudotime files are read.

## Generated rows and crossfit verification

| Topology | Genes | Directed candidates | Inner folds | Edge-fold examples | Temporal missing | Valid temporal targets |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| dyn-LL-500-1 | 18 | 306 | 3 | **918** | 0 | 918 |
| dyn-LI-500-1 | 7 | 42 | 3 | **126** | 0 | 126 |
| dyn-BF-500-1 | 7 | 42 | 3 | **126** | 0 | 126 |
| **Total** | | 390 distinct candidates across datasets | | **1,170** | 0 | 1,170 |

Each dataset: **500 original simulated cells**, with 300 outer fitting, 100 graph-scoring and 100 untouched expression-test cells. Each inner fold uses 200 cells for fitting-only features and ridge fits, and 100 different cells to construct temporal/conditional-expression pseudo-targets. BF uses two branch-specific pseudotime sources; LL and LI use one.

## Observed target distribution

| Topology | Mean temporal reproducibility target | Percentage with temporal target below .02 or above .98 | Mean held-out conditional expression gain | Fraction positive expression gain | Pearson correlation between two pseudo-targets |
| --- | ---: | ---: | ---: | ---: | ---: |
| dyn-LL | .559 | 50.8% | .0315 | 65.1% | .128 |
| dyn-LI | .577 | 54.0% | .0410 | 88.9% | .275 |
| dyn-BF | .574 | 66.7% | .0932 | 94.4% | −.040 |

The temporal reproducibility target measures the **sign stability** of signed lagged-pseudotime evidence across folds; it is not a true-edge label. The conditional expression-gain target measures held-out prediction improvement with a **source-blind alternative TF scaffold** and may be negative. Neither target is causal.

No missing temporal evidence occurs **in these three clean datasets**; missing temporal branches are deliberately represented by `m_present=0`, `z_signed=NaN` and an absent temporal target in regression tests. A missing feature is not a negative regulatory label.

## Validated safeguards

- Correct 3×candidate-universe count; unique `(fold_id,source,target)`.
- Disjoint inner train/validation cells; outer graph-scoring/test cells never used for meta-features/targets.
- No access to `GroundTruthNetwork.csv`, no gold-derived TF labels.
- Exactly 12 allowed feature columns in JSON manifest; validation columns marked `TARGET`/`TARGET_ONLY` and expressly not permitted as model inputs.
- Reverse-edge signed pseudotime scores are antisymmetric, whereas the temporal consistency target is symmetric under direction reversal (by design).
- Original expression/pseudotime hashes and hashed fold cell IDs recorded.

## Before any reliability-model fit

1. **Temperature saturation:** 51–67% of temporal labels near extreme 0 or 1 at predeclared temperature .05; run development-only calibration and sensitivity checks.
2. **Dependence:** do not randomly split edge-fold rows and infer generalization; reverse-pair rows share labels, and folds reuse cells. Use pair/group and independent network-level validation.
3. **Weak surrogate alignment:** temporal consistency vs expression utility correlations range −.040 to .275; do not combine them as if they were the same kind of regulatory truth.
4. **Small-network expression gains:** 89–94% of conditional gains positive in two 7-gene models; assess target discrimination and alternative scaffolds.
5. **Future evidence:** pseudotime is one source and synthetic simulator time may be optimistic; genuine biological motif, ATAC or curated sources require independent provenance and gold-overlap checks.

Complete original CSV and JSON tables remain in the GitHub workflow's artifact. This repository intentionally does not commit raw expression data or potentially large derived training tables.

**Next milestone:** calibrate and audit these pseudo-targets, then compare a small constrained reliability gate against equally informed fixed and nonlearned stability baselines on held-out network contexts.
