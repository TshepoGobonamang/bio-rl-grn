# M2B empirical decision — adaptive reliability versus fixed and nonlearned priors

**Date:** 2026-11-17  
**Decision:** **NO-GO for a claim that learned edge-adaptive temporal reliability is superior**; **GO for M2C controlled redevelopment and independent evaluation**. Retain M2B as a transparent negative/mixed baseline, not a failed software implementation.

## Reproducibility and scope

- Official successful [GitHub Actions run #44091373596](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/44091373596)
- GitHub artifact `m2b-leave-one-network-reliability`, ID **14331543589**.
- Source commit `4276cfaeb9af5199bf0df094c38d5b5131168d65` (PR #12).
- Input: 1,170 cross-fitted gold-blind M2A records from [run #44089516809](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/44089516809).
- BEELINE v4 dyn-LL-500-1, dyn-LI-500-1, dyn-BF-500-1; **one sampled expression realization per topology**.
- Each experiment uses the same 300 fitting / 100 graph-scoring / 100 untouched expression-testing cells, all-measured-gene candidate regulator universe, maximum 20 edits, ridge alpha 1.0, fixed temporal utility weight 0.015, seed 42.
- M2B meta-model: regularized Ridge (alpha 30), imputation and standardized inputs; leave-one-network-topology-out training, equal weight per training network, bounded reliability output; features exclude validation targets and gold labels.
- **Independent verification**: re-computed AUROC and non-interpolated average precision from **all 18** method×network labeled ranking files; each agrees with saved metrics to numerical precision (maximum absolute error about 5.6×10⁻¹⁷).

## 1. Network recovery (average precision, AP; higher is better)

| Method | dyn-LL (18 gold / 306 candidates) | dyn-LI (7 / 42) | dyn-BF (10 / 42) |
| --- | ---: | ---: | ---: |
| Expression-only | 0.486373 | 0.580080 | 0.291758 |
| Fixed temporal rho=1 | **0.656639** | 0.619913 | **0.338515** |
| Global learned rho | 0.629451 | 0.616883 | 0.326896 |
| Nonlearned stability gate | 0.628982 | 0.613251 | 0.333038 |
| **Adaptive learned gate** | 0.634249 | **0.622294** | 0.327954 |
| Shuffled temporal control | 0.521048 | 0.558683 | 0.308020 |

Adaptive versus fixed temporal AP: **−0.022390** (LL), **+0.002381** (LI), **−0.010561** (BF). Fixed wins on 2/3; adaptive on 1/3 by a very small margin. An unweighted across-topology mean is 0.528166 for adaptive and 0.538356 for fixed, **but is descriptive only** because prevalence and candidate universes differ.

Top-20 selected true edges: LL **11 adaptive vs 12 fixed**; LI **5 vs 5**; BF **7 vs 7**. The top-20 overlap between adaptive and fixed is **17/20, 19/20, 20/20**, respectively. Thus the adaptive model generally yields only small graph changes compared with static temporal weighting.

## 2. Orientation and reciprocal edges

| Network | Fixed direction concordance | Adaptive concordance | Fixed reciprocal pairs@20 | Adaptive reciprocal pairs@20 |
| --- | ---: | ---: | ---: | ---: |
| dyn-LL | **0.722222** | 0.666667 | **2** | 4 |
| dyn-LI | 0.571429 | 0.571429 | **5** | 6 |
| dyn-BF | 0.750000 | 0.750000 | 8 | 8 |

The learned reliability gate has **not improved direction concordance** on any of these networks. Pseudotime-based orientation remains an association heuristic, not evidence of causal TF regulation.

## 3. Untouched expression-test loss (lower is better)

| Network | Expression-only | Fixed temporal | Global | Stability | Adaptive |
| --- | ---: | ---: | ---: | ---: | ---: |
| dyn-LL | **0.876705** | 0.936497 | 0.894116 | 0.897279 | 0.897279 |
| dyn-LI | **0.686084** | 0.689238 | 0.688283 | 0.689238 | 0.688283 |
| dyn-BF | **0.569576** | 0.577608 | 0.577608 | 0.577608 | 0.577608 |

Adaptive is better than fixed in LL (−4.19% relative loss) and LI (−0.14%), tied in BF. **But adaptive is worse than expression-only on all three networks:** +2.35% LL, +0.32% LI, +1.41% BF. The strict dual success criterion of improved regulatory ranking *and no material reduction in expression generalization* is not established.

## 4. Meta-surrogate prediction MSE (lower is better)

| Held-out network | Fixed rho=1 | Global learned | Nonlearned stability | Adaptive learned |
| --- | ---: | ---: | ---: | ---: |
| dyn-LL | 0.355161 | **0.160675** | 0.247551 | 0.172974 |
| dyn-LI | 0.348845 | 0.170228 | 0.194031 | **0.154939** |
| dyn-BF | 0.368885 | 0.187158 | 0.273063 | **0.185246** |

The learned model predicts a temporal *sign-reproducibility surrogate*, **not true regulatory-edge probability**. Its mean MSE is only marginally better than the global baseline on these three contexts; this surrogate improvement does not translate into stronger GRN recovery.

## 5. Why the original reliability target may be misaligned

The current pseudo-target is `(1 + tanh(z_train * z_validation / tau²)) / 2`. For reversed edges, both signed temporal margins change sign, so their **product and the reliability target remain unchanged**. Thus the target largely measures *consistency of either direction*, not which direction is correct. Predicting it well is not sufficient to train a direction-aware reliability gate.

This is a mechanistic limitation of the current training target, not evidence that adaptive integration is impossible. A revised target must still be independent of benchmark gold edges: candidate options include cross-fitted directional *incremental held-out lagged gain* and independently sourced perturbation/TF-binding support, with explicit calibration tests. Avoid manufacturing a directional truth label from the evaluation gold network.

## 6. Research decision and next milestone

**Do not advance directly to RL** or claim learned biological prior reliability superiority from M2B. Retain reproducible code and all negative/mixed outputs.

**M2C:** 
1. Establish a **direction-sensitive, fold-held-out, source-blind** evidence usefulness surrogate; test invariance/symmetry explicitly and distinguish predictive direction from causality.
2. Introduce at least one **independent biological evidence modality** where feasible (e.g. TF-binding motif/ATAC/perturbation), with provenance and reference-overlap audit. Pseudotime alone is not heterogeneous evidence.
3. Evaluate source-availability masks, simulated missingness, wrong-sign/conflicting evidence, and topology-specific reliability changes against **fixed, global, and nonlearned stability gates**.
4. Use **uninspected generating networks** or true held-out biological contexts for confirmatory evaluation. Hold out whole network/topology contexts, not random edge-fold rows.
5. Match graph-edit costs/utility queries and keep untouched expression-test losses and AP/direction concordance as co-primary diagnostics.
6. Pre-specify an acceptable noninferiority margin for expression-test loss, and test whether any adaptive benefit survives across independent network seeds; bootstrap uncertainty at independent network/seed level.

**Go/no-go for M3 (RL):** Only after an adaptive gate shows a reproducible incremental advantage over a fully informed fixed/nonlearned gate, or after a compelling separate hypothesis for RL search is validated against matched greedy alternatives.

### Caveat
All three topologies were previously inspected in M1.2 and M2A. This is an engineering/development dataset, not prospective confirmation of an unseen biological context. The small candidate universes in LI/BF give substantial top-K and directional-metric variability.
