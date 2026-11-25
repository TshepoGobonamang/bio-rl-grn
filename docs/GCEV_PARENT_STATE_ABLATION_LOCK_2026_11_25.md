# LOCKED validation — direct parent-state feature ablation

**Protocol:** GCEV-PARENT-LOCK-2026-11-25-001
**Freeze date:** 2026-11-25. Frozen **before inspecting any outcomes from seeds 7201–7232**.
**Research aim:** Does the pair of explicit parent-state features, `parent_count` and `parent_target_correlation`, add reproducible predictive information beyond all existing source-quality and next-action features?

### Scientific scope and distinction
The existing `no_graph_state` model drops **six** features: two explicit parent-state features, plus four source-conditioned action features. Therefore previous GCEV-minus-no-graph-state results did not isolate explicit parent-state information. The new strict comparison:
- `gcev`: source features + parent-state + action features.
- `no_parent_state`: same source features + **identical action features**, but excludes only `parent_count`, `parent_target_correlation`.
- `no_action_state` secondary diagnostic: identical source + parent features, without action features.
- `no_graph_state`, `global_source`, `fixed_sources`, `conditional`, Pearson, Spearman, greedy remain matched comparators.

All learned gates use **the identical source-value target records, target-group folds, imputation/normalization, ridge penalty, output clipping, source-value scale and outer partitions**. All source and action features are computed from the same source-blind training graph scaffold, regardless of the ablated feature set. Caveat: action features themselves are functions of graph state; this isolates the **direct incremental predictive contribution of the explicit parent-state columns**, not all possible causal effects of parent graphs.

### Confirmation generators and primary analysis
32 freshly drawn independent DAGs from RNG seeds **7201–7232**, excluding all earlier pilot/benchmark seeds 11/42/77, 2001–2048 and 6101–6132.

For each DAG reuse the same graph, exogenous Gaussian noise, candidate edge set and truth-derived *synthetic* evidence across four mechanisms: linear, saturation, pairwise interaction and combined. Reuse the prior fixed generator definitions (Bernoulli 0.30, signed coefficients absolute Uniform[0.35,0.80], n=7 or 9, n_cells=140 or 180). Source regimes target_switch/global_stable paired with size (8 DAGs per size×regime stratum). Simulated source files are generated once per DAG, copied unchanged across all four mechanisms; no biological independence claim.

Two preregistered evidence conditions: clean (0 missing, 0 corruption) and stressed (50% source dropout and 25% corruption), with deterministic intervention-independent RNG. Exactly **32 × 4 × 2 = 256** fits, each evaluated on exactly the same source and outcome partitions. Preregistered source boost 0.30, shrinkage 0.05, ridge alpha 1, 3 nested folds, 2 source-blind parent scaffold members, top K=7, graph-edit query budget 1600. Continue existing 60/20/20 outer split and grouping of source-value targets by held-out target gene.

### Primary hypothesis (one predeclared comparison)
Let `dAP_m(DAG)=AP_gcev - AP_no_parent_state` in each clean mechanism m, and
`delta_DAG=(dAP_linear + dAP_saturation + dAP_interaction + dAP_combined)/4`.
One DAG = one independent statistical unit, not four mechanism observations, source targets, folds, edge rows or expression cells. Primary statistic = mean delta over the 32 independent DAGs, equal stratum weight.

Stratified percentile bootstrap (4,000 draws, fixed RNG 20261125) resamples all 8 DAGs within each of 4 size×source-regime strata and retains each DAG's four paired mechanisms in every draw.

**Primary success criterion fixed ex ante:** 95% CI lower bound for mean delta > 0 AND mean delta ≥ +0.01 absolute AP, AND upper 95% CI for the corresponding mean relative untouched-test NMSE difference (GCEV/no_parent_state, clean, averaged within DAG over four mechanisms) ≤ +0.02. If any criterion fails, direct parent-state benefit is **NOT established**.

### Additional preregistered measurements
- `no_parent_state` versus `no_graph_state` to estimate the effect of action features in the absence of parent features, and GCEV versus `no_action_state` to estimate action-context contribution with parent features present.
- Factorial interaction-of-AP with parent-state feature gain, linear vs interactions, and stress vs clean (descriptive, multiple comparisons not corrected).
- Nested target OOF MSE `gcev` vs `no_parent_state`; distinguish source-value-target prediction from true-edge recovery.
- Frequency/variance of `parent_count` and `parent_target_correlation`, including whether `parent_count` is constant; if >95% constant, report limited meaningful parent-state variation.
- Fixed-weight competitiveness, AP/AUROC/P@7/orientation on fully specified synthetic references, untouched expression-test NMSE, source edit-change fraction and source valuation magnitude, model/graph utility evaluation counts, training seconds.
- Seal fitted rankings and source file hashes before reference evaluation, maintain separate evaluation stage, verify 256 complete conditions, 32 DAG identity, shared expression noise, source hashes, partitions and exact reproducibility within numeric tolerance.
- All failures/retries logged. No replacement of unfavorable seeds and no modifying any mechanism strength after outcome inspection.

### Interpretation and stop rule
The goal is to test whether explicit parent-state features provide an **independent benefit**, not to force an accuracy gain. An AP CI overlapping zero or a NMSE noninferiority failure is a negative primary result. Earlier 48-network and 32-DAG mechanism suites remain inspected *development evidence*, not rerun as independent confirmation.

This benchmark uses **truth-derived simulated priors** and a deliberately simple ridge-based learner; it is not independent biological validation, not source-disjoint TF binding evidence and not an RL achievement. Actual parent-state effects mediated by action features cannot be excluded by this feature ablation. A separate biological source-disjoint test remains required for a leading computational biology journal.

**Protocol commit must predate the first evaluation of seeds 7201–7232.**
