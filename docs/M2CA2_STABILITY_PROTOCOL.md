# M2C-A2: Temporal-target stability and incremental-information experiment

## Decision this milestone addresses

[M2C-A informativeness audit](../reports/M2C_A_target_informativeness_audit_2026_11_17.md) found a direction-sensitive yet unstable held-out lagged-gain target (cross-fold rank correlations 0.171/0.243/0.025 on dyn-LL/LI/BF), with weaker post-hoc reference orientation than the original signed pseudotime margin. **No adaptive learner or reinforcement learning policy is trained here.** This stage tests stability, control specificity, and incremental information before M2C-B.

## Predeclared experimental design

**Development structures:** `dyn-LL-500-1`, `dyn-LI-500-1`, `dyn-BF-500-1` (previously inspected).

- Fixed outer seed **42** and 60% fitting / 20% graph-scoring / 20% untouched expression-test partition; only 60% fitting cells are used to construct directional evidence.
- Inner K-fold **3**, with **predeclared** seeds 42, 73, 137, bin caps **8, 12, 16**, lag **1**, minimum 3 cells per bin.
- In each F_b → V_b fold evaluate (i) original **fitting-only** signed temporal margin `z_signed`, (ii) directed held-out source-specific lagged gain `gain`, (iii) pairwise difference `margin = gain(s→t)-gain(t→s)`. The previous target feature `z_signed` is computed on F_b only.
- A **shuffled pseudotime control** with seed 42 and 12 bins, shuffling observed values separately in each F_b and V_b and within each branch; cell IDs and NA masks remain unchanged. Never shuffle the test gold labels.
- Report edge-fold counts, branch support, missing targets, negative/clipped gains, three-pair crossfold **Spearman** correlations and three-fold sign consistency for each method. The folds overlap in their training subsets; these are descriptive stability diagnostics, not independent replicates.

**Locked confirmation structures**, not yet used in the project: `dyn-CY-500-1` (cycle), `dyn-TF-500-1` (trifurcating) and `dyn-BFC-500-1` (bifurcating-converging). These topology families are documented in the official BEELINE release; each has 10 expression realizations in the archive. One expression replicate per topology is **not** a confirmatory claim of general biological validity.

Before reading any confirmation gold edges, generate frozen summaries using **only** inner seed **42** and **12 bins**, plus the branch-preserving shuffled negative control. These settings were fixed before inspection of the three confirmation families.

## Predeclared prospective incremental-information comparison

For every candidate directed edge aggregate scores over the 3 inner folds. Compute:

1. `z_signed` (original signed temporal margin, primary comparator).
2. `gain` (held-out added-regulator relative predictive usefulness).
3. `direction_margin` = gain forward minus gain reverse (mathematically antisymmetric).
4. `shuffled_margin` computed from within-branch shuffled pseudotime.
5. **Fixed blend:** `0.75 * scale95(z_signed) + 0.25 * scale95(direction_margin)`, where `scale95(a) = clip(a / max(quantile(|a|,0.95),1e-10), -1,1)`. The scaling uses candidate-score magnitudes **without gold data** and was fixed prospectively. This blend is a **diagnostic comparator**, not a trained adaptive reliability mechanism.

Write all ranking CSVs **before** importing any gold network. The locked evaluation then computes non-interpolated average precision, AUROC, direction concordance (ties count 1/2), precision@K, and numbers of reference orientations corrected or newly spoiled relative to `z_signed`. Report per structure; do not tune bin count, seed, blend weight or target based on these gold metrics.

**Important:** a positive reference statistic is biological **association to the simulator's gold edges**, not proof of real-world causal TF regulation. An antisymmetric score is not itself evidence.

## Artifacts

Per network:
- `m2ca2_fold_scores.csv`: scores for each directed candidate × inner fold × setting, all **gold-free**
- `m2ca2_stability.csv`: fold stability, availability, sign consistency, clipping and negative-gain rates
- `m2ca2_edge_scores.csv`: fold-averaged unsigned and signed candidate scores
- `m2ca2_manifest.json`: data hashes, training/validation fold hashes, role and frozen settings

For confirmation networks only:
- `ranked_z_signed.csv`, `ranked_gain.csv`, `ranked_direction_margin.csv`, `ranked_fixed_blend.csv`, `ranked_shuffled_margin.csv` (saved **before** gold)
- `m2ca2_reference_audit.csv` and `m2ca2_orientation_changes.json` (evaluation-only)

Aggregate:
- `m2ca2_development_stability_summary.csv` and `m2ca2_confirmation_stability_summary.csv`
- `m2ca2_locked_confirmation_metrics.csv` (evaluation-only)

## Acceptance and go/no-go

Advance to M2C-B **only if** the revised evidence demonstrates informative, reproducible variation beyond the original signed temporal margin. Specific checks:
- New gain/margin crossfold correlation substantially better than the M2C-A baseline under predeclared bin/seed settings, or demonstrate a robust shrinkage/aggregation rationale without gold-based tuning.
- Shuffled pseudotime should not perform comparably to the proper signal on independent networks.
- Locked fixed blend should provide **genuine incremental** AP or orientation compared with `z_signed`, without an excess of new orientation errors.
- Do not select the best setting by reference AP from held-out topologies. Further evaluation on untouched independent network seeds and actual biological modalities remains essential.

A failure is a valid research outcome: **retain original signed temporal evidence and do not learn from a demonstrably unstable target**.

## Execution

GitHub workflow: [M2C-A2 stability / confirmation](../.github/workflows/m2ca2-stability.yml). It downloads BEELINE v4 once, generates all development and locked-confirmation candidate scores first, evaluates gold only in the final step, and uploads the derived tables, not raw expression files.

```bash
python -m pip install -e '.[dev]'
python -m pytest -q
python scripts/run_m2ca2.py --raw-root data/raw/beeline --results-root results/m2ca2 --role development
python scripts/run_m2ca2.py --raw-root data/raw/beeline --results-root results/m2ca2 --role confirmation
python scripts/evaluate_m2ca2.py --results-root results/m2ca2 --raw-root data/raw/beeline
```

## Limitations

- Cells sorted by pseudotime are independent snapshots, not the same cell tracked longitudinally.
- The current lagged model controls only target autoregression and a bin-time trend, not other parent genes. This is a diagnostic baseline, not validated causal inference.
- Multiple inner folds, seed splits and candidate edges are dependent; do not calculate p-values treating them as independent experiments.
- These additional topologies may be unfamiliar to the project but come from the same BoolODE simulation family; statistical generalization requires new network realizations/seeds and independent experimental support.
- In an exceptionally small candidate universe, Precision@20 may include nearly every candidate; compare AP and direction carefully.

**Source for available BEELINE synthetic topology families:** https://murali-group.github.io/Beeline/reproducing-results.html .
