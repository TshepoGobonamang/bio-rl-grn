# GCEV direct parent-state ablation — preregistered computational results
**Experiment ID:** GCEV-PARENT-LOCK-2026-11-25-001  
**Status:** 256/256 completed, independently audited. **Primary superiority claim NOT established**.  
**Date:** 25 November 2026  
**Protocol was committed before any confirmation outcome** to GitHub commit `b77f235a2efd0d10a100413c5bf8ddc59a72cdd1` (protocol narrative), and the frozen executable JSON/code were committed before running seeds 7201–7232. Frozen executable implementation and JSON source were committed before execution: `bf185df34efbd0262c0ff3c2c779dde2265234ca`. The completed full benchmark and independent audit ran locally. The separate GitHub Actions full run is not yet verified. The experiment is not independent biological validation.

## Research question
Does adding **explicit parent-state columns** to a source- and action-aware reliability model improve directed TF-target edge recovery compared with a matched adaptive learner without those two columns? The original prior ablation `no_graph_state` eliminated both parent and action features; this study separates the feature blocks.

**Controlled models.** `GCEV` uses SOURCE + PARENT + ACTION; `no_parent_state` uses SOURCE + ACTION (removes only `parent_count` and `parent_target_correlation`); `no_action_state` uses SOURCE + PARENT; `no_graph_state` uses SOURCE only. All learned weights share identical nested paired graph-edit source-value labels, target-holdout folds and ridge settings. Action features still depend indirectly on the graph's fitted parent scaffold; this is **an explicit-feature block ablation**, not a counterfactual intervention that deletes all graph dependence.

## Methods and independent experimental units
32 completely new simulation DAGs (seeds 7201–7232), each evaluated under four paired gene-expression mechanisms (linear, saturating effects, pairwise regulator interactions, both), and two shared evidence conditions (clean; 50% missing + 25% corrupt): 256 fits. The 32 DAGs, not 128 mechanism outcomes, are the independent units. Each comparison uses the same expression train/internal/test split, candidate edge universe, source inputs and oracle DAG. The sources are **synthetic and derived from simulation truth**. The same DAG-dependent noise and sources are held fixed across the four mechanisms. No human binding, independent cellular replicates or physical TF regulation is inferred from this simulation.

## Main preregistered result
GCEV versus the `no_parent_state` model in **clean** source conditions, averaged over four mechanisms **within each independent DAG**:

- **Mean paired AP change:** **+0.005644**, 95% stratified DAG-bootstrap CI **[-0.000168, +0.011448]**.
- Predeclared success required lower CI > 0 and effect ≥ +0.010 AP. Neither condition was established.
- **Relative untouched expression-test NMSE:** **+0.070%**, 95% CI **[-0.226%, +0.382%]**; the +2% noninferiority criterion was satisfied, but no test-loss superiority claim follows.
- Verdict: **PRIMARY_PARENT_FEATURE_ADVANTAGE_NOT_ESTABLISHED**.

## Comparison by generating mechanism
| mechanism          |   GCEV_AP |   No_parent_state_AP |   No_action_state_AP |   No_graph_state_AP |   Fixed_AP |   Parent_incremental_AP |   Parent_relative_test_NMSE |
|:-------------------|----------:|---------------------:|---------------------:|--------------------:|-----------:|------------------------:|----------------------------:|
| linear             |   0.57026 |              0.57162 |              0.56654 |             0.57234 |    0.64828 |                -0.00136 |                     0.0013  |
| saturation         |   0.5368  |              0.52638 |              0.51597 |             0.5115  |    0.59784 |                 0.01042 |                     0.00411 |
| interaction        |   0.57678 |              0.57725 |              0.57634 |             0.59126 |    0.63579 |                -0.00048 |                    -0.00121 |
| combined_mechanism |   0.53919 |              0.52519 |              0.52034 |             0.51447 |    0.59595 |                 0.014   |                    -0.0014  |

Parent-block AP changes (secondary mechanism-level descriptions, not independent tests) ranged from -0.0014 to +0.0140. Fixed source weighting had the highest mean AP in every mechanism. Across all four clean mechanisms, mean GCEV−fixed AP was -0.063710, 95% independent DAG CI [-0.088927, -0.038403]. In this truth-derived-prior simulation, the fixed comparator is materially stronger on average.

## Why explicit graph features offered limited leverage
- In **all 256 condition-specific source-valuation datasets**, `parent_count` had **only one distinct value: exactly 2**. The scaffold size was set to 2 by design. Therefore, it offered no predictive variation within the experiment.
- The other explicit parent-state column, `parent_target_correlation`, did vary (mean within-condition SD 0.1088), but the full model did not achieve a reproducible AP advantage.
- Mean source-value-label OOF squared errors (descriptive; lower better): GCEV 0.005933, no-parent 0.005901, no-action 0.005500. This does not show a source-utility prediction benefit from adding the parent block.
- A source changed the preferred one-step edit in an average 15.36% of paired comparisons. This is a diagnostic, not a biological effect.
- The action feature block was not consistently beneficial either: `no_parent_state`−`no_graph_state` mean paired AP +0.002722 (95% CI [-0.012993, +0.019274]).

## Independent audit and limitations
The independent evaluator re-opened sealed ranking CSVs, recomputed noninterpolated AP/AUROC for **3,328 rankings** with maximum permitted numeric tolerance 1e-11, checked all 256 conditions, fully matched split/source hashes, reference availability and no gold access in fit. It returned **no errors**. The archived hashes are in `SOURCE_CHECKSUMS.json`.

**Interpretive boundaries:** Negative/non-significant evidence is a genuine outcome, not an invitation to change the test's frozen seeds or success threshold. Our sources are simulated from the oracle truth; no real biological source disjointness or effect is established. The graph model's 2-parent cap constrains structural variation, and the remaining action features can carry indirect parent information. A future algorithmic extension could deliberately introduce variable graph-state depth, but MUST use a fresh frozen holdout and compare with simpler fixed methods. More urgent for a top-tier computational biology submission is independently sourced biological validation.

## Figure files
- `F01_PARENT_STATE_INCREMENTAL_AP`: mechanism-specific incremental AP with network-level 95% intervals.
- `F02_MODEL_COMPARISON_BY_MECHANISM`: GCEV and fixed/ablated controls.
- `F03_PAIRED_DAG_EFFECT_DISTRIBUTION`: network-level benefit heterogeneity.
- `F04_SOURCE_VALUE_OOF_MSE`: source-value-label prediction quality (not regulatory truth).

**Reproduction:** `python scripts/gcev_parent_state_confirm.py --protocol configs/GCEV_PARENT_STATE_PREREG_2026_11_25.json --output results/gcev_parent_state_32_dags --workers 4` followed by `python scripts/verify_gcev_parent_state.py --root results/gcev_parent_state_32_dags`. Use exact source commit and protocol SHA256 as logged; environment Python 3.11+, numpy/pandas/scikit-learn. Report no top-tier biological claims until independently replicated real perturbation/binding comparisons exist.
