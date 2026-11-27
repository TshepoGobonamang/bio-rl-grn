# BIO-V7 — Real K562 matched predictive comparison

**Date:** 27 November 2026. **Status:** all 56 prespecified method–condition comparisons executed on the 131 internal-score cells; the 133 untouched-test cells remain sealed. **Scientific decision:** this developmental experiment does not establish an incremental GCEV advantage. The primary GCEV-versus-graph-independent difference is small and its conditional interval spans zero. Expression-only selection and fixed source weighting have lower primary point errors than GCEV.

This is an actual predictive comparison using the frozen K562 expression cohort, sequence-derived motif scores and experimental DNase context. It does not evaluate regulatory edge accuracy, direct binding or causal effects. No independent edge-truth set or new perturbation outcomes were accessed.

## 1. Prospective execution record

The [locked protocol](../docs/BIOV7_LOCKED_K562_MATCHED_COMPARISON_2026_11_27.md) and [executable configuration](../configs/biov7_k562_comparison_2026_11_27.json) were committed before fitting at `41baca2de083b8b13ed602bcfcbf804e910f71cb`. Implementation was committed at `d136b66dfb59a6cde055f40668d56202a8806b7a`. The freeze manifest and all final graph choices were committed at `ac5c86c820ec2b2178f769639d1e5c0e7de24ade` **before** internal-score evaluation. Development continues in [draft PR #51](https://github.com/TshepoGobonamang/bio-rl-grn/pull/51); no historical PR chain was merged by this experiment.

The [freeze manifest](biov7_k562_2026_11_27/freeze_manifest.json) has SHA256 `7b8865e7cd31b97f7577005ebec96f38f457cb0800c2a305ff70538e21a56826`. All eight transform-by-prior conditions and seven methods were frozen together. There was no hyperparameter search, gene reselection, cell repartitioning or favorable-case substitution after evaluation.

## 2. Cohort, transformations and matched decisions

The [BIO-V3B source](BIOV3B_REAL_K562_ONLY_MATRIX_2026_11_26.md) is GEO GSM9224207, with 658 source-labelled K562 cells and 60,671 genes. Its exported `X` remains **count-like**, rather than independently certified raw UMI counts. Source processing documentation describes UMI counting and quality filtering, but this does not resolve every exported-matrix or biological-QC limitation.

The frozen [BIO-V4B cohort](BIOV4B_REAL_K562_TF_AWARE_COHORT_2026_11_26.md) retains 394 fit, 131 internal-score and 133 untouched-test cells, 16 TFs and 48 additional target-only genes. Its 1,008 nonself TF-to-gene candidates are developmental: this universe followed the earlier 512-gene cohort's inadequate TF coverage. RPLP1 and RPLP2 remain under the original narrow exclusion rule.

The primary scale is `log1p(X)`. The prespecified sensitivity is `log1p(10000 × X / library_total)`, where each cell's denominator uses **all 60,671 original genes**, not only the selected 64. Only permitted fit/internal rows were transformed or summarized. Integrity checks of whole-file bytes did not authorize test-cell analysis.

Each target receives the same two-parent, TF-only scaffold for every method. Parents are selected from fitting expression by absolute Pearson correlation, with lexical ties. Each method then inserts exactly one remaining TF per target: **128 scaffold edges + 64 insertions = 192 graph edges**, selected from 880 residual actions. These predictive graphs may contain cycles; they are not asserted to be causal DAGs.

This constrained adaptation repairs an incompatibility in the historical scaffold interface, which allowed non-TF parents and candidate-dependent scaffolds. It leaves previous synthetic results intact. Nested source valuation uses three fit-cell folds; within each training fold, a 75/25 search split chooses matched source-present/source-ablated actions, evaluated on disjoint fit-cell validation. Gates retain target-group out-of-fold Ridge learning. Regression penalty is 1.0, gate penalty 10.0, source boost 0.3 and utility shrinkage 0.05.

Final expression gains use a common scaffold selected on all 394 fit cells. Their three-fold coefficient-validation scores are **selection statistics**, not unbiased estimates of the complete scaffold-selection procedure. Final contexts can differ from nested training contexts. Removing explicit parent features also leaves indirect scaffold information in action features.

## 3. Primary results

Primary case: original `log1p(X)` and the original mixed-caller BIO-V6 prior. NMSE averages squared prediction error over the indicated genes and the 131 internal-score cells, dividing each target's errors by its variance in the 394 fit cells (population variance; floor `1e-9`). Lower is better. All regression coefficients were fitted on those 394 cells.

| Method | NMSE, 64 genes | NMSE, 48 target-only genes | Insertions changed versus expression-only |
|---|---:|---:|---:|
| Expression-only | 0.401575 | 0.368196 | 0 |
| Fixed sources | 0.402125 | 0.366089 | 33 |
| Global source gate | 0.401575 | 0.368196 | 0 |
| No graph-state features | 0.402822 | 0.368196 | 2 |
| No parent-state features | 0.402598 | 0.367898 | 3 |
| No action-state features | 0.401324 | 0.368196 | 2 |
| Full GCEV | 0.402598 | 0.367898 | 3 |

Source: [complete 56-row results](biov7_k562_2026_11_27/summary.csv). The no-action-state ablation has the lowest primary point error; it was not promoted to a new primary method. Full GCEV and no-parent-state select identical primary graphs. That observation does not establish that parent features never matter in other settings.

| Prespecified contrast | GCEV minus comparator NMSE | Conditional paired-cell 95% interval |
|---|---:|---:|
| **Primary: no graph-state features** | **−0.000223** | **[−0.000636, +0.000175]** |
| Secondary: fixed sources | +0.000473 | [−0.003136, +0.004088] |
| Secondary: expression-only | +0.001023 | [−0.000407, +0.002488] |

Negative differences favor GCEV. All three primary-case intervals span zero. These are 10,000 paired cell-bootstrap resamples of each cell's mean error difference across the same 64 genes. They condition on this fitted model, selected genes and one cohort; **they are not independent biological-replicate intervals** and omit model-training uncertainty. There is no confirmatory population-superiority test or multiplicity correction. Complete results are in [contrasts.csv](biov7_k562_2026_11_27/contrasts.csv).

## 4. What the sources actually did

The [source-action audit](biov7_k562_2026_11_27/source_action_audit.csv) confirms the anticipated accessibility invariant. DNase coverage is constant across candidate TFs at a target. Adding `0.3 × target_weight × target_accessibility` to every action cannot change their ordering. The implementation enforces this algebraic cancellation to prevent floating-point offset artifacts.

Across every two-source condition, DNase changed **0 of 192 nested decisions** and all DNase utility labels were exactly zero. The 192 rows are 64 targets across three folds, not independent experiments. Nevertheless, the shared regression predicted positive GCEV DNase weights for **47 of 64 primary targets**, with mean weight **0.035305**. Those weights reflect regression extrapolation; they are not observed DNase action utility or evidence that accessibility is biologically irrelevant. Constant DNase offsets have no direct action effect, while DNase-dependent source-context features can alter learned motif weights.

Motif inclusion changed **96/192** nested choices on the original scale. Its mean paired utility was **−0.003427**, so changing decisions did not imply useful changes on average. Full GCEV's mean motif weight was 0.039244, with 39/64 positive target weights. On the normalized sensitivity, motif inclusion changed 78/192 choices and mean utility was −0.004915. The motif×accessibility product was not introduced as a third independent source.

## 5. Prespecified sensitivities

Both matched DNase processing versions preserve full GCEV's original-scale selections and the primary GCEV contrasts. Some ablated models change one choice, so complete case-level outputs remain necessary. Fixed-source selections are identical across all four prior variants within each expression transform, as expected from accessibility cancellation.

Removing DNase entirely changes learned contexts: original-scale motif-only GCEV NMSE is **0.403298**, versus **0.401852** for no-graph-state. The contrast is +0.001446, with conditional interval [−0.000177, +0.003121]. This descriptive difference is not evidence of a separate experimentally demonstrated accessibility mechanism.

Under full-library normalization, GCEV NMSE is **0.699570**, no-graph-state/expression-only **0.699506**, and fixed sources **0.706015**. GCEV minus fixed is −0.006444 [−0.011730, −0.001303], while GCEV minus no-graph-state is +0.000065 [−0.000575, +0.000726]. All four prior variants produce identical selections per method on this scale. The favorable fixed-comparator sensitivity does not replace the primary result. NMSE magnitudes across transforms should not be compared as a preprocessing competition because transformed targets and variance denominators differ.

## 6. Verification, cost and limits

The [audit summary](biov7_k562_2026_11_27/audit_summary.json) records identical candidates, TF-only parents, equal graph budgets, frozen-before-evaluation decisions, no reference network, no perturbation outcomes and **zero test cells used**. The local freeze took 58.30 seconds; internal evaluation took 7.09 seconds in the recorded Python 3.12.14 environment. These are run-specific timings, not a matched-compute superiority experiment: learned methods incur extra nested model fitting.

Independent-environment [GitHub reproduction #46004482269](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/46004482269) passed at commit `77d0552c2a885540e5e665d541dbd4d7f7eafd91`. The workflow compares all 3,584 target-level choices, 56 method summaries and 24 contrasts with the locked local reference at absolute tolerance `1e-10`. The downloaded [remote verifier](biov7_k562_2026_11_27/github_independent_verification.json) reports **zero maximum difference** in the 56 summaries and 24 contrasts, and all **3,584 identical decisions**, using identical source bytes. [Artifact provenance](biov7_k562_2026_11_27/github_reproduction_provenance.json) preserves the run and archive fingerprints. Both [PR CI #46004489861](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/46004489861) and [push CI #46004482116](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/46004482116) passed.

The separate [manual NumPy audit](biov7_k562_2026_11_27/independent_biov7_numeric_audit.json) solved all **448 primary-case regressions** without importing scikit-learn or the evaluator. Maximum differences were **2.98×10⁻¹⁴** for individual normalized cell errors and **9.99×10⁻¹⁶** for per-target NMSE. The [saved-output verifier](biov7_k562_2026_11_27/independent_verification.json) also passed 111 checksum entries, 56 error arrays and all bootstrap contrasts, with maximum metric discrepancy **8.88×10⁻¹⁶**. These establish computational reproduction of this experiment, not replication in a new biological cohort.

Existing limitations persist: gene-level rather than active-transcript promoter definitions; no BH-significant motif edges in BIO-V5; incomplete publication-level DNase assay QC; different cultures for expression and accessibility; a selected developmental gene universe; and no independent regulatory truth. Consequently no AUROC, AP/AUPRC, F1 or causal-accuracy values are reported. No neural graph model or trained RL policy was evaluated.

## 7. Research decision and next gate

The earlier [48-DAG benchmark](GCEV_REGISTERED_CONFIRMATION_48_NETWORKS_2026_11_25.md), [32-DAG mechanism study](GCEV_MECHANISM_CONFIRMATION_2026_11_25.md) and [separate 32-DAG parent-state study](GCEV_PARENT_STATE_ABLATION_32_DAGS_2026_11_25.md) also failed their primary superiority criteria. Those oracle-source simulations remain distinct from this real-expression experiment; their negative findings are preserved.

BIO-V7 completes the developmental predictive comparison and identifies a specific source-representation limitation. The next phase should qualify **TF-specific evidence**, or explicitly evaluate accessibility as a target-context interaction in a newly frozen design, and obtain a **fresh, uninspected external perturbation/binding panel** after predictions are frozen. The already inspected six BIO-V2 TF contrasts cannot serve as pristine confirmation. The 133 expression-test cells remain reserved for a later locked evaluation.
