# BIO-V7: locked real K562 matched predictive comparison

**Protocol:** GCEV-BIOV7-K562-MATCHED-COMPARISON-2026-11-27-001  
**Date:** 27 November 2026 (Africa/Gaborone)  
**Parent revision:** `ff50587142c7995c281f156be598815b68682ac7`  
**Status at lock:** No BIO-V7 K562 fit or internal-score outcome has been evaluated.

## Research question

Does the existing GCEV target-group out-of-fold source-weighting policy improve a matched one-TF graph edit over graph-independent source adaptation, fixed source weights, and expression-only selection on the frozen real K562 cohort?

This is a developmental predictive comparison. The 64-gene universe was explicitly redesigned during BIO-V4B after inadequate TF coverage. It is not a fresh biological confirmation, and predictive error does not identify direct regulatory relationships.

## Data boundaries

Use exactly 394 fit, 131 internal-score, and 133 untouched-test cells from BIO-V3B; retain the BIO-V4B order of 16 TFs and 48 additional target-only genes. There are 1,008 nonself TF-to-gene pairs. Verify all source manifests and critical hashes in the executable protocol JSON. Do not select genes again or call the generic `run_fit` routine, which creates fresh partitions.

During freeze, only fitting expression rows may be transformed, modeled, checked against the original 60,671-gene matrix, or summarized. During evaluation, only fitting and internal-score rows may be used. Reading compressed bytes and verifying whole-file checksums are integrity operations; no statistics or predictions for the 133 test cells are allowed.

The original six-TF BIO-V2 perturbation responses remain exploratory and are not inputs. No new regulatory truth or perturbation outcomes are opened in this phase.

## Explicit constrained adaptation

The historical scaffold builder admits all genes and computes some expression gains with a different scaffold per candidate. BIO-V7 instead uses the same two-parent scaffold for all candidate actions at a target, restricted to the frozen TFs. Select parents by absolute fitting-cell Pearson correlation, breaking ties lexically. Insert exactly one additional eligible TF for each of 64 targets: 128 fixed scaffold edges plus 64 selected insertions, or 192 model edges for every method. There are 880 possible residual insertions. These local predictive graphs may contain cycles; no DAG restriction or causal orientation is asserted.

For source valuation, preserve the existing three-fold fit-cell outer split and 75/25 train-only search split. Build each scaffold on that search-fit subset. Select source-present and source-ablated actions on search-score cells and evaluate their paired difference on disjoint fit-cell validation. Retain existing source, parent, action features, group-by-target OOF Ridge gates, and clipping of predicted utility divided by 0.05 to [0,1].

For final choices, use a common scaffold chosen on all 394 fit cells. Average three-fold coefficient validation gains while holding that scaffold fixed, then rank the 880 residual gains globally. These are fit-only selection scores, not unbiased cross-validation estimates of the complete scaffold-selection procedure. Apply the inherited target-averaged OOF weights to those final gains. The final scaffold can differ from the nested training contexts; this limitation is recorded, not hidden.

## Accessibility invariance

DNase coverage is a target-level constant across candidate TFs. For a fixed target, adding an offset `0.3 * weight(target) * accessibility(target)` cannot change the action argmax. Enforce this exact algebraic invariant so floating-point common-offset rounding does not count as a source effect.

The DNase source-present and source-absent action labels must therefore be identical and its observed source-utility labels zero. Shared gate regression can still produce nonzero DNase weights; those are extrapolation of the regression, not empirical evidence of DNase action usefulness. DNase values may affect other learned source weights through source-context features. Report that distinction.

Use the two parent sources, not their product as a third independent source. The three processing versions reuse the same assay/biosamples and are sensitivity conditions, not independent experiments.

## Frozen design

The companion JSON is authoritative for exact constants.

| Element | Locked choice |
|---|---|
| Primary expression transform | log1p(original count-like X) |
| Normalization sensitivity | log1p(10,000 × X / full 60,671-gene per-cell total) |
| Primary prior version | Original mixed-caller BIO-V6 |
| Prior sensitivities | Matched v1; matched v2; motif-only |
| Methods | Expression-only; fixed sources; global source gate; no graph state; no parent state; no action state; full GCEV |
| Ridge coefficient penalty | 1.0, predictor scaling learned from each fitting subset |
| Source-gate Ridge penalty | Existing 10.0 |
| Source boost / shrinkage | 0.3 / 0.05 |
| Seed / folds | 20261127 / 3 |
| Hyperparameter search | None |
| Refit for internal evaluation | 394 fitting cells only |
| Untouched test | Remains sealed |

GEO describes STARsolo UMI counting and barcode QC. The exported X remains classified as count-like rather than independently certified raw UMI. Full-library normalization is therefore an explicitly declared sensitivity, not a silent replacement of the original expression scale.

All eight transform-by-prior conditions and all seven methods must be frozen before the first internal-score readout. No best preprocessing/caller may replace the primary case after seeing outcomes.

## Endpoint and uncertainty

Primary endpoint: mean normalized squared prediction error (NMSE), averaging the 64 target-gene errors and 131 internal-score cells. Each denominator is the target variance in the 394 fitting cells (population variance, floor 1e-9). Primary contrast: GCEV minus no-graph-state NMSE; negative favors GCEV. Secondary contrasts: GCEV minus fixed sources and expression-only. Report the 48 additional target-only genes separately as a prespecified descriptive endpoint.

For each contrast average the paired error difference over the same 64 genes per cell, then perform 10,000 paired cell-bootstrap resamples using seed 20261127 and percentile 95% intervals. These intervals condition on the fitted models, selected genes, and one source cohort; they do not quantify independent biological-replicate or regulatory-network generalization. Do not treat the 64 genes or 192 nested contexts as independent biological replicates. No confirmatory population superiority test is declared.

Report all method scores, selected-action differences, learned source weights, observed source-utility labels and timings. Runtime is descriptive; fixed and learned methods have equal graph budgets but different training costs. Do not report AUROC/AP/F1 because no independent edge-truth set is available.

## Execution and records

1. Commit this protocol and executable JSON before fitting.
2. Validate the implementation with targeted synthetic leakage, invariance, and budget checks.
3. Run `freeze`, which saves source labels, gates, full candidate contexts, residual scores and final selections for every condition, with implementation and input checksums.
4. Commit the freeze manifest and selected graphs before `evaluate`.
5. Evaluate the 131 internal-score cells once under the locked design. Keep all outcomes, including null or unfavorable comparisons.
6. Print the current ten-track project status, update the consolidated research report, and provide reproducible results.

The next scientific phase must address TF-specific evidence and independent biological validation. It is not authorized to reinterpret developmental expression prediction as proof of regulation.

