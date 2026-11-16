# M2 design proposal: cross-fitted adaptive evidence reliability

**Status:** design proposal for review; NOT an implemented reliability model, not experimental results, and not a training prescription for an RL agent.

## 1. Why this is our next experiment

The [M1.2 pilot](../reports/M12_pseudotime_three_topologies_2026_11_16.md) showed:

| BEELINE topology | Greedy AP | Fixed-pseudotime greedy AP | Untouched expression loss, greedy → fixed pseudotime |
|---|---:|---:|---:|
| dyn-LL | 0.4864 | 0.6566 | 0.8767 → 0.9365 |
| dyn-LI | 0.5801 | 0.6199 | 0.6861 → 0.6892 |
| dyn-BF | 0.2918 | 0.3385 | 0.5696 → 0.5776 |

Pseudotime improved ranking/direction, yet worsened untouched expression loss on all three. Temporal-only was unusually effective on dyn-BF (AP 0.7427), but not dyn-LL (AP 0.2886). Therefore test whether **conditional evidence trust**, rather than a single source-wide coefficient, adds value. Do NOT take these data as independent confirmatory evidence: results have already been inspected.

## 2. Essential correction to single-source aggregation

Original weighted mean:
`p_e = sum(m_e,k rho_e,k z_e,k) / (epsilon + sum(m_e,k rho_e,k))`.
With only one nonmissing source and small epsilon, this becomes approximately `z_e` **independently of rho**. A single-source gate would therefore have negligible effect.

**M2 pilot:** use `s_e(G) = rho_e(G) * z_e^time`, with `0 <= rho <= 1` and signed `-1 <= z <= 1`. Missing source evidence has a separate availability mask (never silently map unavailable to negative evidence).

**Future multi-source:** use magnitude-sensitive aggregation, e.g. a fixed-denominator sum of `m rho z`, plus explicit source-coverage features; distinguish total evidence strength from relative source weights. The formula must be unit-tested for one-source attenuation and all-missing cases.

## 3. Staged architecture: prove the gate before making it graph-adaptive

### Stage M2A — fixed context gate (no graph dependence)

For edge e=(i,j), input signed temporal evidence `z_e`, source-availability `m_e`, branch support, the cross-fit variation/sign agreement of temporal margins, number of usable bins, and *source-blind* expression association features computed from fitting cells.

`rho_e = sigmoid(g_phi(h_e^fit, z_e, m_e, context))`.

Start with **small regularized logistic regression or a shallow network**. With only 7–18 genes and 42–306 directed candidates, a large neural reliability network invites overfitting and pseudo-replication.

The expression-only greedy baseline and temporal-greedy baseline must retain exactly the same candidate universe, edit budget, fit/scoring/test split and ridge model.

### Stage M2B — true graph-conditioned gate

After M2A is validated, add local parent-set/context features:
`rho_e(G) = sigmoid(g_phi(h_e^fit, z_e, m_e, h_e^graph(G)))`.

If graph context affects *existing* selected edges, the one-step utility difference is the **full** `U(G + e) - U(G)`, not just `delta expression loss + rho_e z_e`. Recompute (and cache) any affected selected-edge reliabilities during graph edits. Freeze the meta-learner parameters `phi*` before greedy/RL training; update its **outputs** as G changes, not its parameters. This keeps the reward well-defined.

## 4. Training objective: conditional usefulness, not biological truth

The reliability learner cannot see evaluation gold-reference edges, directly or through a TF annotation derived from gold.

Reserve 20% untouched expression test and 20% graph-edit scoring as in M1.1. Within the remaining 60% fitting cells, generate K-fold source-blind graph scaffolds and independent pseudo-targets, keeping all feature construction within the respective fitting fold:

1. For each inner fold b, construct pseudotime evidence on `F_b`; form the same statistic independently on validation cells `V_b` using branch-aware bin counts.
2. Measure a **held-out signed temporal margin** and/or directional reproducibility target. This tests whether the temporal sign is stable out of fitting cells; it does NOT certify correct regulator direction.
3. Independently compute held-out graph-conditioned expression-gain targets on source-blind scaffolds `G_b^(-e)`, using an inner training/validation partition. Exclude edge e from scaffolds.
4. Train the reliability model to predict conditional **out-of-fold temporal reproducibility** and examine expression-gain agreement, with source anchors, weight decay and coverage-aware regularization.
5. Include an expression-only predictor of the same pseudo-target to isolate whether the temporal evidence adds information. Report predictive losses/calibration on fresh inner folds.

Crucial: the expression-reconstruction surrogate and the direction-reproducibility surrogate can disagree. Do not silently treat the former as a ground-truth directional target. The model's `rho` remains a **task-conditional usefulness weight**, not a calibrated probability that the biological edge exists.

**Minimum feature columns:** `source,target,source_id,z_signed,m_present,coverage,branch_n,bin_n,sign_stability,bootstrap_variance,expression_baseline,context`. Include provenance and fold keys; never include gold labels in model input.

## 5. Graph objective and matched baselines

For a simple source-frozen gate use

`U_phi(G) = -L_expr(G) + lambda_t * sum_{e in E(G)} rho_phi(e,G) z_e - lambda_s |E(G)|`.

In stage M2A, `rho_phi(e)` is graph-invariant and the temporal contribution has straightforward incremental gains. In stage M2B, changes in `G` can alter reliability for existing selected edges, so always score the **full utility difference**.

The temporal scale `lambda_t` is chosen on development data only, using a predeclared search range that is identical for fixed and adaptive methods. The prior M1.2 constant `0.015` is a historical baseline, **not a universally optimal coefficient** across network sizes.

**Matched comparisons (minimal):**
- expression only (no temporal contribution)
- fixed temporal coefficient (M1.2)
- cross-fit stability heuristic (nonlearned gate)
- learned **global** reliability gate, shared across edges
- learned context-/edge-dependent reliability gate
- temporally shuffled evidence control
- optional oracle reliability on synthetic truth as an explicitly labeled upper-bound **only**, never as trained model or fair competitor.

## 6. Dataset and fold design

Development: M1.2 topologies dyn-LL, dyn-LI and dyn-BF may be used for exploratory engineering only; they have already been examined. Within each dataset, maintain **fit / graph-scoring / final-expression-test** independence.

Cross-fit within fitting partition, using branch masks and a dynamically chosen viable bin count (branch-level fold sizes may be too small for 16 bins). Require a minimum number of observed cells per bin in both the fitting and validation split.

External confirmation: new simulation seeds/independent generating graphs and a previously uninspected network family, with an explicit train-topology versus held-out-topology split. Real biological validation requires independent motifs, accessibility or curated evidence, with reference-overlap auditing. BoolODE provided time may be optimistic compared with independently estimated pseudotime.

## 7. Pre-specified evaluation and failure conditions

Report **by topology**, with replicate-level intervals, for:
- Noninterpolated average precision (AP), AUROC, Precision@K, orientation accuracy with ties, and reciprocal pairs@K.
- Untouched expression-test loss of the selected graph (never use it for tuning).
- Cross-fitted reliability pseudo-target prediction score/calibration and source dropout/corruption behavior.
- Runtime/number of utility evaluations, feature availability, learned-weight distribution and gate saturation.
- Information matching and a shuffled temporal negative control.

**Success is NOT just higher AP than expression-only.** Evaluate against a fully matched fixed-weight temporal baseline and require the adaptive gate to retain or improve directional ranking **without a material expression-generalization regression** (predefine a noninferiority margin and confidence interval on an untouched development-independent evaluation set).

If the gate collapses to all-ones/all-zeros, fails to outperform cross-fit stability heuristics, or performs comparably with shuffled temporal data, report the negative finding and do not assert adaptive reliability novelty.

## 8. Suggested code organization and checks

Proposed files (not yet created):
- `src/biorlgrn/evidence.py`: evidence table schema, masks, source normalization/provenance.
- `src/biorlgrn/reliability.py`: fold-local feature extraction, cross-fit targets, small gate.
- `src/biorlgrn/utility.py`: explicit frozen-phi graph utility and exact edge-delta functions.
- `scripts/run_m2.py`: matched fixed/global/adaptive/shuffled comparisons.
- `tests/test_reliability.py`: one-source attenuation, all-missing/zero-evidence distinction, no reference use, fixed-phi rewards, bin/fold isolation, label shuffle and deterministic seed.
- `configs/m2_crossfit.yaml` and `.github/workflows/m2-beeline.yml`: reproducible launch and artifact output.

**Implementation gate:** first construct and audit the fold-local reliability meta-dataset, then train the smallest gate and compare against fixed and nonlearned-stability baselines, then (only if supported) add graph-conditioned reliability and eventually reinforcement learning.

## Related scientific context

- Greenfield/Inferelator multi-study prior weighting and Inferelator 3.0 motivated principled prior incorporation rather than assuming priors are always true.
- Review of GRN prior knowledge: https://doi.org/10.1038/s44320-025-00088-3
- Scribe (Qiu et al., 2020): pseudotime-based associations do not constitute temporally coupled causal observations; https://doi.org/10.1016/j.cels.2020.02.003
- PseudoGRN (2025): pseudotime-directed inference is an active existing research topic; https://pmc.ncbi.nlm.nih.gov/articles/PMC12756061/

This document records a **proposed** procedure for discussion. No M2 results or learner implementation are implied.
