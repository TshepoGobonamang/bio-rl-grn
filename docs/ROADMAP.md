# Engineering and experimental roadmap

## M0 — Import and audit the tested E1 starter

- [x] Import the tested Python package, CLI, pytest suite, GitHub Actions workflow, configuration and data preparation scripts.
- [x] Verify CI environment installation, six tests and smoke demonstration on Python 3.10–3.12.
- [x] Exclude downloaded data, caches, experiment outputs and secrets from Git.
- [ ] Record dependency versions, dataset release, SHA256 hashes and license requirements.

**Gate:** full installation and unit tests pass in CI.

## M1 — BEELINE E1: true graph-conditioned baseline

- [ ] Verify official BEELINE dyn-LL-500 replicate, expression orientation, gene identifiers and gold reference.
- [ ] Record an independently defined candidate regulator universe; do not derive the TF list from gold edges.
- [ ] Run Pearson, Spearman and greedy graph-conditioned ridge models with identical candidates.
- [ ] Test contextual marginal gains, add/remove actions, ranked-edge coverage, deterministic seeds and no gold access while fitting.
- [ ] Separate training, internal graph-scoring and untouched final expression evaluation when feasible.
- [ ] Track graph-utility evaluation counts, wall-clock time and edge selection trajectories.

**Gate:** reproducible first genuine benchmark run; results remain exploratory until confirmatory data separation is complete.

## M2 — Cross-fitted prior reliability

- [ ] Define per-edge evidence, presence masks, source-confidence anchors and graph/context features.
- [ ] Construct conditional predictive-gain targets without test gold labels.
- [ ] Fit reliability model with prediction loss, source anchoring and regularization.
- [ ] Freeze reliability parameters within each graph-editing trajectory.
- [ ] Compare no priors, static prior weights, learned context-only weights and graph-adaptive weights.
- [ ] Test missing, erroneous and conflicting priors; audit label provenance and overlap.

**Gate:** isolate measurable value of adaptation without learning biological truth from gold evaluation references.

## M3 — Matched sequential RL comparison

- [ ] Implement masked add/remove/stop policy and a reproducible graph-editing environment.
- [ ] Optimize identical graph utility used by the matched greedy control.
- [ ] Validate the undiscounted telescoping reward against terminal graph utility minus costs.
- [ ] Compare static, greedy and RL under equal utility-query and wall-clock budgets.
- [ ] Log seed sensitivity, convergence, action counts and graph trajectories.

**Gate:** publish accurate positive or negative evidence about the incremental value of RL.

## M4 — Independent benchmark and biological validation

- [ ] Expand to multiple BEELINE/SERGIO networks, DREAM5 and relevant human BEAR-GRN datasets.
- [ ] Ensure modality-matched and prior-matched comparisons.
- [ ] Evaluate prior-absent predictions against independent experimental references.
- [ ] Run interpretability interventions, ablations, prior corruption and runtime scalability experiments.
- [ ] Export complete reproducibility manifests, plots, confidence intervals and manuscript tables.

**Gate:** claims reflect independent evidence and appropriate biological limitations.

## Mandatory safeguards

- No reference gold edges in prior construction, policy reward, model selection or candidate-TF universe.
- Repeated optimization on an internal scoring partition may overfit; retain untouched evaluation contexts.
- Unannotated biological edges are not automatically proven negatives.
- Ridge conditional gain is predictive association, not causal verification.
- Learned evidence weights are not inherently calibrated probabilities of true edges.
- An RL agent is not a scientific contribution unless compared with strong matched non-RL controls.

**E1 baseline completed:** [2026-11-16 E1 real-data audit](../reports/E1_dyn_LL_500_run_43894848848_review.md) reports exploratory dyn-LL-500-1 results (greedy AP 0.4353 versus Pearson AP 0.3430; no independent test partition). These results **do not** establish RL or prior-reliability performance.

**M1.1 implementation complete:** Three-way held-out expression evaluation, evaluation-only gold provenance and directionality scoring, reciprocal-edge heuristic ablation, and a multi-replicate workflow have been merged into main. CI passed on Python 3.10, 3.11 and 3.12.

**M1.1–M1.2 experimental reviews completed:** [M1.1](../reports/M11_3_replica_beeline_2026_11_16.md) confirmed held-out expression improvements for expression-only greedy but weak direction concordance. [M1.2](../reports/M12_pseudotime_three_topologies_2026_11_16.md) showed improved AP from temporal information on three synthetic topologies, with expression-generalization tradeoffs.

**M2A completed:** [Real M2A report](../reports/M2A_crossfit_datasets_2026_11_16.md); 1,170 fold-edge reliability records. Gold-blind feature/target separation and quality audit are implemented.

**M2B implemented and empirically evaluated:** [M2B decision audit](../reports/M2B_adaptive_reliability_decision_2026_11_17.md), [official successful run](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/44091373596). The adaptive gate beats expression-only AP but does **not** consistently outperform a fixed temporal source weight on AP, orientation or untouched expression generalization. **NO-GO on asserting learned reliability superiority and on advancing directly to RL.**

**Next engineering milestone — M2C:** [Issue #13](https://github.com/TshepoGobonamang/bio-rl-grn/issues/13) — gold-blind *direction-sensitive* cross-fitted usefulness targets, independent biological evidence, controlled corruption/missingness, and uninspected network-context validation. The current reproducibility target is symmetric for reversed edges and is not a direction-truth label.


**M2C-A implementation — direction-sensitive held-out pseudotime target:** [merged PR #14](https://github.com/TshepoGobonamang/bio-rl-grn/pull/14) adds independently scored, signed directional lag-gain targets, antisymmetric forward-vs-reverse margins, branch-aware missing masks, fold isolation and an official [dataset workflow](https://github.com/TshepoGobonamang/bio-rl-grn/actions/workflows/m2c-direction.yml). See [protocol](M2C_DIRECTIONAL_TARGET_PROTOCOL.md). This is **target construction only**, not evidence of causal TF direction or adaptive gate superiority. Retain M2B negative/mixed baseline. **Next research gate:** audit M2C target distribution and stability, then evaluate new adaptive gates on held-out network contexts and independent evidence sources.
