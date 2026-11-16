# Engineering and experimental roadmap

## M0 — Import and audit the tested E1 starter

- [ ] Import the tested Python package, CLI, pytest suite, GitHub Actions workflow, configuration and data preparation scripts.
- [ ] Verify fresh checkout, environment installation, tests, smoke demonstration.
- [ ] Exclude downloaded data, caches, experiment outputs and secrets from Git.
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

**First actionable deliverable:** import the existing locally tested E1 scaffold, verify CI, then run the first authentic BEELINE replicate.
