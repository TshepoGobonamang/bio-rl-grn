# Bio-RL-GRN — Reproducible GRN Inference Experiments

**Status:** E1 expression-only baseline has run on BEELINE dyn-LL-500-1. M1.1 held-out expression and directionality improvements have been merged and CI-tested. The new M1.1 real-data workflow has not yet run; RL and adaptive evidence learning remain future milestones.

Bio-RL-GRN investigates whether **graph-dependent regulatory network construction** and **adaptive weighting of heterogeneous biological evidence** can improve directed TF-target network reconstruction compared with fully matched static, greedy, and reinforcement-learning controls.

## Method-to-experiment contract

1. **Graph-conditioned expression model:** fit each target gene with regularized regression using its currently selected TF parents; evaluate on data excluded from fitting.
2. **Adaptive evidence reliability:** train biological evidence weights using cross-fitted, out-of-sample conditional predictive gains, with missingness masks and source anchoring; **do not** interpret learned weights as calibrated biological truth probabilities.
3. **Graph-editing decision process:** choose add/remove/stop operations, score via changes in a fixed graph utility, and match search and compute budgets against greedy controls.
4. **Leakage prevention:** reserve gold-reference edges strictly for final scoring, do not derive candidate TFs or biological priors from held-out reference edges, and retain independently assessed expression data.
5. **Interpretability:** test policy explanations through paired evidence interventions, rather than claiming explanations are faithful from weights alone.

## Experiments

| Stage | Target | Evidence for advancement |
| --- | --- | --- |
| E1 | BEELINE dyn-LL-500 replicate 1, expression-only ridge, Pearson/Spearman and greedy control | Verified data orientation, sample splits, ranked edge universe, deterministic tests, graph-dependent marginal gains |
| E2 | Cross-fitted reliability learner + fixed/static controls | No gold leakage; predictive-gain target validation; robustness to source corruption and missingness |
| E3 | Masked add/remove/stop RL policy | Identical utility, candidate universe, and matched wall-clock / utility-query budgets against greedy |
| E4 | Repeated synthetic/BEELINE and DREAM5 experiments | Replication over networks and random seeds; confidence intervals; prespecified comparisons |
| E5 | Independent biological validation using available human RNA/ATAC and perturbation references | Prior-overlap audit, edge provenance, held-out regulatory evaluation |
| E6 | Explanation faithfulness, robustness and scaling | Intervention response, memory/runtime and ablation analysis |

## Reproducibility principles

- No private patient, student, or unpublished research data in the public repository.
- Store downloaded benchmark archives outside Git and record release versions, SHA256 checksums, and licenses.
- Keep gold-standard edges entirely outside model-fitting and reward-generation paths.
- Record RNG seeds, environment versions, TF-universe rules, hyperparameters, candidate count, and pre/post-evaluation commits.
- Pilot results from toy data or an internal tuning split are not confirmatory scientific results.

**E1 exploratory results:** [Review the first official BEELINE run](reports/E1_dyn_LL_500_run_43894848848_review.md). These results support continued evaluation, not claims of RL superiority.

**M1.1 directionality + untouched-test experiment:** See [protocol](docs/M11_PROTOCOL.md) and [manual workflow](https://github.com/TshepoGobonamang/bio-rl-grn/actions/workflows/m11-beeline.yml). Select 3 replicates for initial validation. The workflow downloads the official archive on the runner and uploads metrics and trajectories, not raw expression data. No M1.1 real-data metrics exist yet.

See [docs/ROADMAP.md](docs/ROADMAP.md) for milestones, [Issue #2](https://github.com/TshepoGobonamang/bio-rl-grn/issues/2) for the E1 record, and [Issue #4](https://github.com/TshepoGobonamang/bio-rl-grn/issues/4) for M1.1 validation.
