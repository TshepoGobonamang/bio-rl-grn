# GCEV mechanism confirmation: paired factorial DAG interventions (2026-11-25)

**Decision:** The predeclared interaction-specific GCEV ranking advantage is **NOT SUPPORTED**. The result is a reproducible computational simulation with oracle-derived sources, **not independent biological validation**. Negative findings remain in the record.

## Frozen research contract and independent re-scoring

- Locked prior-to-evaluation protocol: `docs/GCEV_MECHANISTIC_CONFIRMATION_LOCK_2026_11_25.md` (first GitHub commit `f446a3d8764148dbfc4a885f334004bb74318fba`).
- Frozen protocol JSON source SHA256: `a5b2cd5fc293bc220aec4df519dcb6ccb8e69bac7a390b23b73f10a665cff701`.
- Code and workflow committed before running seeds 6101–6132; GitHub PR #40.
- Full suite of 28 tests passed locally; distinct development-only DAG 9001 exercised the end-to-end runner before the confirmation seeds.
- **32 independent simulated regulatory DAGs** × **4 matched structural mechanisms** × **2 evidence conditions** = **256 fit/evaluation conditions**.
- All 32 DAGs successfully completed with **zero failures or post-hoc exclusions**.
- **2,816** directed-edge ranking files were independently re-scored for AP/AUROC; **2,816/2,816** matched their saved metrics within 1e−11. All 256 conditions passed explicit sealed-ranking, evaluation-only reference, shared-partition, shared-source and target-leakage checks. Audit: `independent_integrity_audit.json`.
- Bootstrap resampling unit: the **whole regulatory DAG**, stratified by 7/9 genes × target-switch/global-stable evidence regime; 4,000 samples, fixed RNG seed 20261125. Repeated mechanism/stress conditions are paired, not independent samples.
- Evaluated source scores were constructed from synthetic reference truth; no genuine ChIP, motif, perturbation or assay-independent prior was used. Zero biological relevance may be inferred from oracle source accuracy alone.

## Primary hypothesis: pairwise interactions should increase the graph-state benefit

Let `dAP(S,I) = AP(GCEV) − AP(no_graph_state)` for matched generator flags saturation `S` and pairwise interaction `I`. The predeclared per-DAG main interaction effect was:

`theta = ( [dAP(0,1)−dAP(0,0)] + [dAP(1,1)−dAP(1,0)] ) / 2`.

| Preregistered statistic | Observed |
|---|---:|
| Mean paired interaction effect `theta` | **−0.002080 AP** |
| Stratified network bootstrap 95% CI | **[−0.016369, +0.011090]** |
| Predeclared superiority threshold | Mean >= +0.020 AP and 95% CI lower > 0 |
| Verdict | **NOT MET** |
| Mean relative untouched-test NMSE (GCEV vs no-state) in interaction-containing settings | **−0.6976%** |
| NMSE 95% CI | **[−1.2342%, −0.1930%]** |
| NMSE noninferiority upper limit | +2%; criterion met, but primary AP criterion failed |

There is **no prospective evidence in this factorial test** that introducing two-parent nonlinear interactions increases the incremental regulatory edge-ranking benefit of the GCEV graph/action-context feature block. The estimate is small and the interval spans no effect. The point result is incompatible with framing the previous post-hoc linear/nonlinear difference as independently replicated.

## Ranking accuracy on clean synthetic sources

The following are descriptive means across 32 paired DAGs per mechanism, not independent 128-DAG averages.

| Generator mechanism | Fixed-source AP | Adaptive without graph-state AP | GCEV AP | GCEV − no-state AP |
|---|---:|---:|---:|---:|
| Linear | **0.6047** | 0.5387 | 0.5411 | +0.0024 |
| Saturation only | **0.5763** | 0.5167 | 0.5097 | −0.0070 |
| Pairwise interaction only | **0.5964** | 0.5370 | 0.5264 | −0.0106 |
| Both mechanisms | **0.5797** | 0.5223 | 0.5241 | +0.0017 |

**Fixed-source weighting is best on mean AP in all four mechanisms.** These priors have artificially high oracle-derived scores, so this is not evidence that fixed biological priors always outperform learning. It does mean the current algorithm cannot be claimed superior even within this generous simulated evidence design.

## Source-induced graph edit behavior

| Mechanism | Source inclusion changes chosen edge, clean | Combined 50% missing + 25% corrupted | Mean absolute held-out source value, clean |
|---|---:|---:|---:|
| Linear | 17.17% | 14.72% | 0.0193 |
| Saturation | 14.30% | 15.53% | 0.0122 |
| Interaction | 16.05% | 15.38% | 0.0186 |
| Combined | 14.53% | 15.07% | 0.0131 |

In approximately 83–86% of source-ablation comparisons, source inclusion **did not change the winning next edge**; such decisions therefore provide zero action-utility difference. Interactions did not materially raise this action-change rate. Source-valuation targets may be too sparse/weak to justify a complicated learned gate. This is an interpretation to test, not a demonstrated causal explanation.

## Secondary design effects and limitations

- Saturation-specific paired AP effect: **+0.00145** [95% CI −0.01328, +0.01509] — no conclusive effect.
- Factorial interaction × saturation AP effect: **+0.02170** [95% CI −0.01277, +0.05632] — uncertain, secondary.
- Mean interaction effect by design stratum: 7-gene/global-stable **+0.0019**; 7-gene/target-switch **−0.0082**; 9-gene/global-stable **−0.0163**; 9-gene/target-switch **+0.0142**. These are eight-DAG, descriptive strata; do not interpret as independent discovered subgroups.
- The `no_graph_state` ablation retains `candidate_count` and removes six extra `GRAPH_FEATURES`, including **action-dependent source-prior contrasts** as well as parent-state signals. The comparison is therefore a *bundled graph/action-context feature ablation*, **not a pure topological parent-state-only effect**. A stricter three-way ablation retaining action contrasts while removing only parent-state features is required before claiming mechanistic isolation.
- The source-decision policy is a **one-step matched graph edit**, not full graph-state RL or a learned search policy. Computing paired source valuations incurs extra nested ridge-loss evaluations; equal one-edit action budgets do **not** imply equal wall-clock cost. No RL advantage is claimed.
- Synthetic DAGs share a simulator and parent-ordering convention; even 32 graph realizations do not constitute external biological validation or establish direct causal regulation from expression.
- The new mechanism-specific networks are now **inspected**. Their seeds 6101–6132 are no longer usable as untouched confirmation for any method modification.

## Recommended research and journal decision

**Do not claim a graph-state-induced AP improvement or a biological method breakthrough.** The realistic contribution at present is a carefully controlled negative/mechanistic benchmark and a rigorous method for measuring paired evidence utility. For a stronger Computer Science-driven computational-biology paper, prioritize: (1) a stricter pure-parent-state ablation against an action-context baseline on a newly frozen test set; (2) quantify effective action-change/utility signal and compare value-supervision learning versus strong nonlearned rules; (3) verify at least one independent biological assay reference and external dataset. Any redesigned algorithm requires fresh independent networks and a preregistered contract.

Figures F01–F05 are downloadable in the companion result archive, including SVG source. Do not reuse a plotted confidence interval as a separate experimental replicate.
