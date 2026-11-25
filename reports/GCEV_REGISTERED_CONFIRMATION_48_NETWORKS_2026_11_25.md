# GCEV — Preregistered controlled synthetic confirmation (25 November 2026)

**Protocol ID:** `GCEV-LCK-2026-11-25-001`  
**Frozen protocol SHA256:** `556c6745948cc89b10c9b7ee92a64d4b5637bcefbe180c0c062082b7d2f0272c`  
**Protocol publicly frozen on GitHub before experimental execution:** [commit cb7b3d2](https://github.com/TshepoGobonamang/bio-rl-grn/commit/cb7b3d2aff6827c98c370ba14fbc19293405a1c3)  
**Analysis code and tests committed before confirmation execution:** [commit c4d98e9](https://github.com/TshepoGobonamang/bio-rl-grn/commit/c4d98e9b48c5fd99dc0d630d39003ac2d2f559e8)  
**Data role:** *Simulation-only, oracle-derived source evidence, not real independent biological validation*.  
**Decision:** **PREREGISTERED PRIMARY CRITERION NOT MET**.

## Design and multiplicity

- **48 independent generated DAG realizations** with uninspected network seeds 2001–2048; previous development seeds 11, 42 and 77 are excluded.
- Exactly **24 linear** and **24 nonlinear** data-generation DAGs; 24 seven-gene/140-cell and 24 nine-gene/180-cell cases.
- Two source regimes, each used on 24 realizations: target-specific source reliability switches, and source reliabilities that are stable across target indices.
- Each DAG was evaluated under **four matched source states**: clean, 50% source dropout, 25% source-score inversion, and combined dropout/inversion. Thus there are **192 fit/evaluate units**, not 192 independent networks.
- All 48 intended networks completed; **zero failed or excluded networks**. Eleven methods produced **2,112 sealed full-universe rankings**, all with an evaluation-only complete synthetic reference.
- The principal ablation holds fixed the **same per-target paired source-value labels, target-out-of-fold groups, source features, source measurements, splits and ridge learner**. Only the graph-state feature block is removed for `no_graph_state`.
- All directed candidate edges are assessed in the same synthetic truth universe. Absence of an oracle-generated edge is a true negative **in this synthetic generator only**, not in real ChIP/curated data.

## Primary test: GCEV versus graph-independent adaptive reliability

On *clean* evidence, the observed **mean per-network AP difference** (GCEV minus `no_graph_state`) was **+0.007715**, with a **95% bootstrap interval [−0.010750, +0.026722]**. Bootstrap draws (4,000, RNG seed 20261125) resample independent networks within each of the eight predeclared graph-family × size × evidence-regime strata, never edges, cells or folds.

The preregistered directional improvement criterion fails: the interval's lower bound is below zero. GCEV exceeded `no_graph_state` on 21/48 networks and tied on two; these counts are descriptive, not a substitute for the locked primary interval.

The companion relative untouched-expression-test NMSE contrast (GCEV minus `no_graph_state`) averaged **−0.000804**, with **95% CI [−0.007940, +0.006036]**. The specified upper limit of +0.02 is satisfied, but because the AP criterion fails, the **joint primary success rule is NOT MET**.

**Conclusion:** We cannot claim a reproducible overall incremental benefit from graph-state features in this locked synthetic benchmark.

## The strongest fixed baseline is better on clean synthetic evidence

| Method | Mean AP, clean | Mean AUROC, clean | Mean untouched test NMSE, clean |
|---|---:|---:|---:|
| Fixed source weighting | **0.620365** | **0.828134** | 0.878056 |
| GCEV | 0.545860 | 0.802644 | 0.859175 |
| Graph-independent adaptive | 0.538145 | 0.806540 | 0.860085 |
| Global source weight | 0.523710 | 0.800305 | 0.861103 |
| Expression-only conditional | 0.497349 | 0.793515 | 0.855729 |
| Greedy expression model | 0.496515 | 0.795693 | **0.851707** |

The clean **GCEV − fixed** AP difference was **−0.074505** (95% network-bootstrap interval **[−0.110828, −0.036043]**). In this simulated setting, fixed weighting has better ranking recovery than GCEV. Neither the fixed model's higher AP nor GCEV's lower NMSE establishes causal regulatory inference.

## Source corruption and missingness — report all four declared conditions

| Scenario | Mean AP: GCEV | Mean AP: no graph state | GCEV − no graph state (95% CI) | Mean AP: fixed |
|---|---:|---:|---|---:|
| Clean | 0.545860 | 0.538145 | +0.007715 [−0.010750, +0.026722] | **0.620365** |
| 50% missing | 0.507068 | 0.505981 | +0.001087 [−0.008541, +0.012163] | **0.548540** |
| 25% score inversion | 0.509399 | 0.519944 | −0.010545 [−0.027204, +0.005771] | **0.522990** |
| Missing + inverted | 0.492404 | 0.505907 | **−0.013503 [−0.025987, −0.001651]** | 0.486952 |

The combined stress condition produces a **negative graph-state ablation effect**, and its descriptive 95% interval excludes zero. Because stress outcomes are secondary and jointly correlated, it should not be interpreted as a standalone unadjusted confirmatory significance claim.

The fixed source-weight method declines from 0.620365 AP in clean conditions to 0.486952 in combined stress. GCEV declines from 0.545860 to 0.492404. This is evidence of **different sensitivity patterns**, not proof that GCEV wins under stress: its combined-condition advantage over fixed weighting is only +0.005452 AP, with 95% CI [−0.021591, +0.033662].

## Source and graph-family heterogeneity (secondary, descriptive)

- Linear DAGs: mean clean GCEV − `no_graph_state` AP = **−0.0253** (24 networks).
- Nonlinear DAGs: corresponding difference = **+0.0407** (24 networks).
- `global_stable` source regime: clean difference = **+0.0277** (24 networks).
- `target_switch` regime: clean difference = **−0.0122** (24 networks).

This heterogeneity is hypothesis-generating. We did not preregister family-specific superiority or an interaction-test threshold; changing the main claim to the favorable nonlinear stratum now would invalidate the original primary test.

## Computational and interpretation boundaries

The model requires **189 source-valuation ridge-loss evaluations** per run for seven-gene networks and **297** for nine-gene networks, before comparison ranking and testing. Mean GCEV target/training elapsed time was around **0.26 s per run** on this particular multi-process local environment. Because fixed weighting does not fit an equivalent source-value learner, a fixed-weight AP comparison is **not matched for training compute**; no computational-efficiency superiority is asserted.

**Validation:** An independent post-fit audit recomputed AP and AUROC from all 2,112 frozen ranking CSVs and checked their file hashes. The maximum absolute recomputation discrepancy was **1.11×10⁻¹⁶**, with zero failures across 192 fit/evaluate units.

The evidence sources in this experiment are generated with reference-oracle knowledge by design. Although fitting code never reads the reference file, reference information enters through the synthetic evidence generator. This is a legitimate *in-silico controlled source experiment*, NOT source-disjoint biological validation. The 48 networks derive from two parameterizations of one ordered-DAG generation framework; they do not represent 48 independent biological experiments or multiple external GRN benchmarks. There are no trained RL policies and no externally verified direct TF-binding labels.

## Final scientific decision

**No-go for a journal claim that graph-conditional evidence valuation consistently beats graph-independent adaptation or strong fixed weighting.**

**Go for reporting a validated computational mechanism and an informative negative/mixed preregistered experiment.** The next genuinely new confirmation should test externally obtained biological evidence and compare graph-state features under predeclared conditions, without reusing the now-inspected seeds 2001–2048 as an untouched test set. Do not alter weights, seed selection or inclusion criteria after seeing these results to transform this study into a positive one.

## Figures and generated artifacts

- `figures/F01_GRAPH_STATE_INCREMENTAL_AP.*` — Primary and stress-condition mean paired AP differences and network-stratified 95% CIs.
- `figures/F02_GCEV_VS_FIXED_AP.*` — GCEV versus fixed-weight baselines under four conditions.
- `figures/F03_GENERATOR_HETEROGENEITY.*` — Descriptive linear/nonlinear subgroup distribution; not a preregistered subgroup superiority test.
- `per_network_method_metrics.csv` — All 2,112 method rows; `paired_network_effects.csv` — 192 network/condition paired contrasts; `network_bootstrap_95.csv` — prespecified uncertainty; `method_means.csv` — method comparison; `network_execution_ledger.json` — all 48 intended network IDs and source hashes.
- `independent_ranking_audit.json` — full independently recomputed metric and file-integrity audit.
- Frozen source: `configs/GCEV_BENCHMARK_PREREG_2026_11_25.json` and its hash.
