# GCEV preregistered synthetic benchmark — locked before evaluation
**Protocol ID:** GCEV-LCK-2026-11-25-001  
**Freeze date:** 2026-11-25  
**Status:** No network from the following confirmation seed set has yet been evaluated in this project at the time of this protocol commit. Previous *development pilot* seeds 11, 42, and 77 are excluded.

## Primary, falsifiable question
Does adding graph-state features to a target-held-out, paired source-value learner produce higher directed-edge **noninterpolated average precision (AP)** than an otherwise identical learned source-value model with all graph-state features removed?

**Primary contrast:** GCEV minus `no_graph_state` on the **clean** source condition, paired per independent generated network, equal weighting across networks. Both use the **same per-fold target labels, source features, source values, outer splits, fitted ridge model and candidate universe**. The only intended difference is the graph-state feature block. Fixed and graph-independent adaptive controls are mandatory.

**Success rule, fixed before testing:** stratified network-level bootstrap 95% percentile CI lower bound for ΔAP must be greater than zero, **and** its 95% upper bound for mean per-network relative untouched-test NMSE difference (GCEV vs no-graph-state) must not exceed +2%. If either fails, the primary claim is NOT supported.

## Dataset / independent units: 48 generated DAG realizations
Exactly 48 distinct RNG seeds **2001 through 2048 inclusive**. Stratum allocation by `i = seed - 2001`:
- `family = linear` if `i % 2 == 0`, else `nonlinear`;
- `n_genes=7; n_cells=140` if `floor(i/2)%2 == 0`, else `n_genes=9; n_cells=180`;
- `evidence_regime = target_switch` if `floor(i/4)%2 == 0`, else `global_stable`.
Exactly 6 networks in each of 8 family × size × evidence strata. No cherry-picking of generated networks; all failures and declared exclusions logged, no posthoc replacement after reference results.

Both graph families are acyclic, with predecessors `i<j`, edge Bernoulli(0.30), sign equally likely ±1, coefficient absolute value Uniform[0.35,0.80], Gaussian noise of SD 1. Graphs with fewer than two genuine regulatory edges are re-drawn **before simulating expression or priors**, with retry count logged. Linear: `x_j = sum_i w_ij*x_i + epsilon_j`. Nonlinear: `x_j=1.1 sum_i w_ij*tanh(x_i) + 0.3*tanh(x_first*x_second)` if at least two parents, plus independent Gaussian noise. Generator code is fixed before evaluating confirmatory seeds.

## Synthetic evidence, explicitly NOT biology
Two manufactured source streams SIM_A and SIM_B, each edge observed with probability 0.77, scores 0.90 for oracle-positive and 0.10 for oracle-negative (Gaussian SD 0.07, clipped 0–1).
- `target_switch`: SIM_A is reliable on even-index targets, SIM_B on odd-index targets; the wrong source inverts the oracle edge label.
- `global_stable`: SIM_A independently flips the oracle label with probability 0.12, SIM_B with probability 0.25, without target-specific source assignment.

**All source scores are derived from the simulation truth, and can NEVER be presented as real independent biological evidence.** The simulator generates source CSVs before fit; the estimation code receives source CSVs and expression, not the oracle edges. Gold is opened for ranking metrics only after ranking outputs/hashes are frozen. This test distinguishes computational ablations, not biological causal validation.

## Fully paired conditions per network — 192 fit-and-evaluate jobs
`clean=(missing0, corrupt0)`, `dropout50=(missing0.50, corrupt0)`, `corrupt25=(missing0, corrupt0.25)`, `combined=(missing0.50, corrupt0.25)`. Exactly the same network realization, expression splits, uncorrupted source file and baseline hyperparameters in all four. Missing evidence remains masked; corruption inverts retained source scores as a declared test manipulation.

Locked settings: fit/internal score/final test 60/20/20%; ridge alpha 1; 3 nested source-valuation folds; 2 source-blind parents; source boost 0.3; shrinkage 0.05; top K=7; graph utility query budget 1600; seed equals network seed for partitions. Do not tune on these outcomes.

Required methods: GCEV, `no_graph_state`, target-held-out `global_source`, `fixed_sources`, expression-only conditional, Pearson, Spearman, greedy, plus all code-produced comparison baselines. All methods use identical candidate universe and source files per network. Greedy uses its existing graph-edit utility; it is NOT a compute-budget matched RL baseline, so RL comparisons are forbidden.

## Statistics and artifacts
Primary average ΔAP with 4,000-replicate percentile **stratified bootstrap by 8 generator×size×evidence strata, sampling entire independent network realizations** (not edge/target/fold/cell as independent samples), bootstrap RNG seed **20261125**. Mean per-network relative untouched NMSE difference with same draws. Report paired ΔAP against fixed and global as secondary, plus AUROC, orientation on assessable directed reference pairs, P@7, method-specific NMSE, changes to next graph edit, nested loss evaluations, and actual GCEV training time. For source corruption, compare each changed condition with the corresponding network's clean condition.

Required files: exact `PREREG` JSON and SHA, per-network-per-condition method metrics, network-paired differences, stratum-level effects, 95% intervals, failure ledger, code commit SHA and frozen ranking hashes. No p-values using correlated edge or cell records. No post-hoc selection of favorable seeds, scenario subsets or baselines.

## Scientific boundary
This is a **moderately broad controlled in-silico methods benchmark** with two synthetic expression families, two source regimes and 48 independently generated DAGs. It is NOT broad independent biological validation, not source-disjoint ChIP evidence and not a demonstration of RL superiority. A top-tier computational biology claim would also require external unseen biological datasets and relevant methods comparison. Negative and mixed outcomes are final results under this protocol, not a trigger to rewrite it.

**Freeze is effective with this GitHub commit BEFORE executing the 48-seed benchmark.** Code from v0.2 has already been evaluated on development seeds 11,42,77; all confirmation seeds listed above remain unused until this lock.
