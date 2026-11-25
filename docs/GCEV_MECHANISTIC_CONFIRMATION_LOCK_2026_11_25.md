# LOCKED mechanistic confirmation: paired graph interactions, not post-hoc slicing

Protocol **GCEV-MECH-LCK-2026-11-25-001**, frozen 2026-11-25 **before running any confirmation seed in this design**. Predecessor exploratory seed sets 11/42/77 and 2001–2048 have already been inspected; neither may be used as independent confirmation.

## Question
Did the positive GCEV-minus-graph-independent AP difference on previously examined nonlinear networks arise because **pairwise regulatory interactions** make graph-state features more useful? Separate saturation from multiregulator interaction by *within-DAG factorial intervention*.

## Independent units and frozen seeds
32 independently drawn directed acyclic graphs with RNG seeds **6101–6132**, none in earlier runs. For index `i=seed−6101`, `n_genes=7` if `i%2==0` else 9; `n_cells=140` if 7 genes else 180; source reliability regime `target_switch` if `(i//2)%2==0` else `global_stable`. Exactly 8 DAGs in each of four size×evidence strata. Graph: predecessor-only edges, Bernoulli(0.30), weights ±Uniform[0.35,0.8], redraw until ≥2 true edges; exactly the same generated graph and sampled exogenous Gaussian noise across all four interventions on that DAG.

## Four paired expression-generating mechanisms
For target j, let P_j be the true parents. For binary saturation S and pairwise interaction I, compute:
```
x_j = sum_{i in P_j} w_ij * ((1-S)*x_i + S*1.1*tanh(x_i))
      + I*0.30*tanh(x_{p1}*x_{p2})*[|P_j|>=2] + eps_j
```
where (S,I) = (0,0) linear; (1,0) saturation only; (0,1) pairwise interaction only; (1,1) both. eps ~ standard Gaussian, sampled ONCE per DAG and reused; all four conditions use identical weights, DAG, source file scores, missingness mask, train/score/test cell IDs, and candidate edges. No new oracle edges are created by the interaction term: it requires pre-existing distinct true parents.

Two **synthetic, oracle-derived**, partially observed sources SIM_A and SIM_B are generated deterministically ONCE per DAG following the v0.2 source rules: 77% observed probability, 0.90/0.10 mean high/low scores with SD 0.07, target-switch reliable source by even/odd target or globally stable per-edge flip probabilities 0.12 and 0.25. They are copied without modification to all four mechanism datasets. These are NOT biological evidence.

## Evidence stress and matched models
For each of the 4 mechanisms run TWO scenarios: `clean` (0 missing, 0 corrupt) and `combined` (50% source dropout, 25% corruption), on the *same source files* with intervention-independent RNG seeds. Total exactly **32 DAGs × 4 mechanisms × 2 scenarios = 256 fits**. Common settings: ridge alpha 1; 3 inner folds; 2 source-blind scaffold parents; source boost 0.30; shrinkage 0.05; top K=7; graph budget=1600; 60/20/20 expression partition using DAG seed. Methods MUST include GCEV, `no_graph_state`, `fixed_sources`, `global_source`, `conditional` and `greedy`; Pearson and Spearman where original code provides them. No tuning on these 32 seeds or reference metrics.

## Preregistered primary comparison and decision
Define `dAP(S,I) = AP_GCEV(S,I) − AP_no_graph_state(S,I)`, for each DAG in the `clean` scenario.

The **primary paired factorial interaction effect** per DAG is
```
theta_DAG = 0.5 * [(dAP(0,1)−dAP(0,0)) + (dAP(1,1)−dAP(1,0))]
```
Estimate its mean giving **one equal weight to each DAG**. Stratified percentile bootstrap of 4,000 draws, RNG 20261125, re-sampling whole DAGs (not four observations/cells/edges) within the 4 design strata.

An *interaction-specific graph-state advantage* is SUPPORTED only if the predeclared theta 95% bootstrap CI lower bound is >0, its point mean ≥0.02 absolute AP, **and** the GCEV-vs-no-state relative untouched-test NMSE (averaged over the two I=1 clean conditions per DAG) has bootstrap 95% upper bound ≤+0.02. Other outcomes, including negative theta, must be reported. This rule is mechanistic simulation evidence, not proof of biology.

## Secondary (not independent confirmation) analyses
- Saturation-only contrast, I×S contrast, interaction-only vs linear and both vs linear. Report intervals **as descriptive; no multiple-testing-free superiority claims**.
- GCEV vs fixed evidence and expression-only, AP/AUROC/orientation/P@7, untouched test NMSE.
- Source-ablation **changed-edit rate**, mean absolute paired utility, nonzero valuation share, source disagreement, parent-target correlation, graph indegree/percentage of multiregulator targets, and nested loss evaluations/timing. These are diagnostic associations; especially source disagreement and graph-state summaries are not independent biological gold.
- Stress `combined` vs `clean` on same DAG. Distinguish the effect of corrupt source quality from the effect of the graph mechanism.
- Record all 256 outcomes and exclusions; 0 permissible cherry-picked replacements.
- Back-check that source file hash, graph reference hash, per-network split, and oracle DAG edge count are identical across paired interventions; gold is passed only to evaluation after method rankings are sealed.
- Full rankings sealed with original project's `fit_manifest.json` hashes; report the number of successes/failures, protocol hash, code commit and verification of aggregation.
- **No independent biology, no trained RL advantage, no direct causal TF-edge inference from observational signal**.

## Stop rule
The results of 6101–6132 are once-only confirmation. Do not change mechanism strengths, source regimes, metrics, reporting or excluded runs based on results. Any follow-up belongs to a new separate experiment with new seeds and a new freeze. A journal-quality biological claim additionally requires genuinely independent biological source data.

Follow the existing [48-network preregistered report](https://github.com/TshepoGobonamang/bio-rl-grn/pull/39) as historical, already-inspected evidence; this design is intentionally different and stricter.
