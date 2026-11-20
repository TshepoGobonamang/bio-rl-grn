# Working manuscript blueprint (claim-neutral, not a completed paper)

**Proposed manuscript title:** *Graph-Conditioned Reliability Signals for Gene Regulatory Network Inference: A Leakage-Aware Experimental Benchmark.*

**Venue:** not fixed; CIBB 2026 dates/submission are past. Select a current 2027+ venue only after reading page/word limit and research scope.

## Abstract contract (write only after final test run)

1. **Problem:** Directed GRN inference from expression alone is confounded by indirect associations; external prior evidence may be unavailable, overlapping or wrong.
2. **Method:** Cross-fitted parent-conditional predictive usefulness, explicit source-availability masks, graph-based construction, static and adaptive controls; optional budget-matched RL *only if executed*.
3. **Evaluation:** Pinned BEELINE/supplementary source benchmark with frozen candidate universe, untouched test expression and evaluation-only references; missingness/corruption analysis where evidence exists.
4. **Findings:** Report the observed AP, orientation, held-out loss and uncertainty against strongest baselines; negative is valid. Do not insert values from one exploratory realization as if independent confirmation.
5. **Conclusion:** State exactly which benefits were observed, which failed, and which claims (causal validation, direct binding, RL superiority) were not established.

## Proposed sections

### 1. Introduction
Directed regulatory graphs support mechanistic hypotheses but observational expression and pseudotime do not uniquely identify regulatory direction. Gene-list/TF eligibility, incomplete biological priors, and benchmark/reference overlap threaten fair evaluation. Motivate a method that tests *incremental usefulness*, not merely source score magnitude.

Contribution claims must be derived from experiments, not from implementation:
- [methods] source-blind cross-fitting for graph-parent-conditional predictive gain;
- [methods] availability-aware source evidence and explicit independent provenance constraints;
- [benchmark] fair graph-budget comparisons with untouched expression tests and evaluation-only directed references;
- [optional] comparative RL search only if actually benchmarked.

### 2. Methodology

Define expression matrix `X` with cells × genes; frozen candidate directed graph universe `E`, source `s`, target `t`, source-blind parent scaffold `S_t^(-s)`. For inner fold b:
```
L0 = MSE_{V_b}( ridge(X_{F_b, S_t^(-s)}, X_{F_b,t}) )
L1 = MSE_{V_b}( ridge(X_{F_b, S_t^(-s) union {s}}, X_{F_b,t}) )
Delta_{s->t,b} = clip((L0-L1)/max(L0,eps),-1,1)
```
State parent-selection method, feature scaling, ridge penalty, fold partition and preprocessing (fit-only). A positive Delta is out-of-fold association/predictive utility, not causal TF binding or regulatory direction. No use of evaluation reference or CollecTRI overlap to construct a held-out truth label.

Set reliability model `r(e,S,m)` from source-specific score `S`, observed missingness mask `m` and training-only features. Compare it with a **fixed weight**, global/source-only and nonlearned-stability heuristic. Distinguish source absent vs source evidence explicitly negative. Retain sign of curated regulation as activation/inhibition, **not TF-target arrow reversal**. Do not claim `r` is P(edge=true).

Define graph objective `J(G)` from fit-only model and internal score split; freeze graph edits/budget. The outer untouched test split is not available to edit selection, hyperparameter selection or policy reward computation.

If RL is evaluated, define Markov state, legal add/remove/stop actions, actual utility increment reward, matching budgets, training environment and comparator greedy/random policies. If no RL evaluation, remove this subsection from final manuscript.

### 3. Datasets and provenance
Table with original source ID / version / checksum, organism/cell context, number of cells/genes/regulators, independently observed ground truth, prior source and provenance overlap, train/internal/test split, gold access stage. Separately report exploratory real hESC CollecTRI integration and GSE216595 *descriptive* TFmap diagnostics; do not count these as independent accuracy validations until qualified.

Known GSE216595: 60,997 × 5,047, 87 TF labels. `batch` has 8 values but no verified culture replicates. TFmap source was seen with likely headerless format; no supplied exact join outputs.

### 4. Results (only populated from final workflow)
R1 — leakage / split integrity and reference eligibility.
R2 — graph-parent-conditional signal vs marginal/temporal signal, sign and rank stability, against fixed baseline.
R3 — AP, AUROC, orientation, early precision and untouched expression test loss by network, graph/realization-level uncertainty.
R4 — source missingness, corruption, conflicting priors and ablation; compute/source coverage.
R5 — graph utility trajectory, search budget, runtime; RL only if actually run.
R6 — external biological plausibility, under carefully disclosed source/reference provenance; no unsupported direct causal edge assertions.

### 5. Discussion
Interpret positive, null and negative tradeoffs. A stability signal can be reproducible but not biologically directional. Parent-conditional prediction can exploit co-expression and indirect effects. A curated prior can share underlying studies with reference edges. Unequal candidate universes or query budgets invalidate headline comparisons. The potential cell-line endpoint shift HUES66→iNP also constrains interpretations.

### 6. Reproducibility and limitations
Include code tag, Python/environment, hardware, source URLs/MD5/SHA256, seeds, fold/hash manifest, frozen candidate selection, reference access audit and exact complete rerun commands. Acknowledge benchmark heterogeneity and limited independent biological validation.

## Required figures/tables and their data contracts

| Figure/table | Minimum underlying artifact | No-artifact fallback |
|---|---|---|
| F1: graph-condition / evidence flow | frozen methods source + code diagram | Method figure only; no empirical claims |
| T1: source and candidate universe | manifest + eligible edge counts | Descriptive source table with unsupported entries marked |
| T2: baseline accuracy and generalization | `per_network.csv`, `results.csv` | Withhold paper-level performance claims |
| F2: AP and orientation across networks | `per_network.csv`, `uncertainty.csv` | No pooled error bars |
| F3: missingness/corruption sensitivity | `ablations.csv` and actual evidence-source masks | Omit modality-specific claims |
| F4: utility vs query budget | `trajectories.csv` and matched budgets | Omit RL/efficiency superiority |
| T3: external source provenance and caveats | original assays, D11/D12/D15B artifacts | Report NO-GO and provenance limitations |

## Minimum review checklist

Every row in a results table is reproducible from a frozen artifact; all metrics have unambiguous evaluation definitions; no candidate gold leakage; no independent-culture claim based on a sequencing batch; no training on held-out reference; no figure uses the unpublished draft branch's synthetic fixtures as external biology. No invented citation/acceptance/accuracy value.

**Sources:** BEELINE Pratapa et al. (Nature Methods 2020) https://doi.org/10.1038/s41592-019-0690-6 ; hESC prior and perturbation provenance from the repo's existing M2D2/M2D14 reports.
