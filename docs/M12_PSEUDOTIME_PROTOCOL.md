# M1.2 — Pseudotime-informed GRN directionality pilot

## Research question

Can adding **pseudotime-derived directional association** to the unchanged graph-conditioned ridge objective improve correct TF-target orientation **without degrading average precision (AP) or untouched-test expression prediction**?

## Scientific caution: pseudotime is NOT true longitudinal observation

BEELINE's synthetic BoolODE generator may use actual simulation time as a pseudotime proxy. Such an ordering is potentially more informative than pseudotime estimated from noisy real scRNA-seq expression. Because different cells are sampled at different positions, time-sorted snapshots do not create time-coupled gene trajectories. Scribe (Qiu et al., Cell Systems 2020, DOI 10.1016/j.cels.2020.02.003) found severe causal-inference losses with pseudo-ordering compared to true temporal coupling. This stage must use the phrase **pseudotime-informed direction heuristic**, never causal discovery.

## Primary computation

For **each branch independently**, use **fitting cells only**:
1. Align the rows of PseudoTime.csv to ExpressionData.csv cell IDs, retaining NA branch masks.
2. Sort fit cells by branch-specific pseudotime; form 16 equal-count bins (at least 3 cells per bin).
3. Aggregate expression for each bin and standardize each gene's binned trajectory.
4. For candidate i → j, compare a ridge autoregression for gene j at bin b+1 (using j at b and a linear bin trend) against the same model also including i at b. The relative in-sample error reduction is D(i → j), clipped to [0,1].
5. Define signed margin **d(i,j) = D(i→j) − D(j→i)**; therefore d(j,i) = −d(i,j). Average margins across available pseudotime branches, weighted by the branch's fit-cell count.
6. During graph construction maximize **incremental expression utility + 0.015 × d(i,j)** using the original ridge scorer and identical 20-edit limit.

The 0.015 temporal coefficient, 16 bins, lag 1 and ridge alpha 1.0 are fixed before any real-run test results are examined. **No claim of optimal tuning is made.** The temporal evidence estimator is in-sample within fitting cells, so a shuffled-label negative control is essential.

## Four direct comparators

- **greedy**: unchanged M1.1 expression-only graph-conditioned greedy.
- **temporal_greedy**: same greedy expression scoring + pseudotime margin.
- **shuffled_pseudotime**: same weight, but fit-cell time labels independently permuted *within each branch*.
- **temporal_only**: ranking using signed temporal evidence alone, without graph-conditional expression gain.

All use the same candidate TF-target universe (all measured genes; no self-loops), fit/scoring/test partitions, maximum edits, seed and reference labels. Gold is loaded only after rankings are written and untouched-test expression losses computed.

## Datasets

1. dyn-LL-500-1 — long linear synthetic network (previously assessed).
2. dyn-LI-500-1 — independent short linear topology.
3. dyn-BF-500-1 — bifurcating topology; handle branch-specific missing values.

Each comprises a different **synthetic generating structure**, but this is still not human biological external validation. Real pseudoTime may be estimated rather than generated and require further fold-local provenance auditing.

## Report the following for each topology

AP (noninterpolated sklearn average precision), trapezoidal PR-AUC separately, AUROC, Precision@K, direction accuracy (ties count 0.5), top-K reciprocal pairs, selected-graph held-out expression MSE, and reference-exclusion audits. Include an explicit shuffled pseudotime comparison; analyze temporal-only rankings as a diagnostic, not a complete GRN method.

Because graph sizes differ, avoid taking unweighted averages of AP across networks without also reporting baseline prevalence and chance-normalized performance. Report raw per-topology results first.

## Important limitations

- Pseudotime **scale, branch structure, root selection, ordering errors** influence evidence; in this first experiment the supplied trajectory is used as-is and **not** independently estimated.
- Bin-averaged expression from separate cells is not a real temporal time series; lagged association may be driven by shared differentiation trends.
- The autoregressive gain and graph search both use development data; only the third split provides untouched expression evaluation.
- The candidate universe permits both directions. A feedback loop is not automatically biologically incorrect.
- AP on a simulated gold network is **not equivalent** to causal edge recovery in real biology.
- The graph-utility penalty is fixed for the pilot and should not be tuned on the same reference labels.
- Follow-up studies should use inferred pseudotime, RNA velocity/perturbation-based support, and multiple independent network seeds, with matched-information comparators.

## Execution

The workflow in [GitHub Actions](../.github/workflows/m12-beeline.yml) is configured to run **once on first merge**, and manually when dispatched later. It retrieves the official BEELINE v4 archive and runs the three synthetic topologies without committing raw data.

Run locally:
```bash
python scripts/extract_beeline_v4.py --archive BEELINE-data.zip --run dyn-LL-500-1 --output data/raw/beeline/dyn-LL-500-1
python scripts/run_m12.py --raw-dir data/raw/beeline/dyn-LL-500-1 --processed-dir data/processed/m12/dyn-LL-500-1 --results-dir results/m12/dyn-LL-500-1 --bins 16 --temporal-weight 0.015
```

**Interpretation gate:** only advance to context-adaptive biological evidence weighting after inspecting both the direction metrics and ranking/independent-test metrics, including potential **negative** findings.
