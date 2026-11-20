# Locked evaluation protocol — research design v0.1

## Hypothesis

Can source-blind, cross-fitted parent-conditional predictive usefulness, combined with explicitly available external regulatory evidence, improve held-out directed edge ranking and expression test error compared with an equal-information fixed evidence weight and nonlearned source-stability heuristic?

A null/negative result is valid. The software does not estimate causal effects, and the test-set reference must never select priors, model weights, candidate edges, genes or search budgets.

## Development and final confirmation

Develop on already inspected older BEELINE topologies and smoke data only. Select **previously uninspected** independent network realizations/cohorts *before* viewing their reference edges. Record source accession/commit, checksum, topology, generator/realization ID, expression split seed, biological vs technical replicate structure and gene eligibility. Split fitting 60%, internal graph scoring 20%, untouched expression test 20%; save split permutation SHA256. Use complete benchmark reference only in evaluation.

Do not claim independent observations from repeated split seeds, edge-fold records, sequencing libraries or pseudo-replicated cells. Compute uncertainty at independently generated **network realization** level and cluster by topology family as sensitivity analysis.

## Method definitions

For proposed edge s→t, select a source-blind scaffold of up to two other potential parents using fitting-only target correlations. Fit baseline ridge on fitting cells and score on independent internal scoring cells. Fit augmented ridge including s on the same fitting cells; let Δ=(L_before − L_after)/max(L_before, ε), clipped to [−1, 1]. This is a *predictive association target*. Inner 3-fold cross-fitting uses only the 60% fitting partition for adaptive target generation.

Adaptive regression uses TF-target feature correlations, prior masks and prior strengths; final reliability predictions are out-of-target via GroupKFold. It has no regulatory gold supervision and is not calibrated P(edge). Source absence must never become an asserted negative. Fixed and nonlearned heuristics use exactly the same prior entries and candidate universe.

Greedy graph editing adds the largest **positive** internal graph utility improvement, subject to utility-query budget. Output ranking fills unaccepted positions using a declared conditional-gain fallback; `selected_edges` records actually accepted edits, and full-ranking AP must not be described as all-greedy edited edges. Match model settings / top-K for untouched test loss comparisons.

## Minimum baseline matrix

Pearson, Spearman, conditional predictive-gain, greedy ridge, fixed weighted evidence, nonlearned evidence heuristic, group-cross-fitted adaptive reliability. **Fixed/adaptive/heuristic are omitted rather than imputed if no real prior source exists**. Prespecify reliability weight and ridge penalty; do not tune on final labels.

## Outcomes

Primary: sklearn **average precision (AP)** of held-out gold rankings on a verified complete candidate reference; orientation accuracy for unique unambiguous one-direction reference pairs; untouched expression test NMSE. Secondary: AUROC, Precision@K, source coverage, accepted greedy edits, utility queries and runtime. Pearson/Spearman tie directions must be counted as unassessable, not wins.

Missingness matrix: 0%,25%,50%,75%; source-strength corruption: 0%,10%,25%,50%, using genuine available source modalities. These perturbations do not convert absent evidence to a negative. Preserve control, shuffled evidence and matched availability experiments in the final study.

## Control of information access

The `fit` command has no reference argument. It records SHA256s for expression, whitelist, provided prior and every ranking. The `evaluate` command checks ranking hashes, candidate coverage, and demands `--reference-complete`. A clinical/real positive-only reference cannot be evaluated as if missing edges were false. Never use final gold labels to choose hyperparameters or candidate genes.

## Completion rule

One source-pinned CI/integration run on frozen independent benchmark families, paired comparisons, per-network metrics, network-level uncertainty and reproduced result tables. Do not claim independent biological validation or RL superiority until their specific experimental requirements pass. A final paper can report negative results.
