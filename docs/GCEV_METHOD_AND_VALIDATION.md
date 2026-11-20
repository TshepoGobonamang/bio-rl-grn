# Graph-Conditional Evidence Valuation (GCEV) — v0.2 experimental contract

## What is new **in this software prototype**

Instead of ranking an edge by its held-out conditional expression gain alone, GCEV **estimates the incremental predictive value of allowing one evidence source to influence a fixed-budget graph-edit choice**, conditioned on the existing source-blind parent state. The paired source-ablation observations are training-only labels for a subsequent, target-held-out evidence-utility model.

This is a **proposed one-step method**, not a validated new class of causal GRN inference, a calibrated probability of a true edge, or established algorithmic novelty over every prior method. The 2024 prior-integration optimization work and 2025 prior-knowledge review must be discussed in any eventual literature review. The 2026 Genes2Reg study likewise addresses adaptive fusion of incomplete biological priors, using a different learning objective.

## Nested split, leakage contract and equations

The full experiment first separates the cell set into 60% outer fit, 20% internal graph scoring and 20% untouched expression test. Reference edges are **not** accepted by the `fit` command. Within the outer-fit cells only:

1. A 3-fold partition divides cells into inner train `T` and validation `V`.
2. The inner train is further divided into a source-blind graph-search fitting subset `A` and internal selection subset `B`.
3. For each target gene `t`, choose the source-blind graph parent set `P_t` from `A` only. Exclude its members from the available next-addition candidate set.
4. Rank eligible next edges from train-only conditional expression gains measured on `B`, normalized across eligible candidate actions. For source `k`, hold all *other* source scores fixed. Select exactly **one** candidate with source k omitted (`e_-`) or included (`e_+`), using the same candidate universe and graph-edit count.
5. Refit graphs `P_t ∪ {e_-}` and `P_t ∪ {e_+}` on the *same* inner train `T`; compare expression loss for target `t` on the *same* unseen inner validation cells `V`:

\[
\Delta_{k,t,P_t} = \operatorname{clip}\left(\frac{L_V(P_t\cup\{e_-\})-L_V(P_t\cup\{e_+\})}{\max(L_V(P_t\cup\{e_-\}),\varepsilon)},-1,1\right).
\]

Both choices have precisely one added parent and identical parent-state construction. A negative value means that including evidence source `k` *harmed* predictive generalization of this chosen graph edit. If the same edit wins both comparisons, the paired difference is zero. This is not a causal intervention or observation of a time-series transition.

6. Across three folds, learn per-(source,target) utility using source coverage and strength plus **graph-state features** (parent count, parent-target correlation, search gain/margin, selected-action evidence contrast). Predicted gates are **out-of-fold by target gene**, so the labels for a target never train the predictor used for that target.
7. Convert predicted evidence value to a bounded heuristic weighting: `w=clip(predicted_value/0.05,0,1)`; this threshold is a **development default, not calibrated**. A source with negative predicted value receives zero weight. Missing evidence is masked rather than assigned a biological negative.
8. Freeze the resulting rankings and hashes. The outer score/test and reference edges have not been accessed during training the GCEV weights. The original model uses internal scoring cells to rank candidate edits, and **only the untouched test cells** are used for post-ranking expression generalization.

## Fair controls and limitations

- `conditional`: graph-conditioned expression score, no priors.
- `fixed_sources`: constant evidence contribution without learned utility.
- `global_source`: target-OFF-fold learned global weight per source (graph independent).
- `no_graph_state`: same GroupKFold model/training labels and sources, but graph-state features removed.
- `gcev`: full source and graph-state feature model, grouped out-of-fold by target.
- Existing Pearson, Spearman, greedy, fixed, heuristic and adaptive baselines are retained.

The paired **graph edits** match one another in candidate universe, one-step action count and validation partition. **Compute budget does not match** between GCEV and fixed weighting: GCEV requires additional inner ridge fits and gate training. The number of nested loss evaluations and elapsed training time are explicitly recorded; the performance comparison must include that cost. A fair RL study would require independently matched whole-trajectory utility-query budgets.

The estimator currently values **one source per target and one next graph edit**, not an arbitrary multi-step graph policy or fine-grained calibrated per-edge source trust. The graph state consists of a selected parent scaffold, not all global interactions; the entire graph dynamic remains future research. Each graph's target-specific reward is a predictive *association* with no direct TF-binding guarantee.

## Current actual execution

The controlled engineering pilot produces **synthetic DAG expression and two gold-derived synthetic evidence streams** with different target regimes, then fits models **without passing the evaluation reference into the estimator**. This is useful for correctness/method stress tests but scientifically cannot establish independent biological evidence, genuine source provenance, regulatory causal direction, or superiority on real data.

Reproduce locally:

```bash
python -m pip install -e '.[dev]'
pytest -q
python scripts/controlled_synthetic_pilot.py --output results/SYNTHETIC_ONLY_gcev_pilot --seeds 11 42 77 --cells 130 --genes 7
```

For real sources, use `grnb fit --gcev --prior PRIOR_A.csv --prior-secondary PRIOR_B.csv ...`. **No modality can be declared independent of the test reference without original data/assay accession provenance review.** Do not use `--reference-complete` for an incomplete positive-only human regulatory catalog.

## Prespecified next genuine scientific evaluation

Select genuinely uninspected network realizations/source configurations, freeze the model, source strength/coverage, grouping, hyperparameters, candidate universe, seeds, and outer train/internal/test partitions before opening gold references. Evaluate AP, orientation, AUROC and untouched test loss with network-realization-level uncertainty, source corruption/dropout stress tests and compute cost. Biological evaluation requires source-disjoint assay/reference provenance and real evidence modalities.

**Prior-art context to review:**
- *Optimizing data integration improves gene regulatory network inference in Arabidopsis thaliana* (Bioinformatics, 2024): https://doi.org/10.1093/bioinformatics/btae415
- *Leveraging prior knowledge to infer gene regulatory networks from single-cell RNA-sequencing data* (2025): https://doi.org/10.1038/s44320-025-00088-3
- *Genes2Reg* (BMC Bioinformatics, 2026): https://doi.org/10.1186/s12859-026-06655-1
