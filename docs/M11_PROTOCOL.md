# M1.1 Protocol — Untouched Expression Test & Directionality Audit

## Status and scientific scope

Exploratory successor to [E1 (run 43894848848)](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/43894848848). Not a confirmation that the learned RL agent works (none is implemented yet), nor a reliability-weighting experiment.

## Why M1.1 exists

The E1 greedy baseline had AP 0.4353 and AUROC 0.7706, versus Pearson AP 0.3430 and AUROC 0.7789. The top-20 greedy predictions included eight reciprocal pairs. Greedy graph search repeatedly queried the same 125-cell internal scoring partition, preventing claims of independent expression generalization. E1 scored 18 eligible gold edges, while the benchmark documentation lists 19 total.

## Protocol

1. Verify the original BEELINE v4 archive and extract specified dyn-LL-500 replicates.
2. Convert expression orientation (genes × cells → cells × genes). The TF universe is derived from independent annotation when provided, otherwise all measured genes (disclosed as a limitation).
3. **Three partitions:** 60% fitting, 20% graph-edit scoring, 20% untouched expression test. Split is deterministic; all partition IDs are disjoint.
4. Compare Pearson, Spearman, unpenalized greedy ridge and an **experimental anti-reciprocity heuristic** with fixed reverse-edge penalty 0.01. These use the same fitting/scoring data; all gold labels are withheld.
5. Save complete candidate rankings and graph-edit trajectories.
6. Only after both greedy graphs have been selected, instantiate an independent expression test evaluator and compute null/selected graph loss. Do not use test results to change the model or select penalty.
7. Only after predictions exist, load reference edges, preserving raw row counts and exclusion categories (self-loops, absent genes, unannotated TFs, duplicates).
8. Evaluate AP (noninterpolated sklearn average precision), trapezoidal PR-AUC, AUROC, precision/recall@K, reciprocal pairs among top K, and direction concordance conditional on uniquely directed gold pairs. Tied orientation scores contribute 0.5.
9. Save evaluation-only labeled rankings and source reference filtering audit. Never feed labels to graph fitting or search.

## Running locally

```bash
python -m pip install -e '.[dev]'
python scripts/extract_beeline_v4.py --archive BEELINE-data.zip --run dyn-LL-500-1 --output data/raw/beeline/dyn-LL-500-1
python scripts/run_m11.py --raw-dir data/raw/beeline/dyn-LL-500-1 --processed-dir data/processed/m11-1 --results-dir results/m11/dyn-LL-500-1 --seed 42 --max-edges 20 --alpha 1.0 --reverse-penalty 0.01
python scripts/aggregate_m11.py --root results/m11
```

For GitHub, choose [Actions → M1.1 BEELINE independent-test and directionality audit](https://github.com/TshepoGobonamang/bio-rl-grn/actions/workflows/m11-beeline.yml) after merge. Trigger **Run workflow** with 3 exploratory replicates initially; 10 may be run once the first results are verified. The workflow downloads data on the runner and uploads results only.

## How to interpret results

- **Independent expression test:** graph search never reads this partition; reported null-to-selected loss difference is new evidence beyond the repeatedly queried scoring partition. The cells are from the *same simulation context*, so this is not independent biological validation.
- **Direction accuracy:** conditional on reference edges where the reverse edge is not also gold; counts ties explicitly. Positive orientation results are not proof of causal direction.
- **Reciprocal-edge penalty:** heuristic search restriction only. Reducing reciprocal predictions may mechanically improve the count while worsening overall accuracy. Compare both, with objective weights fixed in advance.
- **Independent replicates:** different dyn-LL-500 simulation replicates may share network topology. Their variation estimates expression/sample variability, not generalization to unseen network structures.
- **Ref count discrepancy:** inspect `reference_row_audit.csv` to reconcile 19 described versus 18 eligible edges in the E1 run.
- **Gold standard:** a missing simulation reference edge is a negative under that generator, but this assumption must not be imported into human GRN validation.
- **Review needed before M2/M3:** verify original archive schema, repeat over additional independent graph topologies, implement TF annotations from independent sources, and validate policy budgets.

## Mandatory outputs

- `m11_summary.json`: run parameters, partitions, separate expression-test losses, reference audit and SHA256 source hashes
- `m11_metrics.csv`: method-level precision-recall/AUROC and directional scores
- `reference_row_audit.csv`: every raw gold edge plus exclusion reason
- `ranked_*.csv`: inference-stage rankings, **unlabeled**
- `evaluation_annotated_*.csv`: evaluation-stage label annotations
- `greedy_trajectory.csv` and `reverse_trajectory.csv`
- `m11_replicate_metrics.csv` and `m11_descriptive_summary.csv` when using multiple replicates

**Data preservation:** never commit raw or processed transcriptomic matrices into the public GitHub repository. The workflow is manual because the official archive is large.
