# Rerunning the matched-DAG GCEV mechanism experiment

## Exactly frozen input and source lineage

- Precommitted scientific lock: [`GCEV_MECHANISTIC_CONFIRMATION_LOCK_2026_11_25.md`](https://github.com/TshepoGobonamang/bio-rl-grn/blob/research/gcev-graph-mechanism-confirmation-2026-11-25/docs/GCEV_MECHANISTIC_CONFIRMATION_LOCK_2026_11_25.md), first commit `f446a3d8764148dbfc4a885f334004bb74318fba`.
- Protocol file: `configs/GCEV_MECHANISTIC_PREREG_2026_11_25.json`, SHA256 `a5b2cd5fc293bc220aec4df519dcb6ccb8e69bac7a390b23b73f10a665cff701`.
- Python code base: research branch `research/gcev-graph-mechanism-confirmation-2026-11-25`, PR #40 (branch represents v0.2 GCEV with additive v0.3 mechanistic runner). Use a pinned commit SHA from GitHub Actions for independent reproduction.
- Underlying original 48-network benchmark was run on seeds 2001–2048 and inspected before this experiment. **These values are not independent test results for a revised GCEV algorithm**.

## Environment and commands

Python 3.11 recommended. Python environment packages specified in the checked-out `pyproject.toml`.

```bash
python -m pip install -e '.[dev]'
python -m pytest -q
sha256sum configs/GCEV_MECHANISTIC_PREREG_2026_11_25.json
PYTHONPATH=src:. python scripts/gcev_mechanism_factorial.py \
  --protocol configs/GCEV_MECHANISTIC_PREREG_2026_11_25.json \
  --output results/gcev_mechanistic_32_dags --workers 4
PYTHONPATH=src:. python scripts/verify_gcev_mechanistic.py \
  --root results/gcev_mechanistic_32_dags
python -m pip install matplotlib
PYTHONPATH=src:. python scripts/plot_gcev_mechanism.py \
  --root results/gcev_mechanistic_32_dags
```

The run should process **32 DAGs × 4 mechanisms × 2 evidence conditions = 256 conditions**, generating 2,816 method rankings; failed networks remain in an explicit ledger and invalidate the aggregate. The figure output is descriptive: there is no independent biological evidence in the simulated source CSVs.

## Primary numerical checkpoint

- `mean(theta_interaction)` = **-0.002080400791375647**.
- 95% four-stratum bootstrap interval **[-0.0163692746174306, 0.011090167298727343]**.
- Mean GCEV-vs-no-graph relative NMSE under interaction conditions **-0.006975687432682248**; bootstrap 95% **[-0.012341994933931469, -0.001930327573126197]**.
- Primary decision: `PRIMARY_CRITERION_NOT_MET`.
- Original local run wall time **88.31 seconds** on its available environment, **not a normalized method-to-method runtime benchmark**.
- Standalone audit re-scored **2,816 rankings**, all consistent with original AP/AUROC metrics within 1e-11; no 32-DAG failure.

## Source files to preserve

`protocol.json`, `execution_ledger.json`, `mechanism_decision.json`, `all_method_metrics.csv`, `per_dag_factorial_effects.csv`, `per_dag_mechanism_comparisons.csv`, `per_condition_source_decisions.csv`, `independent_integrity_audit.json`, all per-network fit/evaluation manifests and ranked-edge CSVs, and generated plots.

## Scientific interpretation boundary

These results are methodological simulations with artificially accurate, oracle-derived priors. They cannot validate TF binding, independent experimental controls or GRN causal direction. Importantly, the `no_graph_state` baseline retains `candidate_count`, whereas the removed feature group includes source-prior action contrasts. The primary contrast isolates **the six-feature graph/action-context bundle**, not a pure topological-parent-only feature component. A separate stricter ablation and external biological data are needed for a top-tier biological methods paper.
