# GCEV explicit parent-state ablation — reproduction

**Protocol:** `GCEV-PARENT-LOCK-2026-11-25-001`.

## Freeze evidence
- GitHub preregistration narrative commit: `b77f235a2efd0d10a100413c5bf8ddc59a72cdd1`.
- GitHub completed executable implementation and JSON: `bf185df34efbd0262c0ff3c2c779dde2265234ca`; committed before confirmation seed evaluation.
- Protocol JSON SHA256: `04210a82fa9b6cbed7300806c32e244f91b2ce857ac75ad7bcd82c9ba6f56a25`.
- Confirmation seeds 7201–7232. Previous pilot/48-net/32-DAG mechanism seeds excluded.

## Commands
```bash
python -m pip install -e '.[dev]'
PYTHONPATH=src:. pytest -q
PYTHONPATH=src:. python scripts/gcev_parent_state_confirm.py \
  --protocol configs/GCEV_PARENT_STATE_PREREG_2026_11_25.json \
  --output results/gcev_parent_state_32_dags --workers 4
PYTHONPATH=src:. python scripts/verify_gcev_parent_state.py \
  --root results/gcev_parent_state_32_dags
```

The separate GitHub workflow `.github/workflows/gcev-parent-state-ablation.yml` performs this comparison and publishes an Actions artifact. GitHub execution status must be established from its actual run; local results and the independent audit are complete.

## Data contract
This is a deterministic simulator with 32 DAGs × four mechanisms × two source conditions. The two synthetic evidence streams were **generated using oracle regulatory truth**, never actual independent biological data. Training has only fit expression and prior scores, while evaluation loads the full synthetic DAG reference after ranking hashes are saved. Independence is by DAG realization. The mechanism conditions share graphs, noise and priors; do not count them as independent units.

## Outputs
`parent_state_decision.json` gives the preregistered primary outcome; `per_dag_parent_effects.csv` is the 32-DAG statistical unit dataset; `per_dag_mechanism_comparisons.csv` is descriptive; `all_method_metrics.csv` has all 3,328 method comparisons. `independent_integrity_audit.json` confirms 3,328 sealed rankings rescored with no errors. Source file checksum manifest `SOURCE_CHECKSUMS.json` verifies the local result files. The complete independent ranking artifacts are in the provided archive, not automatically committed to Git history.

**Interpretation:** The ablation isolates only the explicit parent-state predictors. Action features still depend on the same scaffold; parent count is fixed at two throughout; the primary superiority criterion is not met. No biological independent validation or causal TF binding claim.
