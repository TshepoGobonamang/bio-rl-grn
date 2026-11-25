# Reproducing the GCEV preregistered 48-network synthetic confirmation

**Scientific status:** Controlled synthetic methods benchmark; primary incremental AP superiority NOT MET. No biologically independent source evidence and no RL experiment.

## Frozen protocol / code

- Protocol GitHub precommit: `cb7b3d2aff6827c98c370ba14fbc19293405a1c3`
- Code/test/workflow publication before evaluation: `c4d98e9b48c5fd99dc0d630d39003ac2d2f559e8`
- SHA256 of exact `configs/GCEV_BENCHMARK_PREREG_2026_11_25.json`:
  `556c6745948cc89b10c9b7ee92a64d4b5637bcefbe180c0c062082b7d2f0272c`
- Graph seeds 2001–2048; development seeds 11,42,77 intentionally excluded.
- No biological, human or proprietary data are required. Every source is explicitly simulated using oracle network truth.

## Commands

From the project root with Python 3.11:

```bash
python -m pip install -e '.[dev]'
python -m pytest -q

# Use exactly the locked file, not an edited protocol.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python scripts/confirmatory_gcev.py \
  --protocol configs/GCEV_BENCHMARK_PREREG_2026_11_25.json \
  --output results/confirmatory_synthetic_only \
  --workers 4

# Recompute all sealed ranking AP/AUROC and verify all manifests.
python scripts/verify_confirmatory_results.py \
  --root results/confirmatory_synthetic_only

# Create the descriptive figures without changing analysis rules.
python scripts/plot_confirmatory_gcev.py \
  --results results/confirmatory_synthetic_only \
  --output results/confirmatory_synthetic_only/figures
```

For Windows PowerShell replace `export` with `$env:OPENBLAS_NUM_THREADS="1"` and likewise for the other two variables.

The workflow `.github/workflows/gcev-preregistered-benchmark.yml` includes the same lock and runs unit tests before computing source-conditional outputs. It uploads compact summaries; the provided full offline ZIP contains every generated source, sealed ranking and evaluation manifest.

## Definition of the independent unit

A simulated graph/network seed is one independent analysis unit. Four evidence conditions are *repeated observations on that same graph*, not independent networks. Each of the eight graph-family×gene-count×source-regime strata contains six independent DAGs. Bootstrap resamples these graph units within strata.

## Expected outputs and integrity contracts

- `primary_decision.json` — **PRIMARY_CRITERION_NOT_MET**.
- `per_network_method_metrics.csv` — 192 × 11 = **2,112** rows.
- `paired_network_effects.csv` — 192 paired graph×condition rows.
- `network_bootstrap_95.csv` — stratified network bootstrap intervals.
- `method_means.csv` — mean per-network AP, AUROC and NMSE by condition.
- `network_execution_ledger.json` — all 48 intended network IDs and source hashes; no missing seed.
- `independent_ranking_audit.json` — 192 complete source evaluations, 2,112 recalculated AP/AUROC and zero errors.
- `network_<seed>/data/` — purely generated expression, TF list, synthetic oracles and fake source CSVs (do not relabel as biological).
- `network_<seed>/conditions/<condition>/fit/` — sealed method rankings and fit manifest (gold-free estimator path).
- `network_<seed>/conditions/<condition>/evaluation/` — reference-access records and scoring metrics.

The numeric values reported in `reports/GCEV_REGISTERED_CONFIRMATION_48_NETWORKS_2026_11_25.md` originate from the frozen protocol and complete output table. Different CPU hardware affects runtime but not intended ranking metrics beyond numerical tolerance. The synthetic benchmark is not evidence of direct causal gene regulation or real-source biological validation.
