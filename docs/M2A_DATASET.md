# M2A: cross-fitted temporal-evidence reliability meta-dataset

**Status: implementation of dataset construction only.** The output does not contain biological edge labels, learned reliability weights or any RL policy.

## Data

Official BEELINE v4 `ExpressionData.csv` (genes × cells), `PseudoTime.csv` (one or more branch-specific columns with NA masks), and optionally an independent TF-list file. `GroundTruthNetwork.csv` is **not opened or required**.

The existing M1.1 seeded outer partition gives **60% fitting / 20% internal graph scoring / 20% untouched expression test**. Only fitting cells enter M2A. A shuffled inner 3-fold split uses `F_b` to generate meta-features and ridge fits, and `V_b` for independent temporal-consistency and predictive-utility targets. Sample ID hashes, source hashes, seed and fold configurations are recorded.

## Edge-by-fold schema

Each row describes one directed TF→target candidate in one fold. The manifest explicitly allowlists the **fit-only feature columns**: `z_signed, m_present, coverage, branch_n, bin_n, abs_z, expression_baseline, target_fit_variance, n_genes, scaffold_n, bootstrap_variance, bootstrap_sign_agreement`. Pseudotime margins are computed branch-by-branch; bins adapt from 16 to at least 8 when cell support is limited, with at least three observations per bin. Bootstrap descriptors use *only F_b*.

**Missing evidence remains missing**: no viable branch gives `m_present=0` and `z_signed=NaN`. Zero directional score is not equivalent to missing evidence.

Out-of-fold targets and diagnostic-only columns:

- `validation_z_signed_TARGET_ONLY` and `validation_available_TARGET_ONLY`: computed from `V_b`, **never model features**.
- `y_temporal_reproducibility_TARGET = (1+tanh(z_F*z_V/tau^2))/2`, when both source scores are available, otherwise NaN. A value ~0.5 indicates weak/ambiguous reproducibility, **not a true negative regulatory edge**.
- `y_expression_conditional_gain_TARGET`: relative held-out reduction in target ridge loss after adding candidate edge to a **source-blind** one-parent scaffold. The scaffold parent is the strongest alternative regulator from absolute training-fold correlation, explicitly excluding the candidate source. This target may be negative and does not certify causal direction.
- `expression_loss_before_TARGET_ONLY, expression_loss_after_TARGET_ONLY`: out-of-fold diagnostics, not input features.

Rows from different folds share samples/genes; reverse-edge pairs share temporal statistics. Do **not** treat these rows as independent biological replicates for confidence intervals.

## Usage

```bash
python scripts/extract_beeline_v4.py --archive BEELINE-data.zip --run dyn-LL-500-1 --output data/raw/beeline/dyn-LL-500-1
python scripts/build_m2_crossfit.py --raw-dir data/raw/beeline/dyn-LL-500-1 --output-dir results/m2a/dyn-LL-500-1 --dataset-id dyn-LL-500-1
```

The [M2A GitHub Actions workflow](../.github/workflows/m2a-crossfit.yml) builds three exploratory topology datasets. Artifacts contain only `m2_crossfit_edges.csv` and `m2_crossfit_manifest.json` in each topology folder; no expression matrices or gold references are redistributed.

## Safeguards and limitations

- No reference edges in feature generation, target generation or candidate TF selection.
- No access to outer graph-scoring/test cell *values* in meta-dataset construction.
- Temporal target is out-of-fold relative to feature estimation, though its lag statistic is computed within the validation fold's binned trajectory. An additional nested temporal holdout could strengthen this in later work.
- BoolODE synthetic pseudotime can reflect simulation time and may overestimate the usefulness of measured pseudotime.
- A reproducible directional association is **not** the biological truth of an edge; conditional expression gain is predictive, not causal.
- Next step is an empirical **meta-dataset quality audit** followed by an appropriately regularized model and matched fixed/nonlearned baselines, then independent topology and external biological testing.

## M2A data contract: separate model inputs from validation targets

The original `m2_crossfit_edges.csv` remains available **for auditing only**. A future reliability learner must NEVER select features by dropping columns from this combined file. The preferred exports are:

- **`m2_features.csv`**: directed-edge/fold identifiers and ONLY the 12 pre-declared fitting-only predictors in the manifest. It cannot contain columns ending in `_TARGET` or `_TARGET_ONLY`.
- **`m2_targets.csv`**: matching identifiers plus validation-only temporal/reconstruction targets and diagnostics. It cannot contain fit predictors.
- **`m2a_quality_summary.csv`** and **`m2a_quality_manifest.json`**: audit of missing masks, target saturation, feature/target key correspondence, row counts, no gold input, and inner-fold separation, generated across all network folders.

Run the audit after data generation:

```bash
python scripts/audit_m2_crossfit.py --root results/m2a
```

The audit fails on unexpected predictor columns, unaligned fold-edge IDs, missing evidence encoded as observed, absent targets incorrectly assigned nonmissing soft labels, or invalid 0–1 pseudo-target values.

**Evaluation requirement:** Reverse-direction pairs, folds and networks are not independent samples. Use grouped or topology-held-out evaluation of future reliability models. The soft reliability target measures **within-fold reproducibility of temporal association**, not the biological truth of a regulatory interaction.
