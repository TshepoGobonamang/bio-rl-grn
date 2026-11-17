# M2C-A — Direction-sensitive cross-fitted temporal usefulness targets

**Status:** target construction **only**. No adaptive reliability model, biological truth classification or RL agent is trained by this milestone.

## Why change the M2A target?

The previous M2A reproducibility target `(1+tanh(z_train*z_valid/tau²))/2` is invariant if the edge `i→j` is reversed, because both signed margins reverse together. The official [M2B audit](../reports/M2B_adaptive_reliability_decision_2026_11_17.md) shows that training on this symmetric target does not consistently beat fixed temporal evidence.

## New fold-held-out directional surrogate

For candidate `s→t`, independently sort fitting-fold **F_b** and validation-fold **V_b** cells by each observed pseudotime branch. Aggregate equal-count bins on each partition using **no shared cells**. Fit two regularized target autoregressions on the *F_b* binned profiles:

- Baseline: predict target at bin `q+lag` using target at `q` and normalized bin position `q`.
- Augmented: use the baseline predictors **plus source s at bin q**.

Apply both models, without refitting, to the *V_b* binned profiles, using fitting-fold means and scaling. The directional held-out improvement is:

`g_{s→t} = clip((MSE_base(V_b)-MSE_aug(V_b)) / max(MSE_base(V_b), epsilon), -1, 1)`.

**Crucial:** the autoregressive target changes with edge direction. Thus `g_{s→t}` need not equal `g_{t→s}`. This target can be negative: adding an edge can worsen held-out predictive loss. Negative predictive usefulness is **not** a label asserting a biologically false interaction.

The comparative direction margin is `d_{s,t}=g_{s→t}-g_{t→s}`, hence `d_{t,s}=-d_{s,t}` by construction. This antisymmetry is a mathematical property, **not proof of causal direction**.

If multiple pseudotime branches are usable, average per-branch directional gains weighted by their *validation-cell counts*. A branch requires at least 8 time bins and 3 cells per bin **on both F_b and V_b**. Unusable branches are excluded; if none remains, evidence availability is `0` and gains are stored as **NaN, never zero**.

## Dataset boundaries and outputs

- Outer split (as M1.1): 60% expression fitting (300 of 500 cells), 20% iterative graph scoring (100), 20% untouched expression test (100).
- Inner 3-fold split of the 300: 200 model-fitting and 100 independent validation cells per fold.
- Original 12 **F_b-only features** are retained exactly from M2A by re-running the original fitting-only feature construction; no validation-derived target appears in that file.
- Training labels are computed from V_b only, using regression coefficients trained on F_b; no `GroundTruthNetwork.csv` access whatsoever.

Artifacts per benchmark network:

| File | Content |
|---|---|
| `m2c_features.csv` | Key columns plus 12 fit-fold-only M2A predictors |
| `m2c_targets.csv` | `y_directional_lag_gain_TARGET`, `y_pairwise_direction_margin_TARGET`, validation baseline/full MSE, availability mask |
| `m2c_manifest.json` | Dataset/fold hashes, source SHA256, branch availability and leakage boundaries |

Across networks, `m2c_quality_summary.csv` and `m2c_quality_manifest.json` verify row counts, exact feature-target key alignment, bounded signed gains, missingness behavior and antisymmetric directional margins.

Expected example row counts for the initial 3-fold exploratory networks:

| Network | Genes | Admissible directed edges | Edge × fold rows |
|---|---:|---:|---:|
| dyn-LL-500-1 | 18 | 306 | 918 |
| dyn-LI-500-1 | 7 | 42 | 126 |
| dyn-BF-500-1 | 7 | 42 | 126 |
| **Total** | | 390 | **1,170** |

The number of **usable directional labels** must be measured, not assumed.

## Commands

```bash
python -m pip install -e '.[dev]'
python -m pytest -q

python scripts/extract_beeline_v4.py \
  --archive BEELINE-data.zip --run dyn-LL-500-1 \
  --output data/raw/beeline/dyn-LL-500-1
python scripts/build_m2_direction.py \
  --raw-dir data/raw/beeline/dyn-LL-500-1 \
  --output-dir results/m2c/dyn-LL-500-1 \
  --dataset-id dyn-LL-500-1 --folds 3 --bins 16 --seed 42
python scripts/audit_m2_direction.py --root results/m2c
```

After merge, [M2C-A GitHub Actions](../.github/workflows/m2c-direction.yml) executes the three official benchmark structures once and uploads only derived features, surrogate targets and provenance; **never raw expression matrices or gold reference edges**.

## The necessary scientific limitations

1. Pseudotime-binned *different cells* do **not** constitute a true longitudinal gene-expression time series. Even when the simulator supplies accurate time ordering, this is a directed **predictive association** proxy, not causal regulation.
2. The target measures **out-of-fold marginal temporal predictive gain** conditional only on the target's own history and a time trend. It does not yet account for other candidate regulators or existing graph parent sets. Graph-conditioned *source-blind* scaffolds and larger evidence modalities remain M2C-B work.
3. Sample/edge-fold rows are dependent, and reversed edges share observations. Future learner fitting must use gene-pair-aware and **network-context held-out evaluation**. Do not randomly split individual edge-fold rows.
4. These same three topologies were already examined during M1/M2. They are **development datasets**, not prospective confirmatory results.
5. No gold-standard edge labels, independent motifs, ATAC-seq or knockouts are introduced here. Those are subsequent experimental stages, with independent provenance audits.
6. Multiple branches can have different temporal axes and missingness. Outputs must reflect available branches explicitly.
7. This is a necessary *target construction and integrity milestone*, not evidence that a future adaptive model outperforms a static prior.

## Acceptance gate

Before training a new reliability gate, inspect: signed gain distribution, non-zero difference between reverse edges, proportion of negative gains, saturation/clipping, branch coverage, target reproducibility under seeds, and agreement/conflict with expression-only conditional gain. Only then do training and held-out GRN comparisons follow.
