# M2C-D — Parent-conditional target-quality audit (17 November 2026)

**Official source:** [Successful GitHub Actions run #44111372508](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/44111372508), artifact **m2cd-parent-conditional-gold-blind**, artifact ID **14340607594**. Code [PR #18](https://github.com/TshepoGobonamang/bio-rl-grn/pull/18), merged commit `f7d228dd51f5c22c34c154b900fa0084f1f0ca85`, tests green Python 3.10–3.12.

## Scientific verdict

**GO for parent-conditional target engineering and prospective independent biological evidence collection. NO-GO to claiming this target is a superior regulatory-direction signal or that independent biological priors have already been evaluated.** Parent-conditioned held-out expression gain is cross-fold reproducible and distinguishes conditional from marginal association; evaluation-only reference recovery on existing BEELINE contexts remains weaker than signed pseudotime for AP and direction.

### 1. Source and leakage checks

- **1,698** edge×fold rows from six networks; **566** unique directed candidates when counted separately by network.
- All training predictors and validation target tables are separate; target value can be reconstructed numerically from the recorded held-out before/after losses.
- Per network: 300 outer fit cells, 100 graph-internal scoring cells **excluded**, 100 untouched expression-test cells **excluded**; 3 inner fitting/validation folds within outer fitting only.
- Source-blind scaffold: maximum 2 potential alternative parents chosen from **inner training expression only**, excluding candidate source and target. No gold edges in parent construction.
- **0 independently sourced biological edges** in all six synthetic networks. Every independent-prior mask is zero and corresponding scores are missing (NaN), **not zero as negative evidence**.
- Any JASPAR/CollecTRI inputs for these abstract synthetic gene IDs would be invalid; this experiment does not include them.

### 2. Real measured cross-fold quality

| BEELINE family | Candidate directed edges | Edge×fold rows | Conditional held-out gain mean | Marginal gain mean | Negative conditional gain | Mean abs conditional-marginal difference | Mean 3-fold Spearman |
|---|---:|---:|---:|---:|---:|---:|---:|
| dyn-LL | 306 | 918 | 0.020739 | 0.049241 | 38.56% | 0.028788 | **0.4379** |
| dyn-LI | 42 | 126 | 0.034489 | 0.055415 | 14.29% | 0.020994 | **0.8222** |
| dyn-BF | 42 | 126 | 0.096345 | 0.123522 | 9.52% | 0.052883 | **0.9150** |
| dyn-CY | 30 | 90 | 0.132839 | 0.196106 | 14.44% | 0.079030 | **0.8229** |
| dyn-TF | 56 | 168 | 0.083417 | 0.087879 | 17.86% | 0.026790 | **0.8409** |
| dyn-BFC | 90 | 270 | 0.056335 | 0.109880 | 11.85% | 0.061230 | **0.7279** |
| **Total** | **566** | **1698** | | | | | |

No target clipped at ±1 in these six datasets. Conditional gain is **not identical** to marginal gain. This establishes added predictor information conditional on other expression predictors, not directional causation.

As a descriptive comparison, previous pairwise lagged-gain scores had low rank stability across folds (mean ~0.171/0.243/0.025 for LL/LI/BF in M2C-A), but the estimands are different: the new target predicts *same-cell expression*, not lagged pseudotime bins. Improved numerical stability does not alone imply better regulatory inference.

### 3. Evaluation-only reference diagnostic (not used in target construction)

**CY, TF and BFC were already inspected in M2C-A2**; their reference networks are reused here *only* to diagnose the new, already frozen per-edge scores. The reference edges were read from the M2C-A2 evaluation-only audit artifact, **never in the M2C-D target generator or GitHub workflow**. Average per candidate score across 3 inner folds. Noninterpolated average precision (AP) and AUROC computed via scikit-learn. The following are **not** prospective independent benchmarks.

| Network | Method | AP | AUROC | Correct unique reference direction |
|---|---|---:|---:|---:|
| CY (6 gold/30 candidates) | Parent-conditional | 0.3810 | 0.7847 | 3/6 |
| | Marginal gain | 0.4089 | 0.7847 | 3/6 |
| | Signed pseudotime `z_signed` | **0.7236** | **0.9097** | **6/6** |
| TF (17/56) | Parent-conditional | **0.4725** | 0.5596 | 5/9 |
| | Marginal gain | 0.4665 | 0.5747 | 4/9 |
| | Signed pseudotime | **0.5700** | **0.6440** | **6/9** |
| BFC (15/90) | Parent-conditional | **0.3671** | 0.7964 | 8/13 |
| | Marginal gain | 0.2902 | 0.6720 | **11/13** |
| | Signed pseudotime | **0.5513** | **0.7982** | **12/13** |

Note the boldface only highlights within-method conditional-vs-marginal gains and the overall stronger `z_signed`. *Different source targets and scoring scales should not be conflated with causal evidence.* Reference prevalence: CY 6/30, TF 17/56, BFC 15/90.

### 4. Next experimental gate

**Stage D2 must use biologically meaningful genes and at least one genuinely independent prior.** A realistic next pair is BEELINE real hESC single-cell expression with human gene IDs and independently derived JASPAR 2026 motif-to-promoter predictions, or carefully provenance-filtered CollecTRI/ChIP/perturbation evidence. JASPAR CORE PFMs are binding models, **not** TF-target edges on their own; they must be scanned against version-matched promoter sequence and linked to genes. Curated CollecTRI edges might overlap evaluation references; audit underlying source-level provenance before calling them independent.

Next code must (a) acquire and hash pinned compatible expression and biological evidence, (b) verify gene ID/organism/assembly matches and source overlap, (c) compare no-prior and true-prior/masked/corrupted controls and (d) use new uninspected validation contexts. A soft parent-conditional gain label alone does **not** train a directional biological reliability model.

### Sources
- [BEELINE original benchmark](https://www.nature.com/articles/s41592-019-0690-6).
- [JASPAR 2026 data](https://jaspar.elixir.no/downloads/).
- [CollecTRI source resource](https://github.com/saezlab/CollecTRI) and [independent underlying source datasets](https://zenodo.org/record/8192729).
