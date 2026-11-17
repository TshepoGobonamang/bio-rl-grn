# M2C-D — Parent-conditional predictive targets and independent biological evidence

## Research premise and boundaries

The M2C-A/A2 *pairwise pseudotime lagged-gain* targets proved unstable and failed to provide consistent incremental signal over `z_signed`. This **new** experiment tests a different, gold-blind target: *conditional expression gain from adding a proposed regulator after a source-blind parent scaffold*. It does not claim causality or fit an adaptive reliability learner.

The BEELINE synthetic `g1...` networks represent abstract simulation genes and are not real human TF targets. **Do not map JASPAR/ENCODE/TRRUST records to `g1...`**. The six-network synthetic run intentionally writes `m_independent_prior=0` and `independent_prior_score=NaN` for every candidate. The absence of external priors is scientific integrity, not a failed import.

## 1. Formal gold-blind conditional target

Let `F_b` be inner fitting cells and `V_b` be disjoint inner validation cells, contained entirely within the original 60%-outer-fitting subset.

For each candidate edge `s→t`:

1. On `F_b` **only**, rank potential other candidate parent genes `p != s,t` by absolute Pearson correlation with the *target t*. Select up to `k=2` top parents, a **source-blind** scaffold `S_t^(-s)`.
2. Fit baseline ridge predicting target `t` from `S_t^(-s)` on `F_b`; compute normalized mean squared error on `V_b`.
3. Fit an augmented ridge including candidate `s`, also on `F_b`; evaluate on the *same* `V_b`.
4. Define signed `y_parent_conditional_gain = clip((L_before-L_after)/max(L_before,eps), -1,1)`. A **negative** gain means predictive deterioration after adding the edge, NOT biologically proven non-regulation.
5. Compute expression-only **marginal gain** on the same partition as a comparator; report their disagreement and crossfold rank Spearman.

These are independently held-out **same-sample expression reconstruction** targets, *not* temporal lags or perturbation effects. The parent scaffold is selected using F_b only, with no evaluation gold or external source labels. This is an intentionally simple source-blind two-parent baseline; broader parent selection and regularization sensitivity remain development tasks.

## 2. Separate evidence modalities and provenance contract

`src/biorlgrn/m2_external.py` accepts **only explicit, gene-mapped independent TF-target evidence files** for `--dataset-kind biological`. CSV schema:

```csv
source,target,score,evidence_id
<measured-TF-symbol>,<measured-target-symbol>,0.85,<unique-source-record>
```

`score` must lie in [0,1]. Every file requires a separate JSON provenance manifest with:

- `database` — JASPAR/ENCODE/other actual provider;
- `release` and `source_url`;
- `organism` and `genome_assembly`;
- `assay` — `motif_promoter_link`, `chip_binding`, `atac_motif_link`, `perturbation_effect`, or `curated_binding`;
- `derivation` — **exact** method linking a TF observation to the target gene, for example motif scanning in a defined promoter region using gene annotations of matching assembly;
- `role: predictor_prior_only`;
- `independent_of_evaluation_reference: true`, subject to posthoc overlap audit.

The loader rejects synthetic dataset kinds and synthetic `gN` identifiers, mismatched measured gene names, missing provenance, unbounded confidences, duplicate records and evaluation-file names. It records SHA256 and excludes unmatched gene records without converting absent source evidence to negative labels. The independent source's assertion is **not proof** that its edges are absent from the benchmark reference. A posthoc provenance and directed-edge overlap audit is required before claims of fair biological validation.

A JASPAR CORE **position-frequency matrix is not itself a TF→target interaction**. It must be combined with promoter regions, gene annotations, an explicitly defined binding-site scan, and preferably cell-context chromatin accessibility; binding-site presence alone does not prove in-vivo transcriptional regulation. JASPAR 2026: https://jaspar.elixir.no/downloads/ ; original database references: https://academic.oup.com/nar/article/54/D1/D184/8343514 .

**Current milestone does not import real gene-matched prior edges or run on real biological single-cell expression.** It provides the provenance-enforced interface to do so in a compatible next experiment. The six-network BEELINE target run will explicitly show zero external biological evidence. Test fixtures with invented scores are marked TEST ONLY, never included as biological results.

## 3. Expected BEELINE output and quality criteria

Each dataset exports:
- `m2cd_features.csv`: fitting-only columns (`z_signed`, `m_temporal`, `m_independent_prior`, `independent_prior_score`, expression correlation, variance, scaffold size/identities). No target labels.
- `m2cd_targets.csv`: source-blind held-out conditional target, marginal target, before/after held-out loss; never mixed into predictor inputs.
- `m2cd_manifest.json`: exact outer split, inner fold hashes, candidate universe and available-modality flags, SHA256 and leakage restrictions.

An aggregate `m2cd_quality_summary.csv` measures available biological sources (must be 0 for synthetic), mean signed gains, fraction negative, target clipping, mean absolute conditional-minus-marginal gain and crossfold Spearman stability.

Expected rows for 3 inner folds of 500-cell BEELINE datasets using every measured gene as regulator:

| Family | Genes | Candidates | Rows |
|---|---:|---:|---:|
| dyn-LL | 18 | 306 | 918 |
| dyn-LI | 7 | 42 | 126 |
| dyn-BF | 7 | 42 | 126 |
| dyn-CY | 6 | 30 | 90 |
| dyn-TF | 8 | 56 | 168 |
| dyn-BFC | 10 | 90 | 270 |
| **Total** | | 566 | **1698** |

These are **expected**, not evidence of a completed workflow until a successful official artifact is inspected.

### How to run

```bash
python -m pip install -e '.[dev]'
python -m pytest -q
python scripts/build_m2cd.py --raw-dir data/raw/beeline/dyn-LL-500-1 --output-dir results/m2cd/dyn-LL-500-1 --dataset-kind synthetic --max-parents 2
python scripts/audit_m2cd.py --root results/m2cd
```

On a properly gene-mapped independent biological expression dataset:

```bash
python scripts/build_m2cd.py \
  --raw-dir data/raw/biological/your-dataset \
  --output-dir results/m2cd/your-dataset \
  --dataset-kind biological \
  --evidence data/priors/verified_binding_links.csv \
  --evidence-manifest data/priors/verified_binding_links.provenance.json
```

The independent evidence must be obtained, mapped, and validated before this command; never invent links.

## 4. Decision gates

**Stage D1 (this implementation):** Confirm gold-blind cross-fitted parent-conditional targets, recompute them from held-out losses, check fold variability, and quantify how much conditioning changes marginal gains. This is **target engineering**, not proven regulatory direction.

**Stage D2 (requires independent data):** Obtain *real measured gene IDs*, independently mapped binding/perturbation evidence from a different source than evaluation labels, demonstrate source-availability and overlap audit, compare prior-masked/fixed/nonlearned/learned reliability under controlled dropout and corruption.

**Stage D3:** Compare true unseen context/topology generalization with no priors, fixed/global and adaptive weights plus independent expression-test loss. Only then evaluate whether graph-conditioning or eventually RL adds incremental value.

Do **not** claim independent biological evidence was used in the synthetic BEELINE run. Do **not** interpret conditional ridge gain as causation.

[BEELINE methodology](https://www.nature.com/articles/s41592-019-0690-6); [JASPAR 2026 downloads](https://jaspar.elixir.no/downloads/).
