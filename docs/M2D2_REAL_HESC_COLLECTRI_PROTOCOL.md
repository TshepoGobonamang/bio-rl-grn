# M2C-D2 — Real human hESC expression and CollecTRI TF-target prior integration

**Experimental scope:** independently sourced human TF–target evidence **relative to hESC expression measurements**; source overlap with any ChIP-seq benchmark evaluation reference is **NOT yet certified**. The intentionally prior-enriched hESC cohort is a schema/data-integration smoke test, **not** an unbiased GRN benchmark, regulatory accuracy claim, or RL experiment.

## Sources, versions, public verification

1. Human hESC single-cell RNA-seq: BEELINE v4 `BEELINE-data.zip` (Zenodo record **19009603**), publisher MD5 `9ea93065004f2f12022adeeaf3728edd`. Exact member paths in the release: `BEELINE-data/inputs/scRNA-Seq/hESC/ExpressionData.csv`, `PseudoTime.csv`, `GeneOrdering.csv`. The source inventory confirms 758 expression cells and gene-symbol labels like A1BG. All cell ordering is validated against hESC pseudotime.
2. CollecTRI **version 2.0** curated TF–target evidence, Zenodo record **8192729**:
   - `CollecTRI_regulons.csv`, publisher MD5 `cee4a3943c059e6dd8796ced6dec44f6`;
   - `CollecTRI_source.tsv`, publisher MD5 `c8622a46a49e21e955182f000cf1c067`.
   - Regulons columns are **source, target, weight, resources, references, sign_decision**. The `weight` sign describes CollecTRI's **regulation sign** (activation/repression), **not** TF→target arrow direction or calibrated probability. Supporting source names and citation counts are retained for matched interactions.
   - CollecTRI's curated and text-mined interactions are assembled from external literature/data. Their provenance is distinct from the hESC single-cell expression measurement, but the literature may overlap BEELINE ChIP-seq references, so **no independence from evaluation labels is asserted**.

## Observed real-source matching results (independently verified)

Pinned [schema run #44112971882](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/44112971882) and [actual gene-match run #44113607598](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/44113607598) succeeded.

| Grounded source quantity | Actual count |
|---|---:|
| Human hESC expression genes | **17,735** |
| Human hESC single cells | **758** |
| CollecTRI original regulon source rows | **43,175** |
| CollecTRI unique directed gene pairs | **42,713** |
| Exact directed CollecTRI pairs with both genes measured in hESC | **36,064** |
| Distinct matched curated regulators | **1,079** |
| Distinct matched curated targets | **5,868** |
| Matched pairs in deliberately enriched 8-TF + 24-target demonstration cohort | **233 of 992 candidates** |

The **36,064** represent *distinct curated associations with both symbols present in the expression matrix*; they are not inferred GRN true positives, and source overlap with benchmark ChIP labels is not certified. The 233/992 cohort ratio is achieved through *prior-based selection* and cannot be used as an unbiased biological accuracy measure.

## Two auditable layers

### D2-A — Gene-matched source coverage (no gold)

Parse all hESC gene identifiers from the original expression matrix and join curated CollecTRI source/target symbols using **exact, case-sensitive equality**. Never silently alias or orthology-convert symbols. Output:

- `m2d2_hesc_prior_match_audit.json`: counts of measured genes, matched TFs, targets and directed interactions; SHA256 of all public inputs; provenance and overlap flags.
- `hesc_collecTRI_matched_prior_edges.csv`: each mapped TF-target candidate; binary `observed_curated_association=1`, distinct `signed_regulatory_effect`, number of source rows, supporting resource names and cited-reference counts.

**An absent CollecTRI row is missing knowledge, not a known true negative or experimentally demonstrated absent regulation.** We do not convert an interaction sign of -1 into a reversed arrow.

### D2-B — Explicitly **prior-enriched** human expression integration smoke test

To demonstrate that the parent-conditional feature interface accepts *real measured human expression with a real curated prior*, select **eight high-coverage CollecTRI TF symbols and 24 measured targets** by a fixed deterministic rule based on matching edges. Extract exactly those 32 gene-expression profiles across available hESC cells, align pseudotime IDs, and compute F-only features and V-fold parent-conditional targets.

- This cohort selection **USES prior edges** to maximize integration coverage. It is **biased** and must NEVER be passed off as neutral target-gene selection, regulatory accuracy validation, general-purpose hESC biology or independent held-out confirmation.
- All candidate source/target pairs among selected genes remain eligible; only observed curated interactions have `m_independent_prior=1`; all other pairs have `m_independent_prior=0` and `independent_prior_score=NaN`.
- `independent_prior_score=1` means **observed curated association** for this smoke test, **not a confidence probability or verified binding strength**.
- Retain separate `hesc_real_features.csv`, `hesc_real_targets.csv`, `hesc_exploratory_manifest.json` and gene-list file. No evaluation gold/reference is opened. Parent loss predicts expression, not causal TF regulation.
- The original 60/20/20 split is retained with inner 3-fold cross-fitting; outer test and graph-scoring cells remain untouched during target construction.

## Distinguish two independence claims

- **Supported**: CollecTRI interactions were curated independently of the BEELINE hESC *expression measurements*, and matching human symbols establish source coverage.
- **Not supported yet**: those interactions are statistically or provenance-independent of BEELINE hESC ChIP-seq reference labels. In fact, the curated source may include overlapping ChIP-based and literature evidence. Source-overlap auditing requires the full underlying BEELINE reference network and independent experiment/source identifiers.

**Do not evaluate CollecTRI-matched edges using BEELINE ChIP as an independent holdout until resource/source overlap is reviewed and data splits established.**

## Reproducibility

New code:
- `scripts/inspect_m2d2_public_sources.py`: schema and pinned archive/record MD5 audits.
- `scripts/match_hesc_collecTRI.py`: all-gene exact source join and provenance rows.
- `scripts/run_m2d2_hesc_pilot.py`: 32-gene prior-enriched real-expression integration smoke.
- `tests/test_m2d2_real_source.py`, `tests/test_m2d2_hesc_pilot.py` and the [M2D2 GitHub workflow](../.github/workflows/m2d2-public-sources.yml).

Caution: reference/source overlapping families include ChIP/ENCODE/DoRothEA components; the regulon `resources` and `references` strings enable *screening* but source-level disjointness must be audited against actual held-out reference evidence before claiming external validity.

## Next independent biology **validation** experiment (not this smoke)

1. Obtain and hash a **non-overlapping** hESC-relevant TF perturbation or TF-binding evidence source, including assay, organism, genomic assembly and experiment accessions.
2. Join using stable gene IDs and document synonym ambiguities; lock the measured gene universe **without** using candidate prior edges or reference labels.
3. Predefine no-prior/fixed prior/corruption/missingness/adaptive conditions with identical graph-search budget, cell partitioning and a new independent evaluation context.
4. Quantify 0/25/50/75% source dropout, 0/10/25/50% sign corruption and positive/negative/missing evidence cases using source provenance; do not treat absent curated edges as negative labels.
5. Report independent-test expression generalization, regulatory AP/AUROC, direction accuracy and edge-wise source overlap, with uncertainty over independent datasets/networks. **No RL escalation solely from source coverage.**

Sources: [BEELINE experimental inputs](https://murali-group.github.io/Beeline/reproducing-results.html), [CollecTRI v2.0 pinned release](https://zenodo.org/records/8192729), [CollecTRI primary study](https://doi.org/10.1093/nar/gkad841).
