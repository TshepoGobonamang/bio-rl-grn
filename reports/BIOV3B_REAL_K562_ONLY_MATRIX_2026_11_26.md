# BIO-V3B — Original K562-only expression, source integrity and split qualification

**Date:** 26 November 2026. **Scope:** expression-source qualification and source-blind dataset preparation; no GCEV training, no biological prior or held-out perturbation scoring.

## Source selection decision
The earlier mixed K562/BT142 GEO sample GSM5534702 has 1,512/1,154 numeric code groups, which do not equal the source publication's reported 1,310 K562/1,356 BT142 after filtering. The original Seurat script defines named labels, but no provenance-backed mapping from the exported numeric codes to these names was established. **Do not label either GEO numeric group K562.**

We therefore selected separately identified **GSM9224207 Fresh_K562.h5ad**, a pure-K562 source from GSE306708 (the GEO metadata describes K562; fresh vs fixed conditions are provided as distinct files). Original URL: https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9224nnn/GSM9224207/suppl/GSM9224207_Fresh_K562.h5ad

## Actual source observation

Original H5AD SHA256: `14eb3129e57769899e0d301f515531a7e410e52dc7488fb1f981c43b2ea99ca6`.

- **658 cells × 60,671 gene symbols**, original identifiers unique; 0 ENSG-style feature names, as opposed to the GENCODE V29 ENSEMBL IDs used by the earlier ENCODE response analysis.
- Original `X`: sparse, **1,878,250 nonzero entries**, all finite and nonnegative; all observed nonzero values integer-like (range 1 to 2458). This is **count-like**, not independently certified raw UMI counts.
- `layers`: `Q-Q_normalized`. Its intended transformation must be inspected in source processing methods before choosing the modelling scale.
- `obs`: `n_counts`, `barcode_seq`, `RTp`, `mito_frac`, `assigned_prc`.
- Original `var`: no metadata columns.
- No ENCODE perturbation response or biological source labels were opened during this inspection.

The initial [original-data GitHub source workflow](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45881393706) produced a source-audit JSON and identity-preserving sparse matrix artifact.

## Fit-only split from original verified artifact (first local reproducibility check)

The original-scale expression object was loaded from the GitHub artifact and the deterministic script `scripts/prepare_biov3b_k562_splits.py` used **only outer fitting cells** to choose 512 genes (minimum 5% fit-cell detection and descending log1p-variance, deterministic gene-name tie-break).

- RNG seed: 20261126.
- Fit: **394 cells**.
- Internal graph scoring: **131 cells**.
- Untouched test: **133 cells**.
- Full matrix: 658 × 60,671; selected source-original-scale matrix: 658 × 512.
- Selection **does not** access regulatory reference edges, ENCODE perturbation outcomes or prior databases. Same cell partition and selection code will be independently rerun on GitHub in follow-on CI.
- The 512-gene set is selected from original source symbols, **not yet** linked to versioned ENSG IDs, TF roles, or a verified regulatory prior.

## Decision and limitations

**PASS:** source cell line is explicitly K562, source file integrity, sparse matrix / gene/cell identity, descriptive original-scale matrix export, and fit-only gene-selection protocol.

**REVIEW:** independent source processing confirmation (whether `X` contains raw counts or normalized values), exact GRCh38/GENCODE feature-ID mapping, cell-level quality filtering and source donor/biological replicate interpretation. A single single-cell study does not supply independently cultured donors.

**NO-GO:** using the already inspected six BIO-V2 TF response outcomes as untouched independent confirmation, interpreting zeros as disproven TF regulation, or claiming GCEV advantage over fixed/graph-independent methods.

**Next:** source-specific normalization/QC choice; lock candidate TF universe with an independent TF annotation that has no held-out assay leakage; validate a genuine K562-compatible motif-to-promoter or other source-disjoint prior; freeze ranking predictions; acquire a **new and uninspected** perturbation reference before accuracy evaluation.
