# BIO-V3B — Source selection and K562-only expression audit

**Research decision (26 Nov 2026):** Do not infer K562 membership from numeric `CellType` in GSM5534702. The original publication reports 1,310 K562 and 1,356 BT142 cells after 10x filtering, while the observed GEO export has 1,512 cells under numeric code 1 and 1,154 under code 2. Source-backed raw matrix orientation and barcode joins passed, but the original label-export transformation has not been established. The original authors' Seurat code `RenameIdents("0"="BT142","1"="K562")` and factor levels alone do **not** reconcile this mismatch.

**Alternative frozen before data inspection:** GSM9224207 from GEO GSE306708, source filename `GSM9224207_Fresh_K562.h5ad`. NCBI GEO identifies the sample as K562 and provides a separate fresh K562 expression H5AD. We must not mix this file with the fixed-cell H5AD or assume they are technical/biological replicates.

## Execution
The GitHub workflow `.github/workflows/gcev-biov3b-pure-k562.yml` downloads original source, prints SHA256, tests unique exact identifiers and nonnegative finite values, reports original `obs`/`var` metadata and any raw/layer fields, and exports **unmodified** sparse expression values with ordered cell/gene IDs. The audit does not infer that X contains raw counts without evidence. It does not select genes/TFs using ENCODE responses and does not fit a GCEV model.

Expected artifact `gcev-biov3b-authentic-k562-only-expression` includes original-scale sparse matrix `k562_expression_original_scale.npz`, `k562_cell_ids.csv`, `k562_gene_ids.csv`, `k562_source_audit.json`. The original source H5AD itself is never committed.

## Scientific next gates
1. Verify actual source H5AD execution and inspect matrix dimensions, gene IDs, expression-scale evidence and sample metadata. If source is normalized or has uncertain transformations, explicitly document its scale and avoid treating values as raw counts.
2. Match gene identifiers to a **fixed** GRCh38/GENCODE annotation. Never guess aliases or mix transcript version numbering.
3. Define a reproducible TF/gene candidate selection using **outer fitting expression only**; freeze fit/internal/test partitions and the candidate universe before use of regulatory references.
4. Add an orthogonal, independently curated **K562-compatible prior** with component-level assay/accession provenance. No motif PFM can be treated as an edge without a sequence scan and promoter/TSS mapping.
5. The six BIO-V2 ENCODE TF responses have already been inspected. They are **exploratory development data**, not a pristine held-out confirmatory reference. Source a new untouched TF panel and lock model rankings before accessing outcomes.

**Current research claims:** original cell-line/source identity selection made; original H5AD inspection results depend on the GitHub execution artifact. No new biological GCEV ranking result is claimed.

Sources:
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9224207
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5534702
- https://pmc.ncbi.nlm.nih.gov/articles/PMC9288977/
- https://github.com/petervangalen/MAESTER-2021/blob/fe8001bfd242d4440dcb72be42a28f346adbc047/2_Cell_line_mixes_decontX/2.2_TenX_CellLineMix_decontX.R
