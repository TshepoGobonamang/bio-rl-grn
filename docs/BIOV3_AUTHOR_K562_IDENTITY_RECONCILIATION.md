# BIO-V3 — Independent K562 source identity decision (26 Nov 2026)

## Source evidence and why we must not relabel GEO numerics

The frozen original GEO source **GSE182685 / GSM5534702** consists of a K562/BT142 **mixed** single-cell experiment. The prior successful original-source audit recorded 2,666 exact cell identifier joins and original `CellType` numeric groups of 1,512 (`1`) and 1,154 (`2`). Source SHA256s are frozen as code defaults in `scripts/biov3_reconcile_k562_labels.py`.

Authors' [original Seurat processing script at commit fe8001bf](https://github.com/petervangalen/MAESTER-2021/blob/fe8001bfd242d4440dcb72be42a28f346adbc047/2_Cell_line_mixes_decontX/2.2_TenX_CellLineMix_decontX.R) assigns cluster 0→BT142 and cluster 1→K562, then creates a factor with named levels `c("BT142","K562")` and exports `TenX_CellLineMix_Seurat_Keep.rds`. The accompanying paper reports **1,310 K562 and 1,356 BT142** cells after 10x decontamination, summing to 2,666. These counts are incompatible with either 1,512/1,154 numeric category being the original cell-line count. Therefore it is **not** defensible to cast factor levels 1/2 into the GEO `CellType` field without showing the original per-barcode mapping.

Sources: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5534702 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC9288977/ ; https://github.com/petervangalen/MAESTER-2021/blob/fe8001bfd242d4440dcb72be42a28f346adbc047/2_Cell_line_mixes_decontX/2.2_TenX_CellLineMix_decontX.R . A subsequent 2024 independent analysis explicitly reused the authors' `TenX_CellLineMix_Seurat_Keep.rds`: https://pmc.ncbi.nlm.nih.gov/articles/PMC10980719/ . This is a source file identity, **not a verified public direct download URL**.

## Safe source identity protocol

1. Inspect the *public publisher resource index* `https://vangalenlab.bwh.harvard.edu/resources/maester-2021/` for an actual linked `TenX_CellLineMix_Seurat_Keep.rds`. A missing, timed-out or authorization-protected link is **not** permission to invent a source URL or substitute marker-based guesses.
2. If legitimately obtained, save the original RDS file SHA256 and run `Rscript scripts/biov3_extract_author_seurat_labels.R original.rds original_labels.csv`. The RDS must be the original author-labeled Seurat object, and its exact `meta.data` rownames and named `CellType` (`K562` / `BT142`) must be used unmodified.
3. Run `scripts/biov3_reconcile_k562_labels.py` with original GEO expression, metadata, original author-derived labels, and original RDS SHA. Require **all 2,666 barcode IDs match exactly**; no suffix trimming, inferential mapping, cluster relabeling or partial rescue. Cross-tabulate GEO numeric codes against original author identities; if mixed, report the disagreement, **but still permit barcode-derived K562 export** because the authoritative original labels do not depend on the GEO code.
4. Export `K562_gene_by_cell_original_counts.tsv.gz` plus SHA256 and a fail-closed manifest **only after a successful source cell identity match**. The exported gene symbols remain unmapped to GENCODE v29; no GCEV model fit is authorized until gene namespace, count transform, TF universe, external prior independence, and model-data partition are audited.

## Alternative if original RDS cannot be retrieved

Do **not** keep rewriting numeric-label guesses. Evaluate a *separately authenticated* single-K562-source dataset (for example, 10x Genomics `10k Human K562-r Cells, Singleplex Sample`, noting K562-r is a **distinct clone** and Flex fixed-RNA differs from the 2021 3′-v3 protocol) or another untreated K562 3′-RNA study with an independently checked unambiguous sample identity. Freeze that alternate cohort and its gene annotation *before* looking at perturbation-based evaluation outcomes. Never silently substitute K562-r expression as if it were the same clone as ENCODE K562 CRISPRi.

## Acceptance and publication limits

**PASS only when** an original-author barcode source has a recorded file hash, 100% exact cell-ID match, fully reconciled cell-type assignments and no source conflicts. A numeric `CellType` field, marker enrichment or the published number of cells alone is **NOT identity proof**.

No ENCODE BIO-V2 perturbation results enter these selection steps. Those six perturbations have already been inspected and can only support exploratory comparisons, not an uninspected biological confirmation. Regulatory priors remain unqualified, and this experiment does not establish direct TF-target causality or GCEV superiority.

## Executable single-cell-line fallback, explicitly not identical K562

10x Genomics published a one-sample **K562-r** (ATCC CRL-3344) Flex/fixed-RNA dataset [10k K562-r singleplex](https://www.10xgenomics.com/datasets/10k-human-k562-r-cells-singleplex-sample-1-standard). **K562-r is a distinct clone/protocol from the 2021 K562–BT142 mixture and ENCODE's K562 CRISPRi panel**. It can only serve as a separately labelled exploratory training candidate subject to context-comparability review.

An independent optional 10x HDF5 source workflow is implemented in `.github/workflows/gcev-biov3-k562r-source-candidate.yml`. It *tests a plausible 10x CDN asset path*, not a previously verified official file URL; failure to retrieve an authentic 10x HDF5 leaves this alternative blocked. After a successful download, the exporter checks gene/cell ID uniqueness and integer counts, records SHA256 and freezes 800-cell 60/20/20 train/internal-score/test splits. Gene selection uses **fit-only variance** and exports 256 Ensembl features without any ENCODE outcome or regulatory reference. It does not fit GCEV or assert equivalent cell-line biology.
