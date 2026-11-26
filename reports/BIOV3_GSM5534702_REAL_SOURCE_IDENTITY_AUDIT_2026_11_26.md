# BIO-V3 — Authentic K562/BT142 independent source inspection (26 Nov 2026)

## Reproducibility record

- Successful original GEO run: [GitHub Actions #45817674955](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45817674955), artifact `gcev-biov3-original-k562-training-source-schema`, ID **15116282713**.
- GEO series GSE182685, sample [GSM5534702](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5534702), a **mixture** of K562 and BT142 human cell lines, independent from our six ENCODE K562 CRISPRi *response* files (no overlap yet found by direct metadata accession; source-level study overlap still requires review).
- Expression `GSM5534702_K562-BT142_10x_scRNAseq_Expression.txt.gz`: source SHA256 `8bf3f175fe8635fec316263901bab18a4201d0199783482122142dd4387ffa75`.
- Cell metadata `GSM5534702_K562-BT142_10x_scRNAseq_MetaData.txt.gz`: source SHA256 `ef87b7605684e66f52a9f53c8feb0863b9bf20a3aed37aa1081f6b18489558d3`.
- Original authors' publicly available source: [MAESTER-2021, 10x cell-line processing script](https://github.com/petervangalen/MAESTER-2021/blob/fe8001bfd242d4440dcb72be42a28f346adbc047/2_Cell_line_mixes_decontX/2.2_TenX_CellLineMix_decontX.R) (source commit `fe8001bfd242d4440dcb72be42a28f346adbc047`).

## Observed, checksum-backed data

| Original-source quantity | Observed |
|:--|--:|
| Gene rows | **19,542** |
| Cell columns | **2,666** |
| Cell metadata records | **2,666** |
| Exact expression-cell to `orig.ident` joins | **2,666 / 2,666** |
| `CellType` original code `1` | 1,512 cells |
| `CellType` original code `2` | 1,154 cells |

The expression matrix is **gene × cell** with the first data field being an *unlabelled gene-row index*. The first line consists of **2,666 actual cell IDs**, not an `id` header; data rows have one additional first field holding gene identity. The original shallow scanner erroneously reported 2,665 joins because it skipped the first header cell ID. The full source parsing and exact join audit verified **all 2,666 cell identities** with no barcode transformations.

## Critical biological identity ambiguity

GEO exported `CellType` as numeric strings `1` and `2`; their direct mapping to K562/BT142 is **not documented in these two TSV files**.

Authors' script renames initial Seurat clusters `0 -> BT142`, `1 -> K562` and stores a factor with levels `c('BT142','K562')` (lines 85–87). A mapping from these initial clusters/factor levels to the exported numeric codes would be **plausible but is not demonstrated** merely by the numeric values. We must not assign `CellType=1` or `CellType=2` to K562 without verifying the exact export pipeline or source cell IDs against the authors' processed Seurat/RDS object.

Exploratory original expression marker values by numeric code:

| Marker | Code 1 (1,512 cells): mean / fraction positive | Code 2 (1,154 cells): mean / fraction positive |
|:--|:--|:--|
| HBG2 | **98.44 / 0.8155** | 86.32 / 0.9064 |
| PTPRZ1 | 4.27 / 0.5331 | **5.51 / 0.5719** |

These marker differences cannot confidently identify code-to-cell-line mapping. HBG2 and PTPRZ1 were referenced as source cell-type markers in the authors' code; the expression summaries here were measured directly from the source matrix. Neither coding direction is promoted to established biological identity.

## BIO-V3 eligibility decision

**PASS:** Original source download/integrity, full cell-barcode and gene-by-cell orientation reconciliation, availability of real single-cell expression without use of held-out ENCODE response labels.

**BLOCKED:** K562-specific cell subset identity. Gene identity namespace and version matching to GRCh38/GENCODE V29 also require inspection; data provenance against future regulatory prior and held-out assay sources is not verified.

**NO-GO:** Declaring an exact set of 1,512 or 1,154 cells as unperturbed K562 training data; fitting GCEV on the mixed expression without source qualification; claiming independent biological validation.

## Next scientific action

Inspect the authors' GEO-exporting code or source Seurat object `TenX_CellLineMix_Seurat_Keep.rds` to verify the numeric-code-to-K562 mapping *using exact barcodes*. If not recoverable, choose a genuinely single-cell-line untreated K562 source. Only then select candidate TFs/targets using expression-only data, source a motif-to-promoter prior with verified sequence/assembly and component provenance, and freeze all model rankings before evaluating an **uninspected** ENCODE perturbation panel.

The inspected six BIO-V2 treatment outcomes cannot count as independent final reference for these new models because they were viewed before prospective model freezing.
