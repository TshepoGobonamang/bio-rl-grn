# M2C-D15 — Original 10x processing-code reconstruction and TFmap reconciliation

Date: 2026-11-19

## Sources and verified rules

Original laboratory source: https://github.com/fengzhanglab/Joung_TFAtlas_Manuscript/blob/main/map_cells_10X.py

The script reads `ORF_lib_barcode_map.csv` (TF barcode to ORF isoform) and `counts_cells.csv` (cell barcode with TF-barcode/read-count pairs), and emits `Cell_ORFs.csv`. In the original script:
- a cell with no TF-count entries is labeled `NA`;
- a single candidate with count > 1 is assigned its mapped ORF, otherwise `AMB`;
- for multiple candidates sorted by count, if top count < second count × 1.25 the cell is labeled `AMB`; otherwise the leading barcode is assigned to its ORF.
These are historical code rules, not proof that the released H5AD was processed with this exact file/version.

The repository README describes `count_barcodes_10X.py` as the prior fastq parsing stage. GEO GSM6681047 states that a custom Python script produced the cell/TF assignment CSV (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6681047).

**Unresolved:** There is no verified public script path yet showing exact H5AD construction of `obs['TF']`, `obs['batch']`, or anomalous `obs['temp']`. GSE216595 H5AD's 87 observed TF labels do not prove which raw barcode entries were filtered, which ambiguous cells were removed, or whether there are experimental controls.

## Implementation

`scripts/reconcile_m2d15_tfmap_h5ad.py` accepts original H5AD and GSM6681047 CSV.gz with **explicit actual column names**. Produces source SHA256s, exact matches, mismatches, missing labels, unjoined H5AD cells, unused TFmap records and exact join coverage.

- Prefer **sample + cell** exact join (specify both fields) when both files expose supported sample IDs.
- A cell-only join is allowed as a diagnostic but is **not source-aware certification** if cell IDs can repeat across samples. Never claim source-disjoint validation from cell-only matching.
- No automatic stripping `-1` barcode suffixes or inventing TF/replicate mappings.
- Conflicting TF values per same join key trigger a hard failure. Identical duplicate keys are reported and deduplicated.
- A successful join is **not** verification of independent cultures, matched controls or direct regulatory edges.

## How to execute after observing real TFmap headers

```bash
python -m pip install anndata pandas pytest
python -m pytest -q tests/test_m2d15_tfmap_reconciliation.py
python scripts/reconcile_m2d15_tfmap_h5ad.py \
  --h5ad local/GSE216595_180124_perturb.h5ad.gz \
  --tfmap local/GSM6681047_180124_TFmap.csv.gz \
  --h5ad-cell-column _h5ad_index_exact \
  --tfmap-cell-column '<OBSERVED_CELL_FIELD>' \
  --h5ad-tf-column TF \
  --tfmap-tf-column '<OBSERVED_TF_FIELD>' \
  --h5ad-sample-column batch \
  --tfmap-sample-column '<OBSERVED_SAMPLE_FIELD>' \
  --output results/m2d15/reconciliation.json
```

Do not substitute placeholder column names, or assert that `batch` values are identical to GEO sample identifiers until confirmed. H5AD `temp` must remain unmodified until reconstructed from the authors' full processing scripts.

## Current gate

**Methodology implemented, real source reconciliation not yet executed.** D15 requires the actual M2D8 source observation JSON or original TFmap CSV to identify real field names and a real paired H5AD download. Experimental control and culture-replicate provenance are still NO-GO.
