# M2C-D6B — Actual GEO TF-barcode matching readiness

## Status
D6 created a public-data TF mapping schema workflow but a completed authentic mapping artifact has **not been verified**. D6B implements a strict **sample-aware exact join** with explicit operator-selected field names. Tests use synthetic fixtures; real match rates remain **unknown**.

## Procedure
1. Run [D6 real source workflow](../.github/workflows/m2d6-geo-schema.yml) and inspect the artifact column names, sample identifiers and hashes.
2. Verify the mapping supplement's true sample/barcode identifiers; obtain an H5AD `obs`-index export with explicit matching sample IDs and no counts. The mapping and H5AD schema may use different barcode naming conventions; **do not silently normalize suffixes**. If they disagree, independently document an accession-specific concordance rule and test it.
3. Run D6B audit with names as observed, rather than assuming `barcode`, `TF`, `sample` are present:

```bash
python scripts/audit_m2d6b_barcode_join.py \
  --tfmap local/real_tfmap.csv.gz --cells local/h5ad_obs_ids.csv \
  --map-barcode VERIFIED_MAPPING_BARCODE_COLUMN \
  --cell-barcode VERIFIED_OBS_BARCODE_COLUMN \
  --tf-column VERIFIED_TF_COLUMN \
  --map-sample VERIFIED_MAPPING_SAMPLE_COLUMN \
  --cell-sample VERIFIED_OBS_SAMPLE_COLUMN \
  --output results/m2d6b
pytest -q tests/test_m2d6b_barcode_join.py
```

4. Examine unmatched cells, unmatched map records and conflicts before attempting expression-effect estimation. Exact duplicate mapping records are deduplicated; distinct TF assignments for one (sample,barcode) raise an error. Sample/library labels are **not biological replicates**.
5. Independently identify negative controls, independent cultures, biological replicate hierarchy, batch structure, and perturbation confidence scores. The audit unconditionally refuses to certify these.

## Gate
This audit yields `barcode_join_audit.csv` and `barcode_join_summary.json` with match counts and file hashes, but never TF response labels, perturbation significance tests, independent regulatory truth or model accuracy. Proceed to replication-aware inference only after authentic metadata and protocol support controls and biological replicate structure.
