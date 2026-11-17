# M2C-D8 — actual external-source observation

This milestone replaces promises of real-data execution with a runnable GitHub Actions source download and metadata audit. The workflow downloads the public GSE216595 TF-to-cell mapping supplement from GSM6681047 to ephemeral runner storage, verifies gzip integrity, records SHA256, and saves actual observed field names, non-missing counts and unique-value counts in `actual_tfmap_observation.json`.

## Validation boundaries

- No parsed source artifact or run is claimed successful until the GitHub Actions job finishes and its artifact is examined.
- The source inspection **does not** infer TF column semantics, map barcodes to H5AD cells, infer experimental controls, or count biological replicates.
- The mapping table's rows or unique barcodes are **not** equivalent to independent experiments or proof of perturbation efficacy.
- The expression H5AD file is not yet downloaded in this workflow. The D6B barcode match and D7B experimental-design audit remain unfulfilled.
- No TF-response effect estimate or binary truth labels are generated.

## Next steps

Inspect the GitHub workflow run `m2d8-observe-geo.yml` and download `m2d8-actual-source-observation`. If the job fails, use its log to correct the real source URL/format rather than replacing results with synthetic data. If successful, pin the observed SHA256 and mapping columns, then acquire the H5AD metadata and run sample-aware barcode matching. Only after control groups and genuine biological replicates are independently verified can D5 effect estimation be applied to the actual assay.

The code and a synthetic test have been committed; observed biological results and external validation remain pending.
