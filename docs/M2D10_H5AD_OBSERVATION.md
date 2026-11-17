# M2C-D10: H5AD metadata observation before perturbation inference

The original GSE216595 TF-to-cell mapping audit has not yet yielded a verified artifact in the connected run lookup. This milestone therefore implements the **other independently required input**: backed inspection of the published H5AD expression object without loading its dense expression matrix or generating biological labels.

## Run on actual downloaded source
Download the original GEO processed file outside Git, verify accession and size, then:

```bash
python -m pip install anndata pandas pytest
python scripts/inspect_m2d10_h5ad.py \
  --h5ad local/GSE216595_180124_perturb.h5ad.gz \
  --output results/m2d10/observed_h5ad_schema.json
pytest -q tests/test_m2d10_h5ad.py
```

Record the download URL, source release, file length and generated SHA256. Do not commit source expression or barcode assignments to the public repository.

The JSON reports source SHA256, cell and gene counts, unique index checks, exact observed metadata field names, missingness and distinct-value counts. It does **not** infer that an `obs` field means biological replication or control status.

## Next decision based on real results
Inspect both actual outputs, `m2d8/actual_tfmap_observation.json` and `m2d10/observed_h5ad_schema.json`. Identify sample and barcode fields from observed names, construct a documented map between any transformed cell IDs, then run D6B exact sample-aware barcode matching. If no compatible keys are present, report NO-GO and examine the authors' TF assignment code; do not silently modify barcode suffixes.

GSE216595 pooled TF overexpression after iNP differentiation does not by itself satisfy the D5 paired culture design. Controls and biologically independent cultures require direct methods and metadata verification. Independent regulatory accuracy, binary truth labels and RL remain prohibited until these prerequisites are satisfied.

**Status:** inspection code and synthetic tests committed; original GEO H5AD not executed and authentic matching/replication results unknown.
