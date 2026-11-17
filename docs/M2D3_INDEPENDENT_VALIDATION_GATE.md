# M2C-D3 — Independent validation provenance eligibility gate

## Scientific status
M2D2 confirmed 17,735 hESC measured genes, 36,064 matched CollecTRI TF→target pairs, and a prior-enriched 32-gene software pilot. It also identified 2,848 edges common to the curated prior and the BEELINE hESC ChIP reference. This does **not** establish independent biological recovery, and the prior-enriched pilot cannot be used for unbiased comparison.

M2D3 introduces a **fail-closed accession-level audit** before any out-of-sample biological claims. It does **not** certify independence. Even disjoint assay accession strings may describe related or reused experiments; expert source-level review remains mandatory.

## Inputs
The CLI `scripts/audit_m2d3_validation_provenance.py` takes:

* `--prior`: CSV with `source,target,evidence_accessions`; semicolon-delimited, documented experiment IDs for **each** interaction. CollecTRI publication references are **not automatically** equivalent to experimental accession IDs. Do not fill missing IDs with invented values.
* `--assay`: external experiment CSV with `source,target,label,assay_accessions`. Labels must be explicitly assayed binary 0/1, with negatives justified by assay design. Missing edges and unrecorded interactions must never be labeled zero.
* `--manifest`: JSON including `candidate_universe_expression_only=true`, `predictions_frozen_before_assay_access=true`, `distinct_assay_origin_documented=true`, `assay_species` and `assay_type`. Boolean flags are asserted by the user/pipeline; audit their supporting artifacts independently.
* `--output`: local output directory; no reference edge files are permitted in fitting/target generation.

Example:

```bash
python scripts/audit_m2d3_validation_provenance.py \
  --prior local/prior_accessions.csv \
  --assay local/heldout_assay_edges.csv \
  --manifest local/frozen_assay_manifest.json \
  --output results/m2d3-provenance

pytest -q tests/test_m2d3_validation_provenance.py
```

## Outputs and interpretation
`m2d3_provenance_edges.csv` has an explicit per-edge status:

- `direct_accession_overlap`: same accession appears in both evidence records; excluded from independent claims.
- `unknown_prior_accession`: prior-record experimental provenance incomplete; excluded.
- `unknown_assay_accession`: validation assay provenance missing; excluded.
- `pending_source_level_review`: recorded accessions differ, but **not** certified independent.
- `pending_global_prior_provenance_review`: no curated edge is recorded; that is not proof of source independence or a true negative.

`m2d3_provenance_summary.json` always sets `independent_validation_certified=false`. Certification requires a separately documented source-family, study/experiment accession and processing-overlap investigation plus a specifically predeclared holdout protocol. Do not flip the flag merely because this audit passes.

## Next gate after D3
1. Locate a hESC-relevant independently generated TF perturbation or assay dataset; preserve experiment IDs and data license, genomic assembly, TF/gene identifier mapping, assay control design, cell context and source provenance.
2. Construct an **expression-only-defined** candidate gene universe, freeze models and edge scores before inspecting assay labels. CollecTRI can then be evaluated as an optional model input but must not define the evaluation universe.
3. Independently document label design, including experimentally assayed non-responders as potential negatives, with statistical power and missing-assay masks.
4. Review overlap among publications, experiment IDs, source families and reused ChIP/ENCODE/DoRothEA datasets; exclude unresolved cases from independent scoring.
5. Only after clearance: compare expression-only, fixed prior, adaptive prior, dropout/sign-corruption controls at matched candidate/compute/search budgets; quantify held-out expression loss, positive-label retrieval and confidence intervals. Avoid causal claims from observational association alone.

This milestone is a **safety and reproducibility infrastructure deliverable**, not an AUROC/AUPRC result and not permission to advance into RL.
