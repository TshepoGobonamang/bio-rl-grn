# M2C-D7B — evidence-backed real assay-readiness gate

The real-source workflow `m2d7-real-source.yml` is committed. **A successful real TF-map run/artifact has not been verified in this milestone**. A general repository “Tests” workflow success does not substitute for the GEO source-inspection output.

## Deliverable

`scripts/decide_m2d7b_evidence_dossier.py` consumes the real D6/D7 schema JSON, D6B barcode-join summary JSON, and an independently curated experiment-design dossier JSON. It refuses advancement unless source hashes match and the dossier provides human-reviewed control/replicate/provenance statements with supporting records.

```bash
python scripts/decide_m2d7b_evidence_dossier.py \
  --schema results/m2d7/observed_geo_schema.json \
  --join results/m2d6b/barcode_join_summary.json \
  --dossier local/verified_experiment_design.json \
  --output results/m2d7b/readiness.json

pytest -q tests/test_m2d7b_evidence_dossier.py
```

**No automatic biological certification.** Even a complete dossier is only `eligible_for_review`; experiment/source independence remains an expert methodological decision. Blank, fabricated or placeholder references must not be provided.

## Next experimental action

1. Execute the real TF-map audit workflow and save its actual artifact, checksum and run URL.
2. Obtain and hash the GSE216595 H5AD source. Inspect its real cell/sample metadata; disk-backed inspection supports `--h5ad`.
3. Perform an exact, sample-aware D6B TF map-to-cell join. Establish the meaning of mismatches and assignment conflicts.
4. Independently examine experimental methods to distinguish cell barcodes, technical libraries, replicate cultures and suitable matched negative controls.
5. Supply auditable evidence records for these decisions. Then run the D7B gate and review whether perturbation effect estimation is methodologically defensible.

Until then, **do not generate binary regulatory truth labels, independence AUROC/AUPRC or adaptive/RL claims**.
