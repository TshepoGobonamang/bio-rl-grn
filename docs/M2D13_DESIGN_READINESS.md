# M2C-D13 — Experimental control and replicate qualification

This milestone consumes a genuine `m2d12-observed-tf-design/observed_tf_design.json` from the [M2D12 workflow](https://github.com/TshepoGobonamang/bio-rl-grn/actions/workflows/m2d12-tf-design.yml). It checks that all TF / batch / temp counts reconcile, reports TF representation across batches, and **requires documented external evidence** before calling an observed label a valid control or an observed batch an independent biological replicate.

## Run

```bash
python -m pytest -q tests/test_m2d13_design_readiness.py
python scripts/decide_m2d13_design_readiness.py \
  --d12 results/m2d12/observed_tf_design.json \
  --evidence configs/m2d13_evidence_UNVERIFIED.json \
  --output results/m2d13/design_readiness.json
```

The supplied manifest intentionally has all flags false and therefore yields NO_GO. Do not change flags to true merely to unblock the pipeline. Every supported control label and design assertion must be cross-checked with the original experimental sample sheet/publication and cite an exact source. A TF appearing in multiple sequencing batches is **not** proof of multiple biological cultures. A GFP-looking label is **not** proof of a matched empty-vector control, and an overexpression response is not direct binding.

## Gating decisions

- Without the **actual** D12 JSON (source hash, full counts and original metadata labels), the gate must not be presented as having been run on real data.
- Without assay-matched controls and independent culture replicate provenance, no confirmatory differential response testing.
- Even if control, replication, barcode assignment and matched cell context pass, the result is only `ELIGIBLE_FOR_FURTHER_DESIGN_REVIEW`. Independent validation additionally needs a source-disjoint regulatory standard and directness assessment.
- GSE216595's pooled HUES66-to-iNP differentiation context differs from BEELINE hESC. Context-matched inference cannot be assumed.

**Status: code committed, synthetic tests authored; real D12 JSON has not been supplied to this milestone.**
