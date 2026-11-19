# M2C-D13 — Real M2D12 artifact review (19 Nov 2026)

## Provenance

User-supplied artifact: `m2d12-observed-tf-design.zip`, containing `observed_tf_design.json`.

- ZIP SHA256: `35b14f96a04f5b9ceae882b76f21b2e0832bca3c4cef1a5f8e5e28992712f1fa`
- Referenced original GEO H5AD SHA256: `a2e06453911473a0b816cf943b6730a718670a4e082a497b578e99ce637ecccd`
- Source name in D12: `GSE216595_180124_perturb.h5ad.gz`
- Cell observations: **60,997**
- TF construct/isoform labels: **87**
- Distinct `batch` values: **8**
- Distinct `temp` values: **6**
- Observed TF × batch × temp groups: **682**

## Integrity checks: PASS

All TF, batch, temp marginal sums and 682 joint-count sums equal **60,997** cells.

TF labels across `batch`: **78** labels appear in 8 batches; **6** in 7; **1** in 6; **2** in 5. These are cell assignments across sequencing batch labels, **NOT independently verified cultures**.

The five largest construct labels are HOPX-NM_139212 (6,645 cells), MAZ-NM_001276276 (4,472), PAX6-NM_001310160 (2,932), CENPA-NM_001042426 (2,844), and ID4-NM_001546 (2,142).

## Critical metadata anomaly: `temp`

The value `temp` (literal string) occurs in **59,449** cells. For the other **1,548** cells, the `temp` value is one of the following **five TF labels** and **exactly equals their `TF` field** in all reported joint groups:

- FEZF2-NM_018008: 703
- SOX9-NM_000346: 378
- EMX2-NM_004098: 230
- ASCL1-NM_004316: 155
- EOMES-NM_005442: 82

Thus `temp` cannot defensibly be described as experimental temperature or six assay conditions. The precise column construction and meaning require authors' processing code or source documentation. Do not silently recode, correct, drop or infer a treatment/control field.

## Evidence qualification: NO-GO

- D12 `possible_control_label_candidates_not_verified` = [].
- Actual identity of appropriately matched experimental negative controls: **not established**.
- Biological culture replicate IDs and independent culture design: **not established**.
- TF-cell source-barcode mapping reconciled: **not established**.
- Independent regulatory edge truth or source-disjoint gold standard: **not established**.
- Cell-context comparability to BEELINE hESC: **not established** (HUES66 perturbation screen is neural-progenitor differentiated).

This review supports **descriptive metadata quality analysis only**; it does not authorize differential expression estimates, parent-conditional ground-truth regulatory edge metrics, causal claims, or reinforcement learning reward fitting.

## Immediate next scientific task

Investigate `temp` provenance and inspect original paper/TF barcode mapping. If there is no documented valid negative control and independently replicated culture design, seek an alternative, better-controlled independent perturbation dataset rather than manufacturing control labels. Retain this GSE216595 dataset for exploratory TF-construct coverage and potential context-shift analysis.

This report is based on an uploaded audit JSON, **not direct reexecution** of the original H5AD on this review branch.
