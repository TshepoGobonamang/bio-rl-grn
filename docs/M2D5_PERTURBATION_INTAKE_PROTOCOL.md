# M2C-D5 — GSE216595 TF perturbation intake and replicate-safe response pipeline

## Verified public-source facts (17 November 2026)
NCBI GEO [GSE216595](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216595) contains HUES66 human embryonic stem cells transduced with 90 TF ORFs and profiled by scRNA-seq 7 days later. GEO lists:
- `GSE216595_180124_perturb.h5ad.gz`, approximately **1.4 GB**, processed Scanpy object.
- `GSE216595_RAW.tar`, approximately **671 MB**, raw 10x per-sample files.
- [GSM6681047](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6681047): `GSM6681047_180124_TFmap.csv.gz`, approximately **532.8 KB**, TF-to-barcode assignment file.
- GSM6681039–GSM6681046: eight scRNA-seq sample accessions; **sample accession count does not equal biological replication count**.
- [GSM6681043](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6681043) documents Cell Ranger 3.1.0, Scanpy 1.7.2, custom TF assignment, reference assembly hg38.

The archive files have **not** been downloaded, their checksums have **not** been computed, and the real H5AD `obs`/`var` and TF map columns have **not** been inspected. Do not invent column mappings or mark replicates/control groups as verified from GEO alone.

## Coding delivered

`scripts/prepare_m2d5_perturbation_responses.py` works only on independently verified normalized CSV **exports of raw count values**, with:
- cells × genes counts, unique cell identifiers in first CSV column, untransformed nonnegative integer UMI counts;
- explicit metadata CSV columns `cell_id,tf,condition,replicate_id,batch`;
- `condition` exactly `control` or `overexpression`, with controls assigned to the intended comparison TF;
- `replicate_id` reflecting **biologically independent cultures** after experimental methods audit, never single-cell barcodes, library lane numbers, or technical sequencing replicates;
- within-batch paired control and TF-overexpression replicate IDs, at least two valid biological pairs for a descriptive contrast.

For each TF/batch/target, it sums per-cell raw counts into replicate-level pseudobulk, computes log2(1 + CPM), calculates treatment minus matched control for each replicate, and emits mean effect, SD and sign agreement. This is **descriptive**, not differential-expression inference. It does not generate p-values, binary labels, edge truth, perturbation-response thresholds or causal TF binding claims. An explicit `--biological-replicates-verified` flag is a recorded operator assertion, not scientific certification; published methods and sample identities need independent documentation.

```bash
python scripts/prepare_m2d5_perturbation_responses.py \
  --counts local/validated_raw_umi_cells_by_genes.csv \
  --metadata local/validated_cell_assignments.csv \
  --output results/m2d5-perturbation \
  --biological-replicates-verified \
  --min-cells 10 --min-pairs 2

pytest -q tests/test_m2d5_perturbation_responses.py
```

## Blocking checks before running on GEO data

1. Fetch and hash GEO source files *outside Git*; inspect full `obs`/`var`/TF map schema and library IDs; validate unique cell/barcode joins and TF assignment confidence.
2. Verify whether the experiment includes a clearly identified, compatible negative control (e.g. empty vector), not merely cells with missing assignments.
3. Document independent biological cultures and perturbation replicate design; **do not interpret eight GSM scRNA samples as eight independent cultures**.
4. Decide whether matching same-`replicate_id` treatment and control is valid based on actual biological design. If unpaired, use a different prespecified model rather than assigning artificial pair IDs.
5. Preserve dataset lineage, checksum, library/batch effects, normalization and provenance; freeze D4 candidate universe and model scores before assay evaluation.
6. Review TF overlap and assay-context differences (HUES66 vs BEELINE hESC), provenance and any CollecTRI evidence overlap.
7. Only then build actual assay-derived response labels, explicit missingness/tested masks and confidence intervals. Without suitable controls or replicate independence, restrict to exploratory single-cell descriptive comparisons.

## Research gate

**D5 status: intake code and synthetic tests written; real GEO file download, assay schema reconciliation, authentic biological replication and response inference pending.** No real biological response targets have yet been generated. Do not proceed to independent GRN AUROC/AUPRC, adaptive-prior superiority or RL training on the basis of this infrastructure alone.
