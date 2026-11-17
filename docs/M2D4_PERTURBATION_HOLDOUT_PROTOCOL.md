# M2C-D4: external TF perturbation assay intake and label-blind candidate universe

## Candidate dataset verified at GEO (17 Nov 2026)
**Primary candidate:** [GSE216595](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216595), subseries of GSE216481, human HUES66 embryonic stem cells; 90-TF overexpression library; 7-day scRNA-seq. [GSM6681043](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6681043) documents RNA-seq, hg38 alignment, a processed h5ad file and TF-to-cell-barcode assignments. BioProject PRJNA894481.

**Secondary exploratory candidate:** [GSE349847](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349847), hESC CRISPRi dose-response of NANOG/OCT4; use only after independently confirming public downloadable processed files and full control design.

**Warning:** TF overexpression response is evidence of downstream perturbation sensitivity, **not proof of direct TF binding or causal edge specificity**. Seven-day responses may be mediated by differentiation and indirect cascades. Distinct hESC lines HUES66 and BEELINE hESC cannot be assumed to represent identical states.

## Status and dataset requirements

M2D4 includes the expression-only candidate selection script, but **does not yet ingest raw perturbation data or produce experimental labels**. Source files must be fetched, hash-pinned, and their metadata inspected separately. No unsupported result or independence claim is made.

Before opening any assay expression or response labels:
1. Lock TF catalog with a documented source independent of CollecTRI and the BEELINE ChIP reference. A prior-derived TF set is prohibited for confirmatory candidate selection. Check if annotation catalogs contain overlapping evaluation-derived design choices.
2. Use only original hESC expression to rank genes by variance and choose a deterministic fixed TF/candidate universe. Record expression and catalog SHA256, freeze candidate file SHA256, record commit and date.
3. Lock the model architecture, prediction scores, expression split and evaluation thresholds; hash outputs.
4. Obtain public perturbation cell-to-TF mapping, nonperturbed/control cells, expression counts, replicate/batch identifiers, perturbation efficiency and assays. Distinguish biological replicates from cells before calculating significance. Never treat single cells as independent experimental biological replicates.
5. Use controls matched for batch and lineage; assess off-targets and indirect changes. No response-difference threshold may be chosen by reference to frozen predicted edges.
6. Check whether any candidate prior/source was trained or curated from the same perturbation experiments. Missing accession information remains unresolved.
7. Predeclare two outcomes: (a) held-out perturbation response retrieval (positive response versus explicitly tested non-response, adjusted for experimental testing coverage); and (b) **separate** direct binding evaluation, only when an independent binding assay is available.

## Code and reproducible command

`scripts/prepare_m2d4_candidate_universe.py` creates label-blind directed candidates and a manifest. Source gene symbols must match exactly; ambiguous symbols and aliases need a separately audited mapping plan. The required TF catalog is a one-column CSV named `gene`.

```bash
python scripts/prepare_m2d4_candidate_universe.py \
  --expression local/hesc_expression_cells_by_genes.csv \
  --tf-catalog local/independently_defined_tfs.csv \
  --tf-catalog-origin "PREDECLARED_SOURCE_NOT_COLLECTRI_OR_REFERENCE" \
  --n-genes 128 --n-tfs 24 --output results/m2d4-candidate-freeze

pytest -q tests/test_m2d4_candidate_universe.py
```

The TF catalog origin argument is **only a logged assertion**, not automatic verification; reviewer confirmation and source ID/provenance are required. The expression input must be obtained independently of perturbation labels. Expression variance selection is a preliminary **pipeline default**, not a claim that 128 genes and 24 TFs are optimal.

## Acceptance criteria for next implementation

- Downloaded source file and recorded SHA256, licensing, sample manifest and stable metadata for assay and controls.
- Demonstrated per-perturbation cell assignments plus experimental replicate structure; a documented negative-label strategy.
- Frozen model outputs and candidate SHA before any assay-label matching.
- No source-overlap certification until experiments, accessions, and source families are reviewed.
- Data preparation creates missing masks and does not label 'not observed' as negative.
- Report statistics only after controls and matched source provenance gate pass.

**No AUROC, AUPRC, RL reward or causal TF binding claim is valid at this stage.**
