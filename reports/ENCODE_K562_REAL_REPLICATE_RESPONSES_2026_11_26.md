# Original ENCODE K562 CRISPRi expression data: source-verification and descriptive responses

**Date:** 26 November 2026  
**Data state:** ORIGINAL ENCODE EXPRESSION PROCESSED; confirmatory GCEV biological validation **NOT YET PERFORMED**

## Authentic execution and artifacts

- [Live completed GitHub Actions run #45797283352](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45797283352) — `source-metadata` job **success**.
- GitHub [derived descriptive response artifact ID 15106010909](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45797283352/artifacts/15106010909).
- GitHub [source metadata audit artifact ID 15106390003](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45797283352/artifacts/15106390003).
- Original data: seven predeclared ENCODE K562 CRISPRi RNA-seq experiment accessions; two non-targeting isogenic control biosamples; 15 targeted TF biosamples total. Exactly **17** original per-biosample RSEM gene-quantification TSVs.
- All file paths, replicate identifiers, MD5s and original file SHA256s recorded in `original_source_response_manifest.json`. Files were downloaded and individually validated against ENCODE published sizes and MD5 checksums.
- Selected compatible analysis version: **GRCh38, GENCODE V29**. The version was selected after the original-source metadata inspection and **before** any expression-response calculations, because the source offers V24 and V29 quantifications. The six-TF panel was frozen 25 November 2026.
- Original response ZIP SHA256: `8c3bd7332b8574db65fea3a7fd4ebe601f3f557e3d8739656369238f7083652d`; source metadata ZIP SHA256: `4f8fccb75c048717686384b42f112ee185fc174bdd6e1fd416b58f937b32746a`.

## Experiment provenance and per-perturbation treatment/control counts

| TF | Original ENCODE accession | Treatment biosamples | Control biosamples |
|---|---|---:|---:|
| TFDP1 | ENCSR552HTC | 2 | 2 |
| KLF2 | ENCSR620DOP | 2 | 2 |
| NR2F2 | ENCSR688UTN | 3 | 2 |
| RELA | ENCSR052ZPX | 2 | 2 |
| SP1 | ENCSR650QAK | 3 | 2 |
| NR2C2 | ENCSR808HQN | 3 | 2 |

Shared non-specific-target control: **ENCSR095PIC**, exactly the same two control biosamples reused across all six comparisons. Distinct ENCODE biosample records were resolved by the canonical Experiment → Replicate → Library → Biosample hierarchy. They are **isogenic cell cultures**, not six independent populations/donors.

## Descriptive response results (actual GEO/ENCODE data, not model predictions)

Each treatment/control comparison preserved **59,526 exact common gene identifiers** in original GRCh38/V29 gene-quantification tables. For each gene, define the descriptive log2 ratio
`log2((mean treatment TPM + 1)/(mean control TPM + 1))`.
There is no DESeq2/edgeR, FDR, posterior probability, direct TF binding, inferred causal regulatory direction, or statistical significance in these values.

| TF | Mean absolute log2 TPM ratio | Median absolute ratio | Genes |\|log2 ratio\| ≥0.5 (unadjusted) | Genes |\|log2 ratio\| ≥1 (unadjusted) |
|---|---:|---:|---:|---:|
| TFDP1 | 0.154361 | 0.011232 | 5,807 | 1,040 |
| KLF2 | 0.102135 | 0.014355 | 2,229 | 457 |
| NR2F2 | 0.128841 | 0.028569 | 3,666 | 689 |
| RELA | 0.124382 | 0.014355 | 3,577 | 653 |
| SP1 | 0.163686 | 0.026358 | 6,324 | 1,240 |
| NR2C2 | 0.134296 | 0.024833 | 4,108 | 669 |

**Important:** Absolute-ratio thresholds are arbitrary descriptive diagnostics, NOT differentially expressed gene counts, complete regulatory labels or significance results. Many genes have low or zero TPM; pseudocount selection affects ratios, especially at low abundance.

## Remaining qualification and scientific limits

1. The ENCODE portal reports audit not_compliant/warning flags for at least some source records, while the original API frame used by this auditor lacks audit metadata. **Manual QC review is still required.** The program explicitly sets `source_qc_review_pending=True` and no analysis assumes QC has passed.
2. Actual knockdown efficacy has **not** been verified by correctly mapping each target TF to the matched ENSEMBL gene ID and examining its replicate-level expression. No overexpression/CRISPRi response is direct binding proof.
3. Results are descriptive TPM mean ratios, not differential-expression inference with dispersion, FDR and effect-size uncertainty. The two shared controls cannot be counted as separate independent experiments for every TF.
4. **There are no source-disjoint K562 model predictions here.** GCEV, graph-independent learned weighting, fixed-prior and expression-only method rankings were NOT trained/frozen on independent K562 expression, so biological AP/AUROC or advantage over fixed weights CANNOT be computed from this response artifact.
5. Unresponsive genes are not proven non-regulatory edges; significant knockdown effects can be indirect.

## Next experimental contract

- Review original quality-control reports and knockdown efficacy; lock exact annotation and ENSEMBL-to-symbol mapping.
- Identify separately sourced K562 gene-expression training data and source-disjoint priors; freeze candidate TF-target universe and all GCEV/fixed/no-state rankings **before looking at these treatment-response outcomes for model selection**. These six TF responses are now inspected, so reserve a separate biological holdout for a truly prospective confirmation.
- Compare frozen ranks with continuous absolute perturbation-response effects (per-TF Spearman / rank-enrichment); no binary negative truth from missing or nonsignificant effects.
- Report independent TF-level counts and note shared controls and QC restrictions. Do not claim direct causal regulation or RL superiority.

**Decision:** Original ENCODE per-replicate expression intake **PASS for exploratory descriptive response analysis**. **NO-GO for independent biological evidence of GCEV superiority** until separate predictions, QC, knockdown and source disjointness are qualified.
