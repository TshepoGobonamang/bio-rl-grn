# BIO-V2 — Authentic ENCODE K562 CRISPRi replicate quality report
Date: 2026-11-26

## Source and execution
- **Successful GitHub QC job:** [workflow run #45800397436](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45800397436), job `original-source-qc`, artifact `gcev-encode-biov2-original-source-qc`, GitHub artifact ID **15106219232**. The R/DESeq2 job is a separate step; source QC success does NOT imply DE statistical completion.
- Frozen original panel: six perturbations TFDP1, KLF2, NR2F2, RELA, SP1 and NR2C2 against shared ENCSR095PIC non-targeting controls.
- Original replicate-level ENCODE RSEM quantification TSVs: 17/17 source-verified against published file sizes and MD5s; uniquely identified biological biosamples.
- Official GENCODE V29 GRCh38 GTF SHA256 **e77a9568d8783f432a03f5e49d8f52d89e8b0e795382cfcb1d93bbd279b0c2b9**.
- Original output: 59,526 common gene-quantification identifiers, including **58,780 common ENSG IDs**, of which **58,721** matched the exact versioned GENCODE v29 gene annotation.
- 2 isogenic control biosamples, 2–3 distinct isogenic treatment biosamples per TF. These are NOT independent donors.

## Original sample QC (descriptive, unfiltered)
- Total per-sample RSEM estimated counts: **11.63–29.36 million**. These are estimated, unnormalized expression counts, not independently assessed read-mapping quality.
- Within-experiment replicate Spearman correlation on log2(TPM+1) across 58,780 common ENSG genes: **13 pairs**, min **0.8403**, median **0.8646**, no automatic replicate exclusion.
- Median/min pair correlations by experiment: control 0.8623; TFDP1 0.8945; KLF2 0.8674; NR2F2 median 0.8410 / min 0.8403; RELA 0.8723; SP1 median 0.8646 / min 0.8450; NR2C2 median 0.8489 / min 0.8453.
- PCA on 2,000 most variable log2(TPM+1) genes: PC1 **26.13%** and PC2 **9.58%** of sample-level variance. PCA loading/batch structure needs inspection.
- The original ENCODE QC/audit object was unavailable/empty in the retrieved response; **manual ENCODE audit, library layout, mapping metrics and batch verification remain pending**.

## Target TF mRNA reduction (TPM ratio with +1 pseudocount; descriptive)
| CRISPRi TF | Versioned GENCODE V29 ENSG ID | Mean control TPM | Mean treated TPM | Approx. reduction* | Treated replicates below control mean |
|---|---|---:|---:|---:|---:|
| TFDP1 | ENSG00000198176.12 | 54.665 | 9.865 | 80.5% | 2/2 |
| KLF2 | ENSG00000127528.5 | 3.135 | 2.025 | 26.8% | 2/2 |
| NR2F2 | ENSG00000185551.14 | 122.940 | 79.357 | 35.2% | 3/3 |
| RELA | ENSG00000173039.18 | 108.355 | 39.265 | 63.2% | 2/2 |
| SP1 | ENSG00000185591.9 | 25.085 | 21.057 | 15.4% | 2/3 |
| NR2C2 | ENSG00000177463.15 | 22.170 | 18.440 | 16.1% | 2/3 |

*Computed as 100×[1−(treated mean TPM+1)/(control mean TPM+1)]. This is NOT verified functional protein knockdown; transcript abundance can be noisy and CRISPRi may alter growth state. SP1 and NR2C2 have one treated replicate with TF TPM above the control mean; retain all original replicates but flag both for manual QC review. KLF2 has low control baseline expression and modest apparent change.

## Consequences for next stage
1. Count-based DESeq2 is being executed in a *separate job* using the original gene-level RSEM estimated counts, TPM and effective-length offsets via `DESeqDataSetFromTximport`, with low-count and zero-length prefilters. **Do not label the earlier TPM difference thresholds statistically significant.**
2. If DESeq2 completes, report BH-FDR, effect sizes, sample-level QC and exact testable genes. Absence of significance is not a negative TF edge; significant downstream RNA change may be indirect.
3. No GCEV predictor rankings were generated or compared. Do not claim GCEV biological validation, direct binding or independent prior/reference disjointness from this QC report.
4. With only six treatment TFs sharing two control biosamples, TF-level validation is constrained by shared-control dependence and should be framed as exploratory, not independent population confirmation.

## Decision
**PASS — Original replicate-level data integrity, gene-ID mapping and descriptive TF-mRNA QC.**
**REVIEW — weak/heterogeneous TF target-mRNA suppression for multiple perturbations, ENCODE QC and batch provenance.**
**NOT YET COMPLETE — DESeq2 statistical results, direct regulatory validation and GCEV model comparison.**
