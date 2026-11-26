# BIO-V2 — Verified original-data DESeq2 statistical findings (26 November 2026)

## Source and outcome integrity

**Successful reproducibility run:** [GitHub Actions #45816301876](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45816301876), job `deseq2-prebuilt`, artifact **gcev-encode-biov2-deseq2-apt-authentic-results** (ID **15114763379**).

- Original QC was independently passed in [run #45800397436](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45800397436) (artifact **15106219232**).
- 17 original replicate-level RSEM quantification TSVs, verified against original ENCODE MD5/size checks and mapped to independently identified isogenic K562 biosamples. Two shared non-targeting control biosamples, 2–3 treatment biosamples per TF.
- Original assay: six frozen CRISPRi perturbations TFDP1, KLF2, NR2F2, RELA, SP1, NR2C2.
- Original source annotation: GRCh38/GENCODE V29; gene identities preserved as **versioned Ensembl IDs**.
- DESeq2 input: RSEM estimated gene counts, same-source TPM and effective lengths via `tximport`; fixed `sum(counts)>=10` and positive-effective-length filters. DESeq2 Wald contrast `~condition`, BH-adjusted FDR 0.05 **separately within each TF**. No six-contrast-wide FDR.
- R 4.3.3; DESeq2 1.42.0; tximport 1.30.0 (apt-packaged route). A previous invalid-connection coding error was corrected by allowing `read.csv(gzfile(path))` to manage its own connection.
- Machine-checked independent verification JSON `verification.json`: **BIOV2_REAL_DESEQ2_ARTIFACT_CONSISTENCY_VERIFIED**. Six tables' gene rows, tested p-values and adjusted-p counts, FDR/up/down totals and original-source QC identity match the reports.
- No GCEV predictions or biological regulatory AP were computed from these outcomes.

## Fixed complete six-TF results (do not exclude weak perturbations post-hoc)

| Target TF | Treated/control biosamples | Genes retained/tested | Adjusted p available | Genes at BH-FDR <0.05 | Up | Down |
|:--|:--|--:|--:|--:|--:|--:|
| TFDP1 | 2 / 2 | 18,388 | 13,397 | **5,425** | 2,651 | 2,774 |
| KLF2 | 2 / 2 | 18,825 | 10,066 | **510** | 190 | 320 |
| NR2F2 | 3 / 2 | 19,788 | 15,184 | **354** | 92 | 262 |
| RELA | 2 / 2 | 18,584 | 12,819 | **2,858** | 1,555 | 1,303 |
| SP1 | 3 / 2 | 20,264 | 14,764 | **4,491** | 2,080 | 2,411 |
| NR2C2 | 3 / 2 | 19,891 | 14,876 | **549** | 208 | 341 |

KLF2 had 18,825 nonmissing unadjusted p-values but only 10,066 nonmissing BH-adjusted p-values. NR2C2 had 19,891 retained gene rows but 19,890 nonmissing unadjusted p-values. Missing adjusted p-values must not be rewritten as `p=1` or biological negatives; DESeq2 independent filtering and statistical diagnostics require review.

## Within-sample-target transcript response (most informative QC)

Original targeted genes matched by **exact GENCODE V29 versioned ENSG ID**, retrieved directly from the 6 original DESeq2 gene-wise result tables:

| TF | GENCODE V29 gene ID | Target log2FC (treated/control) | Target gene adjusted p | Target mRNA DE at FDR <0.05? |
|:--|:--|--:|--:|:--|
| TFDP1 | ENSG00000198176.12 | **−2.4434** | **4.85e−57** | Yes |
| KLF2 | ENSG00000127528.5 | −0.6152 | **NA** (p=0.1763) | No |
| NR2F2 | ENSG00000185551.14 | −0.5654 | 0.5993 | No |
| RELA | ENSG00000173039.18 | **−1.3351** | **1.11e−34** | Yes |
| SP1 | ENSG00000185591.9 | −0.2703 | 0.4478 | No |
| NR2C2 | ENSG00000177463.15 | −0.2556 | 0.7445 | No |

**Conclusion:** Strong statistical evidence for reduced *target transcript* in TFDP1 and RELA; no BH-significant target transcript reduction for the other four. Protein depletion / downstream functional knockdown remains **unverified for all six**. The substantial thousands of FDR-positive downstream responses in SP1 (4,491) and others **do not rescue** lack of significant target suppression. They may reflect indirect effects, cellular state differences or assay batch/QC factors.

## Biological quality and inferential constraints

- This is **exploratory** because all six outcomes were inspected before a GCEV model prediction was frozen. The panel cannot retroactively count as independent confirmatory GCEV validation.
- All six TFs share the same two control biosamples; contrasts are correlated. Isogenic biosample replicates are not six independent human donors.
- Actual original ENCODE audit/QC warnings, complete sample-batch technical metadata and functional protein-level suppression require further examination. Keep all 17 samples and all six predeclared TFs in source records.
- DESeq2 BH-FDR-positive genes are **differential transcriptional responses**, potentially indirect. FDR-negative or filtered genes are **unknown as direct TF-target pairs**, not reference edge negatives.
- The methods are an RSEM estimated-count and length-offset exploratory analysis; confirm compatibility of gene-level RSEM effective-length normalization with the intended model for publication.
- No source-disjoint, independent regulatory priors and no sealed K562 GCEV vs fixed/no-parent-state predictions have been established. **GCEV biological performance remains unevaluated.**

## Next locked milestone (BIO-V3)

1. Acquire a **separate baseline K562 single-cell expression source**; audit cells, gene IDs, library conditions, perturbation status and source independence *before* any training. Candidate: [GSE182685 / GSM5534702](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5534702), whose processed file is labelled K562–BT142 and **must be subset only after actual cell metadata verification**.
2. Source an external biological prior from a provenance-audited TF motif/sequence resource, potentially [JASPAR 2026 CORE](https://jaspar.elixir.no/downloads/), noting that a PFM alone is NOT a TF-target edge and the underlying ChIP sources may include ENCODE.
3. Predeclare candidate genes/TFs and data partitions independent of the perturbation outcomes; freeze fixed weighting, no-parent-state and GCEV rankings before opening a **new, unseen** perturbation evaluation set.
4. Use TF-level association with reproducible response magnitude, and per-source provenance. Do not use the six already-seen DE tables as a purported pristine final holdout.

**BIO-V2 decision: statistical response data VERIFIED; original assay QC and biological efficacy remain conditional; independent GCEV biological validation NO-GO.**
