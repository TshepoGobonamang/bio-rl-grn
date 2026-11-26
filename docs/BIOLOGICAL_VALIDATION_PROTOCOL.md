# GCEV: independent perturbation response validation (pre-outcome freeze)

**Protocol:** GCEV-BIO-ENCODE-K562-2026-11-25-001 (25 November 2026).

## Research question

Does ranking of K562 TF→target candidates by **frozen GCEV**, graph-independent adaptation, or fixed source weighting associate with independently obtained transcriptional responses to CRISPRi targeting that TF? This measures **perturbation-response plausibility**, not direct regulation or causality of each predicted edge.

## Predeclared assay panel and control

The ENCODE3 K562 CRISPRi RNA-seq screen has a shared **non-specific-target control** `ENCSR095PIC` with two published *isogenic* replicate biosamples. Lock six TF knockdowns before examining their expression: TFDP1 `ENCSR552HTC`, KLF2 `ENCSR620DOP`, NR2F2 `ENCSR688UTN`, RELA `ENCSR052ZPX`, SP1 `ENCSR650QAK`, NR2C2 `ENCSR808HQN`.

ENSEMBL IDs and gene-symbol mappings must be pinned to the same GRCh38 annotation before any prediction/responses are joined. These are **bulk replicate-level** experiments, not scRNA-seq, and may differ from any training expression cell line/context. Isogenic biological replicates are not distinct human donors.

Reference links: https://www.encodeproject.org/experiments/ENCSR552HTC/ ; https://www.encodeproject.org/experiments/ENCSR620DOP/ ; https://www.encodeproject.org/experiments/ENCSR688UTN/ ; https://www.encodeproject.org/experiments/ENCSR052ZPX/ ; https://www.encodeproject.org/experiments/ENCSR650QAK/ ; https://www.encodeproject.org/experiments/ENCSR808HQN/ ; https://www.encodeproject.org/experiments/ENCSR095PIC/ .

## Source-qualification gates, in order

1. Fetch **original ENCODE JSON** for all seven accessions with source SHA256 and audit flags, without accepting search snippet or handwritten manifest as sufficient for actual file qualification. Verify released, CRISPRi RNA-seq, K562, expected target/control, exact control accession, and >=2 **distinct biosample IDs** across nontechnical isogenic replicate numbers. Do not substitute technical replicate/library/runs for biological units.
2. Enumerate real `file_format=tsv`, `output_type=gene quantifications`, `status=released`, `assembly=GRCh38` files; require file-specific `biological_replicates`=one replicate and matching gene annotation between treatment/control. If only aggregate files are available, mark FAIL and investigate raw/alternate accessions rather than pseudo-replicating a single table.
3. Inspect ENCODE not_compliant audit and other QC warnings, knockdown efficacy and treatment/control batch provenance. The assay was publicly released despite some audit flags; these must be reviewed before effect estimates are qualified.
4. Download true per-replicate quantifications after metadata clearance with per-file accession/URL, MD5/SHA256, unmodified transcript/gene-ID preservation. Quantification schema, successful TF knockdown and replicate concordance are independent gates.
5. Obtain a **separately sourced**, K562-matched expression dataset to fit GRN models and independently sourced input evidence (e.g., source-specific promoter motif scans rather than labels from the same ENCODE perturbation experiment). Freeze candidate TF/target universe, methods, and predicted rankings **before unsealing** any ENCODE treatment-vs-control response labels.
6. Compare all three methods on exactly the same measured TF target genes: GCEV, no-parent/graph-state adaptation, fixed weights, expression-only. Main biological test: paired per-TF rank correlation between predicted edge ranks and |perturbation log2 fold change|, with an exact gene ID join and TF knockdown verified. Treat RNA changes as possible *indirect* effects; absence of significance is **unknown**, not a proven non-edge.
7. Secondary response enrichment among reproducible BH-FDR-significant effects may be described as assay-supported response enrichment, **not AP against a complete direct-regulation reference**. With only 6 TFs and shared controls, do not infer 6 independent biological populations or strong population-level CIs. Report TF-level effect sizes, control reuse, and sensitivity to QC exclusion.

## No-go conditions

- Experiment controls unlinked or incompatible; less than 2 biological replicate biosamples per comparison; no per-replicate quantification; incompatible genome annotation; ineffective perturbation; uncertain source/readout overlap; insufficient matched TFs/target genes.
- Source dependence: If upstream CollecTRI/ENCODE supporting studies overlap ENCODE reference outcomes, do NOT claim source-disjoint validation; use an independent motif source or explicitly report unverified disjointness.
- Do not reuse GSE216595/GSE133344 sequencing lanes as culture replicates. Norman GSE133344 GEO specifically states its eight GEM groups derive from one pooled experiment.

**Current decision:** The panel is source-identified. The **live original-data audit and downstream biological validation have not yet been executed or qualified**. Keep the prior negative synthetic GCEV results unchanged.
