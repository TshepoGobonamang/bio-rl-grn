# M2C-D14 — Source-backed experimental provenance and control assessment
Date: 2026-11-19
Status: **partial provenance PASS; control/replicate qualification NO-GO**

## Evidence sources
1. GSE216595 series metadata: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216595
2. Published manuscript and STAR Methods (Joung et al., Cell 2023): https://pmc.ncbi.nlm.nih.gov/articles/PMC10344468/
3. Original authors' analysis code: https://github.com/fengzhanglab/Joung_TFAtlas_Manuscript
4. Sample example and processing statement (GSM6681046): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6681046
5. Previously inspected real D12 JSON, summarized in reports/M2D13_REAL_D12_EVIDENCE_REVIEW_2026_11_19.md.

## Source-backed facts

GEO identifies a 90-TF ORF pooled perturbation screen with 7-day cell culture and 8 sequencing samples (S1–S8, GSM6681039–GSM6681046), plus TF barcode mapping (GSM6681047). The published STAR Methods place the targeted 10x screen in HUES66 hESCs transduced at MOI <0.3 followed by differentiation into induced neural progenitors (iNPs). TF and cell barcodes were separately amplified from whole-transcriptome-amplification cDNA and sequenced to pair constructs with single-cell transcriptomes. GEO expressly states that a **custom Python script** maps cells to TF ORFs.

The manuscript describes **two GFP and mCherry control vectors in the much larger MORF library** (3,548 TF isoforms) and **GFP-controlled qPCR and immunostaining validation experiments**. Neither observation verifies that GFP/mCherry cells are included in this targeted 90-ORF 10x screen, that such cells survived its filtering, or that they are matched control cells in the published H5AD.

The manuscript documents biological replicates for **other** experiments (e.g. targeted reporter/flow-FISH or later validation panels). These counts cannot be transferred to the eight 10x sequencing samples as independently randomized culture replicates.

GSE216595 processed H5AD audit: 60,997 observed cells, 87 observed TF labels, eight batch values. The real D12 `temp` field is anomalous: 59,449 `temp` literal strings; other 1,548 values mirror their `TF` label. No evidence found that this encodes a temperature, control arm, biological replicate, or documented treatment assignment.

## Concrete unanswered questions

1. Locate the original code path creating `obs['TF']`, `obs['temp']`, and `obs['batch']`. Preserve code commit SHA and filtering/assignment criteria.
2. Compare the verified GSM6681047 TFmap columns and per-cell assignment with H5AD metadata using exact **sample-aware** join keys. Check doublets, ambiguously assigned TFs and missing barcodes.
3. Locate targeted 90-TF library sample sheet or construct design; verify **specific empty-vector/GFP/mCherry/no-virus controls** in the **same** 10x screen and their presence/absence in the processed cells.
4. Locate **culture-level independent unit IDs** and exact batch/plate/culture assignments. Sequencing S1-S8 and technical barcodes are not sufficient.
5. Verify that assay endpoints and cellular context are appropriate for the intended independent regulatory evidence. Overexpression-induced expression changes are not proof of direct TF binding.

## D14 gate decision

- **PASS:** Publication/GEO assay provenance, host cell line, pool design, seven-day treatment, sequencing and barcode-pairing mechanism established.
- **UNRESOLVED:** Real 10x source-specific controls, original `temp` transformation, independently randomized culture replicates, raw-to-H5AD barcode join.
- **NO-GO:** Treating other-assay GFP validation as a source-specific 10x negative control; treating sequencing batches as independent cultures; identifying direct TF-target causal ground truth; confirmatory AUROC/AUPRC and RL reward supervision.

**Next preferred action:** inspect the original authors' code and actual TFmap JSON/artifact; if no source-specific controls or culture replicates can be established, nominate and verify an alternate perturbation reference with published biological replicates and matched controls. No invented experimental labels.
