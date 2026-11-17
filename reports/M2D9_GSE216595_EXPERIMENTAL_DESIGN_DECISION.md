# M2C-D9 — Experimental design qualification: GSE216595 (2026-11-17)

## Decision: exploratory-source PASS; confirmatory paired-pseudobulk NO-GO

This is an evidence review based on the published methods, **not** a completed download or inspection of the GSE216595 cell-assignment table.

Sources:
- NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216595
- Primary paper and detailed STAR Methods: https://pmc.ncbi.nlm.nih.gov/articles/PMC10344468/
- BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA894481

### What the source does establish

The GSE216595 subseries describes a **pooled 90-TF ORF overexpression experiment**, cultured 7 days, with 10x scRNA-seq and separate TF-barcode amplicon mapping. The associated paper's STAR Methods specify **HUES66 hESCs transduced at MOI <0.3 and differentiated into induced neural progenitor cells (iNPs)** for the 10x screen. Thus endpoint cells should not be described as untreated/pluripotent hESC expression or directly equivalent to BEELINE's hESC measurement context.

Eight 10x sequencing sample accessions (GSM6681039–GSM6681046) plus one TF-barcode map (GSM6681047) do **not** automatically establish eight independent biological culture replicates.

The paper reports selected TF validations in H1, H9 and 11a iPSC contexts and some other experiments with reported replicates; those should **not** be transferred to this specific pooled 90-TF screen as unverified replicate/control labels.

### Locked computational consequences

1. **Do not invoke** the D5 `paired_responses` estimator on the real GSE216595 screen until independently matched TF-specific control cultures and biological replicate IDs are verified. Artificial pairing on sample/GSM/barcode is prohibited.
2. Before model evaluation, download/verify `GSM6681047_180124_TFmap.csv.gz` and `GSE216595_180124_perturb.h5ad.gz`; preserve SHA256, actual column schemas, exact sample-aware cell barcode joins, and the assignment confidence/multiplicity filters in the authors' code.
3. Identify and document whether the actual screen includes GFP/empty vector/non-targeting controls, their TF mapping and culture/batch provenance. **Unknown controls are not unperturbed labels.**
4. Audit real independent culture replicates separately from 10x lanes and sequenced cells; where the required comparisons are unpaired or pseudoreplicated, choose a validated analysis compatible with that design or classify the screen as exploratory.
5. Note **context shift**: BEELINE hESC expression features versus HUES66 differentiated iNP perturbation response. A failure to recover TF responses across cell types is not equivalent to failed regulatory directionality in matched hESCs.
6. Keep TF-induced downstream expression changes separate from direct TF→target binding, and investigate source overlap with CollecTRI before independent confirmation.

### D9 research decision

- **PASS:** published experimental provenance, organism, pooled perturbation design and timing identifiable.
- **BLOCKED:** source artifact actual parsing, valid controls per TF, independent replicate structure, accurate barcode mapping, source-disjoint independent regulatory validation.
- **NO-GO:** binary regulatory labels, confirmatory AUROC/AUPRC, adaptive prior superiority, and RL policy training.

**Next action is execution of the registered [M2D8 launcher](https://github.com/TshepoGobonamang/bio-rl-grn/actions/workflows/m2d8-launch.yml), retrieval of its genuine source-schema artifact, then a real H5AD obs metadata inspection.** If that remains blocked, consider an alternate well-controlled perturbation reference rather than inventing experimental design metadata.
