# Author review notes — manuscript draft v0.1

Date: 28 November 2026

These editorial notes accompany the manuscript and are not part of its scientific results.

## Proposed paper and central claim

**Title:** A Controlled Evaluation of Graph Conditioned Evidence Weighting for Gene Regulatory Network Inference

**Scope:** a controlled evaluation of a ridge-based, one-step source-utility weighting procedure, followed by a developmental real-expression comparison and a diagnosis of evidence influence. The completed studies do not evaluate the originally envisioned reinforcement-learning system. The paper does not claim primary superiority, causal recovery, or independent biological validation.

The principal supported conclusion is that the tested explicit graph-state representation did not establish an incremental advantage under the declared primary criteria. A second, sharply scoped finding is that a target-constant additive accessibility term cannot change within-target regulator selection, even when a shared learned gate assigns it a positive weight. The manuscript reports the favorable secondary comparisons as well as the unfavorable primary and fixed-weight comparisons.

## Author information and declarations to complete

The draft does not invent authors, affiliations or declarations. Before submission, the authors should supply and approve:

1. Author names, order, affiliations, corresponding author and contact information; ORCID identifiers where desired.
2. Contributions, preferably using the selected journal's contributor taxonomy, including responsibility for the scientific design, software, analysis and writing.
3. Funding and grant identifiers, acknowledgements and conflicts of interest. An unverified statement of “none” should not be inserted automatically.
4. The applicable ethics and consent statement for secondary analysis of the cited public data, after reviewing the source studies and the chosen journal's requirements.
5. Data and code availability wording, a durable tagged release and archive identifier if established, and confirmation that redistribution of every included input is permitted. Current repository and accession references identify the sources; workflow artifacts have finite retention.
6. Any disclosure of AI-assisted writing, coding or analysis required by the chosen journal, stated accurately and approved by the responsible human authors.

## Scientific decisions before submission

The full draft can now undergo author review. The following decisions determine its submission scope:

- **Paper category and journal fit.** The present evidence supports a controlled evaluation and diagnostic account. The novelty claim should remain narrower than adaptive prior integration, which is already addressed by related work.
- **External algorithm comparisons and broader simulations.** GENIE3, GRNBoost2 and other published integration approaches are discussed but were not executed here. Larger networks and materially different generators would test generalization beyond the seven- and nine-gene setting.
- **Independent biological evidence.** A qualified binding or perturbation panel has not been used to score real inferred edges. The manuscript therefore contains no real-data edge AP, AUROC, F1 or causal-recovery claim. Metadata acquisition limitations do not establish absence of a relevant assay.
- **Final within-cohort prediction.** The 133 reserved K562 test cells remain unevaluated. Any later use requires a locked final model and analysis plan. This readout would remain within one cohort and would not replace independent biological or cross-cohort validation.
- **Prior and assay qualification.** The exact raw-UMI status, active-transcript TSS interpretation, full assay QC and possible overlap between prior derivation and future reference experiments require further qualification before stronger biological claims.

These items should not be filled with guessed results or silently reframed as completed work. Additional experiments, if selected, should use new locked evaluation outcomes; the manuscript currently describes only completed analyses.

## Review order

Start with the abstract, Tables 2–4 and the conclusion to confirm the intended claim. Then review the Methods distinctions between synthetic global ranking and K562 scaffold-plus-insertion prediction. Finally check source provenance, Supplementary Table S1 and the complete supplementary metrics. Finalize author declarations and journal-specific presentation only after the scientific scope is agreed.
