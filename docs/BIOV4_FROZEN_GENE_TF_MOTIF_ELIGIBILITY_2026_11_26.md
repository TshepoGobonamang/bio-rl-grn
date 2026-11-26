# BIO-V4 — Frozen K562 gene/TF eligibility and motif-source audit
**Protocol:** GCEV-BIOV4-GENE-TF-JASPAR-2026-11-26-001
**Status:** frozen before opening GENCODE/JASPAR matched-gene results, not a GRN accuracy evaluation.

## Unmodified prior inputs
- Authenticated independent K562-only scRNA: GEO GSM9224207 original H5AD SHA256 `14eb3129e57769899e0d301f515531a7e410e52dc7488fb1f981c43b2ea99ca6`.
- Original-fit-only `selected_gene_symbols.csv` has exactly 512 selected names; original split RNG 20261126; fit/internal-score/untouched-test n=394/131/133. Source and selected-gene SHA256s from `split_manifest.json` are recomputed.
- GENCODE **release 29**, human GRCh38 primary-chromosome gene annotations, source `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz`. Match exact gene symbols (case-insensitive only when a single unambiguous GENCODE gene ID matches), retain versioned Ensembl IDs and gene_type.
- JASPAR **2026 CORE Vertebrates non-redundant MEME PFMs**, `https://jaspar2026.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt`. Match only single-factor motif names to exact unambiguous GENCODE gene symbols; exclude TF complexes and ambiguous aliases.

## Required audit and decision rules
1. Check source file hashes against the verified BIO-V3B split manifest, retain 512-gene original list **unchanged**.
2. Map selected gene names to GENCODE V29 gene IDs; reject ambiguous duplicates, unresolved symbols and nonprimary chromosomes. Report protein-coding genes, mitochondrial/ribosomal symbols, and other annotations. GENCODE may contain multi-locus genes with a shared symbol: these are unresolved, not guessed.
3. Parse genuine JASPAR MEME motif IDs and names, reject multi-TF motifs (`::` etc.), deduplicate motif IDs and match only valid single-factor gene symbols. Do not treat a PFM presence as an observed promoter target edge.
4. For each TF with matched JASPAR motif, verify it is one of the 512 expression-selected genes and expression is nonzero in >=5% of **394 fitting cells**. Report the number of potentially eligible regulators. The source matrix may be count-like but unverified raw counts; this is an expression eligibility check, not RNA differential expression.
5. If fewer than 4 eligible regulators remain, mark current 512-gene selection `NO_GO_FOR_GCEV_MODEL`; do not invent regulators. A separate prospective TF-anchored gene selection would then be required, keeping original partitions and excluding held-out ENCODE responses.
6. Preserve `original_gene_symbol`, source-version gene ID, biotype, motif matrix ID and detection fraction with SHA256 provenance. No TF→target scores unless actual GRCh38 promoter sequences are scanned with fixed PWM threshold and a declared promoter/TSS policy. The next independent stage must pin promoter windows and a reference FASTA before scoring.
7. **Never access** the already-inspected ENCODE K562 CRISPRi responses while selecting genes, regulators, promoter priors, or fitting GCEV. Those six TFs are development-consumed, not untouched final validation.

## Expected outputs
`biov4_qualifying_gene_annotations.csv`, `biov4_single_tf_motifs.csv`, `biov4_eligible_regulators.csv`, `biov4_source_eligibility.json`. All biological-prior scores remain uncomputed.

### Sources
- https://www.gencodegenes.org/human/release_29.html
- https://jaspar2026.elixir.no/downloads/
- Original expression source `reports/BIOV3B_REAL_K562_ONLY_MATRIX_2026_11_26.md`.

**This gate tests annotation and motif availability, not biological accuracy, causal direction or source-disjoint reference independence.**
