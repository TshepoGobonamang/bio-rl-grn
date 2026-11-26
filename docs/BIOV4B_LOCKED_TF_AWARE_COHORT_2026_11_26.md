# BIO-V4B — Frozen TF-aware development cohort after BIO-V4 coverage failure
**Decision date:** 2026-11-26
**Protocol:** GCEV-BIOV4B-K562-TF16-TARGET48-2026-11-26-001
**Type:** transparently redesigned, **development-stage** gene selection after original 512-gene source eligibility NO-GO. This is NOT a pristine unseen confirmation gene universe.

## Why a new cohort?
The BIO-V4 source audit, [run 45888943925](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45888943925), found just **3** eligible JASPAR single-factor regulators (TRPS1, PAX5, HES1) among 512 genes selected by log1p expression variance. That gene universe is **not suitable** for a substantive TF→target graph-state benchmark. Retain and report it unchanged.

## New source and selection rule locked before obtaining the TF-aware result
Use the *same original* GEO GSM9224207 K562-only H5AD SHA256 `14eb3129e57769899e0d301f515531a7e410e52dc7488fb1f981c43b2ea99ca6`, same 658 ordered source cells and unchanged RNG seed 20261126 cell partitions: **394 fit / 131 internal score / 133 untouched test**. No refitting or re-splitting.

Use the full original 60,671 gene names and original-scale sparse X. For each gene, calculate `detection_fraction = n_nonzero_fit_cells / 394` and `variance = var_fit(log1p(original_X))`, on 394 fitting cells **only**, with sparse zeros retained. Required eligibility: detection fraction at least **0.05**, positive variance, **unique** primary-chromosome GENCODE V29 gene ID and `protein_coding` type. Exclude symbol prefixes `MT-`, `RPL[0-9]`, `RPS[0-9]` (ribosomal/mitochondrial domination, planned after observing BIO-V4 quality defects).

**Regulator definition:** a measured eligible gene symbol with one or more **single-TF** JASPAR 2026 CORE vertebrate nonredundant MEME motif names matched by unique case-insensitive exact symbol. Exclude complex (`::`) motifs and ambiguous multi-locus genes. Motif availability does NOT establish TF binding in K562 or an edge.

Select exactly **16 candidate TFs** in decreasing fit-only log1p-variance order, lexicographic symbol tie break. Select exactly **48 additional protein-coding target genes**, ranked by the same fitting-only variance, excluding all 16 chosen TFs and, to keep categories distinct, all remaining motif-annotated TF-eligible genes. Keep the TFs and targets in fixed ordered lists. This yields 64 measured genes and a candidate universe of **16 TFs × 64 targets, excluding self edges = 1008 directed pairs**. If eligibility is insufficient, explicitly report NO-GO without data-driven replacement.

No ENCODE CRISPRi expression, differential-expression p-values, regulatory gold edges or prior target connections can enter feature selection. Both old and new K562 gene universes must be preserved, labelled developmental.

## Exact source locks
- GENCODE release 29 GRCh38 primary-chromosome GTF source SHA256: `e77a9568d8783f432a03f5e49d8f52d89e8b0e795382cfcb1d93bbd279b0c2b9`.
- JASPAR 2026 CORE vertebrate nonredundant MEME source SHA256: `ae4f77a17256c7b8b7b73402810c1ecf4fc411e9ad7832cdc87d586dee9a4a15`.
- Match exact symbols, no inferred aliases.
- Download original K562 complete expression from previously authenticated artifact of [run 45881788164](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45881788164).

## Output contract
Export a **64-gene** original-scale sparse matrix with the same 658 cells and a copy of the original split CSV; ordered TF list with motif provenance/GENCODE IDs, 48 target list, all selected genes with observed fitting-only variance and detection, source/selection checksum manifest, and **NO** edge-level prior scores. No GCEV or baseline model is trained here.

Next stage: GRCh38 promoter sequence/reference assembly and explicit motif scanning for the **same** selected TFs/targets, then matched GCEV vs fixed/graph-independent comparison. The already-inspected six BIO-V2 TF responses must not be used as untouched evaluation reference.

Outcome is intentionally unknown when this protocol is committed.
