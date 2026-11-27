# BIO-V6 — Locked K562 DNase promoter-accessibility evidence
**Protocol:** GCEV-BIOV6-K562-DNASE-2026-11-27-001
**Freeze:** 27 November 2026, before inspecting DNase peak coordinates, promoter overlap or any new model results.

## Biological source and locked cohorts
- Primary untreated K562 DNase-seq assay: **ENCODE ENCSR000EOT** (UW K562 wild type, DNase-seq, ENCODE2). Two annotated isogenic biosample replicates must be independently confirmed from the ORIGINAL ENCODE JSON (`replicate -> library -> biosample`). Reference: https://www.encodeproject.org/experiments/ENCSR000EOT/.
- Choose exactly **one** released GRCh38 `bed narrowPeak`, `output_type=peaks`, uniquely numbered one-biological-replicate file **for each of replicates 1 and 2**, lexicographically smallest eligible ENCODE accession within each replicate (avoid post-score file selection). Record all eligible files, selected files, original assay JSON SHA256, biosample IDs, file-accession MD5/size, date and reported audits. If either source replicate cannot be verified, **NO-GO**. No automatic experiment substitution.
- Original BIO-V5 score source: final deterministic GitHub workflow run **45988423247**, artifact **15192626844** (and independent run 45988869627). 64 frozen **1,201-bp** GENCODE V29 gene-oriented GRCh38 promoter windows; 16 TF × 64 genes minus self = **1,008** candidate pairs. The original motif prior CSV source SHA256 is `ea0bc14356fed5d2a8d4b5452e94e2372dc2f6e2da1d990a151ab1283a86b4e0`. Original K562 expression/target selection and 394/131/133 cell splits unchanged.

## Original peak processing and source eligibility
- Use only original checksummed **DNase-seq GRCh38 narrowPeak BED** from the two distinct biological replicates, never expression, TF ChIP, CRISPRi responses, or external regulatory gold labels.
- Parse at least first 3 BED fields (`chrom`, `start0`, `end0`) from gzip, validate every noncomment row `0 <= start < end`, expected primary chromosome, no inconsistent/mixed assembly, and numerical narrowPeak fields where present; do not fabricate unsupported signals. Do not use peak-call scores as calibrated accessibility probabilities. Preserve file md5, sha256, original row count, and per-chromosome interval totals.
- For each originally fixed promoter interval, compute **fraction of its 1,201 genomic bases covered by the union of intersecting peaks** in each replicate independently: `coverage_r=|promoter ∩ union(peaks_r)|/1201`. This is invariant to overlapping peaks and genomic strand. `accessibility_mean=(coverage_1+coverage_2)/2`; `accessibility_min=min(coverage_1,coverage_2)`; `peak_present_both=(coverage_1>0 and coverage_2>0)`; `replicate_absolute_disagreement=abs(coverage_1-coverage_2)`. Do not exclude promoters without calls or count a called zero as missing evidence.
- As a preregistered **descriptive context-gated motif score**, `joint_motif_dnase_product=motif_prior_score * accessibility_mean` for each TF→target. It is not an extra *independent* evidence stream; it is a deterministic product of the motif and accessibility components. `accessibility_mean` is **target-specific and broadcast** to the 16 eligible TFs, not TF-specific evidence of occupancy. The two primary inputs can be presented separately; the product is diagnostic only.

## Required original-data outputs
- Original ENCODE metadata source audit and exactly two checksummed peak files with different biosample IDs, GRCh38 assembly, DNase assay and acceptable file-specific provenance.
- Exactly **64** promoter rows with source coordinates and per-replicate union overlap, and exactly **1,008** edge rows aligned one-to-one to original BIO-V5 TF/gene directed candidate pairs.
- Save `k562_dnase_promoter_accessibility_64.csv`, `k562_tf_target_accessibility_and_motif_context_1008.csv`, `biov6_source_audit.json`, per-file and derived-data SHA256 and accession/biosample/assembly provenance.
- An independent audit recomputes all 1,008 edge values by joining original promoter accessibility to original motif prior and checks that `joint==motif*accessibility`, that no self-edge or source-missingness changes occur, and that the same promoter has the same accessibility across all TFs.

## Interpretation and research decision
A peak-free promoter is **not** proven inaccessible, and positive DNase overlap is **not** confirmed TF binding or direct regulation. Both evidence sources have distinct measurement mechanisms (sequence motif versus empirical K562 chromatin accessibility), but share promoter coordinates, so they are **not statistically independent**. Neither can be used as an independent TF→target gold standard.

This study must remain blinded to already inspected ENCODE CRISPRi responses; no model fitting or performance comparisons occur here. Next would be prospective fixed GCEV/graph-independent/fixed-weight predictions with an independent uninspected perturbation panel and an independent TF-specific source, if feasible.

**Stop rule:** If no two appropriately provenance-backed source files exist for ENCSR000EOT, mark NO-GO. Any fallback must be separately documented and prospectively locked, not selected by attractive overlap scores.
