# BIO-V6 — Authentic K562 DNase promoter accessibility and matched-pipeline sensitivity

**Date:** 27 Nov 2026. **Research status:** empirical K562 **promoter accessibility context** established, not TF-specific binding or biological GRN validation.

## Original source and immutable comparator

- Frozen preregistration `docs/BIOV6_LOCKED_K562_DNASE_ACCESSIBILITY_2026_11_27.md` committed `ea764120c65960d3057f07d4b59d6ad6beb6f725` **before** original ENCODE peak retrieval.
- The exact original two-biosample ENCODE K562 DNase experiment [ENCSR000EOT](https://www.encodeproject.org/experiments/ENCSR000EOT/) was source-qualified in successful [GitHub Action #45991610752](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45991610752).
- Original file artifact `gcev-biov6-authentic-k562-dnase-replicate-promoter-accessibility`, ID **15193389552**.
- Original biological biosamples: replicate 1 `ENCBS545WKJ`, replicate 2 `ENCBS392ZOW`. Original files selected by the locked lexicographic rule: replicate 1 `ENCFF185XRG` (MD5 `04653af177917b3dda96b9454fd8f90e`) and replicate 2 `ENCFF433TIR` (MD5 `90426fbf4bfa4620a61bc751bfed7fe2`). Both original GRCh38 narrowPeak BEDs were downloaded, size- and MD5-verified.
- Unchanged original [BIO-V5](https://github.com/TshepoGobonamang/bio-rl-grn/pull/49) 64 GRCh38 GENCODE v29 1,201 bp gene-level promoters and 1,008 directed, nonself TF-gene motif priors, motif CSV SHA256 `ea0bc14356fed5d2a8d4b5452e94e2372dc2f6e2da1d990a151ab1283a86b4e0`. No read access to the previously inspected BIO-V2 TF knockdown responses, no GCEV training and no changed K562 partitions.

## Original preregistered two-source observations

| Actual metric | Original |
|:--|--:|
| Fixed promoter windows | 64 |
| TF→gene candidates | 1,008 |
| Promoters with DNase peak in replicate 1 | **47** |
| Promoters with DNase peak in replicate 2 | **41** |
| Promoters with peak in **both** | **41** |
| Promoters with no peak in either | **17** |
| Promoters positive in rep 1 only | **6** |
| Replicate Spearman over 64 promoter coverage fractions | **0.7475** |
| Mean promoter-base coverage across genes | **0.171101** |
| Mean absolute difference between replicate coverage fractions | **0.125742** |
| Nonzero descriptive motif×DNase product scores among 1,008 edges | **384** |

Original peak caller **processing mismatch:** `ENCFF185XRG` used `dnase-peak-calling-step-v-1-0`, while `ENCFF433TIR` used `dnase-call-hotspots-se-step-v-2-1`. Original files are genuine and checksummed, but the differences in processed peak definitions limit quantitative reproducibility across replicates. The original result was **not discarded**.

## Secondary matched-pipeline sensitivity, explicitly post-primary

Protocol `docs/BIOV6_SECONDARY_MATCHED_PIPELINE_SENSITIVITY_2026_11_27.md` committed at `fce279af85625fb6a24b4c7273047b4c02275a67` **before** recalculating the matched versions. Successful independent [GitHub Action #45992406882](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45992406882); artifact **15194436273**.

| Descriptor | v1 matched callers | v2 matched callers |
|:--|--:|--:|
| Peak file replicate 1 | ENCFF185XRG | ENCFF821KDJ |
| Peak file replicate 2 | ENCFF784DIM | ENCFF433TIR |
| Source peak-calling version | dnase-peak-calling-step-v-1-0 | dnase-call-hotspots-se-step-v-2-1 |
| Promoters with peak in replicate 1 | 47 | 51 |
| Promoters with peak in replicate 2 | 41 | 41 |
| Promoters with peak in **both** | 41 | 41 |
| No peak in either | 17 | 13 |
| Mean promoter coverage | **0.190095** | **0.146844** |
| Replicate coverage Spearman | **0.73755** | **0.61949** |
| Mean abs. replicate coverage difference | 0.096989 | 0.080350 |
| Positive motif×accessibility products | **384** | **414** |
| Mean descriptive motif×accessibility product | 0.033351 | 0.025844 |

The two version-specific promoter coverage vectors correlate at **Spearman 0.905734** and differ by mean **0.061245** absolute promoter fraction. Selecting one based on greater nominal coverage would be post-result cherry picking. Both are transparently retained. This is **processing-method sensitivity** on the same biological experiment, not two independent biological validations.

## Interpretation
Each promoter score is `union_covered_DNase_peak_bases/1201` per original replicate. Accessibility mean is the average of replicate coverage fractions, in [0,1]. The genomic interval union prevents counting overlapping peaks twice. A peak-free promoter is **not proven inaccessible**; differences may arise from source processing, noise and threshold choices.

The existing sequence motif and experimental DNase evidence differ in assay mechanism, but **share the same promoter coordinates**. DNase accessibility is a **target-level** covariate broadcast across candidate TFs, **not TF-specific binding**, measured TF occupancy or a calibrated edge probability. The `motif_prior × accessibility_mean` product is derived from the other two and must **not** be described as a third independent source.

An exact-source adapter in `scripts/biov6_two_evidence_adapter.py` exports the two qualified source arrays with observed-zero masks into the existing GCEV `evidence` interface, preserving candidate alignment; **no model has been fitted**. Future model comparisons must include an expression-only baseline and a fair fixed/graph-independent baseline and report results before assessing an independent, uninspected perturbation panel.

Original ENCODE experiment QC warnings/peak-calling metrics, full DNAse read signal and sample batch correspondence with the separate fresh-K562 expression cohort have **not** yet been manually reviewed to publication standard; the samples are both K562 but are not from the same exact culture or collection. K562 cell line replication should not be misconstrued as donor-level biological replication.

## Decision
**PASS:** original ENCODE source identity, 2 distinct source biosamples, GRCh38 original peak integrity, reproducible **64-promoter / 1,008-edge** context matrix, source/provenance and matched-caller sensitivity.

**REVIEW:** pipeline-dependent numerical scales, complete ENCODE assay-level QC, whether experimental DNase context adds value beyond motif information and whether source coverage is meaningfully informative with all methods.

**NOT ESTABLISHED:** independently confirmed in-cell TF binding, direct TF→target regulation, biological ranking advantage of GCEV or a pristine held-out K562 knockdown panel.
