# BIO-V8 — Evidence qualification and external-panel feasibility

**Date:** 27 November 2026. **Project:** CIBB Recalibration / Bio-RL-GRN / GCEV.

**Outcome:** the frozen-score export and source-action diagnostics were executed, independently checked and reproduced on GitHub. External biological qualification remains **blocked**: the six original ENCODE metadata searches returned HTTP 403. This is incomplete source access, not evidence that compatible assays do not exist. No new biological scoring, predictive evaluation or RL training was performed. The 133 reserved expression-test cells remain sealed.

## 1. Research scope and execution

BIO-V7 did not establish an incremental GCEV advantage and exposed a target-constant accessibility term that could not distinguish candidate TF insertions. The current 16-TF cohort also has zero overlap with BIO-V2's six previously inspected TF perturbations. BIO-V8 therefore qualifies evidence representation and external-panel feasibility before further biological scoring.

The [design](../docs/BIOV8_LOCKED_EVIDENCE_QUALIFICATION_2026_11_27.md) and [configuration](../configs/biov8_evidence_qualification_2026_11_27.json) were published at [`9169be4`](https://github.com/TshepoGobonamang/bio-rl-grn/commit/9169be40c3658d45178ed170587f5058ad8c9942). Implementation is [`eb3ca37`](https://github.com/TshepoGobonamang/bio-rl-grn/commit/eb3ca377bbbf02e7577c3f0aab6eca31cd61085f), in [draft PR #52](https://github.com/TshepoGobonamang/bio-rl-grn/pull/52). This is a prospective follow-up design after known developmental outcomes, not retroactive independent preregistration.

The workflow reads saved fitting-stage contexts and predictions plus the existing prior tables. It has no expression-loading or model-fitting interface. The external acquisition module permits original ENCODE JSON metadata only and rejects data-download paths. These restrictions are covered by 36 focused tests. The repository's full CI and synthetic engineering smoke also passed; smoke outcomes are not biological results.

## 2. Frozen outputs

| Output | Actual size | Interpretation |
|---|---:|---|
| Existing candidate-score export | 49,280 rows | Seven methods × eight conditions × 880 residual candidate actions; exact original score text retained |
| Original graph choices | 3,584 decisions | 56 conditions × 64 target decisions |
| Source-action diagnostics | 2,688 rows | Two transforms × three prior variants × seven representations × 64 targets |
| Source-action summaries | 42 rows | Every transform/variant/representation retained |
| Candidate eligibility | 1,008 edges | Distinguishes scaffold membership and residual action eligibility by transform |
| TF coverage | 16 TFs | Available candidates, not externally measured or validated targets |

All scaffold insertion scores remain unavailable. The existing `score` is frozen for a future exploratory across-target projection; `action_score` retains the original within-target decision meaning. Neither is a probability. The different ranking tasks and target offsets are explicitly documented before external outcomes.

## 3. Actual source-action findings

Primary case: original log1p expression and the original mixed-caller DNase prior. The expression rank E, motif M and accessibility A are already fixed. The boost is 0.3. Lower or higher action-change counts are not performance outcomes.

| Representation | TF-varying source values across targets | Changed choices versus expression | Choices differing from motif | Choices differing from global attenuation |
|---|---:|---:|---:|---:|
| Expression only | 0/64 | 0/64 | 33/64 | 9/64 |
| Motif | 64/64 | 33/64 | 0/64 | 24/64 |
| Global mean(A) × motif | 64/64 | 9/64 | 24/64 | 0/64 |
| Matched A × motif | 47/64 | 8/64 | 26/64 | 6/64 |
| Permuted A × motif | 47/64 | 8/64 | 26/64 | 10/64 |
| Legacy additive motif + A | 64/64 | 33/64 | 0/64 | 24/64 |
| Accessibility alone | 0/64 | 0/64 | 33/64 | 9/64 |

Source: [all 42 summaries](biov8_evidence_2026_11_27/source_action_summary.csv). Matched and permuted interactions each change eight actions versus expression, but their own selections differ on eight targets. Equal counts do not imply equal graphs or equal usefulness.

Across all six transform/prior cases, pure accessibility preserves the expression winner in **384/384** target/cases; additive motif-plus-accessibility preserves the motif-assisted winner in **384/384**. This confirms the structural limitation under the actual candidate sets.

For a fixed target, adding the same accessibility value to every TF leaves the ranking unchanged. Multiplying motif by positive accessibility also preserves motif-only ordering, but changes the motif's strength relative to expression. Accessibility is zero for **17/64** primary promoters; the product then removes the motif contribution. Its target mean is 0.1711009055, so a global mean-scaled motif is a necessary attenuation control. The product is derived evidence, not an independent third source.

The one permutation uses sorted targets and **NumPy `default_rng` (PCG64), seed 20261128**. Its [explicit mapping](biov8_evidence_2026_11_27/accessibility_permutation.csv) is retained. No seed or coefficient was selected from performance outcomes.

### Processing sensitivities: changed actions versus expression

| Representation | Original / primary | Original / matched v1 | Original / matched v2 | CP10k / primary | CP10k / matched v1 | CP10k / matched v2 |
|---|---:|---:|---:|---:|---:|---:|
| Motif | 33 | 33 | 33 | 28 | 28 | 28 |
| Global attenuation | 9 | 10 | 9 | 9 | 10 | 6 |
| Matched interaction | 8 | 8 | 10 | 6 | 6 | 7 |
| Permuted interaction | 8 | 8 | 8 | 5 | 7 | 8 |
| Accessibility alone | 0 | 0 | 0 | 0 | 0 | 0 |

All counts have denominator 64. These are the same cohort's deterministic representation probes, not independent experiments. No representation is promoted as a winning model.

## 4. External panel: what is known and what remains blocked

The original ENCODE metadata workflow ran on GitHub and preserved all six failed requests and exact query URLs. Each returned **HTTP 403 Forbidden**. The scientific decision is `NO_GO_FOR_BIOLOGICAL_SCORING`; `search_complete=false`. A successful job means the audit ran and retained this decision, not that biological validation passed. See the [metadata summary](biov8_evidence_2026_11_27/encode_metadata_summary.json) and [request ledger](biov8_evidence_2026_11_27/encode_metadata_request_ledger.json).

Official public ENCODE pages identify the following binding leads. They do not replace original JSON, file/annotation selection, full QC, or source-independence review. No compatible fresh perturbation-response panel was established by this phase; the blocked inventory cannot establish exhaustive absence.

| TF | Candidate binding experiment | Observed metadata and remaining gate |
|---|---|---|
| RUNX1 | [ENCSR414TYY](https://www.encodeproject.org/experiments/ENCSR414TYY/); alternative [ENCSR588AKU](https://www.encodeproject.org/experiments/ENCSR588AKU/) | Released K562 TF ChIP-seq leads; warnings and original controls/files/QC need qualification. Motif derivation is unresolved. |
| CUX1 | [ENCSR000EFO](https://www.encodeproject.org/experiments/ENCSR000EFO/) | Listed replicate numbers 1 and 2 share **ENCBS039ENC**. The distinct-biosample gate is not established by this record. This is a metadata limitation, not a finding that the original experiment had no replication. |
| HMBOX1 | [ENCSR757IIU](https://www.encodeproject.org/experiments/ENCSR757IIU/) | Controls ENCSR129NAJ/ENCSR962APM are listed. Original assay/file QC remains required; ENCBS316CKF is shared with the alternative RUNX1 record. |
| ATF6 | [ENCSR125RFR](https://www.encodeproject.org/experiments/ENCSR125RFR/) | Genetically modified, tagged K562 with control ENCSR953NJZ. Tagging and protein-condition compatibility must be retained in interpretation. |
| SOX6 | [ENCSR788RSW](https://www.encodeproject.org/experiments/ENCSR788RSW/) | Two distinct biosamples are listed; original control/file/QC and exact motif lineage remain unresolved. |

The frozen TCF4 identity is not TCF7L2. [ENCSR888XZK](https://www.encodeproject.org/experiments/ENCSR888XZK/) has a TCF4 description alias but canonical target TCF7L2; it cannot validate the frozen TCF4 gene. No such alias substitution is allowed. The [lead registry](biov8_evidence_2026_11_27/external_binding_metadata_leads.csv) is kept separate from machine-qualified metadata.

## 5. Motif derivation versus future reference assays

The exact frozen JASPAR version pages add a useful provenance distinction:

| TF / frozen matrix | Declared derivation | What this supports |
|---|---|---|
| CUX1 / [MA0754.3](https://jaspar.elixir.no/matrix/MA0754.3/) | Human CAP-SELEX; source PMID 30250250 | Declared in-vitro derivation differs from the candidate cellular ChIP assay; complete original lineage is still not retained. |
| HMBOX1 / [MA0895.2](https://jaspar.elixir.no/matrix/MA0895.2/) | Human HT-SELEX; source PMID 23332764 | Metadata support distinct assay derivation; exact sample-level reconciliation remains incomplete. |
| ATF6 / [MA1466.2](https://jaspar.elixir.no/matrix/MA1466.2/) | Human HT-SELEX; source PMID 28473536 | Metadata support distinct assay derivation; condition and full provenance gates remain. |
| RUNX1 / [MA0002.3](https://jaspar.elixir.no/matrix/MA0002.3/) | Mouse ChIP-seq; source field blank | Exact derivation unresolved; cross-species motif transfer must be disclosed. |
| SOX6 / [MA0515.1](https://jaspar.elixir.no/matrix/MA0515.1/) | Mouse ChIP-seq; source PAZAR | Exact original-experiment mapping unresolved; cross-species/context transfer must be disclosed. |

JASPAR's validation PMID field is not automatically the derivation source. No blanket source-disjointness claim is made. See [the derivation audit](biov8_evidence_2026_11_27/jaspar_derivation_metadata_audit.csv). No new binding outcome was retrieved for this audit.

## 6. Candidate coverage and feasibility gates

All 16 TFs have at least 20 available original-scale residual targets, ranging from **27 to 63**. This does not establish that any external assay measures those targets adequately. The locked future response pilot requires at least **six metadata-ready TF perturbations** and at least **20 common measured residual targets per TF**. Six is a descriptive feasibility rule, not a statistical power calculation.

Some declared sensitivities cannot meet the target gate even with complete measurement: HES1 and HMBOX1 each have only 14 residual target-only genes on the original scale; CP10k leaves PAX5 with six residual targets and TRPS1 with 13. These endpoints must be marked ineligible rather than lowering the gate or inventing scaffold scores. The [coverage table](biov8_evidence_2026_11_27/tf_candidate_coverage.csv) retains both transforms and their intersection.

## 7. Reproduction and limits

[GitHub run #46011226484](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/46011226484) completed at the implementation commit. Its eight deterministic output tables are **byte-identical** to the local run, including all frozen legacy scores and all diagnostics. The complete [workflow artifact](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/46011226484/artifacts/15204377155) includes both the numerical outputs and failed source-search audit; [provenance](biov8_evidence_2026_11_27/github_reproduction_provenance.json) records its checksum. GitHub artifacts have a 75-day retention period; this is not a permanent archival release.

An [independent numerical audit](biov8_evidence_2026_11_27/independent_numeric_audit.json), without importing the new implementation, checked 49,280 legacy rows, 3,584 graph decisions, all 2,688 probe decisions and 1,707 finite source-weight thresholds. There were **zero discrepancies**, including maximum absolute difference 0.0 across 48,507 numeric diagnostic fields. This verifies deterministic preservation and recomputation, not the biological validity of the original scores.

New expression matrices read: **0**. New models fitted: **0**. Internal/test expression cells accessed: **0**. New perturbation or binding outcome payloads retrieved: **0**. The original 394/131/133 split is unchanged and all 133 test cells remain reserved.

## 8. Decision and exact next gate

BIO-V8's numerical/source-action work is complete. Original-source metadata qualification is **blocked by access**, and independent biological validation remains pending. The project has no new evidence of GCEV, accessibility-interaction or RL superiority.

The next executable gate is to recover legitimate original ENCODE metadata access or a provenance-preserved original metadata export, rerun the unchanged inventory, and qualify coverage/QC and source separation. Only then freeze the eligible panel and complete endpoint rules before any outcome download. If coverage fails under a completed inventory, record that result and design any cohort expansion prospectively. A binding pilot requires its own endpoint lock and cannot silently replace the perturbation-response question. The reserved expression-test cells remain sealed throughout this decision.

See the [printed ten-track status](TEN_TRACK_RESEARCH_STATUS_2026_11_27.md) for the complete current project position.
