# BIO-V8 — Evidence qualification and external-panel feasibility

**Design lock:** 27 November 2026. **Parent:** BIO-V7 commit `77d0552c2a885540e5e665d541dbd4d7f7eafd91`.
**Stage:** prospective follow-up qualification after known BIO-V7 developmental results; not an independently registered confirmatory experiment.

The executable choices are in [the configuration](../configs/biov8_evidence_qualification_2026_11_27.json). This document fixes the scope before running the full qualification workflow. BIO-V7 results and preliminary public assay listings are already known. Preliminary inspection of source ranges motivated the representation controls. These facts are disclosed rather than retroactively described as preregistration.

## Why this phase is necessary

BIO-V7 did not establish incremental GCEV superiority. Target-level accessibility was constant across each target's candidate TFs and therefore could not change an additive within-target action ranking. The current 16 TFs also do not overlap the six previously examined BIO-V2 perturbations. The immediate scientific task is to establish which evidence representations can affect the constrained graph decision and whether compatible fresh validation assays exist.

## Fixed boundaries

- Retain the 16 TFs, 64 genes, 1,008 candidates, 128 scaffold edges and 880 residual actions from BIO-V7. Preserve all seven original methods and all eight transform/prior conditions.
- Read only the existing fitting-stage candidate contexts, scores, graph selections and their integrity manifests, plus the existing prior tables. This workflow does not open an expression matrix, fit new models, evaluate internal-score cells or access the 133 reserved test cells.
- Freeze exact legacy prediction scores and eligibility before retrieving new experimental outcome files. An omitted scaffold insertion score remains unavailable; it is not filled with a favorable rank.
- Fetch experiment, target, replicate, control and file **metadata only** from ENCODE. Do not follow data download URLs or retrieve expression quantifications, differential-expression effects, peak coordinates or reference network labels.
- A future binding assay is reserved for validation in this phase. It is not simultaneously added as a model input. The motif-by-accessibility product is a derived interaction of two inputs, not a third independent evidence source.
- Preserve all unfavorable/null BIO-V7 and synthetic results. There is no search over cohorts, transformations, coefficients or permutations for a favorable result.

## A. Source-action diagnostic

Use the existing final fitting-stage expression rank E, motif score M and target accessibility A. The coefficient is fixed at 0.3. Run both original expression transformations and all three accessibility processing variants, reporting every case.

| Representation | Candidate score | Interpretation |
|---|---|---|
| Expression | E | Common action baseline |
| Motif | E + 0.3 M | TF-specific sequence potential |
| Global attenuation | E + 0.3 mean(A) M | Controls for ordinary motif shrinkage |
| Matched interaction | E + 0.3 A M | Target context modulates motif strength |
| Permuted interaction | E + 0.3 A_permuted M | Breaks target matching while preserving the accessibility distribution |

The accessibility mean is over all 64 targets, including observed zeros. Generate one target permutation from sorted target identifiers and NumPy seed **20261128**, and retain its exact mapping. Do not rank-normalize the product within each target: positive target scaling would then cancel. Add the legacy additive M+A and pure-A conditions as algebraic checks, removing common target offsets before choosing an action to prevent numerical tie artifacts.

For each target/case, report observed source values, range, distinct values, source constancy, expression and assisted winners, decision margins, and whether the assisted winner differs from expression/motif/global attenuation. Report the threshold at which a source could alter an action when implemented. These are source-representation diagnostics. A changed action is not evidence of an improved action; no new predictive or biological accuracy is measured here.

For A>0, the product A*M retains motif-only TF ordering within a target. Its potential effect is to rescale the motif contribution relative to expression. At A=0, it suppresses that contribution. Positive predicted source weights likewise do not establish observed source utility.

## B. Prediction export before external outcomes

Export every existing candidate-score row from all 56 BIO-V7 conditions: **49,280 rows = 56 × 880 residual actions**. Validate each input against the original freeze manifest and preserve its source hash. The remote reproduction freeze manifest has SHA256 `50b23a11a04927a1ff367a1ca42ef9dd5757a6a9639950f128d6d78da8d27c93`; its biological decisions reproduce the earlier local freeze.

For future external response analysis, the primary condition remains original log1p expression with the original mixed-caller prior. The existing `score` field is the prospectively fixed across-target projection. It includes target-level offsets that cancel within-target decisions but can affect ranking targets within a TF. The `action_score` remains the within-target decision score. These are different ranking tasks. Neither score is a binding probability, causal effect, nor a validated complete regulatory-edge score. Publish both and their exact definitions, without choosing the projection after seeing external outcomes.

The primary descriptive endpoint, if a future panel passes qualification, is per-TF Spearman correlation of the frozen score with absolute perturbation log2 fold change on common measured residual candidates. Report every TF and the unweighted macro mean; compare GCEV with no-graph-state adaptation as the primary descriptive contrast. Expression-only and fixed sources are secondary controls. Repeat on the 48 additional target-only genes as a declared sensitivity. Scaffold edges remain excluded from this insertion-score endpoint.

## C. External metadata inventory and qualification

Inventory released K562 CRISPRi RNA-seq, shRNA/knockdown RNA-seq and TF ChIP-seq metadata. Match the frozen canonical gene identity, including stable Ensembl identifiers when available. Description aliases are insufficient: the TCF4/TCF7L2 ambiguity is an explicit rejection test.

Capture original metadata response bytes, URL, retrieval time and SHA256. Record the experiment and canonical target, exact biosamples and biological/technical replicate numbers, experimental controls, released GRCh38 output file accessions, annotation, processing version, and audit flags. A search failure or truncation is incomplete coverage, never proof that no assay exists. Include all eligible assays rather than choosing one according to its outcome.

For response metadata readiness, require matching treatment/control assay context and gene annotation, at least two distinct biological biosample identifiers and replicate numbers in each group, and real released per-biological-replicate gene quantification files. GENCODE V29 is the fixed annotation for this pilot; another version requires a separately frozen mapping revision before outcomes. Unresolved ERROR/NOT_COMPLIANT audits block readiness. A separate future outcome gate must still verify knockdown efficacy, numeric integrity and replicate concordance.

Binding metadata are alternative leads, not a substitute for the perturbation-response endpoint. Their controls, biological replicates, GRCh38 reproducible outputs, antibody/target specificity, culture or tagging differences, and assay QC require qualification. Promoter peak overlap would support assay-specific binding, not functional regulation. Missing peak overlap cannot automatically be treated as a true negative edge.

The previously inspected BIO-V2 six contrasts and their shared control are excluded from fresh validation. Check original experimental sources behind the exact JASPAR matrix versions before claiming source-disjoint validation. Fresh accessions do not prove independence when their experiments contributed to motif derivation. Unresolved provenance stays unresolved.

## D. Feasibility and stopping rules

At least **six metadata-ready TF perturbations**, and later at least **20 common measured residual targets per TF**, are required to proceed to a descriptive response pilot. Six is a feasibility rule, not a power calculation or proof of biological independence. The target gate cannot be passed from metadata alone. Shared controls, shared cultures and assay dependencies remain visible; no population-superiority inference is authorized.

If coverage or quality fails, publish the concrete no-go, all exclusions and the unfulfilled gate. Do not replace TFs, widen the cohort, substitute an unrelated endpoint or open test cells to obtain a favorable result. Any later cohort redesign or interaction-model training must have its own explicit developmental design and new lock.

## Deliverables and claim boundary

Produce frozen prediction/eligibility tables; source-action diagnostic tables and permutation map; original metadata registry and search ledger; TF-level coverage and readiness report; and the updated existing ten-track status table.

BIO-V8 can establish representation identifiability and the feasibility of a compatible external panel. It cannot establish that accessibility improves GCEV, that inferred edges are true regulatory interactions, that trained RL helps, or that the project is publication-ready.
