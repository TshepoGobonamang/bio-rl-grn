# M2C-D2 — Real human evidence integration and provenance audit (17 November 2026)

## Research decision

**PASS: the external-to-expression biological prior source is correctly downloaded, checksum-pinned, gene-mapped, and ingested in a real human expression parent-conditional feature pipeline.**

**NOT YET VALIDATED: independent regulatory accuracy or source-disjoint ChIP-seq benchmarking.** A post-frozen-prior overlap audit reveals direct shared edges, but we cannot certify underlying experimental provenance disjointness. The small expression pilot is explicitly curated-prior-enriched and cannot support an unbiased GRN recovery claim. **Do not train an RL policy or assert adaptive biological-prior superiority from this result.**

### Reproducibility

- [Successful official GitHub workflow #44114658037](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/44114658037), artifact `m2d2-public-source-schema`, ID **14342636691**.
- Code: [PR #20](https://github.com/TshepoGobonamang/bio-rl-grn/pull/20) (M2C-D2 real data).
- BEELINE v4, Zenodo record [19009603](https://zenodo.org/records/19009603): `BEELINE-data.zip` published MD5 `9ea93065004f2f12022adeeaf3728edd`, source SHA256 `3e67582fb8a4ba5c2189f6a0bfb14abb1d61f21921459bdeb4a2a9f2b652d816`.
- BEELINE-Networks reference archive, same Zenodo record: MD5 `a7dca84d302bb127bad962ae52815662`. It is opened **only** after the original curated prior table, pilot features and held-out pseudo-targets have been saved.
- CollecTRI version **2.0** [Zenodo 8192729](https://zenodo.org/records/8192729) regulons MD5 `cee4a3943c059e6dd8796ced6dec44f6` and source TSV MD5 `c8622a46a49e21e955182f000cf1c067`.
- CollecTRI original columns: `source,target,weight,resources,references,sign_decision`; output preserves regulation sign separately from the TF→target arrow and records supporting resource labels/cited-reference counts. `weight=-1` is **inhibition**, not reverse graph direction.

## 1. Exact gene/edge coverage (all measured human genes)

| Measured/curated quantity | Official result |
|---|---:|
| BEELINE hESC measured human gene symbols | **17,735** |
| BEELINE hESC expression cells | **758** |
| CollecTRI original source rows | **43,175** |
| CollecTRI unique directed TF–target pairs | **42,713** |
| Directed CollecTRI pairs with both source and target directly measured | **36,064** |
| Distinct matched curated regulators | **1,079** |
| Distinct matched curated targets | **5,868** |

Exact, case-sensitive symbol matching; no guessed aliases, ortholog assignments, or synthetic g1-gN substitutions. This is genuine external *regulatory-association provenance* matching, but **presence** of a curated prior record is not calibrated TF binding probability, and lack of a record does not imply non-regulation.

## 2. Prior-enriched real hESC integration proof-of-concept

To prove ingestion without inventing prior links, select 8 high-coverage curated regulators + 24 targets **using the CollecTRI prior itself**. This strongly biases source coverage; it is a software integration pilot, not a fair selection protocol for GRN evaluation.

| Pilot quantity | Actual |
|---|---:|
| Genes | **32** |
| Candidate directed non-self edges | **992** |
| Prior-observed directed edges | **233** (23.49% by deliberate selection) |
| Outer fitting cells | **454** |
| Outer graph-scoring cells withheld | **152** |
| Outer expression-test cells withheld | **152** |
| Inner cross-fitting folds | **3** |
| Separate fit-only feature rows and V-fold target rows | **2,976 each** |
| Records with present curated prior | **699** (233 edges × 3 folds) |
| Records with unavailable curated prior | **2,277**, with actual prior score **NaN**, not zero |
| Reference/ChIP labels accessed during feature or target generation | **0** |

All 2,976 feature/target keys align without duplicates; the source-blind parent scaffold excludes candidate source and target in every record. The original outer 152-cell test subset was not used for target training or model selection.

### Predictive-target quality diagnostic (NOT gold regulatory accuracy)

- Mean held-out parent-conditional expression gain: **−0.007077**; negative gain in **48.35%** of edge-fold records.
- Mean absolute gain **0.01471**; four targets clipped to ±1.
- Across three inner folds, edge-score Spearman correlations **0.2576, 0.3037, 0.2820**, mean **0.2811** — noticeably lower reproducibility than many BoolODE synthetic benchmarks. This is descriptive only, because directed pairs and folds are not independent.
- Prior-present records had mean gain **−0.00745**, prior-unobserved **−0.00696**. **This difference is not a prior-efficacy result**: source-observed interactions are nonrandom and the entire cohort was selected using the prior. Unobserved edges are not true negatives.
- No GRN reference AUROC, average precision, regulatory orientation or independent final expression-test performance was computed.

## 3. Post hoc reference-overlap audit (separate phase)

Source: BEELINE human `Networks/human/hESC-ChIP-seq-network.csv`, extracted only **after** source selection and the entire pilot prediction/feature generation was frozen.

| Audit quantity | Observed |
|---|---:|
| Unique directed BEELINE hESC ChIP reference edges, excluding self-loops | **436,513** |
| Reference edges whose two endpoints occur in the matched CollecTRI gene universe | **142,442** |
| CollecTRI prior / ChIP shared directed edges | **2,848** |
| Shared edges / all 36,064 matched CollecTRI pairs | **7.90%** |
| Shared edges / 142,442 reference edges in eligible gene universe | **2.00%** |

**Important:** edge-set overlap is *descriptive agreement only*. It neither proves both resources are derived from the same experiment nor that their underlying evidence is fully source-disjoint. Literature-based CollecTRI can incorporate ChIP-related source data; source-level experimental accession overlap has **not** been resolved. Consequently neither source can currently serve as a **verified independent validation reference** for a CollecTRI-prior-trained model.

No ChIP edge file was used to select pilot genes, construct priors, fit parent models, tune weights, or generate targets.

## 4. Research conclusion and next scientific gate

- **Completed:** real human gene ID alignment, genuine curated external TF-edge coverage, original source/PMID provenance retention, guarded missingness, real-expression parent-conditional data generation, and reference-access separation.
- **Not completed:** certifying independence of the prior from BEELINE ChIP experiment sources, expression-only vs fixed-prior vs adaptive graph performance, and independent regulatory edge validation.
- **Next data prerequisite:** a source-disjoint, genuinely held-out perturbation or assay-based TF regulation dataset, or independent promoter-sequence JASPAR motif annotations with documented GRCh38 TSS mapping and different evaluation assay origins.
- Before a performance paper claim, construct an **expression-only-selected** gene/candidate universe without using prior/reference labels, pre-specify source-dropout and wrong-sign corruption controls, then evaluate fixed/no/adaptive evidence under matched graph edit budgets and untouched expression-test losses on previously uninspected context/cohorts.
- All numbers here are exploratory ingestion/provenance or correlation diagnostics, **not** evidence of causal TF regulation or M2B/M2C adaptive gate superiority.

Official source documentation: [BEELINE v1.1 experimental data](https://murali-group.github.io/Beeline/reproducing-results.html), [CollecTRI publication](https://doi.org/10.1093/nar/gkad841).
