# Evidence-to-claim ledger — 20 November 2026

**Purpose:** single source of truth for what the final manuscript may say. Evidence is classified into verified source observations, preliminary benchmark findings, unexecuted methods, and currently unsupported claims. This ledger is NOT a new biological experiment.

| Claim/topic | Level | Acceptable current wording | Evidence | Required before upgrading |
|---|---|---|---|---|
| Expression-only greedy graph selection | Exploratory observed result | On one dyn-LL-500-1 run, greedy AP 0.435288 and Pearson AP 0.342964; not a general superiority conclusion | `reports/E1_dyn_LL_500_run_43894848848_review.md` | Untouched test + independent network realizations |
| Fixed vs adaptive reliability | Negative/mixed observed result | Adaptive gate did not consistently exceed fixed temporal weight | `reports/M2B_adaptive_reliability_decision_2026_11_17.md` | Locked independent evaluation of new hypothesis |
| Pairwise directional lag/gain | Negative observed result | Gain/margin less stable and not incrementally superior across inspected networks | `reports/M2CA2_stability_locked_topology_decision_2026_11_17.md` | A distinct prespecified hypothesis/held-out networks |
| Parent-conditional cross-fitting | Technical/target-level observed | Targets can be constructed in source-blind folds; predictive utility ≠ directed causal regulation | `docs/M2CD_PARENT_CONDITIONAL_EVIDENCE.md`, `reports/M2D2_real_hESC_CollecTRI_provenance_2026_11_17.md` | Matched baseline comparison and independent test outcomes |
| Real human prior coverage | Confirmed descriptive integration | 17,735 hESC measured genes, 36,064 matched directed CollecTRI pairs | `reports/M2D2_real_hESC_CollecTRI_provenance_2026_11_17.md` | Unbiased candidate universe plus independently sourced evaluation |
| Independent prior gold | Unsupported | CollecTRI vs BEELINE hESC ChIP disjointness not established; 2,848 shared observed edges | `reports/M2D2_real_hESC_CollecTRI_provenance_2026_11_17.md` | Experiment/accession-level source exclusions and independent assay reference |
| GSE216595 H5AD schema | Verified original-source metadata | 60,997 cells, 5,047 genes, 87 observed TF labels, eight batch labels; cell indices unique | User-provided D11/D15 artifacts, original H5AD SHA256 below | None for *schema only*; biological design claims remain unsupported |
| GSE216595 TFmap provenance | Partial actual source observation | GSM6681047 compressed CSV SHA256 verified; earlier parser treated barcode-like first line as header | User-provided D15 schema artifact | Headerless first-row audit, actual cell barcode reconciliation |
| D12 `temp` meaning | Unresolved with anomaly | `temp` is not defensibly a six-level experimental temperature; 59,449 literal `temp` rows | `reports/M2D13_REAL_D12_EVIDENCE_REVIEW_2026_11_19.md` | Original processing source assigning `temp` |
| TFmap ↔ H5AD concordance | Unsupported, implementation only | Matching code and CI workflow exist; no source-matched count in supplied artifacts | `scripts/audit_m2d15b_exact_concordance.py` on D15 research branch; [workflow](https://github.com/TshepoGobonamang/bio-rl-grn/actions/workflows/m2d15-reconcile-real-sources.yml) | Real `exact_headerless_concordance.json` and sample-identity assessment |
| Independent culture replicates/negative controls | Unsupported | Eight sequencing batches do not prove independent biological cultures; no matched GFP control verified | `reports/M2D14_EXPERIMENTAL_PROVENANCE_CONTROL_EVIDENCE_2026_11_19.md` (draft PR chain) | Original culture/sample sheet and treatment-vs-control mapping |
| Direct TF regulatory effects in GSE216595 | Unsupported | 7-day overexpression may induce indirect cascades | Original GEO study design; D14 report | Assay-specific direct binding/perturbation evidence and source-disjoint gold |
| RL graph policy win | Not executed/unsupported | Do not assert any RL advantage | Repository roadmap, M2B/M2C gate reports | Actual trained RL model, budget-matched greedy/random controls, held-out test |
| Completed publishable benchmark | Not yet | Manuscript is not empirical-final | Completion charter | Frozen B evaluation and artifact-backed report |

## Verified user-supplied source artifacts, not raw expression data

- D11 ZIP SHA256: `a7b2b14409759ae08ce19bf8060e6ea9249ee2b6441ffee72628a6c1a8316418`
- D12 ZIP SHA256: `35b14f96a04f5b9ceae882b76f21b2e0832bca3c4cef1a5f8e5e28992712f1fa`
- D15 ZIP SHA256: `17ee2c44846dafe71f08329e552818384fd59babc448fa19d1783f2954ec6f6b`
- Original GSE216595 processed H5AD source SHA256 (from all reports): `a2e06453911473a0b816cf943b6730a718670a4e082a497b578e99ce637ecccd`
- Original GSM6681047 TFmap gz source SHA256 (from D15 report): `4c5a391113cc248cab12f5d098fd8d17eb96c50e06db589b5908b0b39088b401`.
- The D15 ZIP contained only `observed_h5ad_schema.json` and `observed_tfmap_schema.json`, **NOT** `exact_headerless_concordance.json`.

## Rules for publication text

- Use AP for noninterpolated sklearn average precision unless computing another precisely defined PR-AUC statistic.
- Use direction or orientation accuracy only on identifiable, reference-supported directional pairs.
- Empty TFmap/CollecTRI evidence does not imply a true biological negative.
- Learned weights are evidence-usefulness scores, not calibrated edge truth probabilities.
- Pseudotime association, predictive gains and overexpression response do not establish direct causality.
- Report negative results and sensitivity/selection trade-offs; do not reuse previously examined reference networks as pristine external validation.
- No figure, abstract performance sentence or claimed improvement without a source artifact and code commit in this ledger.
