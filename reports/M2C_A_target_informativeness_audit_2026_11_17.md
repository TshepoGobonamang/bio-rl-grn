# M2C-A target informativeness audit — 17 November 2026

## Executive decision

**NO-GO: do not train M2C-B on the new M2C-A target unchanged.** The signed gain target is technically well formed, nondegenerate and direction-sensitive, but too unstable across independent inner folds and inferior to the existing fitting-only signed pseudotime margin as an *evaluation-only* reference orientation/ranking signal. It shows partial association with gold edges on LL and BF; this is a **diagnostic**, not evidence of causality or of reliable learned biological priors.

**Source provenance:**
- [Official successful M2C-A GitHub Actions run #44105427524](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/44105427524); `m2c-directional-crossfit-targets` artifact ID **14339056837**.
- All M2C-A targets were constructed without reading `GroundTruthNetwork.csv`. The independent *post-hoc* informativeness audit joins labels from the previous [M1.2 run](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/43901680886), which itself withheld the gold reference until evaluation.
- Original datasets: BEELINE v4 dyn-LL-500-1 (18 gold / 306 admissible directed candidates), dyn-LI-500-1 (7 / 42), dyn-BF-500-1 (10 / 42, 8 uniquely oriented reference pairs).
- Three inner folds of 300 outer-fit cells each; source-only inner fitting features and validation-target tables are separate. Each topology has been inspected during previous milestones; none is an untouched prospective network test.

## 1. Schema/missingness/quality check

| Network | Edge-fold target rows | Available | Positive gain % | Negative gain % | Mean signed gain | Clipped | Mean absolute pair margin |
|---|---:|---:|---:|---:|---:|---:|---:|
| dyn-LL | 918 | 918 | 47.3 | 52.7 | −0.0320 | 0.4% | 0.1755 |
| dyn-LI | 126 | 126 | 58.7 | 41.3 | −0.0182 | 1.6% | 0.2293 |
| dyn-BF | 126 | 126 | 61.9 | 38.1 | +0.0231 | 0.0% | 0.1407 |

The automated official audit verified no missing evidence on these clean simulations, exact fitting-feature/validation-target key alignment, direction-margin antisymmetry for all 1,170 directed fold-edge rows, 0 gold-reference training inputs and disjoint folds. Antisymmetry and a 100% nonzero-margin fraction are **mathematical/implementation properties**, *not* independent evidence of biological directionality.

## 2. Main deficiency — out-of-fold reproducibility of the target

Mean Spearman correlation of edge scores across the **three possible pairs** of inner-fold results:

| Network | New directional held-out gain | New forward-reverse margin | Existing fitting-only signed temporal feature |
|---|---:|---:|---:|
| dyn-LL | **0.171** | 0.153 | **0.430** |
| dyn-LI | **0.243** | 0.185 | **0.653** |
| dyn-BF | **0.025** | −0.063 | **0.516** |

Sign of direction margin was the same across **all three folds** for only **37.3% (LL), 42.9% (LI) and 23.8% (BF)** of directed candidates. The fitting-only signed temporal feature had three-fold sign consistency of **49.0%, 71.4%, and 42.9%** respectively. Because each pair has reciprocal directed rows, these fractions do not represent independent observations.

A different score per direction is insufficient if the direction preferred by that score changes substantially between folds. In BF, the target's mean rank correlation (~0.025) is particularly weak.

## 3. Evaluation-only reference association

After averaging each score over its three folds and joining benchmark gold labels **only for this audit**, calculated noninterpolated average precision (AP), AUROC and orientation concordance for gold interactions with a unique orientation:

| Network | Score | AP | AUROC | Direction correct / eligible |
|---|---|---:|---:|---:|
| dyn-LL | **M2C held-out gain** | **0.198** | 0.660 | **12/18 (66.7%)** |
| | Fitting-only signed pseudotime `z_signed` | **0.329** | 0.779 | **15/18 (83.3%)** |
| dyn-LI | **M2C held-out gain** | **0.388** | **0.404** | **2/7 (28.6%)** |
| | Fitting-only signed pseudotime `z_signed` | **0.423** | **0.665** | **5/7 (71.4%)** |
| dyn-BF | **M2C held-out gain** | **0.655** | 0.747 | **6/8 (75.0%)** |
| | Fitting-only signed pseudotime `z_signed` | **0.663** | 0.806 | **7/8 (87.5%)** |

Reference prevalence/random-ranking AP baselines: LL **18/306=.0588**, LI **7/42=.1667**, BF **10/42=.2381**. Hence the new target is *not meaningless*: AP exceeds prevalence on all three networks. But the simpler pre-existing directional feature is consistently stronger. LI's new target AUROC is **below 0.5**, even though AP is above prevalence; a few high-ranked true edges can coexist with poor global ordering.

On uniquely oriented reference edges, the new gain agrees with gold where the original `z_signed` disagrees on **2 of 18 LL cases**, but also reverses **5 LL cases** that `z_signed` gets right. In LI and BF, it fixes **zero** errors made by `z_signed` and introduces 3 and 1 new orientation errors, respectively.

These results are exploratory and may be optimistic from having inspected the same synthetic topologies repeatedly. No experiment on independent, previously unexamined topology has been performed.

## 4. Interpretation

1. **Construction works**: fit-only input features and held-out loss targets are separated; negative predictive gains and no-source masks are handled; reversed edges need not share the same directional score.
2. **Learnability is questionable**: weak across-fold correlation suggests a noisy label; naive regression would risk learning fold-specific temporal fluctuations instead of transferable edge reliability.
3. **No incremental information demonstrated**: reference AP/orientation of the new target do not beat the original signed pseudotime evidence, and the rare corrections it introduces do not offset new errors. This does not prove that no future combination can help.
4. **Biological-causality boundary**: binned pseudotime of independent cells is not true trajectory history. Positive held-out gain is predictive association only.
5. **Caution on sample counts**: the 1,170 fold-edge records do not constitute 1,170 independent tests. LL has 306 directed candidates and the small networks 42 each; reference orientation comparisons rely on only 18, 7 and 8 unique directions.

## 5. Next development gate (M2C-A2, before M2C-B)

**Do not change the published M2C-A artifact or retrain a complex gate yet.** Proposed follow-up:
- Test *fold stability* under more independent network realizations and inner-fold seeds; predeclare 8/12/16 bins and lag 1 rather than selecting on inspected gold.
- Compare source-only signed temporal margins with the signed held-out gain and their agreement across folds. Quantify a **reproducible, direction-sensitive signed evidence confidence** target, preserving negative gain and missingness.
- Use inner time-respecting cross-fit and repeated simulations; analyze variance from only 16 aggregate bins. Add branch-root uncertainty on BF.
- Include shuffled pseudotime as a negative control and evaluate entire new topologies/independent seeds before choosing a target for M2C-B.
- Prefer retaining `z_signed` as the primary directional feature unless a revised score proves **incremental held-out information beyond it** (AP/orientation or downstream independent prediction under comparable inputs).
- Keep the M2B fixed/global/stability baselines and matched expression-generalization checks.
- Predefine an acceptable expression-test loss noninferiority margin and network-level uncertainty assessment.

**Go criterion for M2C-B:** revised signed surrogate is sufficiently stable and demonstrates complementary, prospective usefulness beyond existing fitting-only temporal evidence on unseen network realizations. These criteria have **not** been met by M2C-A.

## Audit deliverables

Reproducible local diagnostic files: `m2c_quality_summary.csv`, `fold_stability.csv`, `per_edge_averages.csv`, and `target_gold_diagnostic.csv`. Statistical calculations used SciPy Spearman and scikit-learn average precision / ROC AUC; no gold labels used in constructing targets. See the official artifact for all raw output rows.
