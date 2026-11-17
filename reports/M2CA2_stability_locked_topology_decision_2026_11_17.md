# M2C-A2 stability and locked-topology incremental-information audit — 17 November 2026

**Decision: NO-GO to M2C-B on the current held-out lagged-gain target.** The method is implemented, tested, and produced reproducible datasets, but across predeclared development seed/bin settings the gain and its direction margin remain less stable than the original signed pseudotime feature. On three locked additional topologies, a fixed 75/25 combination provides only slight AP improvements on CY/BFC, loses more AP on TF, and does not correct a single existing reference-orientation error. This is a **mixed and substantively negative incremental-information result**, not proof that conditional temporal reliability is impossible.

## Official reproducibility record

- [Successful GitHub workflow run #44108340318](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/44108340318), artifact `m2ca2-stability-locked-confirmation`, ID **14339227999**.
- Code: [merged PR #16](https://github.com/TshepoGobonamang/bio-rl-grn/pull/16), source commit `9d8c01050870928f84041e33c955dfc8d6e1b84f`.
- Downloaded BEELINE v4 archive, verified official checksum with the repository extraction utility.
- **Development:** dyn-LL, dyn-LI, dyn-BF — nine predefined seed×bin settings each, plus a shuffled-time negative control.
- **Locked additional structures:** dyn-CY, dyn-TF, dyn-BFC — seed 42, 12 bins only, plus a shuffled control. These were not previously evaluated in our project, but share the BoolODE simulation platform; **one realization per topology** is not broad external biological validation.
- All score generation used **only outer fitting cells**; no benchmark reference edge file was opened until after all confirmation candidate rankings were saved.
- The exact 15 stored method×topology **noninterpolated average precision and AUROC** pairs were independently recomputed from ranking CSVs and evaluation-only gold-audit tables; **15/15 match** at numerical precision.

## 1. Development: does using more folds/seeds/bins improve stability?

Mean of inner-fold pairwise Spearman rank correlations over **nine predeclared configurations per development topology** (seeds 42/73/137 × bin caps 8/12/16; 3 inner folds in each configuration):

| Development network | Original signed `z_signed` | New held-out gain | New direction margin |
|---|---:|---:|---:|
| dyn-LL | **0.380** | 0.212 | 0.126 |
| dyn-LI | **0.563** | 0.188 | 0.122 |
| dyn-BF | **0.633** | 0.204 | 0.149 |

Thus, across all three development structures the new target remains less stable than the original fitting-only signed pseudotime feature. **Nine configurations do not mean nine independent biological replicates.** The split seeds reuse the same 500 cells and reference topology.

Bin-sensitivity highlights:
- dyn-LL, mean gain stability at bins 8/12/16 = **0.203/0.258/0.175**; signed z = **0.424/0.369/0.346**.
- dyn-LI, gain stability = **0.109/0.236/0.219**; signed z = **0.569/0.533/0.586**.
- dyn-BF, gain stability = **0.306/0.261/0.046**; signed z = **0.716/0.583/0.600**. The BF gain target degrades strongly with 16 bins in these settings.

Shuffling time within each branch while preserving NA masks makes z/gain generally less stable. This establishes specificity to time ordering, **not causal regulation**.

## 2. Locked evaluation: does the new score add information beyond z?

**Reference AP (higher is better), evaluated after all rankings had been frozen:**

| Locked topology | Candidates | Eligible reference edges | Chance prevalence | z_signed AP | Held-out gain AP | Direction margin AP | Fixed blend AP | Shuffled margin AP |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CY | 30 | 6 | 0.2000 | **0.723647** | 0.640314 | 0.676389 | **0.732738** | 0.231926 |
| TF | 56 | 17 | 0.3036 | **0.570020** | 0.396895 | 0.308234 | **0.514190** | 0.276002 |
| BFC | 90 | 15 | 0.1667 | **0.551316** | 0.514320 | 0.439104 | **0.562551** | 0.164451 |

Pre-registered `fixed_blend = 0.75*clip(z_signed/q95abs(z_signed)) + 0.25*clip(margin/q95abs(margin))`. All scales computed **without reference labels**.

Within-topology blend minus z average precision:
- CY **+0.009091**;
- TF **−0.055830**;
- BFC **+0.011235**.

The unweighted mean of AP across these topology contexts is lower for the blend (approximately **0.6032**) than z alone (**0.6150**), but the networks differ in prevalence/candidate counts, so such a mean is **descriptive only**.

**Reference orientation (unique single-direction reference edges):**

| Topology | z_signed | Gain | Blend | Blend corrects z mistakes | Blend introduces new mistakes |
|---|---|---|---|---:|---:|
| CY | **6/6 (100%)** | 6/6 | 6/6 | **0** | **0** |
| TF | **6/9 (66.7%)** | 4/9 | 5/9 (55.6%) | **0** | **1** |
| BFC | **12/13 (92.3%)** | 9/13 | 10/13 (76.9%) | **0** | **2** |

The blend does not improve orientation on **any** of these three structures. Adding a direction-specific target does not, on these data, yield additional correct direction labels beyond the original signed score.

Reference provenance: CY raw 6 / eligible 6; TF raw 20 / eligible 17 (three self-loops excluded); BFC raw 18 / eligible 15 (three self-loops excluded). Self-regulation was excluded consistently in all model candidate universes.

## 3. Are the confirmation results just random time ordering?

Mean three-pair fold Spearman for 12 bins/seed 42 (proper vs shuffled branch-order):

| Topology | Original z proper/shuffled | Gain proper/shuffled | Margin proper/shuffled |
|---|---:|---:|---:|
| CY | **0.926 / 0.136** | **0.700 / −0.068** | **0.820 / −0.023** |
| TF | **0.363 / −0.160** | **0.222 / 0.126** | **0.201 / 0.092** |
| BFC | **0.759 / 0.082** | **0.139 / 0.102** | **0.065 / 0.203** |

The true pseudotime has a clear signal on CY and in the original z score on BFC; the target margin remains unstable on BFC and even has **higher** fold correlation under shuffled time there. The shuffled margin AP is approximately chance prevalence on CY/TF/BFC, consistent with limited specific signal in the non-shuffled scores. **No causal claim follows.**

## 4. Important bounds

- The target is a lagged *predictive association* between **different cells binned by pseudotime**, not an observed time series of the same cell or perturbation-based causal TF edge.
- Original signed score and new gain reuse the same inner fitting cells; validation-score associations are correlated. No p-values are computed from candidate edges, folds or split seeds as if independent.
- No graph optimization, learned reliability gate or RL policy is trained in M2C-A2; this is a target-quality/incrementality experiment, not an expression reconstruction comparison.
- Precision@20 is a weak comparison when only 30, 56 or 90 directed candidates exist, and differences in reference prevalence preclude simple cross-network score averaging.
- No parameter was tuned on CY, TF, or BFC gold labels; after inspection these three structures now also become **development-consumed** for future confirmatory designs.

## 5. Next scientific decision

**Do not train M2C-B on the current unstable gain or pairwise-margin labels.** Preserve `z_signed` as the strongest validated temporal-direction feature in the development experiments.

For the next *distinct* hypothesis, prioritize:
1. **Source-blind graph-conditional held-out usefulness** controlling for other prospective parents, not a binned pairwise lagged gain alone.
2. **Independent evidence modality** (e.g. TF motif, perturbation, ATAC-seq), with explicit provenance, source-availability/conflict masks and no overlap with evaluation gold labels.
3. New independent network realizations/seeds and a genuinely uninspected reference context, rather than reusing the now-inspected CY/TF/BFC.
4. Prospective tuning-free comparison with an expression-only baseline, fixed z, a nonlearned stability heuristic, and adaptive evidence weighting; enforce untouched-expression-loss and direction/AP co-evaluation before RL.

This is a **no-go on this particular label design**, not a no-go on the overall Bio-RL-GRN framework.
