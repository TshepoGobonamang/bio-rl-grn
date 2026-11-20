# CIBB recalibration — locked project-completion charter
**Decision date:** 2026-11-20  
**Status:** scientific protocol / delivery scope freeze. Not an assertion of new empirical success.  
**Authority:** prior source-backed reports and user-supplied D11, D12, D15 artifacts.

## One final research question

> When biological evidence is incomplete, unreliable or conflicting, do cross-fitted **graph-parent-conditional reliability signals** improve directed TF–target ranking, held-out expression prediction, or search efficiency relative to carefully matched expression-only, fixed-weight and nonlearned controls?

This must be answered **even if the answer is no**. Model construction is complete only when the comparison is executed on source-pinned benchmark datasets and reproducible, evidence-bounded results are written. Predictive/temporal association is NOT causal regulation.

**Working manuscript title** (claims neutral):
*Graph-Conditioned Reliability Signals for Gene Regulatory Network Inference: A Leakage-Aware Experimental Benchmark.*

The original CIBB 2026 meeting was 2–4 September 2026; do not describe this ongoing project as a pending CIBB 2026 submission. A future venue must be selected using current author guidelines. Existing accepted/extended proceedings rights, if any, require separate verification.

## Nonnegotiable decision: stop serial M2C-D16, D17... gate proliferation

From now on maintain **one integration branch/PR and one issue**, with four finished work packages and a single final experimental report. No additional synthetic or GEO audit branches unless repairing a demonstrated blocker. Existing chained draft PRs (#20–#35) are historical research records; review before merging, and do NOT merge all 16 onto main automatically.

### Work package A — Source audit closeout (GSE216595)

Existing verified observations:
- Original GSE216595 H5AD: 60,997 cells; 5,047 measured genes; 87 TF labels; 8 `batch` values; SHA256 `a2e06453911473a0b816cf943b6730a718670a4e082a497b578e99ce637ecccd`.
- D12 TF/batch/temp marginal sums and 682 joint groups reconcile to 60,997. `temp` field is anomalous; do not interpret as temperature or control.
- Original TFmap source SHA256 `4c5a391113cc248cab12f5d098fd8d17eb96c50e06db589b5908b0b39088b401` observed in D15 schema artifact. First apparent header (`CGTCTACCATGCCCGA`, `MXD3-NM_031300`) looks like a data row. Headerless interpretation **not yet executed/confirmed** from a newer artifact.
- The old schema report listed 76,233 records and 88 nonmissing TF values **under header-as-data parsing**; those are not an authoritative count after correction.
- The D15B source workflow is registered: [M2D15 real-source reconciliation](https://github.com/TshepoGobonamang/bio-rl-grn/actions/workflows/m2d15-reconcile-real-sources.yml). Its actual `exact_headerless_concordance.json` is **not yet available in the supplied artifacts**.

**Finish A** by obtaining a real D15B artifact, verifying first-row treatment and exact barcode match, including missingness/discordance and source hashes. Reconcile sample provenance before interpreting cell-only matches as valid across libraries. If no recoverable culture-level replicate IDs and assay-matched controls exist, mark GSE216595 **exploratory-only**, do not keep it on the main paper's critical path. No invented GFP controls.

**A stop rule:** After one reproducible correct-source reconciliation and one definitive methods/sample-sheet provenance review, unresolved independent control/replicate status is a documented **NO-GO**, not another open-ended sequence of source-inspection workflows.

### Work package B — Publishable comparison protocol (primary paper)

1. **Freeze algorithm/selection rules and code commit before opening new reference labels.** Prior-inspected BEELINE networks (LL/LI/BF/CY/TF/BFC) are development-only even if reused for debugging.
2. **Create a pre-declared held-out evaluation set** from genuinely not-yet-assessed independent network realizations or a different benchmark generator/system. Record generator/source provenance and seeds or accessions, topology identities, candidate edge universe and partition hashes. Replicates of the *same* topology are distinct realizations, not new topology families.
3. Fix expression fitting, internal graph-scoring and **untouched final test** splits before training. Do not use held-out gold edges to select TFs/targets/priors, tune graph search or select weights; save locked ranking outputs before the gold file is accessed by evaluation-only code.
4. Matched baselines: Pearson; Spearman; graph-conditioned ridge/greedy; fixed signed temporal evidence (`z_signed`, where available); **simple fixed reliability weight**; nonlearned evidence-availability/stability heuristic; adaptive graph-conditioned reliability. If a modal source is absent, explicitly mark missing and do not turn missingness into a negative edge.
5. Evidence experiments: report *all* 0/25/50/75% source dropout and 0/10/25/50% source corruption settings only where the modality is real and matched. Include shuffled-source and masked-source controls; do not synthesize external human evidence for abstract g1–gN networks.
6. Primary outcomes are (a) edge average precision [`sklearn.average_precision_score`, labeled AP rather than assuming trapezoidal PR-AUC]; (b) direction/orientation accuracy on assessable unique-reference pairs; (c) untouched expression-test loss. Secondary: AUROC, early precision/P@K, source coverage, graph edit calls and compute runtime. Preserve prevalence/reference coverage. Negative and null results must remain.
7. Report uncertainty at the **independent network realization / graph** level, with topology-cluster sensitivity analysis. No p-values treating correlated cells, edge rows, inner folds or random split seeds as independent biology.
8. Primary success is **not** merely more stable labels or stronger fit on internal scoring cells. Report effect vs strongest fixed/nonlearned control and test-set trade-offs. If AP improves but direction/generalization fails, report that limitation.

**Finish B** with one pinned, end-to-end reproducibility workflow producing `run_manifest.json`, `results.csv`, `per_network.csv`, `uncertainty.csv`, `ablations.csv` and ranked-edge prediction archives, all with exact SHA256/commit/source/seed. Report **no superiority** if the prespecified comparison does not support it.

### Work package C — RL as a falsifiable optional experimental extension

The project software can include a graph-edit RL policy (add, remove, stop) but the *publication conclusion* must not depend on RL wins.

Implement/evaluate RL **only after B's benchmark contract is stable**; use exactly the same candidate universe, graph utility, training data, maximum graph edits, and **utility-query and wall-clock budgets** as greedy or nonlearned search. Report learning curves, utility-vs-budget, generalization, random-policy and no-policy controls. An RL agent that merely uses more reward queries than greedy is not a valid improvement. `stop` must correspond to a learned policy decision, not exhaustion of a fixed edit budget. Include training/execution time.

**C decision:** RL may be reported as a tested extension even when it ties or loses. If implementation is incomplete or its matched controls cannot be executed, omit RL from the paper's title, abstract, and claimed contributions; preserve it as a separate software/research extension.

### Work package D — Final manuscript and release

Deliver exactly:
- a frozen methods section defining target, graph-conditioned scoring, reliability, optional RL, splits and evidence masks;
- a numbered, source-backed experimental dataset table including reference limitations;
- a unified result table and plots (AP, direction, untouched test loss, corruption/dropout, compute budget), without placeholders presented as measurements;
- negative findings and threats to validity, including reference-source overlap and potential indirect effects in perturbation screens;
- a reproducible environment, test suite, command(s), workflow artifacts and code/data availability statement;
- a venue-conforming submission draft only after the empirical results above exist.

**Completion criterion:** Independent rerun from a pinned commit yields matching summary metrics within declared numeric tolerances; all written scientific claims map to source-backed artifacts and no untested model improvement is asserted. Peer-review acceptance is not a software/research completion criterion and cannot be guaranteed.

## Already observed negative/mixed evidence that must be preserved

- E1 dyn-LL-500-1 exploration: greedy ridge AP 0.435288 vs Pearson 0.342964 on 306 candidates / 18 eligible reference edges, but **one replicate**, no untouched expression test; not RL.
- M2B adaptive reliability: did **not consistently beat the fixed temporal weight** on ranking/orientation/untouched expression generalization. Do not claim learned reliability superiority.
- M2C-A2 target: lower rank stability than `z_signed`; preregistered blend AP improvements on CY/BFC but **worse** on TF, with no repaired orientation errors; do not recycle the observed network families as unused confirmation.
- M2C-D2 hESC/CollecTRI: 36,064 measured curated prior directed edges; 2,848 direct edge overlap with BEELINE hESC ChIP reference, and no certified experiment-level source disjointness. The 32-gene prior-enriched pilot is software feasibility, not an unbiased regulatory benchmark.

## Decision tree

**Path 1 (preferred):** B demonstrates a reproducible incremental graph-conditioned or reliability benefit under correctly held-out evaluation; paper reports measured improvements with caveats. RL is evaluated if C runs and comparisons are fair.

**Path 2 (valid publication outcome):** B demonstrates no consistent advantage; finish a transparent benchmark/negative-results methods paper assessing stability, biological prior provenance and graph-edit search. Do not rewrite the hypothesis or tune on test results to find an apparent win.

**Path 3 (true blocker):** Public data / authentic reference provenance prohibits confirmatory biological claims. Complete the reproducible computational benchmark, mark independent biological validation as unverified, and seek separate well-controlled perturbation evidence for an explicitly future extension.

**No new modelling stage is created solely to bypass a negative result.**

## Immediate execution order

1. Collect real D15B `exact_headerless_concordance.json` and `headerless_tfmap_audit.json`, accept or close GSE216595 as exploratory.
2. Freeze B's benchmark configurations, datasets, baseline matrix and commit; choose a genuinely uninspected test set.
3. Execute *one* versioned, end-to-end comparison and independently reproduce the reported metrics.
4. Decide C using B's measured constraints; run RL with matched budgets **only if executed**.
5. Freeze claims and write the complete manuscript from the result artifacts.

## Sources

- Existing repository evidence: `reports/E1_dyn_LL_500_run_43894848848_review.md`, `reports/M2CA2_stability_locked_topology_decision_2026_11_17.md`, `reports/M2D2_real_hESC_CollecTRI_provenance_2026_11_17.md`, `reports/M2D13_REAL_D12_EVIDENCE_REVIEW_2026_11_19.md`.
- BEELINE reference methodology: https://doi.org/10.1038/s41592-019-0690-6
- CIBB 2026 historical call and schedule: https://cibb2026.teralab.ai/calls/
