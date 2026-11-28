---
title: A Controlled Evaluation of Graph Conditioned Evidence Weighting for Gene Regulatory Network Inference
lang: en-GB
date: 28 November 2026
---

# Abstract

Biological priors can assist gene regulatory network inference, but their usefulness depends on how evidence changes an inference decision. We evaluated graph-conditional evidence valuation (GCEV), a one-step procedure that learns source weights from nested, paired comparisons of graph edits with and without an evidence source. Target-held-out learning separates the utility observations for a gene from the model predicting its weights. Three prospectively locked simulation studies used 48, 32 and 32 independently seeded directed acyclic networks to test graph-state information, a nonlinear-interaction hypothesis and an explicit parent-feature ablation. In the 48-network primary comparison, the estimated GCEV-minus-control average-precision difference was 0.007715, with a 95% network-bootstrap interval of −0.010750 to 0.026722. Fixed weighting exceeded GCEV by 0.074505 average precision on clean evidence. Neither subsequent study met its primary success criterion. A developmental K562 comparison used 394 fitting cells, 131 evaluation cells, 16 transcription factors, 64 genes, sequence-motif priors and DNase accessibility. The primary GCEV-minus-source-feature-only prediction-error difference was −0.000223, with a conditional paired-cell interval of −0.000636 to 0.000175. Accessibility represented as an additive target-constant term could not alter within-target transcription-factor selection, despite positive learned accessibility weights for 47 of 64 targets. Multiplicative representations changed decisions, but their predictive usefulness was not evaluated. These results do not establish an advantage for the tested graph-state representation. They identify a distinction between evidence availability, influence on an action and observed predictive utility, while delimiting the conclusions supported by small simulations and a single real-data cohort.

**Keywords:** gene regulatory network inference; biological priors; adaptive evidence weighting; graph-state features; controlled evaluation; source utility; single-cell transcriptomics

# 1 Introduction

Gene regulatory network inference seeks to identify relationships through which transcription factors influence gene expression. Single-cell transcriptomic measurements provide many observations of cellular expression states, yet an expression association alone does not identify direct regulation. Different networks can explain similar observations, and the available samples capture only part of the biological variation relevant to a regulatory system. Comparative evaluations have therefore examined both reconstruction against reference networks and sensitivity to the datasets, assumptions and metrics used for assessment [@beeline].

Prior biological information offers a way to constrain this problem. Sequence motifs indicate potential recognition sites; chromatin accessibility describes genomic regions in which regulatory activity may be possible; binding assays record occupancy under particular experimental conditions; and perturbation measurements quantify responses to interventions. These sources have different meanings and resolutions. A target-level accessibility measurement, for example, does not by itself distinguish which transcription factor regulates that target. Methods that integrate such information must specify both what a source measures and how it affects inference [@scenic; @scenicplus; @celloracle; @encode].

Adaptive weighting appears attractive when prior quality differs across sources or targets. A fixed source weight may be inappropriate when coverage is incomplete or a prior transfers poorly between conditions. However, adding a learned weighting model introduces its own estimation problem. The training signal must distinguish useful evidence from evidence that is merely present, and the learner needs features that can explain variation in source usefulness. An apparent benefit may arise from generic adaptation, extra fitting effort or a change in the candidate universe rather than from information about the current graph. Work on optimising prior integration and adaptive fusion already addresses parts of this problem [@integration2024; @genes2reg].

We study a specific hypothesis: information about a source-blind parent scaffold and the proposed next graph edit may help predict when an evidence source will improve a graph-edit decision. Our procedure, graph-conditional evidence valuation (GCEV), constructs utility observations through paired predictive comparisons. For each target and source, it selects one candidate edge with the source available and one with it ablated, keeping the candidate set and number of inserted edges fixed. It then compares the selected models on disjoint validation cells. A second learner predicts these source utilities for held-out target genes and maps the predictions to bounded weights.

The evaluation asks four questions. First, does the full feature representation improve edge ranking over a matched adaptive learner without the graph-state block and over fixed source weighting? Second, does introducing two-parent nonlinear interaction explain any observed graph-state benefit? Third, do explicit parent features add information beyond source and action features? Fourth, in a real-expression setting, do the priors change the intended decisions and improve prediction? These questions concern the implemented one-step procedure. No reinforcement-learning policy or neural graph model is evaluated.

The experimental sequence permits both prospective tests and later diagnosis. A 48-network study was locked before evaluation. Its descriptive heterogeneity motivated two further, separately locked studies using disjoint network seeds. A K562 analysis subsequently applied the weighting idea to a constrained transcription-factor candidate universe. Its final graphs were frozen before evaluation on an internal-score partition. A final saved-output analysis examined the action-level behavior of the priors. That analysis followed the predictive results and is reported as a diagnosis, not as a retroactively confirmatory test.

The contribution is a controlled empirical assessment of this source-valuation formulation and the interpretation of its learned weights. The primary tests did not establish the proposed advantage. The evidence-use analysis nevertheless exposes a concrete failure of interpretation: a nonzero learned weight can coexist with an evidence term that is structurally unable to change the decision being valued. We report the unfavorable comparisons, the conditions under which the observation applies and the remaining limits on generalization.

# 2 Related work

## 2.1 Expression based network inference

Expression-based methods differ in their statistical representation and in the assumptions connecting that representation to regulatory edges. GENIE3 decomposes inference into target-specific prediction tasks and uses tree ensembles to assign regulator importance. GRNBoost2 provides a gradient-boosting implementation suited to scalable network inference. PIDC uses multivariate information to characterize dependencies among genes [@genie3; @grnboost2; @pidc]. These methods establish relevant expression-based alternatives, but they were not executed as external baselines in the present study. Our principal comparison instead holds the learner, utility labels and evidence fixed while removing selected feature blocks. It therefore estimates an incremental effect within a particular implementation, rather than a ranking against the field's leading algorithms.

BEELINE showed the importance of evaluation across synthetic networks, curated models and experimental datasets, with distinct behavior across these settings [@beeline]. This motivates our separation of synthetic edge recovery from real-expression prediction. It also identifies a limitation of the present evaluation: additional random realizations from one small-network generator do not provide the diversity of several external benchmark families.

## 2.2 Integration of motifs and chromatin evidence

SCENIC combines expression-derived modules with motif-based analysis to infer regulatory programs, while SCENIC+ integrates expression, accessibility and motif information to infer enhancer-driven regulatory networks. CellOracle combines prior network information with single-cell expression to model regulatory configurations and in silico perturbation responses [@scenic; @scenicplus; @celloracle]. PANDA provides an earlier example of network inference through integration of motif, protein-interaction and expression information [@panda]. Inferelator 3.0 extends prior-informed network inference to large single-cell datasets, while LINGER incorporates atlas-scale external information in single-cell multiome inference [@inferelator; @linger]. These methods use richer structures and objectives than our target-level additive accessibility term. Our invariance result should not be interpreted as a limitation of accessibility in those methods.

More directly relevant, DIOgene selects gene-specific prior-integration strength using predictive error and a simulated null comparison in Arabidopsis. Genes2Reg learns adaptive fusion of heterogeneous knowledge into a gene-coupling prior [@integration2024; @genes2reg]. GCEV consequently cannot be presented as the first adaptive use of regulatory evidence. The question examined here is narrower: whether a matched, cross-fitted estimate of the predictive value of a source-informed edit benefits from explicit parent and action features. GLM-Prior separately evaluates sequence-derived prior construction and its use in downstream network inference, emphasizing the influence of prior quality [@glmprior]. Establishing a general advantage over these integration methods would require additional head-to-head experiments.

## 2.3 Evidence contribution and evaluation design

Recent evaluation work has examined whether multimodal GRN methods actually use the information supplied by each modality. BEAR-GRN applies controlled modality perturbations and several reference definitions to assess modality contribution and other aspects of network inference [@beargrn]. The present study shares the need to distinguish availability from contribution, but operates at a different resolution: it measures whether a source changes a particular one-edge choice and whether that choice improves held-out expression prediction. A source may affect a global network or a different downstream endpoint even when it has no effect on this local decision.

Evaluation also depends on the target of inference. DREAM5 assessed network recovery across diverse methods, STREAMLINE examines topological properties, and CausalBench evaluates network inference using large single-cell perturbation datasets [@dream5; @streamline; @causalbench]. Analyses of single-cell transcriptome inference further show why measured RNA need not capture every determinant of regulatory activity [@xue2023]. These studies motivate explicit endpoint definitions; they were not additional benchmarks executed here.

Training and evaluation design are equally important. Varma and Simon demonstrated how using the same cross-validation process for model selection and error estimation can bias performance estimates [@varma2006]. We separate nested source-value construction, target-held-out gate prediction and later performance evaluation, while acknowledging that the final K562 scaffold itself was selected using all fitting cells. The fit-only scores used to rank its residual actions therefore remain selection statistics. Reproducible computation and prospective analysis locks reduce opportunities for outcome-dependent changes, but neither establishes biological validity or guarantees an unbiased estimate of transfer to a new cohort.

# 3 Methods

## 3.1 Study design and scope

We evaluated graph-conditional evidence valuation (GCEV), a ridge-based procedure that predicts the usefulness of an evidence source for selecting a graph edit. Three successive simulation experiments used 48, 32 and 32 distinct generated networks. A subsequent developmental experiment evaluated prediction in a source-labelled K562 expression cohort, followed by a numerical audit of how evidence representations influence graph decisions. Experiment-specific repository protocols fixed designs, seeds, primary contrasts and hyperparameters before the corresponding evaluations. The studies did not train a reinforcement-learning policy or neural graph model. Predictive source utility, simulated edge recovery and biological regulatory accuracy were treated as distinct endpoints.

## 3.2 Predictive loss and evidence-assisted scores

For target $t$ and parent set $P$, we fitted linear ridge regression [@ridge; @sklearn] with an intercept and penalty $\alpha=1$. Predictor columns were standardized using training-set means and population standard deviations, with standard deviations floored at $10^{-9}$; the response was not standardized. On a disjoint evaluation set $V$ containing $n_V$ observations, loss was

$$
L_t(P;T,V)=
\frac{n_V^{-1}\sum_{i\in V}(x_{it}-\widehat f_{t,P,T}(x_{iP}))^2}
{\max\{\operatorname{Var}_{T}(x_t),10^{-9}\}}.
$$

An empty parent set predicted the training response mean. For an eligible insertion $a=(u,t)$, the expression gain was

$$
g_a=\operatorname{clip}\!\left(
\frac{L_t(P)-L_t(P\cup\{u\})}{\max\{L_t(P),10^{-9}\}},-1,1
\right).
$$

Gains were converted to average-tie ranks over the relevant complete candidate universe, $r_a=(\operatorname{rank}(g_a)-1)/\max(N-1,1)$, rather than separately within targets. This normalization is part of the selection rule and does not represent an effect size. Given source score $0\le v_{sa}\le1$, observed-evidence mask $m_{sa}$ and target–source weight $w_{st}$, the composite score was

$$
S_a=r_a+0.3\sum_s w_{st}m_{sa}v_{sa}.
$$

Fixed two-source weighting set all weights to one; expression-only selection omitted the source terms. Missing observations contributed no source term and were not interpreted as negative edges. Scores and weights were not calibrated probabilities. Exact score ties were resolved lexicographically.

## 3.3 Nested construction of source-utility labels

Source-utility labels used only the experiment's outer fitting cells. A shuffled three-fold cell partition defined, for each fold, a training portion and a disjoint validation portion. Each training portion was further randomly divided approximately 75/25 into search-fit and search-score subsets. The parent scaffold was selected on search-fit cells, and candidate gains were estimated using search-score cells. Thus, action selection and source-utility measurement used distinct observations.

For each target, fold and source $s$, two competitors shared the scaffold and inserted exactly one residual candidate. The source-ablated competitor maximized $r_a+0.3\sum_{q\ne s}m_{qa}v_{qa}$; the source-present competitor added $0.3m_{sa}v_{sa}$. Other sources retained fixed unit weights for this comparison. Both selected graphs were then refitted on the same complete training portion and scored on that fold's validation portion. The label was

$$
y_{fts}=\operatorname{clip}\!\left(
\frac{L^-_{fts}-L^+_{fts}}{\max\{L^-_{fts},10^{-9}\}},-1,1
\right).
$$

This label measures the predictive consequence of making one source available under the specified selection rule. It does not estimate source accuracy against regulatory truth or the optimal joint contribution of learned source weights. When both competitors select the same insertion, their utility difference is zero. Fold partitions used the experiment seed $z$, and search partitions used $z+51f+13$, for zero-indexed fold $f$.

## 3.4 Target-held-out gate learning and feature ablations

A second, separate cross-validation layer grouped all label rows by target gene. Three-fold, nonshuffled GroupKFold excluded every label belonging to the held-out target group when predicting its utilities. Features remained available for those targets; this was label exclusion by target, not removal of their expression from every other modelling context.

**Table 1. Features used by the source-utility gate.** The full gate used seven source features, two explicit parent features and four action features, plus one-hot source identity. Feature definitions are operational summaries of the selected scaffold and action; they are not full graph representations.

| Feature block | Variables |
|---|---|
| Source | Observed-source fraction; observed-score mean, maximum and population standard deviation; union coverage of other sources; mean absolute source disagreement over jointly observed candidates; residual candidate count |
| Explicit parent | Number of scaffold parents; mean absolute search-fit Pearson correlation between parents and target |
| Action | Expression gain of the source-ablated winner; source-ablated top-two score margin; tested-source score of the source-present winner; its tested-source score minus that of the source-ablated winner |

Features were standardized on gate-training rows, with standard deviations floored at $10^{-8}$. Ridge regression used penalty 10 and an intercept. For each held-out label row, predicted utility was mapped to a nonnegative weight, then averaged across its three cell-fold contexts:

$$
w_{st}=\frac{1}{3}\sum_{f=1}^{3}
\operatorname{clip}\!\left(\widehat y_{fts}/0.05,0,1\right).
$$

The global comparator substituted the training targets' source-specific mean utility for regression prediction. The historical no_graph_state comparator retained source features and source identity; its name denotes removal of explicit parent/action features, not independence from all graph context. Its source summaries still depend on the scaffold-restricted candidate set. The no_parent_state model removed only the two explicit parent features; no_action_state removed only action features. All used the same labels and target groups. Action features still encode scaffold-dependent information, so the parent ablation does not eliminate all graph dependence. Parent count was fixed at two in these evaluated contexts. Final weights averaged predictions from nested contexts; they were not obtained by applying a refitted gate to final-scaffold features. Figure 1 distinguishes this shared training logic from the final evaluation procedures.

![Figure 1. Source-utility construction and evaluation. The paired competitors share a scaffold, a candidate universe and one inserted edge. Source-utility measurement uses validation cells separate from the action-selection subsets. Target-held-out gate prediction is a separate grouping layer. The synthetic and K562 studies use different final selection and evaluation procedures, as detailed in Sections 3.6 and 3.8.](figures/F01_design.png)

## 3.5 Synthetic networks and controlled evidence

The first experiment used seeds 2001–2048, balanced over linear/nonlinear mechanisms, seven genes/140 cells versus nine genes/180 cells, and two source-reliability regimes. Each topologically ordered graph included predecessor edges independently with probability 0.30, with random signs and magnitudes uniform on [0.35,0.80]. Graphs were redrawn until at least two edges were present, for at most 100 attempts. Linear expression combined weighted parent values and independent standard-normal noise. The nonlinear mechanism replaced each parent contribution with $1.1w_{ut}\tanh(x_u)$ and, where two parents existed, added $0.30\tanh(x_{p_1}x_{p_2})$.

Two simulated sources were explicitly manufactured from the generating edge truth. Each candidate was observed with probability 0.77. Scores combined a value of 0.90 or 0.10 according to an effective truth label with Gaussian noise of standard deviation 0.07, then were clipped to [0,1]. In the target-switch regime, alternating target indices assigned reliability to opposite sources, while the other source inverted the truth label. In the globally stable regime, sources independently flipped labels with probabilities 0.12 and 0.25. Source-generation randomness used seed $z+888$. Four repeated conditions retained original evidence, removed 50% of observed entries, inverted 25% of retained scores, or combined both operations. Corruption used $v\mapsto1-v$, with source-specific seeds $z+500$ and $z+501$. These were oracle-derived simulations, not independent biological priors.

The mechanism experiment used 32 new DAGs, seeds 6101–6132, crossed with saturation absent/present and interaction absent/present. Each DAG retained identical weights, exogenous noise, source scores and cell partitions across its four mechanisms. Clean and combined-stress evidence produced 256 repeated network–condition fits. The explicit parent-feature experiment repeated that factorial design on seeds 7201–7232.

## 3.6 Synthetic selection, endpoints and uncertainty

Synthetic cells were partitioned 60/20/20 by seeded permutation. All genes were eligible regulators, excluding self-edges. For historical final edge scoring, each candidate used a two-parent scaffold selected by absolute training correlation while excluding that candidate regulator; gains were measured on the internal-score cells. This differs from the common-scaffold nested source-valuation contexts. Composite scores ranked all directed candidates globally. Predictive evaluation refitted models on fit plus internal-score cells using the top seven ranked edges, then evaluated the reserved test cells. No acyclicity restriction was imposed on estimated graphs.

Primary edge recovery used noninterpolated average precision (AP) against the complete simulated reference. The evaluation stage accessed that reference after rankings were frozen, although simulation truth had already been used to construct evidence inputs. Secondary endpoints included AUROC, precision at seven, direction recovery with ties retained separately and normalized prediction error. Direction accuracy was the proportion of eligible true-edge pairs for which the true direction received a strictly higher score than its reverse, among untied pairs. We report the mean within-network accuracy over assessable networks and the tie fraction separately; accuracy is undefined if every pair ties. Pearson, Spearman, conditional-gain and greedy procedures supplied additional controls. Greedy search used at most 1,600 candidate-loss queries; its reported ranking appended unaccepted fallback edges, so seven ranked edges need not represent seven accepted improvements.

Paired effects were summarized per independent DAG. Percentile 95% intervals used 4,000 stratified bootstrap resamples, seed 20261125: eight family-by-size-by-source strata of six DAGs for the first study, and four size-by-source strata of eight DAGs for each factorial study. Repeated mechanisms and stresses were not treated as independent networks. The first primary criterion required a positive lower AP-difference bound for GCEV minus no_graph_state and an upper relative test-NMSE difference bound of at most 0.02. The mechanism primary contrast was the interaction effect on that AP difference, averaged over saturation levels; it additionally required a point effect of at least 0.02. The parent experiment averaged GCEV minus no_parent_state across four clean mechanisms within each DAG, requiring an AP effect of at least 0.01, a positive lower bound and the same 0.02 NMSE bound.

## 3.7 K562 cohort and biological evidence representations

The developmental cohort originated from GEO GSM9224207 and contained 658 source-labelled K562 cells and 60,671 measured genes. The exported expression matrix was count-like; its exact raw-UMI status remained unverified. The original split, seed 20261126, retained 394 fitting, 131 internal-score and 133 reserved test cells. After inadequate TF coverage in an earlier gene universe, a declared developmental redesign selected 16 motif-annotated TFs and 48 additional genes by fitting-cell log1p variance. Eligibility required at least 5% detection, positive variance, an unambiguous primary-chromosome GENCODE v29 protein-coding annotation [@gencode] and the locked symbol exclusions. Selection favored motif-available TFs; the target category excluded motif-eligible regulators. All 64 genes were prediction targets, giving 1,008 nonself TF-to-gene candidates.

Motif scores used 17 JASPAR 2026 CORE vertebrate profiles [@jaspar] for the 16 TFs and 1,201-base GRCh38 promoters, −1000 to +200 relative to strand-aware gene-level TSSs. Sequences were retrieved from the UCSC hg38 sequence service [@ucsc]. Both strands were exhaustively scanned using uniform-background log-odds matrices. The best min–max-scaled motif score was contrasted with 128 deterministic exact dinucleotide-preserving promoter shuffles. The prior was the clipped positive excess above the null median, divided by one minus that median; a median within $10^{-10}$ of one yielded zero. The empirical p-value was the fraction of the 129 observed-plus-null values at least as large as the observation, with comparison tolerance $10^{-12}$. Benjamini-Hochberg adjustment was applied across all 1,008 candidates [@bh]. The continuous scores were not thresholded by adjusted significance and remained motif potential, not binding labels; all observed zeros retained observed masks.

DNase context came from two biosamples in ENCODE ENCSR000EOT [@encode]. Each replicate's promoter score was the fraction of promoter bases covered by the union of its narrowPeak intervals; the source value averaged both replicates. The primary file pair used different processing versions. Two matched-version sensitivities and a motif-only condition were retained. Accessibility was target-specific and broadcast unchanged across TF candidates. Expression and accessibility originated from different cultures.

## 3.8 Matched K562 prediction and source-action audit

BIO-V7 used seed 20261127. The primary transformation was $\log(1+X)$; a prespecified sensitivity used $\log(1+10^4X/\text{cell total})$, with totals over all 60,671 original genes. For every target, two eligible TFs with highest absolute fitting correlation formed one scaffold shared by all methods. Residual gains were averaged over three cell folds while holding this full-fit scaffold fixed. They were selection statistics, not unbiased validation of scaffold selection. Seven methods each chose one additional TF per target, yielding 128 scaffold edges and 64 insertions from 880 eligible actions. All eight transformation-by-prior cases were frozen before internal-score evaluation.

Ridge coefficients were fitted on the 394 cells. The primary endpoint averaged normalized squared errors across 64 targets and 131 internal-score cells; the 48 additional targets formed a secondary endpoint. The primary contrast was GCEV minus no_graph_state. Conditional percentile intervals resampled 131 paired cell-level mean error differences 10,000 times, seed 20261127. They omitted training uncertainty and biological-replicate variation; sensitivity contrasts were descriptive and unadjusted for multiplicity. Reserved test-cell outcomes and regulatory edge accuracy were not evaluated.

Target-constant additive accessibility cannot change within-target ordering, since $S_{ut}+c_t$ and $S_{ut}$ have the same maximizer. The implementation explicitly removed such offsets for selection while preserving original scores. BIO-V8 examined frozen contexts without fitting or outcome evaluation. Seven representations compared expression, motif, mean-accessibility-scaled motif, matched accessibility×motif, permuted accessibility×motif, additive motif-plus-accessibility and accessibility alone. A single target permutation used seed 20261128. Product terms were derived interactions, not additional independent sources. These diagnostics measured whether representations could alter decisions, not whether altered decisions improved prediction or regulation.

# 4 Results

## 4.1 The primary synthetic comparison did not establish a graph state advantage

All 48 intended network realizations completed without exclusions. Across the four source conditions and eleven evaluated methods, the study produced 2,112 full-universe rankings. The independent experimental unit was the generated network; repeated source conditions and method outputs were paired observations within that unit.

On clean source evidence, mean average precision was 0.545860 for GCEV and 0.538145 for the matched source-feature-only model. Their mean per-network difference was 0.007715, with a stratified 95% bootstrap interval of −0.010750 to 0.026722 (Figure 2). The lower bound did not exceed zero, so the predeclared AP criterion was not met. GCEV exceeded the comparator on 21 networks and tied on two. These counts describe the observed variation and do not replace the primary interval.

The companion mean relative difference in untouched expression-test NMSE was −0.000804, with an interval of −0.007940 to 0.006036. Its upper bound satisfied the specified 0.02 noninferiority allowance. Because the AP criterion failed, the joint primary rule was not satisfied. This combination illustrates why predictive error and regulatory edge recovery must be reported separately.

Fixed two-source weighting achieved the highest clean AP among the eleven evaluated methods in Table 2. GCEV minus fixed weighting was −0.074505 AP, with a network-bootstrap interval of −0.110828 to −0.036043. Thus, the added adaptation did not improve clean ranking recovery in this controlled source setting. Precision at seven was 0.508929 for GCEV, 0.517857 for the no-parent/action-feature model and 0.559524 for fixed two-source weighting. Their mean untied direction accuracies were 0.578402, 0.579620 and 0.689397, respectively. Pearson and Spearman tied every eligible reference-pair direction, so their untied direction accuracy was undefined (Supplementary Table S2). Pearson correlation had the lowest mean test NMSE (0.844371) and highest mean AUROC (0.858474) among all eleven methods, despite lower AP than fixed weighting. These endpoints therefore favored different models. GCEV exceeded expression-only conditional ranking by 0.048511 AP (95% interval 0.023896 to 0.075128) and greedy ranking by 0.049345 AP (0.012857 to 0.084605). These favorable secondary comparisons do not isolate the contribution of the explicit parent/action feature block and do not replace its primary comparison.

**Table 2. All eleven methods on clean synthetic evidence.** Values are equal-weight means across 48 generated networks. AP, AUROC and precision at seven use the complete synthetic edge reference; NMSE uses untouched expression-test cells. Higher ranking metrics and lower NMSE are preferable. Legacy single-source controls use only the primary prior.

| Method | AP | AUROC | Precision at 7 | Test NMSE |
|:--|--:|--:|--:|--:|
| Fixed two-source weighting | 0.620365 | 0.828134 | 0.559524 | 0.878056 |
| Legacy single-source fixed | 0.596476 | 0.818221 | 0.532738 | 0.875504 |
| Full GCEV | 0.545860 | 0.802644 | 0.508929 | 0.859175 |
| Legacy single-source heuristic | 0.533077 | 0.798141 | 0.497024 | 0.855846 |
| Global source weighting | 0.523710 | 0.800305 | 0.514881 | 0.861103 |
| Legacy single-source adaptive | 0.520019 | 0.819785 | 0.502976 | 0.849593 |
| No parent or action features | 0.538145 | 0.806540 | 0.517857 | 0.860085 |
| Expression-only conditional | 0.497349 | 0.793515 | 0.482143 | 0.855729 |
| Greedy expression model | 0.496515 | 0.795693 | 0.461310 | 0.851707 |
| Absolute Pearson correlation | 0.453455 | 0.858474 | 0.538690 | 0.844371 |
| Absolute Spearman correlation | 0.446882 | 0.852002 | 0.532738 | 0.847077 |

![Figure 2. Synthetic edge-ranking comparisons across the four predeclared source conditions. Panel A shows mean paired GCEV-minus-no-graph-state AP differences with network-stratified 95% bootstrap intervals. The clean comparison is primary; stress comparisons are secondary. Panel B shows mean AP for GCEV, the no-graph-state model and fixed weighting. The four conditions reuse the same 48 networks.](figures/F02_synthetic_evaluation.png)

## 4.2 Source stress altered performance without establishing a robustness advantage

With 50% source dropout, the GCEV-minus-no-graph-state difference was 0.001087 AP, with an interval of −0.008541 to 0.012163. With 25% score inversion, the difference was −0.010545, with an interval of −0.027204 to 0.005771. Under combined dropout and inversion, it was −0.013503, with an interval of −0.025987 to −0.001651. The combined-condition interval excluded zero in the unfavorable direction, but this was a secondary contrast without a separate multiplicity-adjusted confirmatory claim.

Fixed weighting decreased from 0.620365 AP in the clean condition to 0.486952 under combined stress. GCEV decreased from 0.545860 to 0.492404. The smaller decline for GCEV does not establish a better endpoint under stress: its combined-condition difference from fixed weighting was only 0.005452 AP, with an interval of −0.021591 to 0.033662. We therefore interpret the results as different sensitivity patterns rather than demonstrated robustness superiority.

Exploratory stratification suggested a clean GCEV-minus-no-graph-state effect of −0.0253 AP in linear networks and 0.0407 in nonlinear networks. These observations motivated a subsequent mechanism test. They were not predeclared subgroup superiority results and were not substituted for the unsuccessful primary comparison.

## 4.3 Neither follow up study met its primary criterion

The mechanism experiment used 32 new network seeds and paired linear, saturation, interaction and combined mechanisms. Its primary estimand was a difference in the incremental graph-state effect between interaction-present and interaction-absent mechanisms. The observed interaction effect was −0.002080 AP, with a 95% interval of −0.016369 to 0.011090. It failed both the positive-lower-bound criterion and the minimum 0.020 AP point-effect threshold. Consequently, this experiment did not confirm that two-parent interactions explained the earlier favorable nonlinear subgroup.

In the interaction-containing mechanism settings, the mean relative test-NMSE contrast was −0.006976 (95% interval −0.012342 to −0.001930), satisfying the 0.02 noninferiority allowance. This favorable predictive contrast concerns a different endpoint from the unsuccessful interaction AP criterion.

The separate 32-network explicit parent-feature ablation compared full GCEV with a model retaining source and action features but omitting parent count and parent-target correlation. Averaging across the four clean mechanisms within each network gave an AP difference of 0.005644, with an interval of −0.000168 to 0.011448. The positive-lower-bound and minimum 0.010 AP requirements were not met (Table 3). Its relative test-NMSE contrast was 0.000700 (95% interval −0.002264 to 0.003817), also satisfying the noninferiority allowance. Neither later study met its joint primary success rule.

**Table 3. Primary results of the three synthetic studies.** Each interval resamples independently generated networks within the declared strata. The estimands differ and are not pooled. Relative NMSE contrasts use the corresponding locked aggregation rules.

| Study and network count | Primary AP estimand | Estimate | 95% interval | Primary rule |
|:--|:--|--:|:--|:--|
| Main comparison, 48 | GCEV minus no graph state, clean | 0.007715 | −0.010750 to 0.026722 | Not met |
| Mechanism, 32 | Interaction difference in graph-state benefit | −0.002080 | −0.016369 to 0.011090 | Not met |
| Parent ablation, 32 | GCEV minus no parent state, clean | 0.005644 | −0.000168 to 0.011448 | Not met |

The parent-feature diagnosis limits the breadth of this negative finding. Parent count equalled two in all 256 condition-specific datasets, and therefore could not explain variation in source utility. Parent-target correlation varied, with mean within-condition standard deviation 0.1088, but adding the block did not produce the required advantage. Action features still depended on a scaffold. Removing the explicit parent block therefore did not remove every pathway through which graph information entered the learner.

## 4.4 K562 priors had distinct coverage and interpretation limits

All 64 promoter windows and 17 selected motifs were available for the 1,008 nonself TF-to-gene candidates. The motif procedure assigned positive continuous scores to 517 candidates. Forty-seven candidates had nominal empirical p-values below 0.05, but none passed a Benjamini-Hochberg threshold of q below 0.05; the smallest adjusted value was 0.853984. The continuous scores were retained as sequence-derived binding-potential features rather than thresholded discoveries or binding probabilities.

The original DNase files overlapped 47 promoters in one replicate and 41 in the other; 41 promoters overlapped both, and 17 overlapped neither. Mean accessibility across the 64 promoters was 0.171101. The original files used different peak-calling versions. Matched-version sensitivities retained 41 promoters with overlaps in both replicates but changed mean accessibility to 0.190095 and 0.146844. These are processing sensitivities within one DNase experiment, rather than independent biological replications.

## 4.5 The real expression comparison did not establish incremental predictive benefit

Seven methods were evaluated under eight prespecified transformation-by-prior conditions. Within each condition, every method used the same 128 scaffold edges and selected 64 insertions, yielding 192-edge predictive graphs. All 3,584 target-level choices were frozen before the 131-cell readout. Reserved test cells were not used for fitting, selection or performance evaluation.

For the primary original-expression scale and original DNase processing combination, GCEV NMSE was 0.402598, compared with 0.402822 for the no-graph-state model (Table 4). The paired difference was −0.000223, with a conditional 95% cell-bootstrap interval of −0.000636 to 0.000175. GCEV-minus-fixed and GCEV-minus-expression-only differences were 0.000473 and 0.001023, respectively, with intervals spanning zero (Figure 3). These estimates provide no clear primary-case advantage for full GCEV.

**Table 4. K562 primary-case prediction and graph choices.** NMSE averages errors on 131 internal-score cells, normalized by fitting-cell target variance. The additional-target column excludes the 16 candidate TFs from the error average. Insertions are compared with expression-only selection; all graphs retain the same scaffold.

| Method | NMSE across 64 genes | NMSE across 48 additional targets | Changed insertions |
|:--|--:|--:|--:|
| Expression only | 0.401575 | 0.368196 | 0 |
| Fixed sources | 0.402125 | 0.366089 | 33 |
| Global source gate | 0.401575 | 0.368196 | 0 |
| No graph state | 0.402822 | 0.368196 | 2 |
| No parent state | 0.402598 | 0.367898 | 3 |
| No action state | 0.401324 | 0.368196 | 2 |
| Full GCEV | 0.402598 | 0.367898 | 3 |

The no-action-state model had the lowest primary point error, but it was not promoted to the primary method after evaluation. Full GCEV and the no-parent-state model selected identical primary graphs. The target-only error summary also differed from the all-gene summary: fixed weighting had the lowest point NMSE among the seven models on the 48 additional targets. These changes reinforce the importance of retaining the predeclared endpoint and reporting its sensitivities.

![Figure 3. Developmental K562 prediction. Panel A displays primary-case NMSE for all seven methods. Panel B shows GCEV-minus-comparator differences and conditional paired-cell 95% bootstrap intervals. Negative differences favor GCEV. The intervals condition on the selected genes, fitted models and one cohort; they do not quantify between-culture variation or model-training uncertainty.](figures/F03_k562_prediction.png)

The two matched DNase versions preserved full GCEV's primary-scale choices. Removing DNase entirely altered the learned contexts: motif-only GCEV NMSE was 0.403298 versus 0.401852 for the no-graph-state model, a difference of 0.001446 with interval −0.000177 to 0.003121. This does not isolate an independently confirmed accessibility mechanism.

On the full-library normalized sensitivity scale, GCEV NMSE was 0.699570, the no-graph-state value was 0.699506 and the fixed-source value was 0.706015. GCEV-minus-fixed was −0.006444, with interval −0.011730 to −0.001303. In contrast, GCEV-minus-no-graph-state was 0.000065, with interval −0.000575 to 0.000726. The favorable fixed-comparator sensitivity did not replace the primary graph-state comparison. Absolute NMSE values across transformations concern different transformed responses and variance denominators and should not be read as a direct competition between preprocessing schemes.

## 4.6 Evidence influence and learned weight were not equivalent to utility

Accessibility was constant across candidate TFs for each target. Its additive contribution therefore cancelled from every within-target score difference. In each two-source K562 condition, adding accessibility changed zero of the 192 nested source-present versus source-ablated choices, and all observed accessibility utility labels were zero. Nevertheless, the shared ridge gate predicted positive accessibility weights for 47 of 64 primary targets, with mean weight 0.035305. Such predictions can arise from the shared regression fitted across source-target rows. They are not observations of accessibility utility.

Motif inclusion changed 96 of 192 nested choices on the original scale. Its mean paired utility was −0.003427, while the mean learned motif weight was 0.039244 and 39 targets received positive weights. Thus, the ability of a source to change a choice did not establish that its changed choices improved prediction on average. The clipping and averaging used to construct the weights also prevent their interpretation as calibrated utility estimates or probabilities.

The saved-output diagnosis compared seven representations at fixed coefficient 0.3. In the primary case, motif alone changed 33 of 64 decisions relative to expression-only scoring. Multiplying motif scores by global mean accessibility changed nine; multiplying by matched target accessibility changed eight; and multiplying by a frozen permutation of accessibility also changed eight (Figure 4). Equal counts did not imply identical choices: the matched and permuted representations differed from each other at eight targets. Their predictive errors and biological accuracy were not evaluated.

**Table 5. Primary representation diagnosis.** Counts refer to 64 target-specific selections from already saved residual candidates. They quantify action influence at fixed coefficient 0.3, not learned-model performance or a new expression evaluation.

| Representation | Targets with TF-varying source scores | Choices changed from expression only |
|:--|--:|--:|
| Accessibility only | 0 | 0 |
| Motif only | 64 | 33 |
| Motif plus accessibility | 64 | 33 |
| Global mean accessibility times motif | 64 | 9 |
| Matched target accessibility times motif | 47 | 8 |
| Permuted target accessibility times motif | 47 | 8 |

![Figure 4. Action influence of alternative evidence representations. Bars show the number of primary-case target decisions changed relative to expression-only selection, out of 64. Target-constant additive accessibility changes no decisions. Multiplicative accessibility rescales motif relative to expression and can change their combined choice; these counts do not measure usefulness. The matched and permuted representations are fixed saved-score diagnostics.](figures/F04_source_influence.png)

The invariance held across both expression transformations and all three accessibility versions: accessibility-only scoring preserved expression-only choices in all 384 target cases, and additive motif-plus-accessibility scoring preserved motif-only choices in all 384. Positive multiplicative accessibility preserved motif ordering within a target while changing motif's scale relative to expression. A zero accessibility value removed that motif contribution. The product remained a derived representation of two sources, not a third independent assay.

## 4.7 Reproduction supported numerical correctness within the executed settings

The synthetic execution audit recomputed ranking metrics for all 2,112 outputs of the 48-network study, with maximum discrepancy 1.11 × 10⁻¹⁶. Independent-environment K562 reproduction recovered all 3,584 choices, 56 method summaries and 24 contrasts, with zero reported differences in summaries and contrasts. A separate primary-case linear-algebra audit reproduced 448 regressions with maximum individual normalized-error discrepancy 2.98 × 10⁻¹⁴.

The saved-output diagnosis preserved 49,280 legacy residual-score rows and checked 2,688 representation-specific decisions. These checks support numerical reproduction and alignment of the comparisons. They neither add independent experimental units nor establish that the predictions correspond to biological regulation.

GCEV incurred additional fitting work. The 48-network study recorded 189 nested source-valuation loss evaluations for seven-gene networks and 297 for nine-gene networks. Fixed weighting does not perform an equivalent utility-learning stage. Matching candidate sets and inserted-edge budgets therefore did not match total training computation. No claim of computational-efficiency superiority follows from the reported performance comparisons.

# 5 Discussion

## 5.1 What the controlled comparisons establish

The three synthetic studies did not establish the proposed primary advantage of the tested graph-state representation. The real-expression comparison likewise did not show a clear incremental prediction benefit under its primary analysis. Together, these results constrain a specific formulation: a ridge model predicting paired one-edit source utility from the selected source, parent and action features, followed by clipped weighting. They do not establish equivalence to the controls or exclude useful effects under other representations and data regimes.

The fixed-weight result deserves particular attention. The synthetic priors were deliberately informative and derived from complete simulated truth. A fixed contribution could exploit this information without estimating a utility model from a small number of source-target observations. GCEV's utility labels instead rewarded expression prediction after an edit. These objectives need not align: a source that helps rank a true edge can have limited incremental predictive benefit when correlated parents already explain the target. The opposing ordering of methods by AP and NMSE is consistent with this distinction, although the experiments do not identify it as the sole cause of underperformance.

The subsequent mechanism study improved the interpretation of an initially favorable subgroup. The nonlinear-versus-linear difference in the first study could have encouraged a broader claim about interaction-aware inference. Testing a specified interaction estimand on new seeds did not support that explanation. Reporting the failed follow-up preserves the difference between a pattern that generates a hypothesis and evidence that confirms it.

## 5.2 A source must be able to affect the decision being valued

The accessibility diagnosis clarifies three separate questions. Is evidence observed? Can it change the decision? Do the resulting changes improve an evaluation outcome? Accessibility was observed for the chosen promoters, but its additive target-constant representation could not alter a within-target TF ranking. Motif scores could alter that ranking, but their mean nested utility was negative. Neither availability nor influence established usefulness.

The cancellation result follows directly from the scoring equation and is not a new mathematical theorem. Its value here is empirical and methodological: it identifies why a source-use interpretation would be invalid for this executed decision rule. Positive accessibility weights could otherwise be mistaken for learned support for that modality. The shared gate's predictions describe its regression behavior, rather than measured causal effects or calibrated source reliability.

Recent multimodal benchmarking has similarly emphasized testing the contribution of each modality [@beargrn]. Our local analysis complements such evaluations by specifying the decision set on which invariance holds. A target-dependent offset may matter when edges to different targets compete in a global ranking, even though it cancels when candidate TFs compete for one fixed target. The BIO-V8 analysis therefore retained the original across-target ranking score separately from the within-target action score. It did not overwrite a legacy ranking using an equivalence that holds only for a different decision problem.

An interaction between accessibility and motif can avoid the additive invariance, but introduces another question. Positive target accessibility scales all motif scores at that target equally. A reduction in changed choices may therefore result simply from attenuation relative to the expression score. Comparing matched accessibility with both a global-mean attenuation and a fixed permutation helps identify this ambiguity. In the present diagnosis, matched and permuted counts alone provide no evidence of a biological matching benefit. Evaluating that benefit would require newly locked predictive or biological endpoints.

## 5.3 Limits of graph state and source utility estimation

The explicit parent feature set was small, and its parent-count component was constant by design. The action features retained indirect information about the scaffold, so the parent ablation estimated the value of two explicit columns rather than the value of graph information in general. The experiments do not test variable-size parent sets, learned graph embeddings, global topology or a sequential policy. Conclusions about those alternatives would require different experiments.

The utility learner also faces sparse and heterogeneous labels. A label is exactly zero whenever the source-present and source-ablated procedures select the same action. When they select different actions, the difference is measured on a relatively small validation fold and clipped to a bounded range. The model pools rows across sources and predicts for held-out targets; source indicators and features do not guarantee exact zero predictions for a source with zero observed utility. The subsequent weight map removes negative predictions and saturates sufficiently positive values. These properties can generate weights whose interpretation differs from either mean observed utility or an edge's probability of being correct.

Target-held-out learning prevents a target's own utility labels from fitting its predictor, but does not make genes biologically independent. Genes share cells, sources and scaffold construction. Likewise, nested contexts differ from the final graph contexts to which the averaged weights are applied. These are potential sources of error and limits on interpretation; their individual causal contributions were not isolated by the completed experiments.

## 5.4 Scope of the biological evidence

The real-data analysis uses a single source-labelled K562 cohort, selected genes and a constrained one-insertion graph. Its uncertainty intervals condition on the fitted models and observed cells. They do not account for training variability, gene-selection uncertainty or variation across independently collected cultures. The reserved 133 test cells were not analyzed in this manuscript. Their future use could provide a final within-cohort predictive readout, but would still not establish transfer to another biological setting.

The priors have additional limits. Gene-level promoter coordinates are not verified active-transcript starts in these cells. A positive motif score is not evidence of occupied binding, and no edge passed the stated motif multiple-testing threshold. DNase coverage describes promoter accessibility rather than TF-specific occupancy; expression and accessibility came from different K562 cultures. Full publication-level assay QC and the lineage of potentially overlapping prior and reference experiments remain incomplete.

No qualified independent binding or perturbation panel was used to score the inferred real edges. Accordingly, we report no real-data AP, AUROC, F1 or causal-recovery metric. Previously inspected perturbation results cannot retrospectively become untouched confirmation. Binding, perturbation response and expression prediction would also require separate interpretations if future evaluation data become available. The operational difficulty of retrieving external assay metadata is an acquisition limitation, not evidence that an assay or regulatory relation is absent.

## 5.5 Implications for further evaluation

The next computational evaluation should address the limitations that matter for the claim: larger and materially different network-generating systems, meaningful variation in parent state, and relevant external inference baselines. The current simulations isolate selected components but do not constitute broad biological benchmarking. Repeating more seeds within the same seven- and nine-gene framework would improve Monte Carlo precision without addressing that limitation.

If the objective is to improve the source-valuation method, two hypotheses are suggested by the completed analyses: learning source utility with a representation aligned to the actual decision set, and modeling zero-utility or source-specific behavior more explicitly. These are hypotheses for future study, not validated fixes. Their assessment should use new locked outcomes rather than repeatedly tune on the present synthetic seeds or K562 evaluation partition.

# 6 Conclusion

We evaluated whether selected graph-state features improve adaptive weighting of biological evidence in a one-step GRN inference procedure. Three prospectively locked synthetic studies did not meet their primary success criteria, and the developmental K562 comparison did not establish an incremental predictive advantage. Fixed source weighting was stronger on clean synthetic edge recovery. A saved-output diagnosis showed that target-constant additive accessibility could not influence the within-target action being valued, although the learned gate assigned it positive weights. Alternative multiplicative representations changed decisions without establishing that those changes were useful. The findings support a scoped account of the tested formulation and motivate evaluating evidence availability, action influence and predictive or biological utility as separate quantities.

# Data and code availability

The public research repository is https://github.com/TshepoGobonamang/bio-rl-grn [@gcevdata]. The evidence snapshot used for this manuscript is commit `4328d45ec196a1dc05d46a6e36dc8480810597cb`. Stage-specific implementation commits, prospectively locked configurations, source identifiers and reproduction instructions are listed in Supplementary Table S1 and the accompanying evidence ledger. Earlier synthetic experiments require their recorded implementation versions because subsequent development expanded the method grid and adapted the real-data scaffold. Public biological sources comprise GEO GSM9224207 within GSE306708, GENCODE v29, UCSC hg38 sequence, the selected JASPAR 2026 motifs and ENCODE experiment ENCSR000EOT. The source manifests retain exact file identifiers and checksums. No new expression or perturbation experiment was conducted for the manuscript preparation.

# Supplementary information

## S1 Protocol freezes and computational reproduction

**Supplementary Table S1. Protocol and implementation versions.**

Commit links identify exact immutable repository versions. “Locked” denotes a repository precommit, not registration with an independent registry. Auditing saved rankings and rerunning the same source bytes establish computational reproduction, not biological replication.

| Experiment | Protocol record | Implementation or prediction freeze | Reproduction evidence |
|---|---|---|---|
| 48-DAG comparison; seeds 2001–2048 | Narrative/JSON: [cb7b3d2](https://github.com/TshepoGobonamang/bio-rl-grn/commit/cb7b3d2aff6827c98c370ba14fbc19293405a1c3) | Code/tests/workflow: [c4d98e9](https://github.com/TshepoGobonamang/bio-rl-grn/commit/c4d98e9b48c5fd99dc0d630d39003ac2d2f559e8) | Independent rescore of 2,112 rankings; zero failures; maximum AP/AUROC discrepancy $1.11\times10^{-16}$. |
| 32-DAG mechanism factorial; seeds 6101–6132 | Narrative: [f446a3d](https://github.com/TshepoGobonamang/bio-rl-grn/commit/f446a3d8764148dbfc4a885f334004bb74318fba); executable JSON: [0178b37](https://github.com/TshepoGobonamang/bio-rl-grn/commit/0178b372b7b532f25b0aa187a360e256774c95ca) | Runner: [eeb7c12](https://github.com/TshepoGobonamang/bio-rl-grn/commit/eeb7c126638cef61af8ec3baa39bb523c1873e8c); complete pre-evaluation code/configuration/workflow checkout: [7857ad6](https://github.com/TshepoGobonamang/bio-rl-grn/commit/7857ad6b020621cfdf45ca744ddc852df04a8a32) | Independent local audit of 2,816 rankings; all AP/AUROC reconciled within $10^{-11}$; no DAG failure. |
| 32-DAG parent-feature ablation; seeds 7201–7232 | Narrative: [b77f235](https://github.com/TshepoGobonamang/bio-rl-grn/commit/b77f235a2efd0d10a100413c5bf8ddc59a72cdd1) | Executable JSON/code: [bf185df](https://github.com/TshepoGobonamang/bio-rl-grn/commit/bf185df34efbd0262c0ff3c2c779dde2265234ca) | Completed local evaluation and independent rescore of 3,328 rankings. The cited record does not establish the separate GitHub full-run outcome. |
| K562 matched prediction, BIO-V7 | Protocol: [41baca2](https://github.com/TshepoGobonamang/bio-rl-grn/commit/41baca2de083b8b13ed602bcfcbf804e910f71cb) | Implementation: [d136b66](https://github.com/TshepoGobonamang/bio-rl-grn/commit/d136b66dfb59a6cde055f40668d56202a8806b7a); all prediction choices: [ac5c86c](https://github.com/TshepoGobonamang/bio-rl-grn/commit/ac5c86c820ec2b2178f769639d1e5c0e7de24ade) | Run [46004482269](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/46004482269) at [77d0552](https://github.com/TshepoGobonamang/bio-rl-grn/commit/77d0552c2a885540e5e665d541dbd4d7f7eafd91): identical 3,584 choices; zero differences in 56 summaries and 24 contrasts. |
| Saved-score evidence qualification, BIO-V8 | Follow-up design: [9169be4](https://github.com/TshepoGobonamang/bio-rl-grn/commit/9169be40c3658d45178ed170587f5058ad8c9942) | Implementation: [eb3ca37](https://github.com/TshepoGobonamang/bio-rl-grn/commit/eb3ca377bbbf02e7577c3f0aab6eca31cd61085f) | Run [46011226484](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/46011226484) at the implementation commit: eight output tables byte-identical. No new prediction or biological outcome evaluation. |

Later code expanded the synthetic comparator grid; reproducing the historical 11-method outputs therefore requires the corresponding historical version. BIO-V7 introduced a separately declared TF-constrained adaptation. BIO-V8 followed known developmental outcomes and cannot be interpreted as an independent retrospective confirmation. Its external metadata inventory remained incomplete after six HTTP 403 responses. Workflow artifacts have finite retention; they are not a permanent archival release.

## S2 Exact estimands and implementation details

For network $j$, let $E_{Gj}$ and $E_{Bj}$ be test NMSE for GCEV and its specified comparator. Relative prediction error is

$$
\rho_j=(E_{Gj}-E_{Bj})/\max(E_{Bj},\epsilon),
$$

with $\epsilon=10^{-12}$ in the 48-network script and $10^{-10}$ in the factorial scripts. The statistic averages these paired relative differences; it is not the relative difference of grand means.

For saturation $a\in\{0,1\}$, interaction $b\in\{0,1\}$, and $d_{jab}=AP_{G,jab}-AP_{\mathrm{no\_graph\_state},jab}$, the clean mechanism estimand is

$$
\theta_{\mathrm{mean}}=\frac1{32}\sum_j\frac12
\left((d_{j01}-d_{j00})+(d_{j11}-d_{j10})\right).
$$

Its accompanying NMSE contrast averages $(\rho_{j01}+\rho_{j11})/2$ within each DAG. The parent study instead averages GCEV-minus-no-parent AP and relative NMSE across all four clean mechanisms within each DAG, then across the 32 DAGs. Bootstrap resampling preserves these DAG-level paired units.

The motif null seed is the big-endian integer represented by the first eight SHA256 bytes of the UTF-8 string formed by the protocol identifier, a pipe separator and the **target gene symbol**, not its Ensembl ID. Canonical A/C/G/T iteration fixes shuffle reproducibility. With 128 null scores,

$$
p=\frac{1+\#\{v_{\mathrm{null}}\ge v_{\mathrm{observed}}-10^{-12}\}}{129}.
$$

The prior uses positive relative excess above the null median $q$: $\operatorname{clip}(\max\{0,(v_{\mathrm{observed}}-q)/(1-q)\},0,1)$ when $q<1-10^{-10}$, and zero otherwise. BH adjustment covers 1,008 candidates; continuous priors are not BH-thresholded.

## S3 Frozen biological inputs

The 16 TFs, in frozen selection order, were **TRPS1, PAX5, HES1, RUNX1, CUX1, ZIC1, HMBOX1, NFIA, BACH2, MEF2C, TCF4, ZNF423, EBF4, EBF1, ATF6 and SOX6**. TCF4 is not TCF7L2. Motif availability guided this developmental selection and does not demonstrate K562 binding.

Expression used GSM9224207_Fresh_K562.h5ad from GSE306708. Promoters used GENCODE v29 GRCh38 coordinates and UCSC hg38 sequence; motifs used JASPAR 2026 CORE vertebrates. ENCODE DNase experiment ENCSR000EOT supplied biosamples ENCBS545WKJ and ENCBS392ZOW. The primary peak pair was ENCFF185XRG/ENCFF433TIR; matched-v1 used ENCFF185XRG/ENCFF784DIM and matched-v2 used ENCFF821KDJ/ENCFF433TIR. These processed files represent one experiment, not independent validation panels.

The exported expression scale remains count-like rather than certified raw UMI. Gene-level promoters are not verified active-transcript starts; DNase is target context rather than TF-specific occupancy. Expression and accessibility cultures differ. Incomplete assay QC and unresolved prior-to-reference assay lineage limit biological claims.


## S4 Additional comparator definitions

The historical fixed, heuristic and adaptive controls use only the primary source, unlike fixed_sources and the GCEV family, which use both sources. The legacy fixed score adds 0.3 times the observed primary-source score to the global expression-gain rank. The heuristic additionally multiplies that evidence term by the positive part of the conditional gain. The legacy adaptive score adds 0.3 times the rank of a target-held-out ridge prediction of cross-fitted conditional gain, using absolute Pearson correlation, the masked primary prior, its observation mask and a rank-correlation-by-mask feature. These controls do not share the GCEV paired source-utility labels. They are contextual comparisons, whereas no_graph_state is the matched label-and-learner ablation. Greedy search begins from empty parent sets, accepts positive-loss-reduction insertions within its query allowance and appends unaccepted fallback-ranked edges to complete the ranking.

## S5 Direction recovery and complete supplementary data

**Supplementary Table S2. Clean synthetic direction recovery across 48 networks.** Accuracy first divides correct directions by untied assessable pairs within a network, then averages across assessable networks. The tie fraction is averaged separately across all networks. Undefined values indicate that no untied direction was available, not zero accuracy.

| Method | Mean untied direction accuracy | Mean tie fraction |
|:--|--:|--:|
| Fixed two-source weighting | 0.689397 | 0.000000 |
| Legacy single-source fixed | 0.639484 | 0.000000 |
| Full GCEV | 0.578402 | 0.000000 |
| Legacy single-source heuristic | 0.561922 | 0.000000 |
| Global source weighting | 0.547216 | 0.000000 |
| Legacy single-source adaptive | 0.518251 | 0.001894 |
| No parent or action features | 0.579620 | 0.000000 |
| Expression-only conditional | 0.509885 | 0.000000 |
| Greedy expression model | 0.486853 | 0.000000 |
| Absolute Pearson correlation | Undefined | 1.000000 |
| Absolute Spearman correlation | Undefined | 1.000000 |

The supplementary data include all 8,256 saved synthetic method rows, 236 method-by-source-by-mechanism summaries, 92 summaries averaging mechanisms within networks, the original 48-network bootstrap contrasts, all 56 K562 method-condition summaries and 24 paired contrasts, and all 42 source-representation summaries. These are summaries or copies of completed analyses; their assembly introduced no new hypothesis tests or model fits. The original confidence intervals were retained unchanged. File identities, source paths and checksums are recorded in the evidence ledger.

The shared GCEV training-time and nested-loss columns are duplicated across method rows in the historical exports. Their cost summary deduplicates by network, source condition and mechanism. Nested source valuation required 189 loss evaluations for seven genes and 297 for nine. The separate greedy utility-query counter recorded 273 and 483, respectively. A zero in the greedy counter for GCEV does not mean zero training work. These counters concern different stages and do not establish matched total computation.

**Supplementary Table S3. Companion relative test-NMSE results.** Negative contrasts favor GCEV. These are the original locked companion estimands, not ratios calculated from the method means in Table 2. All satisfy the specified upper noninferiority allowance of 0.02; failed AP criteria remain failed.

| Study | Relative NMSE contrast | 95% interval |
|:--|--:|:--|
| Main clean comparison | −0.000804 | −0.007940 to 0.006036 |
| Mechanism interaction-containing settings | −0.006976 | −0.012342 to −0.001930 |
| Parent ablation across clean mechanisms | 0.000700 | −0.002264 to 0.003817 |

# References

::: {#refs}
:::
