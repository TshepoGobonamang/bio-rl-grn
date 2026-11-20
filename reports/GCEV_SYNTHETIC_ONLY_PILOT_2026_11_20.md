# GCEV v0.2 controlled synthetic engineering pilot — 20 November 2026

**Scientific verdict:** The proposed graph-conditioned source-ablation training objective was implemented, exercised on three synthetic networks, and compared with fixed and graph-independent controls. **It does not establish superiority, causal regulation, or independent biological source reliability.** Two evidence streams were deliberately generated from each simulation's gold edges, with complementary target-specific reliability and missing records. These tests are software-engineering pilots only.

## Source and protocol

- Implementation: `src/grnbench/gcev.py`; script: `scripts/controlled_synthetic_pilot.py`.
- Three generated DAG-expression networks, seed = **11, 42, 77**; each uses **130 cells, 7 synthetic genes**; same simulation family.
- Two simulated evidence streams per network (`SIM_A`/`SIM_B`), each retaining about 77% of candidate entries; target-regime-dependent reliable versus inverted evidence. **Sources derive from simulation truth**; no genuine TF motifs or perturbation assays.
- Partition: **60% outer fit, 20% internal graph scoring, 20% untouched expression test**. Paired evidence-value targets use only nested folds within outer fit. Gold reference is opened by the separate evaluator after rankings are sealed.
- Controlled methods: `conditional` (expression-only), `fixed_sources`, `global_source`, `no_graph_state`, `gcev`.
- Three source-valuation folds, one next graph edit per counterfactual; ridge alpha 1.0, parent scaffold up to two, source boost 0.3, gate shrinkage 0.05; GroupKFold by target.
- GCEV requires **42** matched source-value rows and **189** nested ridge-loss computations per network, plus group-regression fits. GCEV and fixed weights are **not training-compute matched**; paired ablations control only the one-step edit count.

## Observed noninterpolated AP (higher is better)

| Seed | Conditional | Fixed sources | Global source | No graph state | GCEV | GCEV − fixed |
|---:|---:|---:|---:|---:|---:|---:|
| 11 | 0.555411 | 0.568979 | 0.526648 | 0.688240 | **0.744589** | +0.175610 |
| 42 | 0.611111 | **0.691667** | 0.538194 | 0.652778 | 0.625000 | −0.066667 |
| 77 | 0.438889 | 0.473037 | 0.473090 | 0.460605 | **0.478027** | +0.004990 |

GCEV beats fixed source weighting on two synthetic realizations and loses on one. Against the graph-state-ablated learned gate, it also wins twice and loses once. These three development-generated networks are **not** independent confirmatory evidence or a basis for p-values.

## Observed untouched expression-test NMSE (lower is better)

| Seed | Conditional | Fixed sources | Global source | No graph state | GCEV |
|---:|---:|---:|---:|---:|---:|
| 11 | 0.952145 | 0.940303 | 1.047803 | 0.940303 | 0.940303 |
| 42 | **0.803825** | 0.806988 | 0.848253 | 0.815891 | 0.815891 |
| 77 | **0.863238** | 0.889073 | 0.889073 | 0.881032 | 0.889073 |

GCEV did not consistently improve generalization. Graph-state feature inclusion improved AP on seed 11 but not seed 42, and did not translate into universal expression-test improvements.

## Source-ablation action diagnostics

| Seed | Labels | Source changes chosen edge | Positive source benefit | Negative source benefit |
|---:|---:|---:|---:|---:|
| 11 | 42 | 4 | 3 | 1 |
| 42 | 42 | 2 | 1 | 1 |
| 77 | 42 | 6 | 4 | 2 |

When source inclusion does not change the selected edge, source-value is exactly zero. This signals **action insensitivity**, not proven evidence invalidity. The learned gate is not a calibrated TF-edge probability.

## Reproducibility and integrity

- **15 local tests passed**, including the test that mutates all outer scoring and final-test cell values while verifying GCEV label/weight invariance to floating-point tolerance.
- Sealed ranking and GCEV training artifact hashes are verified before evaluation references are opened.
- Per-(target,source) weights are trained by leaving that *target gene group* out of the gate model; both paired graph losses are measured on the same inner validation partition.
- `fit_manifest.json` records source SHA256s, partition hash, gated artifact hashes, nested loss-evaluation counts and training elapsed seconds.
- The pilot generator uses the synthetic ground truth to create controlled pseudo-sources but **does not pass that reference to the GRN fitting function**.

## Decision for the paper

The software demonstrates a **distinct training objective**, not a verified biological improvement. The estimand values sources for one next-addition decision in a source-blind parent scaffold, not arbitrary full graph trajectories. New preregistered independent networks and independent biological priors are required to establish relevance; computational cost must be compared fairly. Do not use these engineering outputs in a results section as independent biological evidence.
