# Graph-Conditioned Reliability Signals for Gene Regulatory Network Inference: A Leakage-Aware Benchmark

**DRAFT — results not yet complete, not a submitted or accepted paper**

## Abstract

To be written from the final frozen, source-backed benchmark outputs. No predictive gain or accuracy claims are made here.

## 1. Introduction

Directed GRN reconstruction is challenging because transcriptional correlations are not direct causal regulatory effects and biological evidence sources may be incomplete or overlap with benchmarking references. This study tests whether parent-conditional, source-blind reliability features improve prediction under matched information and compute budgets.

## 2. Methods

Use the definitions and leakage protocol in `docs/EXPERIMENT_PROTOCOL.md`. Keep predictive association distinct from causal regulation. Predefine candidate sets, source masks, baseline parameters, and all held-out splits before accessing reference edges.

## 3. Experiments

Report the final dataset table (source accession, hashes, cell/gene coverage, reference completeness, TF whitelist and independent realization IDs), computational resources, software commit, seeds, and fixed/corruption/dropout scenario matrix.

## 4. Results

**Pending actual independent benchmark execution.** Tables: per-network AP/AUROC, orientation accuracy with ties, untouched expression NMSE, source coverage and sensitivity, graph-edit utility-query and runtime. Synthetic smoke results are software tests, not figures for biological validation.

## 5. Discussion and limitations

Consider negative and null comparisons, direction ambiguity, source/reference overlap, reproducibility across independent graph realizations, and lack of confirmed controls/replicates in GSE216595.

## 6. Code and data availability

Fill in immutable tagged release, public dataset accessions, exact SHA256 and reexecution commands at release time.
