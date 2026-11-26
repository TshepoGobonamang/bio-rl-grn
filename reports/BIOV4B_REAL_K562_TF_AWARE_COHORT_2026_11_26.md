# BIO-V4B — Real K562 TF-aware fit-only development cohort (26 Nov 2026)

## Execution and reproducibility
- Successful real-original-source GitHub Action: https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45889615265
- Artifact: `gcev-biov4b-real-k562-tfaware-64-gene-development-cohort`, ID **15151091546**.
- Locked protocol: `docs/BIOV4B_LOCKED_TF_AWARE_COHORT_2026_11_26.md`, commit `816e4a427be52809aae5096776eda9e4ef396753`, before any new target-selection results.
- Original K562 H5AD SHA256: `14eb3129e57769899e0d301f515531a7e410e52dc7488fb1f981c43b2ea99ca6`; same original 658 cells and sealed fit=394, score=131, test=133 split.
- GENCODE v29 GRCh38 SHA256: `e77a9568d8783f432a03f5e49d8f52d89e8b0e795382cfcb1d93bbd279b0c2b9`; JASPAR 2026 CORE vertebrate nonredundant MEME SHA256: `ae4f77a17256c7b8b7b73402810c1ecf4fc411e9ad7832cdc87d586dee9a4a15`.

## Measured fit-only candidate universe
| Actual quantity | Count |
|---|---:|
| Source genes with K562 expression IDs | 60,671 |
| Eligible unambiguous, protein-coding, >=5%-detected single-TF motif regulators before selection | **394** |
| Eligible nonmotif targets before selection | **9,808** |
| Frozen selected TFs | **16** |
| Frozen selected target-only genes | **48** |
| Frozen measured genes | **64** |
| Candidate directed TF→gene pairs, excluding self loops | **1,008** |
| Original fit/internal/test cells retained | **394 / 131 / 133** |
| Nonzero entries in selected 658x64 matrix | **30,152** |

Selected 16 TFs in locked fit-only variance order:
`TRPS1, PAX5, HES1, RUNX1, CUX1, ZIC1, HMBOX1, NFIA, BACH2, MEF2C, TCF4, ZNF423, EBF4, EBF1, ATF6, SOX6`.

These are genes matching **vertebrate** JASPAR single-factor motifs. This is not evidence of context-specific TF activity, binding to candidate promoters, or direct K562 regulation.

**Important source-eligibility constraint:** gene/TF selection is based on motif-profile availability plus expression, which favors known motif TFs. This is a **development** cohort generated after the earlier 512-gene set failed regulatory coverage (three motif-eligible TFs); it must not be described as an independent unbiased confirmation universe.

## Target-gene exclusion limitation (source-documented, not silently corrected)
The locked exclusion rule rejects `RPL[0-9]` and `RPS[0-9]` but **does not exclude** `RPLP1` or `RPLP2`; both appear among 48 targets. This is a literal application of the declared rule, not evidence that the full target set is devoid of ribosomal products. A broader `^RP[LS]` filter, if desired, requires a transparently versioned sensitivity analysis or prospectively frozen different selection. No labels or outcomes were used to select these targets.

## Evidence and model gates
- All 64 gene names have unambiguous GENCODE IDs and fit-only quantified original-scale expression.
- All 16 TFs have published JASPAR 2026 CORE single-factor motif profiles.
- **No** TF-target binding matrix has been constructed: we have **not** scanned GRCh38 promoter sequences.
- **No** ENCODE CRISPRi outcome labels or regulatory gold edges have been used.
- **No** GCEV, graph-independent or fixed prior models fitted, and no biological validation outcomes asserted.

**Next stage:** BIO-V5 source-disjoint motif-to-promoter scan. Pin exact GRCh38 reference sequence FASTA, GENCODE v29 TSS definitions, promoter window coordinates and threshold/GC controls *before* computing TF→target scores. Since no source-controlled biological prior exists yet, current 1,008 edges are **candidates**, not validated links. Then freeze matched predictions before any *new* uninspected perturbation panel is opened.

## Results files (authentic GitHub artifact)
- `biov4b_fit_only_manifest.json`
- `selected_16_TFs.csv`
- `selected_48_targets.csv`
- `selected_gene_universe.csv`
- `cell_partitions.csv`
- `k562_tfaware64_original_scale.npz`

The previous initial 512-gene set remains archived in the earlier BIO-V3B artifact; do not overwrite it.
