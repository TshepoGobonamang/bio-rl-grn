# BIO-V4 — Authentic GENCODE v29 / JASPAR 2026 TF-eligibility audit (26 Nov 2026)

## Original-source result
- GitHub Actions real-source run: https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45888943925 (successful).
- Artifact: `gcev-biov4-real-gencode-jaspar-eligibility`, ID **15150851690**, checksums saved in `biov4_source_eligibility.json`.
- Original K562-only GSM9224207 H5AD SHA256: `14eb3129e57769899e0d301f515531a7e410e52dc7488fb1f981c43b2ea99ca6`.
- Original fit-only 512 genes and 394/131/133 partitions unchanged; no ENCODE response data read.
- GENCODE v29 GRCh38 gzip SHA256: `e77a9568d8783f432a03f5e49d8f52d89e8b0e795382cfcb1d93bbd279b0c2b9`.
- JASPAR 2026 vertebrate CORE nonredundant MEME SHA256: `ae4f77a17256c7b8b7b73402810c1ecf4fc411e9ad7832cdc87d586dee9a4a15`.

## Observed original 512-gene selection
| Metric | Actual |
| --- | ---: |
| Selected source genes | 512 |
| Unique primary-chromosome GENCODE mapped symbols | 476 |
| Not mapped under fixed primary-chromosome annotation | 36 |
| Ambiguous mapped IDs | 0 |
| Unambiguous protein-coding genes | 410 |
| Mitochondrial prefix | 15 |
| Ribosomal prefix | 76 |
| Genes with JASPAR single-factor motif | **3** |
| Eligible observed TF regulators | **3** |

Exactly:
- TRPS1 (ENSG00000104447.12, MA1970.2), detected in 67.0051% of outer-fit cells
- PAX5 (ENSG00000196092.12, MA0014.4), 51.7766%
- HES1 (ENSG00000114315.3, MA1099.3), 46.7005%

**Decision: NO_GO_INSUFFICIENT_TF_COVERAGE_REVISE_FIT_ONLY_SELECTION.**

## Scientific interpretation
The initial 512-gene selection ranked fitting-only log1p variance, without TF eligibility. It retained many mitochondrial and ribosomal genes and only three genes with unambiguous, measured, single-factor JASPAR 2026 motifs. This is a genuine source-qualification failure for a graph-conditioned TF→target model, not a model outcome. The 512-gene selection must remain archived and its inadequacy reported, not silently overwritten.

## Next authorized development redesign
Freeze an **explicit new** training gene universe selected from the full 60,671-feature K562 expression matrix and the same original 394 fitting cells: (a) 16 unambiguous protein-coding TFs carrying JASPAR 2026 single-factor profiles, (b) 48 additional unambiguous protein-coding target genes not overlapping the TFs; both ranked by fit-only log1p expression variance above 5% detection, lexicographic tie break, with no ENCODE response or regulatory gold. Reuse the exact same 394/131/133 cells. If insufficient eligible TFs exist, fail closed without replacing assay sources. This would be a source-availability-aware proof-of-concept candidate universe and must not be described as unbiased random TF selection.

No new motif-to-target prior scores, promoter sequence scans, GCEV training or independent biological regulatory validation are claimed by this audit.

Source pages: https://www.gencodegenes.org/human/release_29.html ; https://jaspar2026.elixir.no/downloads/ .
