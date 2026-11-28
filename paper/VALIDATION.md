# Manuscript preparation and verification

Date: 28 November 2026. Draft: v0.1.

## Scope of this preparation

The manuscript was assembled from the completed experiment reports, historical saved metrics and verified primary literature. No models were fitted, new hypothesis tests performed, confidence intervals replaced, or reserved K562 test observations opened during manuscript preparation. Recovery of saved synthetic metric CSVs restored the full prespecified endpoint presentation, including precision at seven, direction recovery, ties and cost metadata.

## Content verification

- Independent numerical review reconciled the primary estimates and intervals, all eleven clean synthetic comparator rows, the three original relative-NMSE estimands, the K562 summaries and the source-influence counts with saved results.
- Methods review checked the paired utility construction, order of clipping and averaging, target-held-out gate learning, historical code versions and the distinction between synthetic global ranking and K562 common-scaffold selection.
- Literature review retained 27 cited records, including the closest adaptive-prior and modality-contribution work. The paper makes no first-use claim for adaptive evidence integration and distinguishes cited methods from executed baselines.
- All 27 citation keys resolve; no unresolved insertion markers or scientific placeholders remain.
- The evidence ledger's 94 recorded SHA256 values were reconciled with the exact local input bytes. This check reads only the explicitly recorded documentary and saved-summary files.

## Document verification

The final DOCX was rendered to a 22-page PDF using the selected runtime's Word-compatible renderer. Every page was visually inspected. The final file contains four figures, five main tables and three supplementary tables, with editable native Word equations. Checks covered cell contents, repeated table headers, row boundaries, captions, axis labels, references, equations, page numbers and page bounds.

The review corrected missing Pandoc table styles and renderer discrepancies for non-round equation delimiters and an overbar. Equivalent round grouping, an explicitly defined evaluation-set size and a mean subscript preserve the mathematics. The final four affected pages were inspected again; the other eighteen page images were byte-identical to the previously reviewed version.

Pandoc version: 3.1.11.1. The reference document preserves Pandoc's native styles before applying manuscript formatting. The generated LaTeX source has the same resolved manuscript content, but was not compilation-validated in this environment because TeX format files were unavailable. The reviewed PDF is rendered from the final DOCX.

## Final review files

| File | Bytes | SHA256 |
|---|---:|---|
| GCEV_Journal_Manuscript_Draft_v0_1_2026_11_28.docx | 538940 | eac93c94b792b8f0ad7a609b74624769f7dafe8b6cc06ea514fa3225ad9c5c02 |
| GCEV_Journal_Manuscript_Draft_v0_1_2026_11_28.pdf | 881353 | 20dfed1be2c1d2132e40b152d72c83e80ba7d69542d164ba938ce1b11c1b5e75 |

The review establishes consistency and readable presentation of the completed work. It does not close the independent biological validation gate, establish primary superiority, or constitute journal submission approval. Author information and declarations remain for author review in `AUTHOR_REVIEW_NOTES.md`.
