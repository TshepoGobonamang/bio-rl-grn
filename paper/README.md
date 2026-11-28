# GCEV journal manuscript

This directory contains the complete first manuscript draft, dated 28 November 2026:

**A Controlled Evaluation of Graph Conditioned Evidence Weighting for Gene Regulatory Network Inference**

The draft reports the completed 48-, 32- and 32-network synthetic experiments, the developmental K562 prediction study, and the subsequent saved-score evidence-use diagnosis. It does not report a trained reinforcement-learning policy, a demonstrated primary graph-state advantage, or a validated biological network. All analyses described here were completed before manuscript preparation.

## Reading and editing

- `MANUSCRIPT.md` is the readable manuscript with resolved numerical citations and bibliography.
- `manuscript_source.md` is the editable Pandoc source, including equations, figure captions and supplementary information.
- `manuscript.tex` is the generated standalone LaTeX source with resolved citations. Its content is generated from the same source; the DOCX-derived PDF is the visually checked review edition. A complete TeX installation is needed to compile the LaTeX source.
- `references.csl.json`, `references.bib` and `numeric.csl` provide the reference data and numerical citation style. `reference_verification.json` records the literature verification scope.
- `figures/` contains four reproducible scientific figures in PNG, SVG and PDF formats.
- `data/` contains saved summary and metric tables, including all eleven original synthetic comparators and the original confidence intervals. `data/evidence_ledger.json` and `data/recovered_metrics_manifest.json` record input provenance and checksums.
- `AUTHOR_REVIEW_NOTES.md` identifies decisions and declarations that require the authors' review before journal submission.

The scientific evidence snapshot is repository commit `4328d45ec196a1dc05d46a6e36dc8480810597cb`. Supplementary Table S1 identifies the earlier implementation versions required to reproduce each historical experiment. The present branch adds the manuscript, figures, saved metric summaries and documentation; it does not rerun experiments or inspect the reserved K562 test outcomes.

## Rebuild the manuscript

Use Python 3 with `python-docx`, `matplotlib`, `pandas` and `numpy`, plus Pandoc with citeproc. From the repository root:

```bash
python paper/make_figures.py
python paper/build_manuscript.py --out /absolute/path/to/output
```

The first command reads saved summary tables. The second creates the editable DOCX, refreshes `MANUSCRIPT.md`, and regenerates the standalone LaTeX source. Neither command fits a model, evaluates reserved outcomes, or changes the original experiment results. Render the DOCX to PDF with a current Word-compatible renderer and inspect every page after layout changes. Figure and table captions are part of the manuscript source.

## Interpretation and review scope

Three prospectively locked synthetic studies did not meet their joint primary success rules. The primary K562 contrast also did not establish an incremental predictive benefit. Favorable secondary comparisons remain reported with their scope and uncertainty. The source-use analysis distinguishes evidence availability, influence on a selected action, and measured utility.

The manuscript is ready for substantive author review. Submission preparation still requires author and affiliation details, declarations, journal formatting, a durable release, and a decision about any further evaluation needed for the intended scientific claim. Draft completion does not close the project's independent biological validation gate.
