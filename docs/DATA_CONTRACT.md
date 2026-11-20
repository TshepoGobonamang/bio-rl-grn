# Dataset contract and provenance

## Required input files

- `expression.csv`: first column unique cell ID, remaining columns exactly matched gene identifiers, numeric and finite; cells are rows.
- `tfs.csv`: one column `tf`, containing distinct measured TF gene identifiers. All-gene candidate regulators are permitted only as explicitly labeled toy smoke mode.
- `prior.csv` (optional): `source,target,score` directed edge records with unique pairs and strengths in [0,1]. **Missing row = unobserved evidence, not a negative TF edge**. Provider, exact release, assay, code and underlying study accessions must be recorded outside the scored CSV.
- `reference.csv` (evaluation only): `source,target` with complete eligible directed-edge reference for valid AP/AUROC. Incomplete source-specific positive lists cannot legitimately be used as complete gold.

## Source provenance manifest for real experiments

Maintain a separate human-reviewed manifest of source URL/accession, organism/genome assembly/cell context, collection timepoint, TF and gene naming scheme, SHA256, sample/culture replicate IDs, reference completeness, which constituent sources overlap the biological priors, exact programmatic extraction procedure, and any manipulation of barcodes/identifiers. Without this, classify the result exploratory.

## External perturbations

The historical GSE216595 analysis established H5AD 60,997 cells × 5,047 genes with 87 assigned TF labels, but matched controls and independent culture replicates remain unverified. The TFmap was likely headerless and no final concordance artifact was supplied with the latest project state. Thus it is **not a confirmatory direct TF→target ground-truth reference** in this clean repository. Do not claim eight 10x sequencing samples are eight independent cultures; the odd `temp` metadata values are not experimental temperature.

## Prior-reference overlap

The prior project found 2,848 overlapping directed edges between CollecTRI and a BEELINE hESC ChIP reference, without proving experiment-level source independence. CollecTRI-present edges are not calibration probabilities, and CollecTRI-unobserved edges are not verified biological negatives.
