# M2C-D8B — Source-derived experimental design assessment (2026-11-17)

## Observed evidence, distinguished from executed source audit

The **official NCBI GEO metadata** for [GSE216595](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216595) describes HUES66 human embryonic stem cells transduced with a library of 90 TF ORFs, cultured for seven days and profiled by 10x scRNA-seq. It lists eight scRNA-seq sample accessions GSM6681039 through GSM6681046, and one TF-amplicon mapping accession [GSM6681047](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6681047).

GSM6681047 states that TFs were linked to cells by PCR amplification of transcription-factor constructs and cell barcodes from cDNA; raw count matrices were processed with Cell Ranger, annotated cells in Scanpy, and construct assignments from custom Python. Assembly hg38. The available TF mapping supplement is `GSM6681047_180124_TFmap.csv.gz`. The series includes a processed `GSE216595_180124_perturb.h5ad.gz` file.

Relevant paper: [Joung et al., A Transcription Factor Atlas of Directed Differentiation](https://pmc.ncbi.nlm.nih.gov/articles/PMC10344468/). It reports independent targeted validation of selected TFs and cross-line reproducibility, **not** proof that the eight GEO libraries constitute eight independent biological replicate cultures for all TFs.

## Execution evidence observed in GitHub

The PR-head commit for M2D8 `dd5a7cfbdbc0363ee5f827214eb859a6c3090758` returned a successful general `Tests` workflow run `44120316582` when checked; this **does not** establish execution of `m2d8-observe-geo.yml`. No completed original GEO mapping-audit artifact was available to review via the connected run lookup. This project's available tools did not offer a workflow-dispatch action. A direct attempt to download the GEO supplement in the execution container failed due to blocked DNS; **the mapping CSV was not parsed**.

## Decision table

| Question | Evidence available | Decision |
|---|---|---|
| Human TF-overexpression assay exists | Official GEO series and sample records | PASS (source identity only) |
| TF barcode mapping data listed | GSM6681047 supplement and mapping protocol | PASS (public availability metadata) |
| Actual mapping CSV checksum, columns and TF assignments verified | No executed source artifact reviewed | PENDING |
| Actual H5AD obs/var schema verified | No H5AD download/inspection | PENDING |
| Exact sample-aware cell barcode join counts verified | No join output | PENDING |
| Experimental negative controls usable for every TF | Not verified at cell/construct assignment level | PENDING |
| Independent biological replicate cultures established | Eight GSM sequencing samples are not proof | PENDING |
| CollecTRI/assay experimental provenance certified disjoint | No provenance resolution | PENDING |
| Independent biological regulatory accuracy | No eligible experimental labels | NOT AUTHORIZED |

## Immediate next steps — **do not add another abstract gate**

1. Trigger the exact [M2D8 workflow](../.github/workflows/m2d8-observe-geo.yml) on `feat/m2d8-source-observation`; preserve its **real** run ID and `m2d8-actual-source-observation` artifact, or diagnose an actual failure log.
2. Obtain the H5AD archive (or header/obs-only export from a trusted runner), checksum it, inspect backed metadata and export observed cell IDs.
3. Using the actual source column names, run `scripts/audit_m2d6b_barcode_join.py` and report exact TF assignment match/mismatch/duplicate counts.
4. Read published experimental methods for actual biological independence, appropriate unperturbed controls and pseudobulk pairing. If unpaired, redesign rather than invent paired replicate IDs.
5. If pooling and replication cannot be resolved, classify GSE216595 as exploratory and compare alternative experimentally replicated perturbation references. Only afterwards define TF-response labels or compute independent retrieval statistics.

**Locked scientific status:** source identity documented, real execution unverified, assay-level confirmatory validation blocked. No AUROC/AUPRC or RL claims.
