# BIO-V5 — Real hg38 JASPAR motif-to-promoter prior: authenticated results

**Final completed execution:** 27 November 2026, [GitHub Actions #45988423247](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45988423247), successful. The artifact `gcev-biov5-actual-hg38-jaspar-motif-promoter-prior` is **15192626844**. Main protocol `docs/BIOV5_LOCKED_SEQUENCE_MOTIF_PROMOTER_PRIOR_2026_11_27.md` was committed at `7c5c16ac5cbd08e8eb6baf80a2c974c4810c2fc5` before obtaining any promoter sequence scores. Fixed numerical reproducibility correction in `d1f4d43ed73b3b234dea917d74e12a90504aaca7`.

## Frozen biological inputs
- Original unmodified fit-only K562 cohort from [BIO-V4B Action #45889615265](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45889615265): 16 motif-annotated TFs, 48 additional protein-coding target genes, 64 total genes; original 394 fitting, 131 internal-score, 133 untouched test K562 cells. No cell or gene resplitting or reference-based selection.
- Original GENCODE v29 GRCh38 gene features (SHA256 `e77a9568d8783f432a03f5e49d8f52d89e8b0e795382cfcb1d93bbd279b0c2b9`). **Gene-level** strand-aware TSS, not empirically verified or transcript-specific K562 active TSS. Plus and minus genes were oriented into 5'→3' transcription direction.
- **64 distinct 1,201 bp** promoter DNA regions (−1000 to +200 bp relative to approximate gene-level TSS) obtained through the [original UCSC hg38 sequence API](https://genome.ucsc.edu/goldenPath/help/api.html). Original genomic and gene-oriented DNA SHA256, chr/start/end and source URLs are retained for each gene.
- **17** selected JASPAR 2026 CORE vertebrates nonredundant single-factor MEME motifs representing the exact 16 source-identified TFs (BACH2 has two), original SHA256 `ae4f77a17256c7b8b7b73402810c1ecf4fc411e9ad7832cdc87d586dee9a4a15`. Matrix probability hashes saved per motif.
- No previously inspected ENCODE K562 CRISPRi downstream response or differential-expression outcome accessed; no gold direct edge labels, GCEV predictor fitting or expression-based score tuning.

## Real scanning method
A canonical uniform `b=(.25,.25,.25,.25)` background and JASPAR letter frequencies define `PWM_i(b)=log2(max(p_i(b),1e−6)/.25)`. Each true gene-oriented promoter was scanned exhaustively on **both** DNA strands; the best match over all registered motifs for each TF generated `observed_best_pwm_relative`, normalized by the motif's min/max attainable log-odds score.

For each **target gene**, generate 128 deterministic **exact dinucleotide-preserving** shuffled promoter controls through randomized Euler trails with canonical A/C/G/T iteration order, initialized by SHA256(protocol + gene ID). All TFs share the same control set per promoter. `empirical_p=(1+#null_best≥observed_best)/129`; BH adjustment is across all 1,008 TF→target candidate edges. Continuous prior score is `max(0,(observed_best−null_median)/(1−null_median))`, clipped to [0,1]; a saturated null median safely yields zero. This is **not** a calibrated biological binding probability, and exact zero is *observed* sequence evidence, not missing data.

## Original-source observed results — final deterministic run

| Result | Confirmed |
|:--|--:|
| hg38 gene-level promoters retrieved and SHA256 recorded | **64 / 64** |
| Length per promoter | **1,201 bp** |
| Motifs used | **17 from 16 TFs** |
| TF→gene directed candidates (self edges excluded) | **1,008** |
| Candidate scores with complete original source | **1,008 / 1,008** |
| Positive continuous prior scores | **517** |
| Mean prior score | **0.1711523** |
| Median prior score | **0.0126136** |
| Nominal empirical `p < .05` | **47** |
| BH adjusted `q < .05` | **0** |
| Smallest BH `q` | **0.8539840** |
| Scores `>=.8` | **37** |
| Of those >=.8 scores, empirical `p >= .05` | **25** |

The last row shows why a high PWM-based score must **not** be interpreted as confirmed TF binding: even high relative-score matches can occur often under the dinucleotide controls. In this 1,008-test developmental analysis **no TF→promoter edge passes BH q<.05**. About 5% of tests would fall below unadjusted p<.05 under an ideal null, and the 47 observed nominal hits alone do not demonstrate widespread motif enrichment beyond chance. Moreover, q-values from shared promoters and motif families are dependent; treat them descriptively, not as a complete binding reference.

## Reproducibility correction and audit
Two pre-final scans used an unintended Python set-iteration order for constructing the nucleotide transition dictionary; although both retrieved exactly the **same 64 promoter DNA sequences**, they obtained slightly different shuffled controls due to Python hash randomization. This was identified by an *across-run* audit, corrected by fixed `"ACGT"` iteration (commit `2aea3ddb46ed3c2cd9d0c45b69fbb7c35cf0a4b7`), and guarded by subprocess tests using distinct `PYTHONHASHSEED` values. The successful final run includes the source audit and independently reconciles edge counts, bounded prior scores, p-values with floating-precision tolerance, promoter lengths, per-file checksums, and absence of regulatory labels. **Only the final deterministic artifact is authoritative.**

Real TF→gene scores are archived under:
- `k562_tf_to_promoter_motif_prior_1008.csv` (SHA256 **`ea0bc14356fed5d2a8d4b5452e94e2372dc2f6e2da1d990a151ab1283a86b4e0`**)
- `gencode29_promoter_coordinates_and_sources.csv` (exact start/end and DNA source sha)
- `gene_oriented_hg38_promoters.fa` (64 real sequences)
- `jaspar_selected_motifs_sha256.csv` (17 original matrices)
- `biov5_motif_prior_manifest.json` (all metadata and source/derived hashes).

The strict adapter `scripts/biov5_prior_adapter.py` maps all 1,008 source-observed values into an exact candidate-ordered `{SOURCE_NAME: (np.ndarray scores, np.ndarray masks)}` compatible with the existing GCEV evidence interface; it rejects altered source hashes, changed edge universe or misinterpreted zeros. No actual model fitted in this stage.

## Interpretation, experimental limits and next decision
**BIO-V5 passes as a legitimate, source-provenanced *sequence motif binding-potential prior*.** It does not establish functional TF regulation, chromatin accessibility, in vivo occupancy, experiment-independent TF-target truth, causal direction, calibrated evidence probability or GCEV superiority. Selection remains developmental after the original 512-HVG three-TF NO-GO, and two ribosomal targets RPLP1/RPLP2 remain in the originally frozen 48 targets.

For biological testing, next (a) assess motif prior information against other independent prior streams, e.g., an ATAC accessibility component with K562-specific provenance, carefully excluding any held-out response experiments; (b) freeze identical K562 TF→target predictions for GCEV, graph-independent and fixed methods; and (c) choose a **new TF perturbation panel not already inspected** before outcome retrieval. The existing 6 ENCODE TF response comparisons have already been consumed for method development and cannot certify untouched confirmation.

### Source citations
- UCSC DNA API: https://genome.ucsc.edu/goldenPath/help/api.html
- GENCODE v29: https://www.gencodegenes.org/human/release_29.html
- JASPAR 2026: https://jaspar2026.elixir.no/downloads/
- JASPAR 2026 paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC12807658/
