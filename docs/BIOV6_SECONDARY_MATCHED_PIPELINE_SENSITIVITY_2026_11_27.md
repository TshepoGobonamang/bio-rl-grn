# BIO-V6S — Pipeline-matched DNase sensitivity, declared after primary source-audit result
**Status:** secondary robustness analysis, frozen before examining any pipeline-matched promoter coverage. It is **not** preregistered independent confirmation.

The original BIO-V6 [GitHub run 45991610752](https://github.com/TshepoGobonamang/bio-rl-grn/actions/runs/45991610752) matched both GRCh38 K562 DNase isogenic biosamples, but the fixed lexicographic file-selection rule selected:
- Replicate 1: ENCFF185XRG, ENCODE peak-calling-step-v-1-0
- Replicate 2: ENCFF433TIR, ENCODE call-hotspots-se-step-v-2-1.

This heterogeneous processing makes replicate coverage potentially less comparable. Do **not** discard or replace the first observations (41/64 both-replicate promoters, mean 0.1711009, 384 positive joint motif×DNase scores).

## Fixed secondary sensitivity design
Use the same **ENCSR000EOT**, exact same 64 BIO-V5 promoter coordinates, same 1008 directed motif scores and no ENCODE CRISPRi expression outcomes. Re-query authentic ENCODE metadata and verify original files, same two distinct biosamples, assays/assembly and original MD5/size.

Compare **both** of the source-defined, version-specific, one-file-per-replicate pipelines (never choose whichever scores are favorable):
- v2: `/analysis-step-versions/dnase-call-hotspots-se-step-v-2-1/`
- v1: `/analysis-step-versions/dnase-peak-calling-step-v-1-0/`

For each version and biological replicate, choose lexicographically first released GRCh38 narrowPeak `peaks` file from its authentic eligible source-file records. If any needed file is missing, mark the version unavailable, not guessed or replaced.

Compute the **same** 1,201 bp union peak overlap per original promoter, per replicate coverage, mean, min and joint descriptive motif×mean-DNase product. Compare presence concordance, continuous Spearman correlation between two original biological replicates, differences between v1/v2 mean coverage for each gene, observed peak count, positive motif×accessibility products and source SHA256. Preserve original heterogeneous pipeline baseline.

Any method-selection decision based on these results is **exploratory**; the v1/v2 comparison diagnoses implementation/source sensitivity, not a new independent biological experiment or TF binding validation. A promoter accessibility score is target-only, not a distinct TF-specific direct-binding source; no fitted GCEV superiority conclusions.
