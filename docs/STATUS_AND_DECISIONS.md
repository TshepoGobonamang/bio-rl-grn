# Project status and gated delivery

## Current record: 27 November 2026

The [ten-track report](../reports/TEN_TRACK_RESEARCH_STATUS_2026_11_27.md) is the current status. [BIO-V8](../reports/BIOV8_EVIDENCE_QUALIFICATION_2026_11_27.md) exported and verified frozen scores, executed representation diagnostics, and reproduced the numerical outputs on GitHub. External ENCODE metadata searches failed with HTTP 403, so source inventory remains incomplete and biological scoring is blocked. All 133 expression-test cells remain sealed. No GCEV or RL superiority is established.

The older stage table below is preserved as a **20 November historical checkpoint**, not a current assertion that all later source and comparison work is unexecuted.

## Historical checkpoint: 20 November 2026

**2026-11-20:** clean-room software scaffold built and locally tested; no real independent run is claimed.

| Stage | Gate | State |
|---|---|---|
| Engineering | Package, CLI, unit tests, deterministic smoke data | Implemented; local tests run |
| Evidence | Genuine TF-map/H5AD barcode concordance | Pending; retain old project as historical record |
| Benchmark freeze | Independent unseen network source/checksums/candidate set | Required before confirmatory run |
| Benchmark execution | Per-network AP, orientation, untouched test loss and controls | Not run on new independent datasets |
| Reliability | Fixed vs adaptive vs heuristic under source dropout/corruption | Implemented small-scale; independent validation pending |
| RL | Budget-matched add/remove/stop trained policy | NOT implemented |
| Manuscript | Completed methods/results/figures from measured values | Draft structure only |

The full old experiment history is in `TshepoGobonamang/bio-rl-grn` (issue #36 and PR #37). Do not automatically merge or transfer unreviewed experimental scripts into this clean project. Future code changes should follow this program's data contract and test requirements.
