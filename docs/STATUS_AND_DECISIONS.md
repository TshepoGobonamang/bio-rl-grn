# Project status and gated delivery

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
