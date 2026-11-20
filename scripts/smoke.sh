#!/usr/bin/env bash
set -euo pipefail
python -m pytest -q
python -m grnbench smoke --output results/smoke --seed 42
