"""Generate development sensitivity and locked confirmation M2C-A2 datasets."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd

from biorlgrn.m2_stability import (
    DEVELOPMENT, CONFIRMATION, run_network,
)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw-root", type=Path, required=True)
    p.add_argument("--results-root", type=Path, required=True)
    p.add_argument("--role", choices=["development", "confirmation"], required=True)
    args = p.parse_args()
    names = DEVELOPMENT if args.role == "development" else CONFIRMATION
    output = []
    for name in names:
        record = run_network(
            args.raw_root / name, args.results_root / name, role=args.role)
        output.append(record)
        print(json.dumps(record, indent=2))
    stats = []
    for name in names:
        stats.append(pd.read_csv(
            args.results_root / name / "m2ca2_stability.csv"
        ))
    pd.concat(stats, ignore_index=True).to_csv(
        args.results_root / f"m2ca2_{args.role}_stability_summary.csv",
        index=False
    )


if __name__ == "__main__":
    main()
