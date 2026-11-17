"""Execute leave-one-network-out M2B meta-target validation."""
from __future__ import annotations
import argparse
from pathlib import Path
from biorlgrn.m2_reliability import export_m2b


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--m2a-root", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    scores, meta = export_m2b(args.m2a_root, args.output_dir)
    print(scores.to_string(index=False))
    print(f"Evaluated {len(meta['held_out_networks'])} held-out topologies.")


if __name__ == "__main__":
    main()
