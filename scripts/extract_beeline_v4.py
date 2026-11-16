"""Extract one BEELINE v4 replicate from the official separately downloaded ZIP.

Official Zenodo record: 19009603. MD5 belongs to the publisher's archive.
"""
import argparse
import hashlib
from pathlib import Path
from zipfile import ZipFile

EXPECTED_MD5 = "9ea93065004f2f12022adeeaf3728edd"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True, help="BEELINE-data.zip from Zenodo v4")
    parser.add_argument("--run", default="dyn-LL-500-1")
    parser.add_argument("--output", default="data/raw/beeline/dyn-LL-500-1")
    parser.add_argument("--skip-checksum", action="store_true")
    args = parser.parse_args()
    archive = Path(args.archive)
    h = hashlib.md5()  # Dataset-provided integrity checksum, not a security primitive.
    with archive.open("rb") as src:
        for chunk in iter(lambda: src.read(1 << 20), b""):
            h.update(chunk)
    if h.hexdigest().lower() != EXPECTED_MD5 and not args.skip_checksum:
        raise ValueError(f"Expected v4 MD5 {EXPECTED_MD5}, found {h.hexdigest()}")
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    with ZipFile(archive) as zf:
        files = [n for n in zf.namelist() if not n.endswith("/")]
        matches = [n for n in files if f"/{args.run}/" in f"/{n}" and
                   Path(n).name in {"ExpressionData.csv", "GroundTruthNetwork.csv", "PseudoTime.csv"}]
        if not any(Path(n).name == "ExpressionData.csv" for n in matches):
            raise FileNotFoundError(f"Could not locate {args.run}/ExpressionData.csv in archive")
        for name in matches:
            (out / Path(name).name).write_bytes(zf.read(name))
        if not (out / "GroundTruthNetwork.csv").exists():
            raise FileNotFoundError(f"No GroundTruthNetwork.csv for run {args.run}")
    print("Extracted:", sorted(p.name for p in out.iterdir()))
    print("Destination:", out)


if __name__ == "__main__":
    main()
