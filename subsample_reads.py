#!/usr/bin/env python3
"""
Subsample paired-end FASTQ reads using reservoir sampling.

Finds R1/R2 pairs in raw_data/, randomly selects N read pairs,
and writes them to a subsampled_reads/ output directory.
"""

import gzip
import random
import re
import sys
from itertools import islice
from pathlib import Path

N_PAIRS = 1000
RAW_DATA = Path("raw_data")
OUT_DIR = Path("subsampled_reads")
SEED = 42


def open_fastq(path: Path):
    """Open a plain or gzipped FASTQ file."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


def read_records(fh):
    """Yield one FASTQ record (4 lines) at a time as a tuple."""
    while True:
        lines = list(islice(fh, 4))
        if not lines:
            break
        if len(lines) != 4:
            raise ValueError(f"Truncated FASTQ record: {lines}")
        yield tuple(line.rstrip("\n") for line in lines)


def find_pairs(directory: Path) -> list[tuple[Path, Path]]:
    """Return sorted list of (R1, R2) path pairs found in directory."""
    r1_files = sorted(directory.glob("*_R1_*.fastq*"))
    pairs = []
    for r1 in r1_files:
        r2_name = re.sub(r"_R1_", "_R2_", r1.name)
        r2 = r1.parent / r2_name
        if r2.exists():
            pairs.append((r1, r2))
        else:
            print(f"WARNING: no R2 found for {r1.name}, skipping.", file=sys.stderr)
    return pairs


def reservoir_sample_pairs(r1_path: Path, r2_path: Path, n: int, seed: int):
    """
    Reservoir sampling (Algorithm R) over paired FASTQ records.
    Returns a list of (r1_record, r2_record) tuples of length <= n.
    """
    rng = random.Random(seed)
    reservoir = []

    with open_fastq(r1_path) as fh1, open_fastq(r2_path) as fh2:
        for i, (rec1, rec2) in enumerate(zip(read_records(fh1), read_records(fh2))):
            if i < n:
                reservoir.append((rec1, rec2))
            else:
                j = rng.randint(0, i)
                if j < n:
                    reservoir[j] = (rec1, rec2)

    return reservoir


def write_records(records, path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt") as fh:
        for rec in records:
            fh.write("\n".join(rec) + "\n")


def main():
    OUT_DIR.mkdir(exist_ok=True)
    pairs = find_pairs(RAW_DATA)

    if not pairs:
        sys.exit(f"No R1/R2 pairs found in {RAW_DATA}/")

    for r1_path, r2_path in pairs:
        print(f"Sampling {N_PAIRS} pairs from:\n  {r1_path.name}\n  {r2_path.name}")

        sample = reservoir_sample_pairs(r1_path, r2_path, N_PAIRS, SEED)
        actual = len(sample)
        print(f"  -> {actual} pairs written (requested {N_PAIRS})")

        # Derive output names by inserting 'subsample' before the extension(s)
        suffix = "".join(r1_path.suffixes)  # e.g. .fastq.gz
        stem = r1_path.name[: -len(suffix)]
        out_r1 = OUT_DIR / f"{stem}.subsample{suffix}"
        out_r2 = OUT_DIR / re.sub(r"_R1_", "_R2_", f"{stem}.subsample{suffix}")

        write_records([r for r, _ in sample], out_r1)
        write_records([r for _, r in sample], out_r2)
        print(f"  -> written to {out_r1.name} / {out_r2.name}")


if __name__ == "__main__":
    main()
