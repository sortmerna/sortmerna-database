#!/usr/bin/env python3

"""
fair_share_rfam.py

Compute per-family sequence allocations for Rfam sampling.

Given a target total and a dict of {family: available_sequences}, distributes
the target evenly. Families smaller than the per-family quota contribute all
their sequences; the remainder is redistributed to larger families until the
target is met or all sequences are exhausted.

CLI usage:
  python3 fair_share_rfam.py <target> <rfam_dir>

  Reads all RF*.fa files in rfam_dir, computes allocations, and prints
  tab-separated (stem, count) pairs - one per line - for use by the caller.

  Returns exit code 1 if total available < target (warning printed to stderr).
"""

import sys
from pathlib import Path


def fair_share(sizes: dict, target: int) -> dict:
    """Return {name: n_to_sample} allocating target sequences across families.

    Names are processed in sorted order for reproducibility. If total available
    sequences across all families is less than target, all sequences are taken
    and the returned total will be less than target.

    Args:
        sizes:  {name: available_count}
        target: desired total sequences across all families

    Returns:
        {name: allocated_count}
    """
    if target < 0:
        raise ValueError("target must be >= 0")

    remaining = target
    names = sorted(sizes)
    allocs = {}

    while names:
        per = remaining // len(names)
        small = [n for n in names if sizes[n] <= per]
        if not small:
            for i, n in enumerate(names):
                allocs[n] = per + (1 if i < remaining % len(names) else 0)
            break
        for n in small:
            allocs[n] = sizes[n]
            remaining -= sizes[n]
        names = [n for n in names if n not in allocs]

    return allocs


def count_sequences(path: Path) -> int:
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <target> <rfam_dir>", file=sys.stderr)
        sys.exit(1)

    target = int(sys.argv[1])
    rfam_dir = Path(sys.argv[2])

    files = sorted(rfam_dir.glob("RF*.fa"))
    if not files:
        print(f"Error: no RF*.fa files found in {rfam_dir}", file=sys.stderr)
        sys.exit(1)

    sizes = {f.stem: count_sequences(f) for f in files}
    allocs = fair_share(sizes, target)

    total = sum(allocs.values())
    if total < target:
        print(
            f"Warning: total available ({total}) < target ({target}); "
            f"all sequences will be used",
            file=sys.stderr,
        )

    for stem, n in sorted(allocs.items()):
        print(f"{stem}\t{n}")


if __name__ == "__main__":
    main()
