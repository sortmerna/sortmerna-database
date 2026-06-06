#!/usr/bin/env python3
"""
blast_family_breakdown.py

Parse a SortMeRNA BLAST-format file and count aligned reads per rRNA family
using a pre-built family_map.tsv. Appends one row per family to a cumulative
family_counts.tsv used by the PacBio sweep report.

Usage:
    python3 blast_family_breakdown.py \\
        --blast   aligned.blast[.gz] \\
        --map     family_map.tsv \\
        --seeds   50 \\
        --type    rrna|nonrrna \\
        --out     family_counts.tsv
"""

import argparse
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path


def open_file(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_family_map(map_path: str) -> dict:
    """Return {seq_id: rna_family} from family_map.tsv."""
    family_map = {}
    with open(map_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            family_map[row["seq_id"]] = row["rna_family"]
    return family_map


def count_families(blast_path: str, family_map: dict) -> Counter:
    """Count hits per rna_family from BLAST column 2."""
    counts = Counter()
    with open_file(blast_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            ref_id = cols[1]
            family = family_map.get(ref_id, "Unknown")
            counts[family] += 1
    return counts


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--blast",  required=True)
    ap.add_argument("--map",    required=True)
    ap.add_argument("--seeds",  required=True, type=int)
    ap.add_argument("--type",   required=True, choices=["rrna", "nonrrna"])
    ap.add_argument("--out",    required=True)
    args = ap.parse_args()

    family_map = load_family_map(args.map)
    counts = count_families(args.blast, family_map)

    write_header = not Path(args.out).exists()
    with open(args.out, "a") as f:
        if write_header:
            f.write("num_seeds\trna_type\tfamily\tcount\n")
        for family, count in sorted(counts.items()):
            f.write(f"{args.seeds}\t{args.type}\t{family}\t{count}\n")

    print(f"  {args.type} family breakdown ({sum(counts.values())} hits): "
          + ", ".join(f"{v} {k}" for k, v in counts.most_common(3)))


if __name__ == "__main__":
    main()
