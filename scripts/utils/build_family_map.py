#!/usr/bin/env python3
"""
build_family_map.py

Build a family_map.tsv mapping each sequence ID in a SortMeRNA database to its
rRNA family and SILVA taxonomy levels.

Usage (called from build_sortmerna_index.sh):
    python3 build_family_map.py <output.tsv> <fasta1:rna_family1> [<fasta2:rna_family2> ...]

Columns written:
    seq_id, original, rna_family, domain, phylum, class, order, tax_family, genus, species
"""

import gzip
import re
import sys
from pathlib import Path

HEADER = ["seq_id", "original", "rna_family", "domain", "phylum", "class", "order", "tax_family", "genus", "species"]
N_TAX_LEVELS = len(HEADER) - 3  # 7


def parse_header(line: str) -> tuple[str, str, list[str]]:
    """Parse a FASTA header line into (seq_id, original, taxonomy_levels).

    Handles SILVA format: '>ID taxonomy;levels;here;size=N'
    and bare headers:      '>ID'

    'original' is the full header line (without the leading '>') preserved as-is.
    """
    line = line.lstrip(">").rstrip()
    original = line
    parts = line.split(" ", 1)
    seq_id = parts[0]
    tax_str = parts[1] if len(parts) > 1 else ""
    tax_str = re.sub(r";size=\d+$", "", tax_str)
    levels = tax_str.split(";") if tax_str else []
    levels += [""] * (N_TAX_LEVELS - len(levels))
    return seq_id, original, levels[:N_TAX_LEVELS]


def write_family_map(output_path: str, fasta_family_pairs: list[tuple[str, str]]) -> int:
    """Write family_map.tsv from a list of (fasta_path, rna_family) pairs.

    Returns total number of sequences written.
    """
    total = 0
    with open(output_path, "w") as out:
        out.write("\t".join(HEADER) + "\n")
        for fasta_path, rna_family in fasta_family_pairs:
            opener = gzip.open(fasta_path, "rt") if fasta_path.endswith(".gz") else open(fasta_path)
            with opener as f:
                for line in f:
                    if not line.startswith(">"):
                        continue
                    seq_id, original, levels = parse_header(line)
                    out.write("\t".join([seq_id, original, rna_family] + levels) + "\n")
                    total += 1
    return total


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <output.tsv> <fasta:rna_family> ...", file=sys.stderr)
        sys.exit(1)

    output_path = sys.argv[1]
    pairs = []
    for entry in sys.argv[2:]:
        fasta_path, rna_family = entry.rsplit(":", 1)
        p = Path(fasta_path)
        if not p.exists() and Path(fasta_path + ".gz").exists():
            fasta_path += ".gz"
        elif not p.exists():
            print(f"WARNING: {fasta_path}[.gz] not found - skipping", file=sys.stderr)
            continue
        pairs.append((fasta_path, rna_family))

    total = write_family_map(output_path, pairs)
    print(f"  family_map.tsv: {total:,} sequences -> {output_path}")


if __name__ == "__main__":
    main()
