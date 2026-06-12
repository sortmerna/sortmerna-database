#!/usr/bin/env python3
"""
compute_masking_stats.py

Count sequences with soft-masked (lowercase) bases in a FASTA file.
Soft masking comes from VSEARCH's DUST algorithm applied during clustering
(--qmask dust, the default), which lowercases low-complexity regions in centroid output.

Outputs one TSV row per FASTA:
    subunit  domain  total_seqs  silva_masked  silva_pct

Usage:
    python3 compute_masking_stats.py --fasta FILE --subunit STR --domain STR \\
        --threshold STR --out TSV
"""

import argparse
import gzip
import sys
from pathlib import Path


def open_fasta(path: str):
    p = Path(path)
    return gzip.open(p, "rt") if p.suffix == ".gz" else open(p)


def count_softmasked(fasta_path: str) -> tuple[int, int]:
    """Return (total_seqs, seqs_with_lowercase) from a FASTA."""
    total = masked = 0
    has_lower = False
    with open_fasta(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if total > 0 and has_lower:
                    masked += 1
                total += 1
                has_lower = False
            else:
                if not has_lower and any(c.islower() for c in line):
                    has_lower = True
    if total > 0 and has_lower:
        masked += 1
    return total, masked


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta",          required=True)
    ap.add_argument("--subunit",        required=True, help="e.g. SSU, LSU, 5S, 5.8S")
    ap.add_argument("--domain",         required=True, help="e.g. bacteria, archaea, eukaryota")
    ap.add_argument("--threshold",      required=True, help="e.g. 97, 95, 90, 85")
    ap.add_argument("--out",            required=True, help="Append-mode TSV output")
    args = ap.parse_args()

    if not Path(args.fasta).exists() and not Path(args.fasta + ".gz").exists():
        print(f"WARNING: {args.fasta} not found - skipping", file=sys.stderr)
        return

    fasta = args.fasta if Path(args.fasta).exists() else args.fasta + ".gz"

    total, vsearch_masked = count_softmasked(fasta)
    vsearch_pct = f"{vsearch_masked / total * 100:.2f}" if total else "0.00"

    write_header = not Path(args.out).exists()
    with open(args.out, "a") as f:
        if write_header:
            f.write("threshold\tsubunit\tdomain\ttotal_seqs\tvsearch_masked\tvsearch_pct\n")
        f.write(f"{args.threshold}\t{args.subunit}\t{args.domain}\t{total}\t{vsearch_masked}\t{vsearch_pct}\n")

    print(f"  {args.subunit} {args.domain} {args.threshold}%: {vsearch_masked}/{total} VSEARCH DUST soft-masked ({vsearch_pct}%)")


if __name__ == "__main__":
    main()
