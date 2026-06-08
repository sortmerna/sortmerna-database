#!/usr/bin/env python3
"""
hard_mask_fasta.py

Convert soft-masked (lowercase) bases in a FASTA file to hard-masked (N).
Headers are preserved unchanged. Reports masking statistics to stdout.

Usage:
    python3 hard_mask_fasta.py --input in.fasta[.gz] --output out.fasta[.gz]
"""

import argparse
import gzip
import re
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class MaskingStats:
    total_sequences: int = 0
    masked_sequences: int = 0
    total_bases: int = 0
    masked_bases: int = 0

    @property
    def masked_sequences_pct(self) -> float:
        return self.masked_sequences / self.total_sequences * 100 if self.total_sequences else 0.0

    @property
    def masked_bases_pct(self) -> float:
        return self.masked_bases / self.total_bases * 100 if self.total_bases else 0.0


def hard_mask_sequence(seq: str) -> tuple[str, int]:
    """Replace all lowercase bases with N. Returns (masked_seq, n_masked)."""
    masked = re.sub(r'[a-z]', 'N', seq)
    n_masked = sum(1 for a, b in zip(seq, masked) if a != b)
    return masked, n_masked


def open_fasta(path: str, mode: str = "r"):
    """Open a FASTA file, handling .gz transparently."""
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, mode + "t")
    return open(p, mode)


def hard_mask_fasta(input_path: str, output_path: str) -> MaskingStats:
    """
    Read a FASTA, hard-mask all lowercase bases to N, write output.
    Returns MaskingStats with counts of sequences and bases affected.
    """
    stats = MaskingStats()
    current_header = None
    current_seq_parts: list[str] = []

    def flush_sequence(out):
        nonlocal current_header, current_seq_parts
        if current_header is None:
            return
        seq = "".join(current_seq_parts)
        masked_seq, n_masked = hard_mask_sequence(seq)
        stats.total_sequences += 1
        stats.total_bases += len(seq)
        stats.masked_bases += n_masked
        if n_masked > 0:
            stats.masked_sequences += 1
        out.write(current_header + "\n")
        # Preserve original line width (60 chars)
        for i in range(0, len(masked_seq), 60):
            out.write(masked_seq[i:i+60] + "\n")
        current_header = None
        current_seq_parts = []

    with open_fasta(input_path) as inp, open_fasta(output_path, "w") as out:
        for line in inp:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                flush_sequence(out)
                current_header = line
                current_seq_parts = []
            elif line:
                current_seq_parts.append(line.rstrip())
        flush_sequence(out)

    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input",  required=True, help="Input FASTA (plain or .gz)")
    ap.add_argument("--output", required=True, help="Output FASTA (plain or .gz)")
    args = ap.parse_args()

    if not Path(args.input).exists():
        print(f"ERROR: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    stats = hard_mask_fasta(args.input, args.output)

    print(f"Hard masking complete: {args.input} → {args.output}")
    print(f"  Sequences:      {stats.total_sequences:,}")
    print(f"  Masked seqs:    {stats.masked_sequences:,} ({stats.masked_sequences_pct:.2f}%)")
    print(f"  Total bases:    {stats.total_bases:,}")
    print(f"  Masked bases:   {stats.masked_bases:,} ({stats.masked_bases_pct:.4f}%)")


if __name__ == "__main__":
    main()
