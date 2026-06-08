#!/usr/bin/env python3
"""
compute_masking_stats.py

Count sequences with soft-masked (lowercase) bases in a FASTA file,
and optionally run dustmasker to produce a second independent count.

Outputs one TSV row per FASTA:
    subunit  domain  total_seqs  silva_masked  silva_pct  dust_masked  dust_pct

Usage:
    python3 compute_masking_stats.py --fasta FILE --subunit STR --domain STR \\
        --out TSV [--run-dustmasker] [--dustmasker-bin dustmasker]
"""

import argparse
import gzip
import re
import subprocess
import sys
from pathlib import Path


def open_fasta(path: str):
    p = Path(path)
    return gzip.open(p, "rt") if p.suffix == ".gz" else open(p)


def count_silva_softmasked(fasta_path: str) -> tuple[int, int]:
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
                if not has_lower and re.search(r'[a-z]', line):
                    has_lower = True
    if total > 0 and has_lower:
        masked += 1
    return total, masked


def count_dustmasker_softmasked(fasta_path: str, dustmasker_bin: str) -> int:
    """Run dustmasker and return the number of sequences it would soft-mask."""
    try:
        result = subprocess.run(
            [dustmasker_bin, "-in", fasta_path, "-outfmt", "fasta"],
            capture_output=True, text=True, check=True
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"WARNING: dustmasker failed: {e}", file=sys.stderr)
        return -1  # -1 = not available

    masked = 0
    has_lower = False
    in_seq = False
    for line in result.stdout.splitlines():
        if line.startswith(">"):
            if in_seq and has_lower:
                masked += 1
            in_seq = True
            has_lower = False
        else:
            if not has_lower and re.search(r'[a-z]', line):
                has_lower = True
    if in_seq and has_lower:
        masked += 1
    return masked


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta",          required=True)
    ap.add_argument("--subunit",        required=True, help="e.g. SSU, LSU, 5S, 5.8S")
    ap.add_argument("--domain",         required=True, help="e.g. bacteria, archaea, eukaryota")
    ap.add_argument("--out",            required=True, help="Append-mode TSV output")
    ap.add_argument("--run-dustmasker", action="store_true")
    ap.add_argument("--dustmasker-bin", default="dustmasker")
    args = ap.parse_args()

    if not Path(args.fasta).exists() and not Path(args.fasta + ".gz").exists():
        print(f"WARNING: {args.fasta} not found - skipping", file=sys.stderr)
        return

    fasta = args.fasta if Path(args.fasta).exists() else args.fasta + ".gz"

    total, silva_masked = count_silva_softmasked(fasta)
    silva_pct = f"{silva_masked / total * 100:.2f}" if total else "0.00"

    if args.run_dustmasker:
        dust_masked = count_dustmasker_softmasked(fasta, args.dustmasker_bin)
        dust_pct = f"{dust_masked / total * 100:.2f}" if total and dust_masked >= 0 else "NA"
        dust_masked_str = str(dust_masked) if dust_masked >= 0 else "NA"
    else:
        dust_masked_str = "NA"
        dust_pct = "NA"

    write_header = not Path(args.out).exists()
    with open(args.out, "a") as f:
        if write_header:
            f.write("subunit\tdomain\ttotal_seqs\tsilva_masked\tsilva_pct\tdust_masked\tdust_pct\n")
        f.write(f"{args.subunit}\t{args.domain}\t{total}\t{silva_masked}\t{silva_pct}\t"
                f"{dust_masked_str}\t{dust_pct}\n")

    print(f"  {args.subunit} {args.domain}: {silva_masked}/{total} SILVA-masked "
          f"({silva_pct}%)" +
          (f", {dust_masked_str}/{total} dustmasker ({dust_pct}%)" if args.run_dustmasker else ""))


if __name__ == "__main__":
    main()
