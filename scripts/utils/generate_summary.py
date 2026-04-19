#!/usr/bin/env python3
"""
Generate a markdown clustering summary table from a TSV results file.

Input TSV columns (no header):
  ref_db, db, kingdom, clustering, num_seqs, total_nt

log10(#sequences) is computed here, not in the shell script.

Usage:
    python3 generate_summary.py <results.tsv> --output <summary.md>

Exit codes:
    0 — success
    1 — input file error
"""

import argparse
import math
import sys


def generate_summary(tsv_path, output_path):
  rows = []
  with open(tsv_path) as f:
    for line in f:
      ref_db, db, kingdom, clustering, num_seqs, total_nt = line.rstrip("\n").split("\t")
      num_seqs = int(num_seqs)
      log10_seq = f"{math.log10(num_seqs):.2f}" if num_seqs > 0 else ""
      rows.append((ref_db, db, kingdom, clustering, num_seqs, total_nt, log10_seq))

  with open(output_path, "w") as f:
    f.write("# Clustering Summary\n\n")
    f.write("| # | Ref DB | DB | kingdom | clustering | #sequences | total size (nt) | log10(#seq) |\n")
    f.write("|---|--------|-----|---------|------------|------------|-----------------|-------------|\n")
    for i, (ref_db, db, kingdom, clustering, num_seqs, total_nt, log10_seq) in enumerate(rows, 1):
      f.write(f"| {i} | {ref_db} | {db} | {kingdom} | {clustering} | {num_seqs} | {total_nt} | {log10_seq} |\n")


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("tsv", help="TSV results file written by cluster_sequences.sh")
  parser.add_argument("--output", required=True, metavar="FILE", help="Output markdown file")
  args = parser.parse_args()

  try:
    generate_summary(args.tsv, args.output)
  except OSError as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
  main()
