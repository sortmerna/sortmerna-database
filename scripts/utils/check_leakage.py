#!/usr/bin/env python3
"""
Check that no seed sequences from the clustered database appear in the test
members file. Overlap means data leakage: reads simulated from test members
would be derived from sequences already in the database.

Usage:
    python3 check_leakage.py <seeds.fasta> <test_members.fasta>

Exit codes:
    0 — no leakage
    1 — leakage detected (or input file error)
"""

import argparse
import sys

import skbio
import skbio.io


def load_ids(fasta_path):
  return {seq.metadata["id"] for seq in skbio.io.read(fasta_path, format="fasta", constructor=skbio.DNA)}


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("seeds", help="Clustered database FASTA (seed/centroid sequences)")
  parser.add_argument("test_members", help="Test members FASTA (non-seed sequences)")
  args = parser.parse_args()

  try:
    seed_ids = load_ids(args.seeds)
    member_ids = load_ids(args.test_members)
  except OSError as e:
    print(f"  ERROR: {e}", file=sys.stderr)
    sys.exit(1)

  leaked = seed_ids & member_ids
  if leaked:
    print(f"  ERROR: {len(leaked)} seed sequence(s) found in test members — data leakage detected")
    for seq_id in sorted(leaked):
      print(f"    {seq_id}")
    sys.exit(1)

  print(f"  Leakage check OK: 0 of {len(member_ids)} test members are seeds")


if __name__ == "__main__":
  main()
