#!/usr/bin/env python3
"""
Parse a VSEARCH .uc file and write member IDs and cluster mapping in a single pass.

.uc column layout (tab-delimited):
  1  type          S=seed, H=hit/member, C=cluster summary
  2  cluster#
  3  length
  4  %identity     (* for seeds)
  5  strand
  6  unused
  7  unused
  8  CIGAR
  9  query_label   full header e.g. "ACC.1.1234 Archaea;Phylum;..."
  10 target_label  centroid label (same format), or * for seeds

Only H records are processed; S and C records are skipped.
query_label and target_label include taxonomy after the first space,
so the ID is extracted as the first whitespace-delimited token.

Usage:
    python3 parse_uc.py <file.uc> --member-ids <ids.txt> --mapping <mapping.txt>

Exit codes:
    0 — success
    1 — input file error
"""

import argparse
import sys


def parse_uc(uc_path, member_ids_path, mapping_path):
  with open(uc_path) as uc, \
       open(member_ids_path, "w") as ids_out, \
       open(mapping_path, "w") as map_out:

    map_out.write("member_id\tseed_id\n")

    for line in uc:
      if not line.startswith("H"):
        continue
      cols = line.rstrip("\n").split("\t")
      member_id = cols[8].split()[0]
      seed_id = cols[9].split()[0]
      ids_out.write(f"{member_id}\n")
      map_out.write(f"{member_id}\t{seed_id}\n")


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("uc", help="VSEARCH .uc file")
  parser.add_argument("--member-ids", required=True, metavar="FILE",
                      help="Output file for member IDs (one per line)")
  parser.add_argument("--mapping", required=True, metavar="FILE",
                      help="Output TSV: member_id, seed_id")
  args = parser.parse_args()

  try:
    parse_uc(args.uc, args.member_ids, args.mapping)
  except OSError as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
  main()
