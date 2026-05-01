#!/usr/bin/env python3
"""
parse_cmsearch.py - Filter FASTA sequences based on cmsearch --tblout results.

Sequences with at least one above-threshold hit (inc == '!') are written to
the clean output, trimmed to the hit coordinates [seq_from, seq_to].
Sequences with no qualifying hit go to the flagged output.
A TSV log records the best hit details for every kept sequence.
"""

import argparse
from skbio import DNA, read


def parse_tblout(tblout_fp):
  """Return dict: seq_id -> best (score, evalue, seq_from, seq_to, strand) for '!' hits."""
  hits = {}
  with open(tblout_fp) as f:
    for line in f:
      if line.startswith('#'):
        continue
      fields = line.split()
      if len(fields) < 17:
        continue
      seq_id = fields[0]
      seq_from = int(fields[7])
      seq_to = int(fields[8])
      strand = fields[9]
      score = float(fields[14])
      evalue = fields[15]
      inc = fields[16]
      if inc == '!' and (seq_id not in hits or score > hits[seq_id][0]):
        hits[seq_id] = (score, evalue, seq_from, seq_to, strand)
  return hits


def trim_to_hit(seq, seq_from, seq_to):
  """Return sequence trimmed to cmsearch hit coordinates (1-based, inclusive)."""
  return seq[seq_from - 1:seq_to]


def read_fasta(fp):
  """Yield (full_header, sequence) for each record."""
  for seq in read(fp, format='fasta', constructor=DNA):
    header = seq.metadata['id']
    description = seq.metadata.get('description')
    if description:
      header = f"{header} {description}"
    yield header, seq


def main():
  ap = argparse.ArgumentParser(description=__doc__)
  ap.add_argument('--fasta',   required=True,  help='Input FASTA (DNA)')
  ap.add_argument('--tblout',  required=True,  help='cmsearch --tblout file')
  ap.add_argument('--clean',   required=True,  help='Output: sequences with a qualifying hit')
  ap.add_argument('--flagged', required=True,  help='Output: sequences with no qualifying hit')
  ap.add_argument('--log',     required=True,  help='Output: TSV hit details for kept sequences')
  ap.add_argument('--gene',    default='',     help='Gene label for stats row (e.g. ssu, lsu)')
  ap.add_argument('--domain',  default='',     help='Domain label for stats row (e.g. bacteria)')
  ap.add_argument('--stats',   default=None,   help='Append a summary row to this TSV file')
  args = ap.parse_args()

  hits = parse_tblout(args.tblout)
  print(f"  Sequences with qualifying hit (inc='!'): {len(hits)}")

  kept = flagged = 0
  with open(args.clean, 'w') as clean_f, \
     open(args.flagged, 'w') as flag_f,  \
     open(args.log, 'w') as log_f:

    log_f.write("seq_id\tscore\tevalue\tseq_from\tseq_to\tstrand\n")

    for header, seq in read_fasta(args.fasta):
      seq_id = header.split()[0]
      if seq_id in hits:
        score, evalue, sf, st, strand = hits[seq_id]
        clean_f.write(f">{header}\n{trim_to_hit(seq, sf, st)}\n")
        log_f.write(f"{seq_id}\t{score}\t{evalue}\t{sf}\t{st}\t{strand}\n")
        kept += 1
      else:
        flag_f.write(f">{header}\n{seq}\n")
        flagged += 1

  total = kept + flagged
  pct = f"{flagged / total * 100:.1f}" if total else "0.0"
  print(f"  Kept: {kept}  Flagged: {flagged} ({pct}% flagged)")

  if args.stats:
    with open(args.stats, 'a') as sf:
      sf.write(f"{args.gene}\t{args.domain}\t{total}\t{kept}\t{flagged}\n")


if __name__ == '__main__':
  main()
