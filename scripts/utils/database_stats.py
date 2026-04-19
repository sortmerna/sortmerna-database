#!/usr/bin/env python3
"""
database_stats.py - Generate statistics for rRNA databases

Analyzes FASTA databases and generates comprehensive statistics including
sequence counts, length distributions, GC content, and taxonomic coverage.
"""

import argparse
import json
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import re


def parse_fasta(filepath: str) -> List[Tuple[str, str]]:
  """Parse a FASTA file and return list of (header, sequence) tuples."""
  sequences = []
  current_header = None
  current_seq = []

  with open(filepath, 'r') as f:
    for line in f:
      line = line.strip()
      if line.startswith('>'):
        if current_header is not None:
          sequences.append((current_header, ''.join(current_seq)))
        current_header = line[1:]
        current_seq = []
      else:
        current_seq.append(line.upper())

    if current_header is not None:
      sequences.append((current_header, ''.join(current_seq)))

  return sequences


def calculate_gc_content(sequence: str) -> float:
  """Calculate GC content of a sequence."""
  sequence = sequence.upper()
  gc_count = sequence.count('G') + sequence.count('C')
  total = len(sequence)
  if total == 0:
    return 0.0
  return (gc_count / total) * 100


def calculate_length_stats(lengths: List[int]) -> Dict:
  """Calculate length statistics."""
  if not lengths:
    return {'min': 0, 'max': 0, 'mean': 0, 'median': 0, 'total_bp': 0}

  sorted_lengths = sorted(lengths)
  n = len(sorted_lengths)

  return {
    'min': sorted_lengths[0],
    'max': sorted_lengths[-1],
    'mean': round(sum(sorted_lengths) / n, 2),
    'median': sorted_lengths[n // 2] if n % 2 == 1
              else (sorted_lengths[n // 2 - 1] + sorted_lengths[n // 2]) / 2,
    'total_bp': sum(sorted_lengths)
  }


def extract_taxonomy(header: str) -> Optional[str]:
  """Extract taxonomy from SILVA-style header."""
  # SILVA format: >ACCESSION.start.stop taxonomy
  # Example: >AB000001.1.1234 Bacteria;Proteobacteria;...
  parts = header.split(' ', 1)
  if len(parts) > 1:
    return parts[1]
  return None


def analyze_taxonomy(sequences: List[Tuple[str, str]]) -> Dict:
  """Analyze taxonomic distribution of sequences."""
  domain_counts = defaultdict(int)
  phylum_counts = defaultdict(int)

  for header, _ in sequences:
    taxonomy = extract_taxonomy(header)
    if taxonomy:
      levels = taxonomy.split(';')
      if levels:
        domain = levels[0].strip()
        domain_counts[domain] += 1
        if len(levels) > 1:
          phylum = f"{domain};{levels[1].strip()}"
          phylum_counts[phylum] += 1

  return {
    'domains': dict(domain_counts),
    'top_phyla': dict(sorted(phylum_counts.items(),
                             key=lambda x: x[1],
                             reverse=True)[:20])
  }


def calculate_length_distribution(lengths: List[int], bins: int = 10) -> Dict:
  """Calculate length distribution histogram."""
  if not lengths:
    return {'bins': [], 'counts': []}

  min_len = min(lengths)
  max_len = max(lengths)
  bin_width = (max_len - min_len) / bins if max_len > min_len else 1

  bin_edges = [min_len + i * bin_width for i in range(bins + 1)]
  counts = [0] * bins

  for length in lengths:
    bin_idx = min(int((length - min_len) / bin_width), bins - 1)
    counts[bin_idx] += 1

  return {
    'bin_edges': [round(e, 0) for e in bin_edges],
    'counts': counts
  }


def analyze_database(filepath: str) -> Dict:
  """Analyze a FASTA database and return statistics."""
  print(f"Analyzing: {filepath}")

  sequences = parse_fasta(filepath)

  if not sequences:
    return {
      'file': filepath,
      'error': 'No sequences found',
      'sequence_count': 0
    }

  lengths = [len(seq) for _, seq in sequences]
  gc_contents = [calculate_gc_content(seq) for _, seq in sequences]
  taxonomy = analyze_taxonomy(sequences)
  file_size = os.path.getsize(filepath)

  return {
    'file': filepath,
    'filename': os.path.basename(filepath),
    'file_size_bytes': file_size,
    'file_size_mb': round(file_size / (1024 * 1024), 2),
    'sequence_count': len(sequences),
    'length_stats': calculate_length_stats(lengths),
    'gc_content': {
      'min': round(min(gc_contents), 2) if gc_contents else 0,
      'max': round(max(gc_contents), 2) if gc_contents else 0,
      'mean': round(sum(gc_contents) / len(gc_contents), 2) if gc_contents else 0
    },
    'length_distribution': calculate_length_distribution(lengths),
    'taxonomy': taxonomy
  }


def format_number(n: int) -> str:
  """Format number with thousands separator."""
  return f"{n:,}"


def print_stats(stats: Dict) -> None:
  """Print statistics in a human-readable format."""
  print()
  print("=" * 60)
  print(f"Database: {stats['filename']}")
  print("=" * 60)
  print()
  print(f"File size: {stats['file_size_mb']} MB")
  print(f"Sequence count: {format_number(stats['sequence_count'])}")
  print()

  length_stats = stats['length_stats']
  print("Sequence Length Statistics:")
  print(f"  Min:    {format_number(length_stats['min'])} bp")
  print(f"  Max:    {format_number(length_stats['max'])} bp")
  print(f"  Mean:   {format_number(int(length_stats['mean']))} bp")
  print(f"  Median: {format_number(int(length_stats['median']))} bp")
  print(f"  Total:  {format_number(length_stats['total_bp'])} bp")
  print()

  gc = stats['gc_content']
  print("GC Content:")
  print(f"  Min:  {gc['min']}%")
  print(f"  Max:  {gc['max']}%")
  print(f"  Mean: {gc['mean']}%")
  print()

  taxonomy = stats['taxonomy']
  if taxonomy['domains']:
    print("Taxonomic Distribution (Domains):")
    for domain, count in sorted(taxonomy['domains'].items(),
                                key=lambda x: x[1], reverse=True):
      pct = (count / stats['sequence_count']) * 100
      print(f"  {domain}: {format_number(count)} ({pct:.1f}%)")
    print()


def compare_clustering_levels(stats_list: List[Dict]) -> None:
  """Compare statistics across clustering levels."""
  print()
  print("=" * 80)
  print("Clustering Level Comparison")
  print("=" * 80)
  print()

  sorted_stats = sorted(stats_list, key=lambda x: x['sequence_count'], reverse=True)

  if not sorted_stats:
    print("No statistics to compare.")
    return

  original_count = sorted_stats[0]['sequence_count']

  print(f"{'Database':<40} {'Sequences':>12} {'Reduction':>10} {'Size (MB)':>10}")
  print("-" * 80)

  for stats in sorted_stats:
    reduction = ((original_count - stats['sequence_count']) / original_count) * 100
    print(f"{stats['filename']:<40} "
          f"{format_number(stats['sequence_count']):>12} "
          f"{reduction:>9.1f}% "
          f"{stats['file_size_mb']:>10.2f}")


def main():
  parser = argparse.ArgumentParser(description='Generate statistics for rRNA databases')
  parser.add_argument('input', nargs='+',
                      help='Input FASTA file(s) or directory containing FASTA files')
  parser.add_argument('-o', '--output', help='Output JSON file for statistics')
  parser.add_argument('--compare', action='store_true',
                      help='Compare statistics across input files (for clustering analysis)')
  parser.add_argument('-q', '--quiet', action='store_true',
                      help='Suppress detailed output, only write JSON')
  args = parser.parse_args()

  fasta_files = []
  for path in args.input:
    if os.path.isdir(path):
      for ext in ['*.fasta', '*.fa', '*.fna']:
        fasta_files.extend(Path(path).glob(ext))
    elif os.path.isfile(path):
      fasta_files.append(Path(path))
    else:
      print(f"Warning: Path not found: {path}", file=sys.stderr)

  if not fasta_files:
    print("Error: No FASTA files found.", file=sys.stderr)
    sys.exit(1)

  print(f"Found {len(fasta_files)} FASTA file(s) to analyze")

  all_stats = []
  for fasta in sorted(fasta_files):
    stats = analyze_database(str(fasta))
    all_stats.append(stats)
    if not args.quiet:
      print_stats(stats)

  if args.compare and len(all_stats) > 1:
    compare_clustering_levels(all_stats)

  if args.output:
    output_data = {
      'databases': all_stats,
      'summary': {
        'total_files': len(all_stats),
        'total_sequences': sum(s['sequence_count'] for s in all_stats),
        'total_size_mb': sum(s['file_size_mb'] for s in all_stats)
      }
    }
    with open(args.output, 'w') as f:
      json.dump(output_data, f, indent=2)
    print(f"\nStatistics written to: {args.output}")

  print("\nAnalysis complete!")


if __name__ == '__main__':
  main()
