#!/usr/bin/env python3
"""Extract rRNA feature coordinates from a GFF3 annotation file into BED format.

Parses GFF3 line-by-line and writes one BED line per rRNA feature, converting
GFF3 1-based inclusive coordinates to 0-based half-open BED coordinates.

Usage: python extract_rrna_loci.py input.gff[.gz] output.bed [--name-map assembly_report.txt]
"""

import argparse
import gzip
from pathlib import Path


def load_name_map(report_path):
    """Build RefSeq-accn -> UCSC-name map from an NCBI assembly report.

    Reads tab-delimited rows (skipping # comment lines) and maps column 7
    (RefSeq-Accn, e.g. NC_060925.1) to column 10 (UCSC-style-name, e.g. chr1).
    """
    mapping = {}
    opener = gzip.open if Path(report_path).suffix == ".gz" else open
    with opener(report_path, "rt") as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 10:
                refseq = parts[6].strip()
                ucsc   = parts[9].strip()
                if refseq not in ('na', '') and ucsc not in ('na', ''):
                    mapping[refseq] = ucsc
    return mapping


def extract_rrna_loci(in_path, out_path, margin=0, name_map_path=None):
    """Parse GFF3 and write rRNA loci to BED (0-based, half-open).

    Args:
        in_path:       path to .gff or .gff.gz (GFF3 format)
        out_path:      path to output BED file
        margin:        bp to add on each side of every locus (default 0);
                       start is clamped to 0
        name_map_path: optional NCBI assembly report for RefSeq -> chr name mapping

    Returns:
        int: number of loci written
    """
    if margin < 0:
        raise ValueError("margin must be >= 0")

    name_map = load_name_map(name_map_path) if name_map_path is not None else {}

    count = 0
    opener = gzip.open if Path(in_path).suffix == ".gz" else open
    with opener(in_path, "rt") as f, open(out_path, "w") as out:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'rRNA':
                continue
            chrom = parts[0]
            chrom = name_map.get(chrom, chrom)
            start = max(0, int(parts[3]) - 1 - margin)
            end   = int(parts[4]) + margin
            out.write(f"{chrom}\t{start}\t{end}\n")
            count += 1

    return count


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("gff",      help="Input .gff or .gff.gz file (GFF3 format)")
    ap.add_argument("bed",      help="Output BED file")
    ap.add_argument("--margin", type=int, default=0,
                    help="Extend each locus by this many bp on each side (default: 0)")
    ap.add_argument("--name-map", dest="name_map",
                    help="NCBI assembly report for RefSeq accn -> chr name mapping")
    args = ap.parse_args()
    count = extract_rrna_loci(args.gff, args.bed, margin=args.margin,
                               name_map_path=args.name_map)
    print(f"  rRNA loci: {count} regions -> {args.bed}")


if __name__ == "__main__":
    main()
