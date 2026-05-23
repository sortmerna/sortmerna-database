#!/usr/bin/env python3
"""Extract rRNA feature coordinates from a GenBank flat file (GBFF) into BED format.

Parses GBFF line-by-line (no sequence loading) and writes one BED line per
coordinate range found in every rRNA feature, converting GenBank 1-based
coordinates to 0-based half-open BED coordinates.

Usage: python extract_rrna_loci.py input.gbff[.gz] output.bed
"""

import argparse
import gzip
import re
from pathlib import Path


def extract_rrna_loci(in_path, out_path, margin=0):
    """Parse GBFF and write rRNA loci to BED (0-based, half-open).

    Args:
        in_path: path to .gbff or .gbff.gz
        out_path: path to output BED file
        margin:   bp to add on each side of every locus (default 0);
                  start is clamped to 0, end is clamped to chromosome length

    Returns:
        int: number of loci written
    """
    if margin < 0:
        raise ValueError("margin must be >= 0")

    count = 0
    chrom = None
    chrom_len = None
    in_rrna = False
    loc_buf = ""

    def emit(loc):
        # Each coordinate pair in the location string becomes one BED entry.
        # join() and complement(join()) features produce multiple entries - one
        # per segment - which is intentional: we mask exactly the annotated
        # rRNA nucleotides and leave non-rRNA gaps (e.g. introns) unmasked.
        nonlocal count
        for s, e in re.findall(r"<? *(\d+)\.\. *>?(\d+)", loc):
            start = max(0, int(s) - 1 - margin)
            end   = int(e) + margin
            if chrom_len is not None:
                end = min(end, chrom_len)
            out.write(f"{chrom}\t{start}\t{end}\n")
            count += 1

    opener = gzip.open if Path(in_path).suffix == ".gz" else open
    with opener(in_path, "rt") as f, open(out_path, "w") as out:
        for line in f:
            if line.startswith("LOCUS "):
                if in_rrna and loc_buf:
                    emit(loc_buf)
                parts = line.split()
                chrom = parts[1]
                # LOCUS line: LOCUS <name> <length> bp ...
                try:
                    chrom_len = int(parts[2])
                except (IndexError, ValueError):
                    chrom_len = None
                in_rrna = False
                loc_buf = ""
            elif in_rrna:
                stripped = line.rstrip()
                # Location continuation: exactly 21 spaces, not a qualifier (no /)
                if len(stripped) > 21 and stripped[:21] == " " * 21 and stripped[21] != "/":
                    loc_buf += stripped[21:]
                else:
                    emit(loc_buf)
                    in_rrna = False
                    loc_buf = ""
            # New rRNA feature: 5-space indent + "rRNA"
            if not in_rrna and line.startswith("     rRNA"):
                parts = line[5:].split()
                in_rrna = True
                loc_buf = parts[1] if len(parts) > 1 else ""
        if in_rrna and loc_buf:
            emit(loc_buf)

    return count


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("gbff", help="Input .gbff or .gbff.gz file")
    ap.add_argument("bed", help="Output BED file")
    ap.add_argument("--margin", type=int, default=0,
                    help="Extend each locus by this many bp on each side (default: 0)")
    args = ap.parse_args()
    count = extract_rrna_loci(args.gbff, args.bed, margin=args.margin)
    print(f"  rRNA loci: {count} regions -> {args.bed}")


if __name__ == "__main__":
    main()
