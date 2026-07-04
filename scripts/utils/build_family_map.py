#!/usr/bin/env python3
"""
build_family_map.py

Build a family_map.tsv mapping each sequence ID in a SortMeRNA database to its
rRNA family, source database version, sequence length and SILVA taxonomy levels.

Usage (called from build_sortmerna_index.sh):
    python3 build_family_map.py <output.tsv> <fasta1:rrna_family1> [<fasta2:rrna_family2> ...] \\
        [--silva-ssu-version V] [--silva-lsu-version V] [--rfam-version V]

Columns written:
    seq_id, rrna_family, database_version, seq_length, original,
    domain, phylum, class, order, tax_family, genus, species

The 7 taxonomy columns are populated only for SILVA sequences (rrna_family starting
with "silva"); for Rfam sequences they are left empty, since Rfam headers are
free-text descriptions rather than a ';'-delimited SILVA lineage.
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

HEADER = ["seq_id", "rrna_family", "database_version", "seq_length", "original",
          "domain", "phylum", "class", "order", "tax_family", "genus", "species"]
N_TAX_LEVELS = 7  # domain, phylum, class, order, tax_family, genus, species


def parse_header(line: str, parse_tax: bool = True) -> tuple[str, str, list[str]]:
    """Parse a FASTA header line into (seq_id, original, taxonomy_levels).

    Handles SILVA format: '>ID taxonomy;levels;here;size=N'
    and bare headers:      '>ID'

    'original' is the full header line (without the leading '>') preserved as-is.
    When parse_tax is False the taxonomy levels are all left empty - Rfam headers
    are free-text descriptions, not a ';'-delimited SILVA lineage.
    """
    line = line.lstrip(">").rstrip()
    original = line
    parts = line.split(" ", 1)
    seq_id = parts[0]
    if not parse_tax:
        return seq_id, original, [""] * N_TAX_LEVELS
    tax_str = parts[1] if len(parts) > 1 else ""
    tax_str = re.sub(r";size=\d+$", "", tax_str)
    levels = tax_str.split(";") if tax_str else []
    levels += [""] * (N_TAX_LEVELS - len(levels))
    return seq_id, original, levels[:N_TAX_LEVELS]


def resolve_db_version(rna_family: str, silva_ssu_version: str,
                       silva_lsu_version: str, rfam_version: str) -> str:
    """Return the source database version string for an rRNA family label,
    e.g. 'SILVA 138.2' or 'Rfam 15.1'. Empty when the version is unknown."""
    fam = rna_family.lower()
    if fam.startswith("silva"):
        ver = silva_lsu_version if "lsu" in fam else silva_ssu_version
        return f"SILVA {ver}" if ver else ""
    if fam.startswith("rfam"):
        return f"Rfam {rfam_version}" if rfam_version else ""
    return ""


def write_family_map(output_path: str, fasta_family_pairs: list[tuple[str, str]],
                     silva_ssu_version: str = "", silva_lsu_version: str = "",
                     rfam_version: str = "") -> int:
    """Write family_map.tsv from a list of (fasta_path, rna_family) pairs.

    Returns total number of sequences written.
    """
    total = 0
    with open(output_path, "w") as out:
        out.write("\t".join(HEADER) + "\n")
        for fasta_path, rna_family in fasta_family_pairs:
            is_silva = rna_family.lower().startswith("silva")
            db_version = resolve_db_version(rna_family, silva_ssu_version,
                                            silva_lsu_version, rfam_version)
            opener = gzip.open(fasta_path, "rt") if fasta_path.endswith(".gz") else open(fasta_path)
            with opener as f:
                seq_id = original = None
                levels = []
                seq_len = 0

                def emit():
                    out.write("\t".join(
                        [seq_id, rna_family, db_version, str(seq_len), original] + levels
                    ) + "\n")

                for line in f:
                    if line.startswith(">"):
                        if seq_id is not None:
                            emit()
                            total += 1
                        seq_id, original, levels = parse_header(line, parse_tax=is_silva)
                        seq_len = 0
                    else:
                        seq_len += len(line.strip())
                if seq_id is not None:
                    emit()
                    total += 1
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("output", help="Output family_map.tsv path")
    ap.add_argument("pairs", nargs="+", metavar="fasta:rrna_family",
                    help="FASTA path and rRNA family label separated by ':'")
    ap.add_argument("--silva-ssu-version", default="", help="SILVA SSU version (e.g. 138.2)")
    ap.add_argument("--silva-lsu-version", default="", help="SILVA LSU version (e.g. 138.2)")
    ap.add_argument("--rfam-version", default="", help="Rfam version (e.g. 15.1)")
    args = ap.parse_args()

    pairs = []
    for entry in args.pairs:
        fasta_path, rna_family = entry.rsplit(":", 1)
        p = Path(fasta_path)
        if not p.exists() and Path(fasta_path + ".gz").exists():
            fasta_path += ".gz"
        elif not p.exists():
            print(f"WARNING: {fasta_path}[.gz] not found - skipping", file=sys.stderr)
            continue
        pairs.append((fasta_path, rna_family))

    total = write_family_map(args.output, pairs,
                             silva_ssu_version=args.silva_ssu_version,
                             silva_lsu_version=args.silva_lsu_version,
                             rfam_version=args.rfam_version)
    print(f"  family_map.tsv: {total:,} sequences -> {args.output}")


if __name__ == "__main__":
    main()
