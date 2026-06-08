#!/usr/bin/env python3
"""
compute_masking_stats.py

Count sequences with soft-masked (lowercase) bases in a FASTA file,
and optionally run RepeatMasker (-noint -low, simple/low-complexity only)
to produce a second independent count.

Outputs one TSV row per FASTA:
    subunit  domain  total_seqs  silva_masked  silva_pct  rm_masked  rm_pct

Usage:
    python3 compute_masking_stats.py --fasta FILE --subunit STR --domain STR \\
        --out TSV [--run-repeatmasker] [--repeatmasker-bin RepeatMasker] \\
        [--threads 4]
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


def count_repeatmasker_softmasked(fasta_path: str, rm_bin: str, threads: int = 4) -> int:
    """Run RepeatMasker (-noint -xsmall) to mask simple/low-complexity only and count soft-masked seqs.

    -noint  : only simple repeats and low-complexity (skips transposable elements)
    -xsmall : soft-mask output (lowercase), not hard-mask with N
    -norna  : don't mask small RNA genes (important - input IS rRNA sequences)
    -pa     : parallel batch jobs (not CPU threads)
    """
    import tempfile, shutil
    tmpdir = tempfile.mkdtemp()
    try:
        subprocess.run(
            [rm_bin, "-noint", "-xsmall", "-norna", "-pa", str(threads),
             "-dir", tmpdir, fasta_path],
            capture_output=True, check=True
        )
        # RepeatMasker writes <basename>.masked with soft-masked output
        masked_file = Path(tmpdir) / (Path(fasta_path).name + ".masked")
        if not masked_file.exists():
            # No repeats found - all sequences unmasked
            return 0
        _, rm_masked = count_silva_softmasked(str(masked_file))
        return rm_masked
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"WARNING: RepeatMasker failed: {e}", file=sys.stderr)
        return -1
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta",          required=True)
    ap.add_argument("--subunit",        required=True, help="e.g. SSU, LSU, 5S, 5.8S")
    ap.add_argument("--domain",         required=True, help="e.g. bacteria, archaea, eukaryota")
    ap.add_argument("--out",            required=True, help="Append-mode TSV output")
    ap.add_argument("--run-repeatmasker", action="store_true")
    ap.add_argument("--repeatmasker-bin", default="RepeatMasker")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    if not Path(args.fasta).exists() and not Path(args.fasta + ".gz").exists():
        print(f"WARNING: {args.fasta} not found - skipping", file=sys.stderr)
        return

    fasta = args.fasta if Path(args.fasta).exists() else args.fasta + ".gz"

    total, silva_masked = count_silva_softmasked(fasta)
    silva_pct = f"{silva_masked / total * 100:.2f}" if total else "0.00"

    if args.run_repeatmasker:
        rm_masked = count_repeatmasker_softmasked(fasta, args.repeatmasker_bin, args.threads)
        rm_pct = f"{rm_masked / total * 100:.2f}" if total and rm_masked >= 0 else "NA"
        rm_masked_str = str(rm_masked) if rm_masked >= 0 else "NA"
    else:
        rm_masked_str = "NA"
        rm_pct = "NA"

    write_header = not Path(args.out).exists()
    with open(args.out, "a") as f:
        if write_header:
            f.write("subunit\tdomain\ttotal_seqs\tsilva_masked\tsilva_pct\trm_masked\trm_pct\n")
        f.write(f"{args.subunit}\t{args.domain}\t{total}\t{silva_masked}\t{silva_pct}\t"
                f"{rm_masked_str}\t{rm_pct}\n")

    print(f"  {args.subunit} {args.domain}: {silva_masked}/{total} SILVA-masked "
          f"({silva_pct}%)" +
          (f", {rm_masked_str}/{total} RepeatMasker ({rm_pct}%)" if args.run_repeatmasker else ""))


if __name__ == "__main__":
    main()
