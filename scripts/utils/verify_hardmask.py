#!/usr/bin/env python3
"""
verify_hardmask.py

Sanity check that a hard-masked FASTA matches a soft-masked FASTA with the sole
difference that every lowercase base in the soft-masked file is N in the hard-masked file.

Usage:
    python3 verify_hardmask.py --soft soft.fasta --hard hard_masked.fasta.gz
"""

import argparse
import gzip
import sys
from pathlib import Path


def iter_records(path: str):
    """Yield (id, sequence_str) from a FASTA file (plain or .gz), preserving case."""
    opener = gzip.open if Path(path).suffix == ".gz" else open
    header = None
    parts = []
    with opener(path, "rt") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header.split()[0][1:], "".join(parts)
                header = line
                parts = []
            elif line:
                parts.append(line)
    if header is not None:
        yield header.split()[0][1:], "".join(parts)


def check(soft_path: str, hard_path: str) -> int:
    errors = 0
    soft_iter = iter_records(soft_path)
    hard_iter = iter_records(hard_path)

    i = 0
    while True:
        s_rec = next(soft_iter, None)
        h_rec = next(hard_iter, None)
        if s_rec is None and h_rec is None:
            break
        if s_rec is None or h_rec is None:
            print(f"ERROR: sequence count mismatch (soft exhausted={s_rec is None}, hard exhausted={h_rec is None})")
            errors += 1
            break

        i += 1
        sid, ss = s_rec
        hid, hs = h_rec

        if sid != hid:
            print(f"ERROR seq {i}: id mismatch (soft='{sid}', hard='{hid}')")
            errors += 1
            continue

        if len(ss) != len(hs):
            print(f"ERROR >{sid}: length mismatch (soft={len(ss)}, hard={len(hs)})")
            errors += 1
            continue

        for pos, (s, h) in enumerate(zip(ss, hs)):
            if s.islower():
                if h != "N":
                    print(f"ERROR >{sid} pos {pos}: soft '{s}' should be 'N' in hard, got '{h}'")
                    errors += 1
            else:
                if s != h:
                    print(f"ERROR >{sid} pos {pos}: soft '{s}' != hard '{h}' (no masking expected)")
                    errors += 1

    return errors


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--soft", required=True, help="Soft-masked FASTA (plain or .gz)")
    ap.add_argument("--hard", required=True, help="Hard-masked FASTA (plain or .gz)")
    args = ap.parse_args()

    errors = check(args.soft, args.hard)
    if errors:
        print(f"\nFAILED: {errors} error(s) found.")
        sys.exit(1)
    else:
        print("OK: hard-masked file matches soft-masked file (lowercase -> N only).")


if __name__ == "__main__":
    main()
