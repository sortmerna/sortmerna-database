#!/usr/bin/env python3
"""
verify_hardmask.py

Sanity check that a hard-masked FASTA matches a soft-masked FASTA with the sole
difference that every lowercase base in the soft-masked file is N in the hard-masked file.

Usage:
    python3 verify_hardmask.py --soft soft.fasta --hard hard_masked.fasta.gz
"""

import argparse
import sys

import skbio


def iter_records(path: str):
    """Yield (id, sequence_str) from a FASTA file (plain or .gz)."""
    for seq in skbio.io.read(path, format="fasta", constructor=skbio.DNA, lowercase=True):
        yield str(seq.metadata["id"]), str(seq)


def check(soft_path: str, hard_path: str) -> int:
    errors = 0
    soft_iter = iter_records(soft_path)
    hard_iter = iter_records(hard_path)

    i = 0
    for (sid, ss), (hid, hs) in zip(soft_iter, hard_iter):
        i += 1
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

    # Check both iterators are exhausted
    soft_remaining = sum(1 for _ in soft_iter)
    hard_remaining = sum(1 for _ in hard_iter)
    if soft_remaining or hard_remaining:
        print(f"ERROR: sequence count mismatch (extra soft={soft_remaining}, extra hard={hard_remaining})")
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
