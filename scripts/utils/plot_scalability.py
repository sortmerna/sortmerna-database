#!/usr/bin/env python3
"""Generate scalability plots from SortMeRNA runs at multiple read volumes.

Parses aligned.log and aligned.blast from each scale-point directory produced
by run_scalability.sh, then writes three figures to --output-dir:

  <label>_summary.png      % aligned, runtime, S_min vs. read count
  <label>_evalue_dist.png  E-value distributions (one panel per scale point)
  <label>_identity_dist.png  % identity distributions (overlaid, all scale points)

Usage:
  python plot_scalability.py --output-dir DIR --scale-dirs DIR [DIR ...] --label STR
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_log(log_path):
    """Return dict of key metrics from an aligned.log file."""
    text = Path(log_path).read_text()

    def _int(pattern):
        m = re.search(pattern, text)
        return int(m.group(1)) if m else None

    def _float(pattern):
        m = re.search(pattern, text)
        return float(m.group(1)) if m else None

    return {
        'total_reads': _int(r'Total reads = (\d+)'),
        'aligned':     _int(r'Total reads passing E-value threshold = (\d+)'),
        'min_len':     _int(r'Minimum read length = (\d+)'),
        'max_len':     _int(r'Maximum read length = (\d+)'),
        'mean_len':    _int(r'Mean read length\s+=\s+(\d+)'),
        's_min':       _float(r'Minimal SW score based on E-value = ([\d.]+)'),
        'lambda_':     _float(r'Gumbel lambda = ([\d.eE+\-]+)'),
        'gumbel_k':    _float(r'Gumbel K = ([\d.eE+\-]+)'),
    }


def parse_blast(blast_path):
    """Return DataFrame with pident and evalue columns from a BLAST tabular file."""
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    chunks = []
    try:
        for chunk in pd.read_csv(
            blast_path, sep='\t', header=None, names=cols,
            chunksize=500_000, usecols=['pident', 'evalue']
        ):
            chunks.append(chunk)
    except Exception as exc:
        print(f'  WARNING: could not parse {blast_path}: {exc}', file=sys.stderr)
        return pd.DataFrame(columns=['pident', 'evalue'])
    return pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=['pident', 'evalue'])


def write_top_hits(blast_path, out_path, top_n=20):
    """Count aligned reference hits and write the top_n most frequent to out_path.

    Reads only the sseqid column to keep memory use low for large BLAST files.
    """
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    counts: dict = {}
    try:
        for chunk in pd.read_csv(
            blast_path, sep='\t', header=None, names=cols,
            chunksize=500_000, usecols=['sseqid']
        ):
            for ref, cnt in chunk['sseqid'].value_counts().items():
                counts[ref] = counts.get(ref, 0) + cnt
    except Exception as exc:
        print(f'  WARNING: could not parse {blast_path}: {exc}', file=sys.stderr)
        return

    if not counts:
        return

    top = sorted(counts.items(), key=lambda x: x[1], reverse=True)[:top_n]
    with open(out_path, 'w') as f:
        f.write(f'{"Count":>8}  Reference\n')
        f.write('-' * 60 + '\n')
        for ref, cnt in top:
            f.write(f'{cnt:>8}  {ref}\n')
    print(f'  Saved: {out_path}')


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

_BLUE   = '#2c7bb6'
_ORANGE = '#d7191c'
_GREEN  = '#1a9641'


def _plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _np():
    import numpy as np
    return np


def plot_summary(stats, label, out_dir):
    """3-panel figure: % aligned, runtime, S_min vs. read count."""
    plt = _plt()
    ns       = [s['n'] for s in stats]
    pct      = [s['log']['aligned'] / s['log']['total_reads'] * 100 for s in stats]
    runtimes = [s['runtime'] for s in stats]
    s_mins   = [s['log']['s_min'] for s in stats]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.suptitle(label, fontsize=12, y=1.01)

    ax = axes[0]
    ax.plot(ns, pct, 'o-', color=_BLUE)
    ax.set_xscale('log')
    ax.set_xlabel('Read count')
    ax.set_ylabel('Reads aligned (%)')
    ax.set_title('Aligned rate vs. scale')
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(ns, runtimes, 'o-', color=_ORANGE)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Read count')
    ax.set_ylabel('Runtime (s)')
    ax.set_title('Runtime vs. scale (log-log)')
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    ax.plot(ns, s_mins, 'o-', color=_GREEN)
    ax.set_xscale('log')
    ax.set_xlabel('Read count')
    ax.set_ylabel('S_min (SW score units)')
    ax.set_title('Score threshold vs. scale')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out_path = out_dir / f'{label}_summary.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {out_path}')


def plot_evalue_dist(stats, label, out_dir):
    """E-value distribution for aligned reads, one panel per scale point."""
    plt = _plt()
    np  = _np()
    n_panels = len(stats)
    fig, axes = plt.subplots(1, n_panels, figsize=(4 * n_panels, 4))
    if n_panels == 1:
        axes = [axes]

    for ax, s in zip(axes, stats):
        df = s['blast']
        if df.empty:
            ax.text(0.5, 0.5, 'no alignments', ha='center', va='center',
                    transform=ax.transAxes)
            ax.set_title(f'N = {s["n"]:,}')
            continue
        log_ev = np.log10(df['evalue'].clip(lower=1e-300))
        ax.hist(log_ev, bins=50, color=_BLUE, edgecolor='none', alpha=0.85)
        ax.set_xlabel('log$_{10}$(E-value)')
        ax.set_ylabel('Count')
        ax.set_title(f'N = {s["n"]:,}  ({len(df):,} aligned)')
        ax.grid(True, alpha=0.3)

    fig.suptitle(f'{label} - E-value distribution (aligned reads)', fontsize=11)
    plt.tight_layout()
    out_path = out_dir / f'{label}_evalue_dist.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {out_path}')


def plot_identity_dist(stats, label, out_dir):
    """% identity distribution for aligned reads, all scale points overlaid."""
    plt = _plt()
    np  = _np()
    fig, ax = plt.subplots(figsize=(8, 4))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(stats)))

    for s, color in zip(stats, colors):
        df = s['blast']
        if df.empty:
            continue
        ax.hist(df['pident'], bins=50, alpha=0.55, color=color,
                label=f'N={s["n"]:,}', edgecolor='none', density=True)

    ax.set_xlabel('% Identity')
    ax.set_ylabel('Density')
    ax.set_title(f'{label} - Identity distribution (aligned reads)')
    ax.legend(fontsize=8, title='Scale point')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out_path = out_dir / f'{label}_identity_dist.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {out_path}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--output-dir', required=True, type=Path,
                    help='Directory for output PNG files')
    ap.add_argument('--scale-dirs', nargs='+', required=True, type=Path,
                    help='Scale-point directories (each produced by run_scalability.sh)')
    ap.add_argument('--label', default='sortmerna_scalability',
                    help='Label for plot titles and output filenames')
    ap.add_argument('--n-reads', nargs='+', type=int,
                    help='Read counts matching --scale-dirs order; inferred from '
                         'directory name (scale_<N>) if omitted')
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    stats = []
    for i, d in enumerate(args.scale_dirs):
        log_path   = d / 'smr_out' / 'out' / 'aligned.log'
        blast_path = d / 'smr_out' / 'out' / 'aligned.blast'
        rt_path    = d / 'runtime_seconds.txt'

        if not log_path.exists():
            print(f'  WARNING: missing {log_path}, skipping', file=sys.stderr)
            continue

        if args.n_reads:
            n = args.n_reads[i]
        else:
            m = re.search(r'scale_(\d+)', d.name)
            if not m:
                print(f'  WARNING: cannot infer read count from dir name "{d.name}", '
                      'pass --n-reads explicitly', file=sys.stderr)
                continue
            n = int(m.group(1))

        runtime = int(rt_path.read_text().strip()) if rt_path.exists() else 0
        log     = parse_log(log_path)
        blast   = parse_blast(blast_path) if blast_path.exists() else \
                  pd.DataFrame(columns=['pident', 'evalue'])
        if blast_path.exists():
            write_top_hits(blast_path,
                           args.output_dir / f'{args.label}_top_hits_{n}.txt')

        aligned = log['aligned'] or 0
        total   = log['total_reads'] or 1
        print(f'  {n:>12,} reads: {aligned:,} aligned '
              f'({aligned / total * 100:.2f}%), {runtime}s')
        stats.append({'n': n, 'log': log, 'runtime': runtime, 'blast': blast})

    if not stats:
        print('No valid scale directories found.', file=sys.stderr)
        sys.exit(1)

    stats.sort(key=lambda s: s['n'])

    plot_summary(stats, args.label, args.output_dir)
    plot_evalue_dist(stats, args.label, args.output_dir)
    plot_identity_dist(stats, args.label, args.output_dir)

    print(f'\nPlots written to: {args.output_dir}')


if __name__ == '__main__':
    main()
