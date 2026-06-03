#!/usr/bin/env python3
"""Plot ROC curve sweeping E-value thresholds across SortMeRNA runs.

For each E-value threshold, reads aligned.log from the largest scale-point
directory inside the corresponding --rrna-dirs and --nonrrna-dirs entries,
computes TPR (sensitivity) and FPR (1 - specificity), and plots one (FPR, TPR)
point per E-value threshold.

Multiple series (e.g. v5.0.1 vs v6.0.1, or different databases) can be
overlaid on the same plot by repeating --rrna-dirs, --nonrrna-dirs, and
--series-labels flags in the same order.

Usage:
  python plot_roc_evalue.py \\
      --output-dir plots/ \\
      --evalues 1 0.1 0.05 0.01 \\
      --rrna-dirs    scalability_rrna_ev1 scalability_rrna_ev0.1 scalability_rrna_ev0.05 scalability_rrna_ev0.01 \\
      --nonrrna-dirs scalability_t2t_ev1  scalability_t2t_ev0.1  scalability_t2t_ev0.05  scalability_t2t_ev0.01 \\
      --series-labels "v6.0.1 default" \\
      --rrna-dirs    scalability_rrna_v5_ev1 ... \\
      --nonrrna-dirs scalability_t2t_v5_ev1 ... \\
      --series-labels "v5.0.1 default"
"""

import argparse
import base64
import io
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

_COLORS = ['#2c7bb6', '#d7191c', '#1a9641', '#7b2d8b', '#f07d00']

_FAMILY_ORDER = [
    'silva_ssu_bacteria', 'silva_ssu_archaea', 'silva_ssu_eukaryota',
    'silva_lsu_bacteria', 'silva_lsu_archaea', 'silva_lsu_eukaryota',
    'rfam_5_8s', 'rfam_5s',
]


def find_largest_scale_dir(top_dir):
    """Return the scale_N subdirectory with the largest N, or None."""
    top = Path(top_dir)
    candidates = []
    for d in top.iterdir():
        m = re.match(r'scale_(\d+)$', d.name)
        if m and d.is_dir():
            candidates.append((int(m.group(1)), d))
    if not candidates:
        return None
    return max(candidates, key=lambda x: x[0])[1]


def read_counts(scale_dir):
    """Return (aligned, total) from aligned.log inside a scale directory."""
    log_path = Path(scale_dir) / 'smr_out' / 'out' / 'aligned.log'
    if not log_path.exists():
        print(f'  WARNING: missing {log_path}', file=sys.stderr)
        return None, None
    text = log_path.read_text()
    m_aligned = re.search(r'Total reads passing E-value threshold = (\d+)', text)
    m_total   = re.search(r'Total reads = (\d+)', text)
    if not m_aligned or not m_total:
        print(f'  WARNING: could not parse counts from {log_path}', file=sys.stderr)
        return None, None
    return int(m_aligned.group(1)), int(m_total.group(1))


def compute_series(rrna_dirs, nonrrna_dirs, evalues):
    """Return list of (evalue, fpr, tpr) for one series, sorted by evalue descending."""
    points = []
    for ev, rd, nd in zip(evalues, rrna_dirs, nonrrna_dirs):
        r_scale  = find_largest_scale_dir(rd)
        nr_scale = find_largest_scale_dir(nd)
        if r_scale is None or nr_scale is None:
            print(f'  WARNING: no scale dir found in {rd} or {nd}, skipping e={ev}',
                  file=sys.stderr)
            continue
        r_aligned,  r_total  = read_counts(r_scale)
        nr_aligned, nr_total = read_counts(nr_scale)
        if None in (r_aligned, r_total, nr_aligned, nr_total):
            continue
        tpr = r_aligned  / r_total
        fpr = nr_aligned / nr_total
        n   = r_total
        print(f'  E={ev}: TPR={tpr:.4f} ({r_aligned}/{r_total} rRNA)  '
              f'FPR={fpr:.6f} ({nr_aligned}/{nr_total} non-rRNA, N={n:,})')
        points.append((float(ev), fpr, tpr))
    # sort highest evalue first (loosest threshold = top-right of ROC)
    points.sort(key=lambda x: x[0], reverse=True)
    return points


def plot_roc(all_series, output_dir, label):
    """Plot ROC curve with one point per E-value, one line per series."""
    fig, ax = plt.subplots(figsize=(7, 6))

    # Label offsets cycle per series so annotations from different series don't overlap
    _label_offsets = [(6, 3), (6, -10), (6, 16), (6, -17)]

    all_fprs, all_tprs = [], []
    for si, ((series_label, points), color) in enumerate(zip(all_series, _COLORS)):
        if not points:
            continue
        fprs = [p[1] for p in points]
        tprs = [p[2] for p in points]
        evalues = [p[0] for p in points]
        all_fprs.extend(fprs)
        all_tprs.extend(tprs)
        ax.plot(fprs, tprs, 'o-', color=color, label=series_label, linewidth=1.5)
        # Only label the loosest and strictest E-value to avoid overlap
        label_idx = {0, len(points) - 1}
        base_xy = _label_offsets[si % len(_label_offsets)]
        for i, (fpr, tpr, ev) in enumerate(zip(fprs, tprs, evalues)):
            if i in label_idx:
                ax.annotate(f'E={ev:g}', (fpr, tpr),
                            textcoords='offset points', xytext=base_xy, fontsize=8,
                            color=color)

    x_max = max(all_fprs) * 1.3 if all_fprs else 0.05
    y_min = max(0.0, min(all_tprs) - 0.02) if all_tprs else 0.0
    ax.set_xlim(-x_max * 0.02, x_max)
    ax.set_ylim(y_min, 1.01)
    ax.set_xlabel('False positive rate  (non-rRNA reads classified as rRNA)')
    ax.set_ylabel('Sensitivity / TPR  (rRNA reads correctly classified)')
    ax.set_title(f'{label} - ROC by E-value threshold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    plt.close()
    out_path = Path(output_dir) / f'{label}_roc_evalue.png'
    out_path.write_bytes(base64.b64decode(b64))
    print(f'\n  Saved: {out_path}')
    return b64


def load_family_map(tsv_path):
    """Load source_seq_id -> family_name from rRNA_test_10M_family.tsv."""
    mapping = {}
    with open(tsv_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if '\t' in line:
                src_id, fam = line.split('\t', 1)
                mapping[src_id] = fam
    return mapping


def resolve_family(read_id, family_map):
    """Map a read_id to its rRNA family.

    Rfam direct reads keep their source_seq_id as-is.
    SILVA ISS reads are named {source_seq_id}_{read_number}/{mate}.
    """
    if read_id in family_map:
        return family_map[read_id]
    m = re.match(r'^(.+)_\d+_\d+/[12]$', read_id)
    if m and m.group(1) in family_map:
        return family_map[m.group(1)]
    return 'unknown'


def _ids_from_fasta(path):
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                yield line[1:].split()[0]


def _ids_from_blast(path):
    with open(path) as f:
        for line in f:
            if line.strip():
                yield line.split('\t', 1)[0]


def _build_family_tables(rrna_dirs_list, evalues, family_map):
    """Build per-(ev, scale_N) family counts. Uses first rrna_dirs for totals."""
    first_dir = Path(rrna_dirs_list[0][0])

    # Compute totals once from the first E-value dir's reads files
    family_totals = {}  # scale_N -> {family: int}
    for d in sorted(first_dir.iterdir()):
        m = re.match(r'scale_(\d+)$', d.name)
        if not m or not d.is_dir():
            continue
        n = int(m.group(1))
        reads_fa = d / f'reads_{n}.fasta'
        if not reads_fa.exists():
            continue
        counts = defaultdict(int)
        for rid in _ids_from_fasta(reads_fa):
            counts[resolve_family(rid, family_map)] += 1
        family_totals[n] = dict(counts)
        print(f'    Totals at scale {n:,}: {sum(counts.values()):,} reads across '
              f'{len(counts)} families')

    # For each E-value compute assigned from blast
    result = {}  # ev -> {scale_N: {family: {total, assigned}}}
    for ev, top_dir in zip(evalues, rrna_dirs_list[0]):
        top = Path(top_dir)
        result[ev] = {}
        for d in sorted(top.iterdir()):
            m = re.match(r'scale_(\d+)$', d.name)
            if not m or not d.is_dir():
                continue
            n = int(m.group(1))
            blast_f = d / 'smr_out' / 'out' / 'aligned.blast'
            totals = family_totals.get(n, {})
            if not totals and not blast_f.exists():
                continue
            assigned = defaultdict(int)
            if blast_f.exists():
                for rid in _ids_from_blast(blast_f):
                    assigned[resolve_family(rid, family_map)] += 1
            all_families = set(totals) | set(assigned)
            result[ev][n] = {
                fam: {'total': totals.get(fam, 0), 'assigned': int(assigned.get(fam, 0))}
                for fam in all_families
            }
    return result


def _read_perf_data(rrna_dirs_list, evalues):
    """Read runtime_seconds.txt and peak_rss_mb.txt for each (series, evalue, scale_N).

    Returns dict: series_idx -> evalue -> scale_N -> (runtime_s, ram_mb)
    """
    result = {}
    for si, rrna_dirs in enumerate(rrna_dirs_list):
        result[si] = {}
        for ev, top_dir in zip(evalues, rrna_dirs):
            result[si][ev] = {}
            for d in sorted(Path(top_dir).iterdir()):
                m = re.match(r'scale_(\d+)$', d.name)
                if not m or not d.is_dir():
                    continue
                n = int(m.group(1))
                rt_f = d / 'runtime_seconds.txt'
                ram_f = d / 'peak_rss_mb.txt'
                if rt_f.exists() and ram_f.exists():
                    result[si][ev][n] = (
                        float(rt_f.read_text().strip()),
                        float(ram_f.read_text().strip()),
                    )
    return result


def _make_perf_plots(perf_data, evalues, series_labels):
    """Create runtime and RAM charts. Returns (runtime_b64, ram_b64) - either may be None."""
    plots = []
    for metric_idx, ylabel, title in [
        (0, 'Runtime (seconds)', 'Runtime vs Number of Reads'),
        (1, 'Peak RAM (MB)', 'Peak RAM vs Number of Reads'),
    ]:
        fig, ax = plt.subplots(figsize=(8, 5))
        has_data = False
        color_idx = 0
        for si, slabel in enumerate(series_labels):
            ev_data = perf_data.get(si, {})
            for ev in sorted(evalues, reverse=True):
                scale_data = ev_data.get(ev, {})
                if not scale_data:
                    continue
                ns = sorted(scale_data)
                vals = [scale_data[n][metric_idx] for n in ns]
                color = _COLORS[color_idx % len(_COLORS)]
                color_idx += 1
                lbl = f'{slabel} E={ev:g}' if len(series_labels) > 1 else f'E={ev:g}'
                ax.plot(ns, vals, 'o-', label=lbl, color=color, linewidth=1.5)
                has_data = True
        if not has_data:
            plt.close()
            plots.append(None)
            continue
        ax.set_xscale('log')
        ax.set_xlabel('Number of reads')
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
        buf.seek(0)
        plots.append(base64.b64encode(buf.read()).decode())
        plt.close()
    return plots[0], plots[1]


def render_family_html(tables, evalues, label, output_dir,
                       roc_b64=None, runtime_b64=None, ram_b64=None):
    """Write <label>_family_breakdown.html with per-family sensitivity per E-value."""
    sections = []
    for ev in sorted(evalues, reverse=True):
        ev_data = tables.get(ev, {})
        scale_blocks = []
        for n in sorted(ev_data):
            fam_data = ev_data[n]
            ordered = [f for f in _FAMILY_ORDER if f in fam_data]
            ordered += sorted(f for f in fam_data if f not in _FAMILY_ORDER)
            rows = []
            grand_total = grand_assigned = 0
            for fam in ordered:
                d = fam_data[fam]
                t, a = d['total'], d['assigned']
                grand_total += t
                grand_assigned += a
                pct = f'{a / t * 100:.5f}' if t else '-'
                rows.append(
                    f'      <tr><td>{fam.replace("_", " ")}</td>'
                    f'<td>{t:,}</td><td>{a:,}</td><td>{pct}%</td></tr>'
                )
            if not rows:
                continue
            gp = f'{grand_assigned / grand_total * 100:.5f}' if grand_total else '-'
            rows.append(
                f'      <tr style="font-weight:bold;border-top:2px solid #2c3e50">'
                f'<td>Total</td><td>{grand_total:,}</td>'
                f'<td>{grand_assigned:,}</td><td>{gp}%</td></tr>'
            )
            scale_blocks.append(
                f'    <h3>{n:,} reads</h3>\n'
                f'    <div class="table-wrap"><table>\n'
                f'      <thead><tr><th>rRNA type</th><th>Total reads</th>'
                f'<th>Reads assigned to rRNA</th><th>%</th></tr></thead>\n'
                f'      <tbody>\n' + '\n'.join(rows) + '\n      </tbody>\n    </table></div>'
            )
        if scale_blocks:
            sections.append(
                f'  <section>\n  <h2>E-value = {ev:g}</h2>\n'
                + '\n'.join(scale_blocks)
                + '\n  </section>'
            )

    def _img_tag(b64):
        return f'<img src="data:image/png;base64,{b64}" style="max-width:100%;height:auto">'

    plots_section = ''
    if roc_b64 or runtime_b64 or ram_b64:
        parts = ['<section>']
        if roc_b64:
            parts.append(f'<h2>ROC Curve</h2>\n{_img_tag(roc_b64)}')
        if runtime_b64 or ram_b64:
            parts.append('<h2>Performance Summary</h2>')
            parts.append('<div style="display:flex;gap:1em;flex-wrap:wrap">')
            if runtime_b64:
                parts.append(f'<div style="flex:1;min-width:300px">{_img_tag(runtime_b64)}</div>')
            if ram_b64:
                parts.append(f'<div style="flex:1;min-width:300px">{_img_tag(ram_b64)}</div>')
            parts.append('</div>')
        parts.append('</section>')
        plots_section = '\n'.join(parts)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>rRNA Family Breakdown - {label}</title>
  <style>
    body {{ font-family: sans-serif; padding: 2em; max-width: 960px; margin: auto; }}
    h1 {{ color: #2c3e50; }}
    h2 {{ color: #2c3e50; margin-top: 2em; border-bottom: 2px solid #2c3e50; padding-bottom: 4px; }}
    h3 {{ color: #555; margin-top: 1.2em; }}
    p  {{ color: #444; }}
    .table-wrap {{ overflow-x: auto; margin-bottom: 1em; }}
    table {{ border-collapse: collapse; min-width: 480px; }}
    th, td {{ border: 1px solid #ccc; padding: 6px 14px; text-align: right; white-space: nowrap; }}
    th {{ background: #2c3e50; color: white; text-align: center; }}
    td:first-child, th:first-child {{ text-align: left; }}
  </style>
</head>
<body>
<h1>rRNA Family Breakdown - {label}</h1>
{plots_section}
<p>Reads assigned = reads in <code>aligned.blast</code>.
   Totals and family assignment are derived from the subsampled reads and
   <code>rRNA_test_10M_family.tsv</code> (source sequence ID to family mapping).</p>
{''.join(sections)}
</body>
</html>"""

    out_path = Path(output_dir) / f'{label}_family_breakdown.html'
    out_path.write_text(html)
    print(f'\n  Saved: {out_path}')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--output-dir', required=True, type=Path)
    ap.add_argument('--label', default='sortmerna', help='Plot title prefix and filename stem')
    ap.add_argument('--evalues', nargs='+', type=float, required=True,
                    help='E-value thresholds in the same order as --rrna-dirs / --nonrrna-dirs')
    ap.add_argument('--rrna-dirs', nargs='+', type=Path, action='append', required=True,
                    help='Top-level run_scalability.sh output dirs for rRNA reads at each E-value '
                         '(repeat flag once per series)')
    ap.add_argument('--nonrrna-dirs', nargs='+', type=Path, action='append', required=True,
                    help='Top-level run_scalability.sh output dirs for non-rRNA reads at each E-value '
                         '(repeat flag once per series)')
    ap.add_argument('--series-labels', nargs='+', action='append',
                    help='Label for each series (repeat flag once per series); '
                         'defaults to "series 1", "series 2", ...')
    ap.add_argument('--rrna-family-tsv', type=Path, default=None,
                    help='TSV mapping source_seq_id to rRNA family, produced by '
                         'simulate_rrna_reads.sh (rRNA_test_10M_family.tsv). '
                         'When provided, generates a per-family sensitivity breakdown HTML.')
    args = ap.parse_args()

    n_series = len(args.rrna_dirs)
    if len(args.nonrrna_dirs) != n_series:
        ap.error('--rrna-dirs and --nonrrna-dirs must be repeated the same number of times')

    series_labels = []
    if args.series_labels:
        for grp in args.series_labels:
            series_labels.append(' '.join(grp))
    while len(series_labels) < n_series:
        series_labels.append(f'series {len(series_labels) + 1}')

    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_series = []
    for i, (rrna_dirs, nonrrna_dirs, slabel) in enumerate(
            zip(args.rrna_dirs, args.nonrrna_dirs, series_labels)):
        if len(rrna_dirs) != len(args.evalues) or len(nonrrna_dirs) != len(args.evalues):
            ap.error(f'series {i+1}: number of dirs must match number of --evalues '
                     f'({len(args.evalues)})')
        print(f'\nSeries: {slabel}')
        points = compute_series(rrna_dirs, nonrrna_dirs, args.evalues)
        all_series.append((slabel, points))

    roc_b64 = plot_roc(all_series, args.output_dir, args.label)

    if args.rrna_family_tsv:
        if not args.rrna_family_tsv.exists():
            print(f'\n  WARNING: --rrna-family-tsv not found: {args.rrna_family_tsv}',
                  file=sys.stderr)
        else:
            print('\nComputing per-family sensitivity breakdown...')
            family_map = load_family_map(args.rrna_family_tsv)
            print(f'  Loaded {len(family_map):,} source sequence IDs')
            tables = _build_family_tables(args.rrna_dirs, args.evalues, family_map)
            print('\nCollecting performance data...')
            perf_data = _read_perf_data(args.rrna_dirs, args.evalues)
            runtime_b64, ram_b64 = _make_perf_plots(perf_data, args.evalues, series_labels)
            render_family_html(tables, args.evalues, args.label, args.output_dir,
                               roc_b64=roc_b64, runtime_b64=runtime_b64, ram_b64=ram_b64)


if __name__ == '__main__':
    main()
