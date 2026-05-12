#!/usr/bin/env python3
"""
generate_verification_summary.py - HTML summary of cmsearch verification results.

Reads the stats TSV produced by parse_cmsearch.py (--stats) and writes
a self-contained HTML page showing input / kept / flagged counts per domain,
plus a histogram of seq_from (5' alignment offset) distribution.
"""

import argparse
import sys
from pathlib import Path


MODEL_MAP = {
  ('SSU', 'bacteria'):          'RF00177',
  ('SSU', 'archaea'):           'RF01959',
  ('SSU', 'eukaryota'):         'RF01960',
  ('SSU', 'eukaryota_nuclear'): 'RF01960',
  ('SSU', 'eukaryota_mito'):    'RF02545',
  ('SSU', 'eukaryota_chloro'):  'RF01960',
  ('LSU', 'bacteria'):          'RF02541',
  ('LSU', 'archaea'):           'RF02540',
  ('LSU', 'eukaryota'):         'RF02543',
  ('LSU', 'eukaryota_nuclear'): 'RF02543',
  ('LSU', 'eukaryota_mito'):    'RF02546',
  ('LSU', 'eukaryota_chloro'):  'RF02543',
  ('5S',   'full'):             'RF00001',
  ('5S',   'seed'):             'RF00001',
  ('5.8S', 'full'):             'RF00002',
  ('5.8S', 'seed'):             'RF00002',
}

BIN_LABELS = ['= 1', '2–5', '6–10', '11–50', '51–200', '201–500', '> 500']
BIN_COLORS = ['#27ae60', '#82c91e', '#fab005', '#fd7e14', '#e67e22', '#c0392b', '#7b0000']


def classify_bin(sf):
  if sf == 1:   return 0
  if sf <= 5:   return 1
  if sf <= 10:  return 2
  if sf <= 50:  return 3
  if sf <= 200: return 4
  if sf <= 500: return 5
  return 6


def load_stats(tsv_fp):
  rows = []
  with open(tsv_fp) as f:
    for line in f:
      line = line.strip()
      if not line or line.startswith('#'):
        continue
      parts = line.split('\t')
      if len(parts) != 6:
        continue
      gene, domain, n_input, n_cmsearch_kept, n_cmsearch_dropped, n_trimming_passed = parts
      gene_upper = gene.upper()
      rows.append({
        'gene':              gene_upper,
        'domain':            domain,
        'n_input':           int(n_input),
        'n_cmsearch_kept':   int(n_cmsearch_kept),
        'n_cmsearch_dropped': int(n_cmsearch_dropped),
        'n_trimming_passed': int(n_trimming_passed),
        'rfam_model':        MODEL_MAP.get((gene_upper, domain), '-'),
      })
  return rows


def load_seq_from_dist(log_dir, gene, domain):
  """Return list of bin counts, or None if log file is missing."""
  fp = Path(log_dir) / f"cmsearch_log_{gene.lower()}_{domain}.tsv"
  if not fp.exists():
    return None
  counts = [0] * len(BIN_LABELS)
  with open(fp) as f:
    next(f)  # skip header
    for line in f:
      parts = line.split('\t')
      if len(parts) < 4:
        continue
      try:
        counts[classify_bin(int(parts[3]))] += 1
      except (ValueError, IndexError):
        continue
  return counts


def pct(n, total):
  return f"{n / total * 100:.3f}%" if total else "-"


def render_histogram_svg(rows_with_dist):
  label_w  = 170
  bar_w    = 450
  row_h    = 26
  gap      = 6
  legend_h = 32
  pad_top  = 10
  pad_bot  = 10
  svg_w    = label_w + bar_w + 10

  n     = len(rows_with_dist)
  svg_h = pad_top + legend_h + n * (row_h + gap) + pad_bot

  parts = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w}" height="{svg_h}" '
    f'style="font-family:system-ui,sans-serif;font-size:12px">'
  ]

  # Legend - evenly spaced across bar_w
  item_w = bar_w // len(BIN_LABELS)
  for i, (label, color) in enumerate(zip(BIN_LABELS, BIN_COLORS)):
    rx = label_w + i * item_w
    parts.append(f'<rect x="{rx}" y="{pad_top}" width="12" height="12" fill="{color}"/>')
    parts.append(f'<text x="{rx + 15}" y="{pad_top + 11}" fill="#333">{label}</text>')

  # One stacked bar per row
  for row_i, (label, counts) in enumerate(rows_with_dist):
    y     = pad_top + legend_h + row_i * (row_h + gap)
    total = sum(counts) or 1

    parts.append(
      f'<text x="{label_w - 6}" y="{y + row_h // 2 + 4}" '
      f'text-anchor="end" fill="#333">{label}</text>'
    )

    x = label_w
    for i, (cnt, color) in enumerate(zip(counts, BIN_COLORS)):
      w = cnt / total * bar_w
      if w < 0.5:
        continue
      pct_str = f"{cnt / total * 100:.1f}%"
      parts.append(
        f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="{row_h}" fill="{color}">'
        f'<title>{BIN_LABELS[i]}: {cnt:,} ({pct_str})</title>'
        f'</rect>'
      )
      x += w

  parts.append('</svg>')
  return '\n'.join(parts)


def render_html(rows, version, output_fp, log_dir, title="SILVA Verification Summary", tool="cmsearch --hmmonly --cut_ga"):
  total_input           = sum(r['n_input']           for r in rows)
  total_cmsearch_kept   = sum(r['n_cmsearch_kept']   for r in rows)
  total_cmsearch_dropped = sum(r['n_cmsearch_dropped'] for r in rows)
  total_trimming_passed = sum(r['n_trimming_passed'] for r in rows)

  gene_order   = ['SSU', 'LSU']
  domain_order = [
    'bacteria', 'archaea',
    'eukaryota_nuclear', 'eukaryota_mito', 'eukaryota_chloro',
    'eukaryota',
  ]

  def sort_key(r):
    g = gene_order.index(r['gene'])     if r['gene']   in gene_order   else 99
    d = domain_order.index(r['domain']) if r['domain'] in domain_order else 99
    return (g, d)

  rows = sorted(rows, key=sort_key)

  def badge(domain):
    colours = {
      'bacteria':          '#d4edda',
      'archaea':           '#fff3cd',
      'eukaryota_nuclear': '#cce5ff',
      'eukaryota_mito':    '#f8d7da',
      'eukaryota_chloro':  '#d1ecf1',
      'eukaryota':         '#cce5ff',
    }
    bg    = colours.get(domain, '#e2e3e5')
    label = domain.replace('_', ' ')
    return f'<span style="background:{bg};padding:2px 8px;border-radius:4px;font-size:0.85em">{label}</span>'

  def dropped_cell(n_dropped, n_input):
    p      = n_dropped / n_input * 100 if n_input else 0
    colour = '#c0392b' if p > 5 else ('#e67e22' if p > 1 else '#27ae60')
    return (f'<td style="text-align:right;color:{colour};font-weight:600">'
        f'{n_dropped:,} ({p:.1f}%)</td>')

  body_rows = ''
  for r in rows:
    body_rows += (
      f'<tr>'
      f'<td><b>{r["gene"]}</b></td>'
      f'<td>{badge(r["domain"])}</td>'
      f'<td style="text-align:right;font-family:monospace">{r["rfam_model"]}</td>'
      f'<td style="text-align:right">{r["n_input"]:,}</td>'
      f'<td style="text-align:right;color:#27ae60;font-weight:600">{r["n_cmsearch_kept"]:,}</td>'
      + dropped_cell(r['n_cmsearch_dropped'], r['n_input']) +
      f'<td style="text-align:right;color:#1971c2;font-weight:600">'
      f'{r["n_trimming_passed"]:,} ({pct(r["n_trimming_passed"], r["n_input"])})</td>'
      f'</tr>\n'
    )

  dropped_pct   = total_cmsearch_dropped / total_input * 100 if total_input else 0
  totals_colour = '#c0392b' if dropped_pct > 5 else ('#e67e22' if dropped_pct > 1 else '#27ae60')

  # Histogram section
  rows_with_dist = []
  for r in rows:
    counts = load_seq_from_dist(log_dir, r['gene'], r['domain'])
    if counts is not None:
      label = f"{r['gene']} {r['domain'].replace('_', ' ')}"
      rows_with_dist.append((label, counts))

  hist_section = ''
  if rows_with_dist:
    svg = render_histogram_svg(rows_with_dist)
    hist_section = f"""
<h3 style="margin-top:2em">5&#x2032; Alignment Offset Distribution (seq_from) - Kept sequences only</h3>
<p class="meta">
  Each kept sequence is trimmed to [seq_from,&nbsp;seq_to]. When multiple hits are merged
  into a single span, seq_from&nbsp;=&nbsp;min(seq_from) and seq_to&nbsp;=&nbsp;max(seq_to)
  across all hits. When hits fail the coverage threshold or there is only one hit,
  seq_from and seq_to are taken from the best-scoring hit.
  seq_from&nbsp;=&nbsp;1 means no 5&#x2032; trimming was needed.
  Hover over a bar segment to see exact counts.
</p>
{svg}
"""

  html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{title}</title>
<style>
  body {{ font-family: system-ui, sans-serif; margin: 2em; color: #333; }}
  h2   {{ margin-bottom: 0.25em; }}
  h3   {{ margin-bottom: 0.25em; }}
  p.meta {{ color: #666; font-size: 0.9em; margin-top: 0; }}
  table {{ border-collapse: collapse; min-width: 600px; }}
  th, td {{ border: 1px solid #dee2e6; padding: 8px 14px; }}
  th {{ background: #f8f9fa; text-align: left; }}
  tr:hover {{ background: #f1f3f5; }}
  tfoot td {{ font-weight: 700; background: #f8f9fa; }}
</style>
</head>
<body>
<h2>{title}</h2>
<p class="meta">
  Version: {version} &nbsp;|&nbsp;
  Tool: Infernal {tool}
</p>
<table>
<thead>
  <tr>
  <th>Gene</th>
  <th>Domain</th>
  <th style="text-align:right">Rfam model</th>
  <th style="text-align:right">Input seqs</th>
  <th style="text-align:right">cmsearch kept</th>
  <th style="text-align:right">cmsearch dropped</th>
  <th style="text-align:right">trimming passed</th>
  </tr>
</thead>
<tbody>
{body_rows}</tbody>
<tfoot>
  <tr>
  <td colspan="3">Total</td>
  <td style="text-align:right">{total_input:,}</td>
  <td style="text-align:right;color:#27ae60">{total_cmsearch_kept:,}</td>
  <td style="text-align:right;color:{totals_colour}">{total_cmsearch_dropped:,} ({dropped_pct:.1f}%)</td>
  <td style="text-align:right;color:#1971c2">{total_trimming_passed:,} ({total_trimming_passed / total_input * 100:.1f}%)</td>
  </tr>
</tfoot>
</table>
<p style="font-size:0.85em;color:#666;margin-top:1em">
  <b>cmsearch dropped</b>: no above-threshold hit (inc&nbsp;&#x2260;&nbsp;'!') - excluded from clustering.<br>
  <b>trimming passed</b>: passed cmsearch <em>and</em> trimmed length &ge; 50% of input sequence length - written to verified_*.fasta.
</p>
{hist_section}
</body>
</html>
"""
  Path(output_fp).write_text(html)
  print(f"Verification summary written to: {output_fp}")


def main():
  ap = argparse.ArgumentParser(description=__doc__)
  ap.add_argument('stats_tsv',               help='Stats TSV produced by parse_cmsearch.py --stats')
  ap.add_argument('--output',  required=True, help='Output HTML file')
  ap.add_argument('--version', default='', help='Database version string shown in the HTML header')
  ap.add_argument('--title', default='SILVA Verification Summary', help='HTML page title and heading')
  ap.add_argument('--tool', default='cmsearch --hmmonly --cut_ga', help='Tool invocation shown in the HTML header')
  args = ap.parse_args()

  rows = load_stats(args.stats_tsv)
  if not rows:
    print("No stats rows found - nothing to summarise.", file=sys.stderr)
    sys.exit(1)

  log_dir = Path(args.stats_tsv).parent
  render_html(rows, args.version, args.output, log_dir, args.title, args.tool)


if __name__ == '__main__':
  main()
