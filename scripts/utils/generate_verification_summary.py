#!/usr/bin/env python3
"""
generate_verification_summary.py - HTML summary of cmsearch verification results.

Reads the stats TSV produced by parse_cmsearch.py (--stats) and writes
a self-contained HTML table showing input / kept / flagged counts per domain.
"""

import argparse
import sys
from pathlib import Path


def load_stats(tsv_fp):
  rows = []
  with open(tsv_fp) as f:
    for line in f:
      line = line.strip()
      if not line or line.startswith('#'):
        continue
      parts = line.split('\t')
      if len(parts) != 5:
        continue
      gene, domain, n_input, n_kept, n_flagged = parts
      rows.append({
        'gene':      gene.upper(),
        'domain':    domain,
        'n_input':   int(n_input),
        'n_kept':    int(n_kept),
        'n_flagged': int(n_flagged),
      })
  return rows


def pct(n, total):
  return f"{n / total * 100:.1f}%" if total else "—"


def render_html(rows, silva_ssu_version, silva_lsu_version, output_fp):
  total_input   = sum(r['n_input']   for r in rows)
  total_kept    = sum(r['n_kept']    for r in rows)
  total_flagged = sum(r['n_flagged'] for r in rows)

  gene_order   = ['SSU', 'LSU']
  domain_order = [
    'bacteria', 'archaea',
    'eukaryota_nuclear', 'eukaryota_mito', 'eukaryota_chloro',
    'eukaryota',
  ]

  def sort_key(r):
    g = gene_order.index(r['gene'])   if r['gene']   in gene_order   else 99
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
    bg = colours.get(domain, '#e2e3e5')
    label = domain.replace('_', ' ')
    return f'<span style="background:{bg};padding:2px 8px;border-radius:4px;font-size:0.85em">{label}</span>'

  def flagged_cell(n_flagged, n_input):
    p = n_flagged / n_input * 100 if n_input else 0
    colour = '#c0392b' if p > 5 else ('#e67e22' if p > 1 else '#27ae60')
    return (f'<td style="text-align:right;color:{colour};font-weight:600">'
        f'{n_flagged:,} ({p:.1f}%)</td>')

  body_rows = ''
  for r in rows:
    body_rows += (
      f'<tr>'
      f'<td><b>{r["gene"]}</b></td>'
      f'<td>{badge(r["domain"])}</td>'
      f'<td style="text-align:right">{r["n_input"]:,}</td>'
      f'<td style="text-align:right;color:#27ae60;font-weight:600">{r["n_kept"]:,}</td>'
      + flagged_cell(r['n_flagged'], r['n_input']) +
      f'</tr>\n'
    )

  flagged_pct = total_flagged / total_input * 100 if total_input else 0
  totals_colour = '#c0392b' if flagged_pct > 5 else ('#e67e22' if flagged_pct > 1 else '#27ae60')

  html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>SILVA Verification Summary</title>
<style>
  body {{ font-family: system-ui, sans-serif; margin: 2em; color: #333; }}
  h2   {{ margin-bottom: 0.25em; }}
  p.meta {{ color: #666; font-size: 0.9em; margin-top: 0; }}
  table {{ border-collapse: collapse; min-width: 600px; }}
  th, td {{ border: 1px solid #dee2e6; padding: 8px 14px; }}
  th {{ background: #f8f9fa; text-align: left; }}
  tr:hover {{ background: #f1f3f5; }}
  tfoot td {{ font-weight: 700; background: #f8f9fa; }}
</style>
</head>
<body>
<h2>SILVA Verification Summary</h2>
<p class="meta">
  SSU version: {silva_ssu_version} &nbsp;|&nbsp;
  LSU version: {silva_lsu_version} &nbsp;|&nbsp;
  Tool: Infernal cmsearch --cut_ga vs Rfam CMs
</p>
<table>
<thead>
  <tr>
  <th>Gene</th>
  <th>Domain</th>
  <th style="text-align:right">Input seqs</th>
  <th style="text-align:right">Kept</th>
  <th style="text-align:right">Flagged (dropped)</th>
  </tr>
</thead>
<tbody>
{body_rows}</tbody>
<tfoot>
  <tr>
  <td colspan="2">Total</td>
  <td style="text-align:right">{total_input:,}</td>
  <td style="text-align:right;color:#27ae60">{total_kept:,}</td>
  <td style="text-align:right;color:{totals_colour}">{total_flagged:,} ({flagged_pct:.1f}%)</td>
  </tr>
</tfoot>
</table>
<p style="font-size:0.85em;color:#666;margin-top:1em">
  Flagged sequences had no above-threshold Rfam hit (inc&nbsp;≠&nbsp;'!')
  and are excluded from clustering.
</p>
</body>
</html>
"""
  Path(output_fp).write_text(html)
  print(f"Verification summary written to: {output_fp}")


def main():
  ap = argparse.ArgumentParser(description=__doc__)
  ap.add_argument('stats_tsv',              help='Stats TSV produced by parse_cmsearch.py --stats')
  ap.add_argument('--output',  required=True, help='Output HTML file')
  ap.add_argument('--silva-ssu-version', default='', dest='ssu_version')
  ap.add_argument('--silva-lsu-version', default='', dest='lsu_version')
  args = ap.parse_args()

  rows = load_stats(args.stats_tsv)
  if not rows:
    print("No stats rows found — nothing to summarise.", file=sys.stderr)
    sys.exit(1)

  render_html(rows, args.ssu_version, args.lsu_version, args.output)


if __name__ == '__main__':
  main()
