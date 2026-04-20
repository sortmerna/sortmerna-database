#!/usr/bin/env python3
"""
Generate a pivot-style HTML clustering summary table from a TSV results file.

Input TSV columns (no header):
  ref_db, db, kingdom, clustering, num_seqs, total_nt

Output is a markdown file containing an HTML table with three header levels:
  Level 1 — source database (e.g. Silva 138.2 SSURef NR99, RFAM 15.1)
  Level 2 — kingdom or family (Archaea, Bacteria, Eukaryota, 5S, 5.8S)
  Level 3 — metric (%clust, #seqs)

Rows correspond to clustering configurations (Complete DB, SMR v4.7 sensitive, etc.)
derived from the clustering threshold column.

Usage:
    python3 generate_summary.py <results.tsv> --output <summary.md> --thresholds 97 95 90 85

Exit codes:
    0 — success
    1 — input file error
"""

import argparse
import sys
from collections import defaultdict


# Ordered column definitions: (ref_db_prefix, db_key, kingdom, header_label)
# ref_db_prefix is matched with startswith() to handle version strings
COLUMNS = [
  ("SILVA", "SSU Ref NR 99", "archaea",   "Archaea"),
  ("SILVA", "SSU Ref NR 99", "bacteria",  "Bacteria"),
  ("SILVA", "SSU Ref NR 99", "eukaryota", "Eukaryota"),
  ("SILVA", "LSU Ref NR 99", "archaea",   "Archaea"),
  ("SILVA", "LSU Ref NR 99", "bacteria",  "Bacteria"),
  ("SILVA", "LSU Ref NR 99", "eukaryota", "Eukaryota"),
  ("RFAM",  "5S",            "root",      "5S"),
  ("RFAM",  "5.8S",          "eukaryota", "5.8S"),
]

# Labels for known thresholds; any unlisted threshold falls back to its value
THRESHOLD_LABELS = {
  "99%":  "SILVA NR99 (unclustered)",
  "100%": "RFAM (unclustered)",
  "97%":  "SMR v4.7 sensitive db",
  "95%":  "SMR v4.7 default db",
  "90%":  "SMR v4.7 fast db",
  "85%":  "SMR v4.7 fast db (85%)",
}


def load_tsv(tsv_path):
  data = {}
  with open(tsv_path) as f:
    for line in f:
      ref_db, db, kingdom, clustering, num_seqs, total_nt = line.rstrip("\n").split("\t")
      data[(ref_db, db, kingdom, clustering)] = int(num_seqs)
  return data


def find_value(data, ref_db_prefix, db_key, kingdom, threshold):
  for (ref_db, db, kgm, clust), num_seqs in data.items():
    if ref_db.startswith(ref_db_prefix) and db == db_key and kgm == kingdom and clust == threshold:
      return ref_db, num_seqs
  return None, None


def get_baseline_threshold(ref_db_prefix):
  """Return the unclustered baseline threshold for a given database."""
  if ref_db_prefix == "SILVA":
    return "99%"
  return "100%"


def build_top_headers(data):
  """Derive level-1 header labels and their column spans from actual data."""
  groups = {}
  for ref_db_prefix, db_key, kingdom, _ in COLUMNS:
    baseline = get_baseline_threshold(ref_db_prefix)
    ref_db, _ = find_value(data, ref_db_prefix, db_key, kingdom, baseline)
    if ref_db is None:
      # Fall back to prefix for label
      ref_db = ref_db_prefix
    label = format_source_label(ref_db, db_key)
    groups.setdefault(label, 0)
    groups[label] += 1
  return list(groups.items())  # [(label, colspan), ...]


def format_source_label(ref_db, db_key):
  if "SSU" in db_key:
    return f"{ref_db} SSURef NR99"
  elif "LSU" in db_key:
    return f"{ref_db} LSURef NR99"
  else:
    return ref_db


def generate_html_table(data, thresholds):
  top_headers = build_top_headers(data)

  # Build ordered row list: baselines first (99% for SILVA, 100% for RFAM), then clustering thresholds
  threshold_rows = [("99%", THRESHOLD_LABELS.get("99%", "99%")),
                    ("100%", THRESHOLD_LABELS.get("100%", "100%"))]
  for t in thresholds:
    pct = f"{t}%"
    threshold_rows.append((pct, THRESHOLD_LABELS.get(pct, pct)))

  rows = []

  # Header row 1: source databases
  cells = ['<th rowspan="3">Configuration</th>']
  for label, colspan in top_headers:
    cells.append(f'<th colspan="{colspan * 2}">{label}</th>')
  rows.append(f'  <tr>{"".join(cells)}</tr>')

  # Header row 2: kingdoms / families
  cells = []
  for _, _, _, label in COLUMNS:
    cells.append(f'<th colspan="2">{label}</th>')
  rows.append(f'  <tr>{"".join(cells)}</tr>')

  # Header row 3: metrics
  cells = []
  for _ in COLUMNS:
    cells.append("<th>%clust</th><th>#seqs</th>")
  rows.append(f'  <tr>{"".join(cells)}</tr>')

  # Data rows
  for threshold, row_label in threshold_rows:
    has_data = any(
      find_value(data, ref_db_prefix, db_key, kingdom, threshold)[1] is not None
      for ref_db_prefix, db_key, kingdom, _ in COLUMNS
    )
    if not has_data:
      continue

    cells = [f"<td><b>{row_label}</b></td>"]
    for ref_db_prefix, db_key, kingdom, _ in COLUMNS:
      _, num_seqs = find_value(data, ref_db_prefix, db_key, kingdom, threshold)
      if num_seqs is not None:
        cells.append(f"<td>{threshold}</td><td>{num_seqs:,}</td>")
      else:
        cells.append("<td>—</td><td>—</td>")
    rows.append(f'  <tr>{"".join(cells)}</tr>')

  thead_rows = rows[:3]
  tbody_rows = rows[3:]

  return (
    "<table>\n"
    "  <thead>\n"
    + "\n".join(f"    {r}" for r in thead_rows) + "\n"
    "  </thead>\n"
    "  <tbody>\n"
    + "\n".join(f"    {r}" for r in tbody_rows) + "\n"
    "  </tbody>\n"
    "</table>"
  )


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("tsv", help="TSV results file written by cluster_sequences.sh")
  parser.add_argument("--output", required=True, metavar="FILE", help="Output HTML file")
  parser.add_argument("--thresholds", nargs="+", required=True, metavar="N",
                      help="Clustering thresholds as integers (e.g. 97 95 90 85)")
  args = parser.parse_args()

  try:
    data = load_tsv(args.tsv)
  except OSError as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)

  table = generate_html_table(data, args.thresholds)

  html = (
    "<!DOCTYPE html>\n"
    "<html lang=\"en\">\n"
    "<head>\n"
    "  <meta charset=\"UTF-8\">\n"
    "  <title>SortMeRNA Database Clustering Summary</title>\n"
    "  <style>\n"
    "    *, *::before, *::after { box-sizing: border-box; }\n"
    "    body {\n"
    "      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;\n"
    "      background: #f5f7fa;\n"
    "      color: #222;\n"
    "      margin: 0;\n"
    "      padding: 2.5em 3em;\n"
    "    }\n"
    "    .container { max-width: 1400px; margin: 0 auto; }\n"
    "    h1 {\n"
    "      font-size: 1.7em;\n"
    "      font-weight: 700;\n"
    "      color: #1a3a5c;\n"
    "      margin-bottom: 0.2em;\n"
    "    }\n"
    "    .subtitle {\n"
    "      font-size: 1em;\n"
    "      color: #555;\n"
    "      margin-bottom: 1.2em;\n"
    "    }\n"
    "    .description {\n"
    "      background: #eaf1fb;\n"
    "      border-left: 4px solid #3a7bd5;\n"
    "      border-radius: 4px;\n"
    "      padding: 0.9em 1.2em;\n"
    "      margin-bottom: 1.8em;\n"
    "      font-size: 0.92em;\n"
    "      line-height: 1.6;\n"
    "      color: #333;\n"
    "    }\n"
    "    .description p { margin: 0 0 0.5em; }\n"
    "    .description p:last-child { margin-bottom: 0; }\n"
    "    .description ul { margin: 0.3em 0 0 1.2em; padding: 0; }\n"
    "    .description li { margin-bottom: 0.2em; }\n"
    "    .table-wrap { overflow-x: auto; border-radius: 6px; box-shadow: 0 1px 6px rgba(0,0,0,0.1); }\n"
    "    table { border-collapse: collapse; width: 100%; background: #fff; font-size: 0.88em; }\n"
    "    th, td {\n"
    "      border: 1px solid #d0d7e3;\n"
    "      padding: 7px 13px;\n"
    "      text-align: center;\n"
    "      white-space: nowrap;\n"
    "    }\n"
    "    thead tr:first-child th { background: #1a3a5c; color: #fff; font-size: 0.95em; }\n"
    "    thead tr:nth-child(2) th { background: #2e5f9e; color: #fff; }\n"
    "    thead tr:nth-child(3) th { background: #dce6f5; color: #1a3a5c; font-weight: 600; }\n"
    "    tbody tr:nth-child(odd)  { background: #fff; }\n"
    "    tbody tr:nth-child(even) { background: #f4f7fc; }\n"
    "    tbody tr:hover { background: #ddeeff; }\n"
    "    td:first-child { text-align: left; font-weight: 500; min-width: 200px; }\n"
    "    td { color: #333; }\n"
    "    .footer { margin-top: 1.2em; font-size: 0.8em; color: #888; }\n"
    "  </style>\n"
    "</head>\n"
    "<body>\n"
    "  <div class=\"container\">\n"
    "    <h1>SortMeRNA Database Clustering Summary</h1>\n"
    "    <p class=\"subtitle\">Sequence counts after vsearch clustering at multiple identity thresholds</p>\n"
    "    <div class=\"description\">\n"
    "      <p>This table summarises how many representative sequences remain in each reference database\n"
    "      after clustering at decreasing identity thresholds. Each column group covers a source database\n"
    "      and taxonomic kingdom; each row corresponds to a SortMeRNA database configuration:</p>\n"
    "      <ul>\n"
    "        <li><b>%clust</b> — identity threshold used for vsearch clustering</li>\n"
    "        <li><b>#seqs</b> — number of representative (centroid) sequences retained</li>\n"
    "      </ul>\n"
    "      <p>Baseline rows (SILVA NR99 / RFAM unclustered) show the full unfiltered sequence counts\n"
    "      against which the clustered databases can be compared.\n"
    "      A dash (—) indicates that no sequences were available for that kingdom/database combination.</p>\n"
    "    </div>\n"
    "    <div class=\"table-wrap\">\n"
    f"      {table}\n"
    "    </div>\n"
    "    <p class=\"footer\">Generated by cluster_sequences.sh &mdash; SortMeRNA database pipeline</p>\n"
    "  </div>\n"
    "</body>\n"
    "</html>\n"
  )

  with open(args.output, "w") as f:
    f.write(html)

  print(f"Summary table written to: {args.output}")


if __name__ == "__main__":
  main()
