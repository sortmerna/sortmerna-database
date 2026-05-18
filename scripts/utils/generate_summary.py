#!/usr/bin/env python3
"""
Generate two pivot-style HTML clustering summary tables from a TSV results file.

Input TSV columns (no header):
  ref_db, db, kingdom, clustering, num_seqs, total_nt

Table 1 - sequence counts per database source / kingdom across all thresholds
Table 2 - recommended database configurations with per-domain thresholds and counts

Usage:
  python3 generate_summary.py <results.tsv> --output <summary.html> --thresholds 97 95 90 85

Exit codes:
  0 - success
  1 - input file error
"""

import argparse
import sys


# Ordered column definitions for Table 1: (ref_db_prefix, db_key, kingdom, header_label)
# ref_db_prefix is matched with startswith() to handle version strings
COLUMNS = [
  ("SILVA", "SSU Ref NR 99", "archaea",   "Archaea"),
  ("SILVA", "SSU Ref NR 99", "bacteria",  "Bacteria"),
  ("SILVA", "SSU Ref NR 99", "eukaryota", "Eukaryota"),
  ("SILVA", "LSU Ref NR 99", "archaea",   "Archaea"),
  ("SILVA", "LSU Ref NR 99", "bacteria",  "Bacteria"),
  ("SILVA", "LSU Ref NR 99", "eukaryota", "Eukaryota"),
  ("Rfam",  "5S",            "root",      "5S"),
  ("Rfam",  "5.8S",          "eukaryota", "5.8S"),
]

# Labels for known clustering keys
THRESHOLD_LABELS = {
  "verified": "original (trimmed, filtered)",
  "97%":      "SMR v4.7 sensitive db",
  "95%":      "SMR v4.7 default db",
  "90%":      "SMR v4.7 fast db",
  "85%":      "SMR v4.7 fast db (85%)",
}

# Table 2: flat column list (no source-db grouping)
TABLE2_COLS = [
  ("SILVA", "SSU Ref NR 99", "archaea",   "Archaea SSU"),
  ("SILVA", "SSU Ref NR 99", "bacteria",  "Bacteria SSU"),
  ("SILVA", "SSU Ref NR 99", "eukaryota", "Eukaryota SSU"),
  ("SILVA", "LSU Ref NR 99", "archaea",   "Archaea LSU"),
  ("SILVA", "LSU Ref NR 99", "bacteria",  "Bacteria LSU"),
  ("SILVA", "LSU Ref NR 99", "eukaryota", "Eukaryota LSU"),
  ("Rfam",  "5S",            "root",      "5S"),
  ("Rfam",  "5.8S",          "eukaryota", "5.8S"),
]

# Table 2: per-configuration, per-column lookup spec.
# Each inner list has one entry per TABLE2_COLS column:
#   (ref_db_prefix, db_key, kingdom, threshold) or None = not applicable ('-')
TABLE2_CONFIGS = [
  ("SMR v4.7 sensitive db", [
    ("SILVA", "SSU Ref NR 99", "archaea",   "97%"),
    ("SILVA", "SSU Ref NR 99", "bacteria",  "97%"),
    ("SILVA", "SSU Ref NR 99", "eukaryota", "97%"),
    ("SILVA", "LSU Ref NR 99", "archaea",   "97%"),
    ("SILVA", "LSU Ref NR 99", "bacteria",  "97%"),
    ("SILVA", "LSU Ref NR 99", "eukaryota", "97%"),
    ("Rfam",  "5S",            "root",      "97%"),
    ("Rfam",  "5.8S",          "eukaryota", "97%"),
  ]),
  ("SMR v4.7 default db", [
    ("SILVA", "SSU Ref NR 99", "archaea",   "95%"),
    ("SILVA", "SSU Ref NR 99", "bacteria",  "90%"),
    ("SILVA", "SSU Ref NR 99", "eukaryota", "95%"),
    ("SILVA", "LSU Ref NR 99", "archaea",   "95%"),
    ("SILVA", "LSU Ref NR 99", "bacteria",  "95%"),
    ("SILVA", "LSU Ref NR 99", "eukaryota", "95%"),
    ("Rfam",  "5S seed",       "root",      "100%"),
    ("Rfam",  "5.8S seed",     "eukaryota", "100%"),
  ]),
  ("SMR v4.7 fast db", [
    ("SILVA", "SSU Ref NR 99", "archaea",   "90%"),
    ("SILVA", "SSU Ref NR 99", "bacteria",  "85%"),
    ("SILVA", "SSU Ref NR 99", "eukaryota", "90%"),
    ("SILVA", "LSU Ref NR 99", "archaea",   "90%"),
    ("SILVA", "LSU Ref NR 99", "bacteria",  "90%"),
    ("SILVA", "LSU Ref NR 99", "eukaryota", "90%"),
    ("Rfam",  "5S seed",       "root",      "100%"),
    ("Rfam",  "5.8S seed",     "eukaryota", "100%"),
  ]),
]


def load_tsv(tsv_path):
  data = {}
  with open(tsv_path) as f:
    for line in f:
      ref_db, db, kingdom, clustering, num_seqs, _ = line.rstrip("\n").split("\t")
      data[(ref_db, db, kingdom, clustering)] = int(num_seqs)
  return data


def find_value(data, ref_db_prefix, db_key, kingdom, threshold):
  for (ref_db, db, kgm, clust), num_seqs in data.items():
    if ref_db.startswith(ref_db_prefix) and db == db_key and kgm == kingdom and clust == threshold:
      return ref_db, num_seqs
  return None, None


def get_original_count(data, ref_db_prefix, db_key, kingdom):
  """Look up the pre-clustering baseline: 'verified' for SILVA, '100%' for Rfam."""
  threshold = "verified" if ref_db_prefix == "SILVA" else "100%"
  return find_value(data, ref_db_prefix, db_key, kingdom, threshold)


def get_active_columns(data, thresholds):
  """Return COLUMNS filtered to those with at least one data point."""
  all_thresholds = ["verified", "100%"] + [f"{t}%" for t in thresholds]
  result = []
  for col in COLUMNS:
    ref_db_prefix, db_key, kingdom, _ = col
    has_data = any(
      find_value(data, ref_db_prefix, db_key, kingdom, thr)[1] is not None
      for thr in all_thresholds
    )
    if has_data:
      result.append(col)
  return result


def build_top_headers(data, active_cols):
  """Derive level-1 header labels and column spans from active columns."""
  groups = {}
  for ref_db_prefix, db_key, kingdom, _ in active_cols:
    ref_db, _ = get_original_count(data, ref_db_prefix, db_key, kingdom)
    if ref_db is None:
      ref_db = ref_db_prefix
    label = format_source_label(ref_db, db_key)
    groups.setdefault(label, 0)
    groups[label] += 1
  return list(groups.items())


def format_source_label(ref_db, db_key):
  if "SSU" in db_key:
    return f"{ref_db} SSURef NR99"
  elif "LSU" in db_key:
    return f"{ref_db} LSURef NR99"
  else:
    return ref_db


def generate_html_table(data, thresholds):
  """Table 1: #seqs per source-db / kingdom / threshold. Empty columns are dropped."""
  active_cols = get_active_columns(data, thresholds)
  top_headers = build_top_headers(data, active_cols)

  rows = []

  # Header row 1: Configuration | source database groups
  cells = ['<th rowspan="2">Configuration</th>']
  for label, colspan in top_headers:
    cells.append(f'<th colspan="{colspan}">{label}</th>')
  rows.append(f'  <tr>{"".join(cells)}</tr>')

  # Header row 2: kingdom / domain labels
  cells = []
  for _, _, _, label in active_cols:
    cells.append(f'<th>{label}</th>')
  rows.append(f'  <tr>{"".join(cells)}</tr>')

  # Original (pre-clustering) row
  orig_cells = ['<td><b>original (trimmed, filtered)</b></td>']
  has_orig = False
  for ref_db_prefix, db_key, kingdom, _ in active_cols:
    _, num_seqs = get_original_count(data, ref_db_prefix, db_key, kingdom)
    if num_seqs is not None:
      orig_cells.append(f'<td>{num_seqs:,}</td>')
      has_orig = True
    else:
      orig_cells.append('<td>-</td>')
  if has_orig:
    rows.append(f'  <tr>{"".join(orig_cells)}</tr>')

  # Clustered threshold rows
  for t in thresholds:
    threshold = f"{t}%"
    has_data = any(
      find_value(data, ref_db_prefix, db_key, kingdom, threshold)[1] is not None
      for ref_db_prefix, db_key, kingdom, _ in active_cols
    )
    if not has_data:
      continue
    cells = [f'<td><b>{threshold}</b></td>']
    for ref_db_prefix, db_key, kingdom, _ in active_cols:
      _, num_seqs = find_value(data, ref_db_prefix, db_key, kingdom, threshold)
      cells.append(f'<td>{num_seqs:,}</td>' if num_seqs is not None else '<td>-</td>')
    rows.append(f'  <tr>{"".join(cells)}</tr>')

  thead_rows = rows[:2]
  tbody_rows = rows[2:]

  return (
    "<table>\n"
    "  <thead>\n"
    + "\n".join(f"    {r}" for r in thead_rows) + "\n"
    + "  </thead>\n"
    "  <tbody>\n"
    + "\n".join(f"    {r}" for r in tbody_rows) + "\n"
    + "  </tbody>\n"
    "</table>"
  )


def generate_config_table(data):
  """Table 2: recommended per-domain thresholds and counts for each database configuration."""

  # Drop columns where every config row has no data
  active_indices = []
  for i, _ in enumerate(TABLE2_COLS):
    has_data = any(
      col_specs[i] is not None and
      find_value(data, col_specs[i][0], col_specs[i][1], col_specs[i][2], col_specs[i][3])[1] is not None
      for _, col_specs in TABLE2_CONFIGS
    )
    if has_data:
      active_indices.append(i)

  active_cols2  = [TABLE2_COLS[i]   for i in active_indices]
  active_specs  = [[col_specs[i] for i in active_indices] for _, col_specs in TABLE2_CONFIGS]
  config_labels = [label for label, _ in TABLE2_CONFIGS]

  rows = []

  # Single header row
  cells = ['<th>Configuration</th>']
  for _, _, _, label in active_cols2:
    cells.append(f'<th>{label}</th>')
  rows.append(f'  <tr>{"".join(cells)}</tr>')

  for config_label, col_specs in zip(config_labels, active_specs):
    cells = [f'<td><b>{config_label}</b></td>']
    for spec in col_specs:
      if spec is None:
        cells.append('<td style="text-align:center">-</td>')
        continue
      ref_db_prefix, db_key, kingdom, threshold = spec
      _, num_seqs = find_value(data, ref_db_prefix, db_key, kingdom, threshold)
      if num_seqs is None:
        cells.append('<td style="text-align:center">-</td>')
      else:
        if "seed" in db_key:
          disp = "seed"
        elif threshold == "verified":
          disp = "original"
        else:
          disp = threshold
        cells.append(
          f'<td style="text-align:right">'
          f'<span style="font-size:0.8em;color:#666">{disp}</span><br>'
          f'{num_seqs:,}</td>'
        )
    rows.append(f'  <tr>{"".join(cells)}</tr>')

  thead_rows = rows[:1]
  tbody_rows = rows[1:]

  return (
    "<table>\n"
    "  <thead>\n"
    + "\n".join(f"    {r}" for r in thead_rows) + "\n"
    + "  </thead>\n"
    "  <tbody>\n"
    + "\n".join(f"    {r}" for r in tbody_rows) + "\n"
    + "  </tbody>\n"
    "</table>"
  )


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("tsv", help="TSV results file written by cluster_sequences.sh")
  parser.add_argument("--output", required=True, metavar="FILE", help="Output HTML file")
  parser.add_argument("--thresholds", nargs="+", required=True, metavar="N",
                      help="Clustering thresholds as integers (e.g. 97 95 90 85)")
  parser.add_argument("--silva-version", default=None, metavar="VER",
                      help="SILVA release version (e.g. 138.2)")
  parser.add_argument("--rfam-version", default=None, metavar="VER",
                      help="Rfam release version (e.g. 15.1)")
  args = parser.parse_args()

  try:
    data = load_tsv(args.tsv)
  except OSError as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)

  table1 = generate_html_table(data, args.thresholds)
  table2 = generate_config_table(data)

  db_versions = []
  if args.silva_version:
    db_versions.append(f"SILVA release <b>{args.silva_version}</b> (SSURef NR99 and LSURef NR99)")
  if args.rfam_version:
    db_versions.append(f"Rfam release <b>{args.rfam_version}</b> (5S / RF00001, 5.8S / RF00002)")
  versions_html = (
    "<p><b>Source database versions:</b></p><ul>"
    + "".join(f"<li>{v}</li>" for v in db_versions)
    + "</ul>"
  ) if db_versions else ""

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
    "    .container { max-width: 1600px; margin: 0 auto; }\n"
    "    h1 { font-size: 1.7em; font-weight: 700; color: #1a3a5c; margin-bottom: 0.2em; }\n"
    "    h2 { font-size: 1.3em; font-weight: 700; color: #1a3a5c; margin-top: 2.5em; margin-bottom: 0.4em; }\n"
    "    .subtitle { font-size: 1em; color: #555; margin-bottom: 1.2em; }\n"
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
    "    .table-wrap { overflow-x: auto; border-radius: 6px; box-shadow: 0 1px 6px rgba(0,0,0,0.1); margin-bottom: 1em; }\n"
    "    table { border-collapse: collapse; width: 100%; background: #fff; font-size: 0.88em; }\n"
    "    th, td { border: 1px solid #d0d7e3; padding: 7px 13px; text-align: center; white-space: nowrap; }\n"
    "    thead tr:first-child th { background: #1a3a5c; color: #fff; font-size: 0.95em; }\n"
    "    thead tr:nth-child(2) th { background: #dce6f5; color: #1a3a5c; font-weight: 600; }\n"
    "    tbody tr:nth-child(odd)  { background: #fff; }\n"
    "    tbody tr:nth-child(even) { background: #f4f7fc; }\n"
    "    tbody tr:hover { background: #ddeeff; }\n"
    "    td:first-child { text-align: left; font-weight: 500; min-width: 220px; }\n"
    "    td { color: #333; }\n"
    "    .footer { margin-top: 1.2em; font-size: 0.8em; color: #888; }\n"
    "  </style>\n"
    "</head>\n"
    "<body>\n"
    "  <div class=\"container\">\n"
    "    <h1>SortMeRNA Database Clustering Summary</h1>\n"
    "    <p class=\"subtitle\">Sequence counts after vsearch clustering at multiple identity thresholds</p>\n"
    "    <div class=\"description\">\n"
    "      <p>Sequence counts per source database and taxonomic domain across all clustering thresholds.\n"
    "      Columns with no data across all thresholds are omitted.\n"
    "      A dash (-) indicates no sequences were available for that combination.</p>\n"
    f"      {versions_html}\n"
    "    </div>\n"
    "    <div class=\"table-wrap\">\n"
    f"      {table1}\n"
    "    </div>\n"
    "    <h2>Recommended Database Configurations</h2>\n"
    "    <div class=\"description\">\n"
    "      <p>Recommended clustering threshold per gene/domain for each SortMeRNA database configuration,\n"
    "      with the resulting sequence count shown below each threshold label.\n"
    "      Bacteria SSU is clustered one threshold step lower than other SILVA domains.\n"
    "      Rfam 5S and 5.8S use seed sequences for the default and fast configurations.\n"
    "      A dash (-) indicates the configuration does not apply to that source database.</p>\n"
    "    </div>\n"
    "    <div class=\"table-wrap\">\n"
    f"      {table2}\n"
    "    </div>\n"
    "    <p class=\"footer\">Generated by cluster_sequences.sh - SortMeRNA database pipeline</p>\n"
    "  </div>\n"
    "</body>\n"
    "</html>\n"
  )

  with open(args.output, "w") as f:
    f.write(html)

  print(f"Summary table written to: {args.output}")


if __name__ == "__main__":
  main()
