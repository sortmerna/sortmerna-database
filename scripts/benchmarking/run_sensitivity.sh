#!/usr/bin/env bash
# run_sensitivity.sh  Experiment 2: sensitivity across database configurations.
#
# Runs SortMeRNA for each of the three Set/db pairs and computes per-set sensitivity.
# Use the e-value identified as optimal from the Experiment 1 ROC curve.
#
# Usage:
#   bash run_sensitivity.sh <output_dir> <threads> [OPTIONS]
#
# Arguments:
#   output_dir     Directory for all outputs
#   threads        Number of SortMeRNA threads (default: 4)
#
# Options:
#   --evalue FLOAT      E-value threshold (required - select from Experiment 1 ROC curve)
#   --rrna-sim-dir DIR  Directory containing set{1,2,3}_rrna_reads.fasta
#                       (default: RRNA_SIM_DIR env var)
#   --index-dir DIR     Directory containing SortMeRNA index subdirectories
#                       (default: INDEX_DIR env var)
#   --force             Re-run SortMeRNA even if output already exists
#
# Required env vars:
#   SMR_BIN          Full path to SortMeRNA binary
#   SMR_VERSION      SortMeRNA version string (e.g. 6.0.2)
#   INDEX_DIR        Directory containing SortMeRNA index subdirectories
#   RRNA_SIM_DIR     Directory containing set{1,2,3}_rrna_reads.fasta
#   SMR_DB_ROOT_DIR  Root of the sortmerna-database repository
#
# Outputs (all under <output_dir>):
#   set<N>/smr_out/out/aligned.fasta   SortMeRNA-detected rRNA reads
#   set<N>/smr_out/out/aligned.log     SortMeRNA summary log
#   set<N>/runtime_seconds.txt         Wall-clock runtime in seconds
#   set<N>/peak_rss_mb.txt             Peak resident set size in MB
#   set<N>/sensitivity.txt             Detected / total reads (fraction)
#   sensitivity_summary.html           HTML summary across all three Sets

set -euo pipefail

OUTPUT_DIR="$1"
THREADS="${2:-4}"

EVALUE=""
RRNA_SIM_DIR_OPT=""
INDEX_DIR_OPT=""
FORCE=false

shift 2 || true
while [[ $# -gt 0 ]]; do
    case "$1" in
        --evalue)       EVALUE="$2";           shift 2 ;;
        --rrna-sim-dir) RRNA_SIM_DIR_OPT="$2"; shift 2 ;;
        --index-dir)    INDEX_DIR_OPT="$2";    shift 2 ;;
        --force)        FORCE=true;            shift   ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

[[ -n "${RRNA_SIM_DIR_OPT}" ]] && RRNA_SIM_DIR="${RRNA_SIM_DIR_OPT}"
[[ -n "${INDEX_DIR_OPT}" ]]    && INDEX_DIR="${INDEX_DIR_OPT}"

: "${SMR_BIN:?SMR_BIN env var not set}"
: "${SMR_VERSION:?SMR_VERSION env var not set}"
: "${INDEX_DIR:?INDEX_DIR not set - pass --index-dir or export INDEX_DIR}"
: "${RRNA_SIM_DIR:?RRNA_SIM_DIR not set - pass --rrna-sim-dir or export RRNA_SIM_DIR}"
: "${SMR_DB_ROOT_DIR:?SMR_DB_ROOT_DIR env var not set}"
: "${EVALUE:?--evalue is required (select the optimal value from the Experiment 1 ROC curve)}"

mkdir -p "${OUTPUT_DIR}"

# Set number -> (reads file, database config, label)
declare -A SET_READS SET_DB SET_LABEL
SET_READS[1]="${RRNA_SIM_DIR}/set1_rrna_reads.fasta"
SET_READS[2]="${RRNA_SIM_DIR}/set2_rrna_reads.fasta"
SET_READS[3]="${RRNA_SIM_DIR}/set3_rrna_reads.fasta"
SET_DB[1]="smr_v${SMR_VERSION}_sensitive_db"
SET_DB[2]="smr_v${SMR_VERSION}_default_db"
SET_DB[3]="smr_v${SMR_VERSION}_fast_db"
SET_LABEL[1]="Set 1 - sensitive (97%)"
SET_LABEL[2]="Set 2 - default (90-95%)"
SET_LABEL[3]="Set 3 - fast (85-90%)"

echo "============================================"
echo "Sensitivity benchmark (Experiment 2)"
echo "E-value:     ${EVALUE}"
echo "Output:      ${OUTPUT_DIR}"
echo "Threads:     ${THREADS}"
echo "rRNA sim dir: ${RRNA_SIM_DIR}"
echo "Index dir:   ${INDEX_DIR}"
echo "============================================"

for set_num in 1 2 3; do
    reads="${SET_READS[$set_num]}"
    db_config="${SET_DB[$set_num]}"
    label="${SET_LABEL[$set_num]}"
    set_dir="${OUTPUT_DIR}/set${set_num}"
    smr_workdir="${set_dir}/smr_out"
    aligned_log="${smr_workdir}/out/aligned.log"
    aligned_fa="${smr_workdir}/out/aligned.fasta"
    sensitivity_file="${set_dir}/sensitivity.txt"

    echo ""
    echo "--------------------------------------------"
    echo "${label}"
    echo "--------------------------------------------"

    if [[ ! -f "${reads}" ]]; then
        echo "  ERROR: reads file not found: ${reads}" >&2; exit 1
    fi

    db_fasta="${INDEX_DIR}/${db_config}/${db_config}.fasta"
    db_idx="${INDEX_DIR}/${db_config}/idx"

    if [[ ! -f "${db_fasta}" ]]; then
        echo "  ERROR: database not found: ${db_fasta}" >&2; exit 1
    fi

    mkdir -p "${set_dir}"

    n_total=$(seqkit stats -T "${reads}" | tail -1 | cut -f4)
    echo "  Reads:    $(basename "${reads}") (${n_total} total)"
    echo "  Database: ${db_config}"

    if [[ -f "${aligned_log}" ]] && [[ "${FORCE}" == false ]]; then
        echo "  Already exists: aligned.log (use --force to re-run)"
    else
        [[ "${FORCE}" == true ]] && rm -rf "${smr_workdir}"
        echo "  Running SortMeRNA..."
        start=$(date +%s)
        "${SMR_BIN}" \
            --ref     "${db_fasta}" \
            --reads   "${reads}" \
            --idx-dir "${db_idx}" \
            --workdir "${smr_workdir}" \
            --fastx --blast 1 \
            --threads "${THREADS}" \
            -e "${EVALUE}" &
        smr_pid=$!
        peak_rss_mb=0
        while kill -0 "${smr_pid}" 2>/dev/null; do
            rss=$(ps -p "${smr_pid}" -o rss --no-headers 2>/dev/null || echo 0)
            rss_mb=$(( (rss + 0) / 1024 ))
            (( rss_mb > peak_rss_mb )) && peak_rss_mb=${rss_mb}
            sleep 5
        done
        wait "${smr_pid}" || { echo "  ERROR: sortmerna failed"; exit 1; }
        end=$(date +%s)
        runtime=$(( end - start ))
        echo "${runtime}"      > "${set_dir}/runtime_seconds.txt"
        echo "${peak_rss_mb}"  > "${set_dir}/peak_rss_mb.txt"
        echo "  Done: ${runtime}s - peak RSS: ${peak_rss_mb} MB"
    fi

    n_detected=$(seqkit stats -T "${aligned_fa}" | tail -1 | cut -f4)
    sensitivity=$(awk "BEGIN {printf \"%.4f\", ${n_detected}/${n_total}}")
    echo "${sensitivity}" > "${sensitivity_file}"
    echo "  Sensitivity: ${n_detected}/${n_total} = $(awk "BEGIN {printf \"%.2f\", ${sensitivity}*100}")%"
done

echo ""
echo "============================================"
echo "Writing HTML summary..."
echo "============================================"

OUTPUT_HTML="${OUTPUT_DIR}/sensitivity_summary.html"
EVALUE="${EVALUE}" OUTPUT_HTML="${OUTPUT_HTML}" \
python3 - "${OUTPUT_DIR}" "${SMR_VERSION}" <<'PYEOF'
import sys, os, datetime
from pathlib import Path

output_dir  = Path(sys.argv[1])
smr_version = sys.argv[2]
evalue      = os.environ['EVALUE']
date_str    = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

sets = [
    ("1", "Set 1 - sensitive (97%)",  f"smr_v{smr_version}_sensitive_db", "set1_rrna_reads.fasta"),
    ("2", "Set 2 - default (90-95%)", f"smr_v{smr_version}_default_db",   "set2_rrna_reads.fasta"),
    ("3", "Set 3 - fast (85-90%)",    f"smr_v{smr_version}_fast_db",      "set3_rrna_reads.fasta"),
]

rows = ""
for sn, label, db_config, reads_file in sets:
    set_dir = output_dir / f"set{sn}"
    sens_f   = set_dir / "sensitivity.txt"
    rt_f     = set_dir / "runtime_seconds.txt"
    rss_f    = set_dir / "peak_rss_mb.txt"
    aligned  = set_dir / "smr_out" / "out" / "aligned.fasta"

    if not sens_f.exists():
        continue

    sensitivity = float(sens_f.read_text().strip())
    runtime     = rt_f.read_text().strip()  if rt_f.exists()  else "-"
    peak_rss    = rss_f.read_text().strip() if rss_f.exists() else "-"
    pct         = f"{sensitivity * 100:.2f}%"

    rows += (
        f"      <tr>"
        f"<td>{label}</td>"
        f"<td>{db_config}</td>"
        f"<td>{reads_file}</td>"
        f"<td>{pct}</td>"
        f"<td>{runtime}</td>"
        f"<td>{peak_rss}</td>"
        f"</tr>\n"
    )

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Sensitivity Summary - Experiment 2</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 40px; color: #333; }}
  h1   {{ color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 8px; }}
  h2   {{ color: #34495e; margin-top: 32px; }}
  .description {{ background: #f8f9fa; border-left: 4px solid #3498db; padding: 10px 16px; margin: 12px 0; }}
  .table-wrap  {{ overflow-x: auto; margin: 16px 0; }}
  table  {{ border-collapse: collapse; min-width: 600px; }}
  th, td {{ border: 1px solid #ddd; padding: 8px 14px; text-align: left; }}
  th     {{ background: #2c3e50; color: #fff; }}
  tr:nth-child(even) {{ background: #f2f2f2; }}
  .meta   {{ font-size: 0.9em; color: #666; margin-top: 4px; }}
  .footer {{ font-size: 0.8em; color: #999; margin-top: 40px; border-top: 1px solid #ddd; padding-top: 8px; }}
</style>
</head>
<body>
<h1>Sensitivity Summary - Experiment 2</h1>
<p class="meta">Generated: {date_str}</p>
<section>
<h2>Sensitivity across database configurations</h2>
<div class="description">
<p>Each Set is tested against its matched SortMeRNA database configuration using reads simulated
from non-seed cluster members at the matching clustering threshold. E-value: {evalue}
(optimal threshold selected from Experiment 1 ROC curve). Sensitivity = detected rRNA reads / total reads.</p>
</div>
<div class="table-wrap"><table>
  <thead>
    <tr>
      <th>Set</th><th>Database</th><th>Reads file</th>
      <th>Sensitivity</th><th>Runtime (s)</th><th>Peak RSS (MB)</th>
    </tr>
  </thead>
  <tbody>
{rows}  </tbody>
</table></div>
</section>
<p class="footer">Generated by run_sensitivity.sh - SortMeRNA database pipeline</p>
</body>
</html>"""

Path(os.environ['OUTPUT_HTML']).write_text(html)
print(f"  Done - {os.environ['OUTPUT_HTML']}")
PYEOF

echo ""
echo "============================================"
echo "Sensitivity analysis complete"
echo "============================================"
echo "  Summary: ${OUTPUT_HTML}"
echo ""
echo "Next step: Review sensitivity_summary.html and compare across database configurations."
