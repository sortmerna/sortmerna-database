#!/usr/bin/env bash

################################################################################
# simulate_rrna_reads.sh
#
# Simulate Illumina rRNA reads from non-seed cluster members for three database
# configurations and produce one pooled FASTA per Set, used in Experiment 1
# (scalability, rRNA sensitivity) and Experiment 2 (sensitivity across db configs).
#
# Each Set uses non-seed members from the matching database configuration so that
# simulated reads are real rRNA sequences absent from the database being tested:
#
#   Set 1 (sensitive, 97%):   all types at 97% non-seed members
#   Set 2 (default, 90-95%):  bacteria SSU at 90%, all others at 95%; Rfam at 97%
#   Set 3 (fast, 85-90%):     bacteria SSU at 85%, all others at 90%; Rfam at 97%
#
# Rfam 5S and 5.8S use 97% members for all three Sets because the default and fast
# databases use seed-only Rfam (no threshold-based clustering), so 97% members are
# absent from all three databases and serve as a consistent rRNA specificity signal.
#
# Outputs (all under <output_dir>):
#   set<N>/<type>/iss_R1.fastq   - InSilicoSeq R1 output
#   set<N>/<type>/iss_R2.fastq   - InSilicoSeq R2 output
#   set<N>/<type>/reads.fasta    - per-type simulated reads (FASTA, R1+R2 concatenated)
#   set<N>_rrna_reads.fasta      - pooled reads across all rRNA types for Set N
#   rrna_simulation_summary.html - HTML summary report
#
# Usage: bash simulate_rrna_reads.sh [output_dir [threads]] [OPTIONS]
#
# Positional:
#   output_dir        Where to write simulated reads (default: $RRNA_SIM_DIR or data/rrna_sim)
#   threads           Threads for InSilicoSeq (default: 4)
#
# Options:
#   --reads-per-type INT  Target reads per rRNA type per Set (default: 12500)
#   --model STR           InSilicoSeq error model: HiSeq, NovaSeq, MiSeq (default: NovaSeq)
#   --seed INT            Random seed (default: 42)
#   --clustered-dir DIR   Directory containing *_test_members.fasta files
#                         (overrides CLUSTERED_DIR env var)
#   --force               Re-run InSilicoSeq even if output already exists
#   -h, --help            Show this help
#
# Required env vars (can be overridden with options above):
#   CLUSTERED_DIR    Directory with *_test_members.fasta files from cluster_sequences.sh
#
# Environment variables (optional):
#   RRNA_SIM_DIR     Default output directory (overridden by positional arg)
#
# Requires: iss (InSilicoSeq), seqkit
#
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

POSITIONAL=()
N_READS_PER_TYPE=12500
ISS_MODEL=NovaSeq
RAND_SEED=42
CLUSTERED_DIR_OPT=""
FORCE=false

show_help() {
    grep '^#' "$0" | grep -v '#!/usr/bin/env bash' | sed 's/^# \?//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --reads-per-type) N_READS_PER_TYPE="$2"; shift 2 ;;
        --model)          ISS_MODEL="$2";         shift 2 ;;
        --seed)           RAND_SEED="$2";         shift 2 ;;
        --clustered-dir)  CLUSTERED_DIR_OPT="$2"; shift 2 ;;
        --force)          FORCE=true;             shift   ;;
        -h|--help)        show_help ;;
        *)                POSITIONAL+=("$1");     shift   ;;
    esac
done

OUTPUT_DIR="${POSITIONAL[0]:-${RRNA_SIM_DIR:-data/rrna_sim}}"
THREADS="${POSITIONAL[1]:-4}"

[[ -n "${CLUSTERED_DIR_OPT}" ]] && CLUSTERED_DIR="${CLUSTERED_DIR_OPT}"
CLUSTERED_DIR="${CLUSTERED_DIR:?CLUSTERED_DIR not set - pass --clustered-dir or export CLUSTERED_DIR}"

OUTPUT_DIR="$(mkdir -p "${OUTPUT_DIR}" && cd "${OUTPUT_DIR}" && pwd)"

echo "============================================"
echo "rRNA Read Simulation Script"
echo "Clustered dir:    ${CLUSTERED_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Reads per type:   ${N_READS_PER_TYPE}"
echo "ISS model:        ${ISS_MODEL}"
echo "Random seed:      ${RAND_SEED}"
echo "Threads:          ${THREADS}"
echo "Force re-run:     ${FORCE}"
echo "============================================"
echo ""

for tool in iss seqkit; do
    if ! command -v "${tool}" &>/dev/null; then
        echo "Error: ${tool} not found. Install via: conda env create -f environment.yml"
        exit 1
    fi
done

# Ordered list of rRNA types - consistent across all Sets
RRNA_TYPES=(
    silva_ssu_bacteria
    silva_ssu_archaea
    silva_ssu_eukaryota
    silva_lsu_bacteria
    silva_lsu_archaea
    silva_lsu_eukaryota
    rfam_5s
    rfam_5_8s
)

# Minimum source sequences required to run ISS; types below this are skipped
MIN_SOURCE_SEQS=10

# ---------------------------------------------------------------------------
# Source file lookup: set_source <set_num> <type_name> -> prints path
# Bacteria SSU uses a one-step lower threshold in default/fast to match
# build_sortmerna_index.sh configuration.
# Rfam uses 97% members for all Sets (see header comment).
# ---------------------------------------------------------------------------
set_source() {
    local set_num="$1"
    local type="$2"
    local C="${CLUSTERED_DIR}"
    case "${set_num}_${type}" in
        # Set 1 - sensitive (97% for all SILVA and Rfam)
        1_silva_ssu_bacteria)  echo "${C}/silva_ssu_bacteria_97_test_members.fasta" ;;
        1_silva_ssu_archaea)   echo "${C}/silva_ssu_archaea_97_test_members.fasta"  ;;
        1_silva_ssu_eukaryota) echo "${C}/silva_ssu_eukaryota_97_test_members.fasta" ;;
        1_silva_lsu_bacteria)  echo "${C}/silva_lsu_bacteria_97_test_members.fasta" ;;
        1_silva_lsu_archaea)   echo "${C}/silva_lsu_archaea_97_test_members.fasta"  ;;
        1_silva_lsu_eukaryota) echo "${C}/silva_lsu_eukaryota_97_test_members.fasta" ;;
        1_rfam_5s)             echo "${C}/rfam_5s_97_test_members.fasta"    ;;
        1_rfam_5_8s)           echo "${C}/rfam_5_8s_97_test_members.fasta"  ;;
        # Set 2 - default (bacteria SSU at 90%, others at 95%, Rfam at 97%)
        2_silva_ssu_bacteria)  echo "${C}/silva_ssu_bacteria_90_test_members.fasta" ;;
        2_silva_ssu_archaea)   echo "${C}/silva_ssu_archaea_95_test_members.fasta"  ;;
        2_silva_ssu_eukaryota) echo "${C}/silva_ssu_eukaryota_95_test_members.fasta" ;;
        2_silva_lsu_bacteria)  echo "${C}/silva_lsu_bacteria_95_test_members.fasta" ;;
        2_silva_lsu_archaea)   echo "${C}/silva_lsu_archaea_95_test_members.fasta"  ;;
        2_silva_lsu_eukaryota) echo "${C}/silva_lsu_eukaryota_95_test_members.fasta" ;;
        2_rfam_5s)             echo "${C}/rfam_5s_97_test_members.fasta"    ;;
        2_rfam_5_8s)           echo "${C}/rfam_5_8s_97_test_members.fasta"  ;;
        # Set 3 - fast (bacteria SSU at 85%, others at 90%, Rfam at 97%)
        3_silva_ssu_bacteria)  echo "${C}/silva_ssu_bacteria_85_test_members.fasta" ;;
        3_silva_ssu_archaea)   echo "${C}/silva_ssu_archaea_90_test_members.fasta"  ;;
        3_silva_ssu_eukaryota) echo "${C}/silva_ssu_eukaryota_90_test_members.fasta" ;;
        3_silva_lsu_bacteria)  echo "${C}/silva_lsu_bacteria_90_test_members.fasta" ;;
        3_silva_lsu_archaea)   echo "${C}/silva_lsu_archaea_90_test_members.fasta"  ;;
        3_silva_lsu_eukaryota) echo "${C}/silva_lsu_eukaryota_90_test_members.fasta" ;;
        3_rfam_5s)             echo "${C}/rfam_5s_97_test_members.fasta"    ;;
        3_rfam_5_8s)           echo "${C}/rfam_5_8s_97_test_members.fasta"  ;;
        *) echo "" ;;
    esac
}

# ---------------------------------------------------------------------------
# simulate_type: simulate reads for one rRNA type in one Set
# Writes per-type reads.fasta; prints n_source and n_reads to stdout as
# two tab-separated values for the caller to capture.
# ---------------------------------------------------------------------------
simulate_type() {
    local type_dir="$1"
    local type_name="$2"
    local source="$3"

    local out_fasta="${type_dir}/reads.fasta"

    # Source absent or empty
    if [[ ! -f "${source}" ]] || [[ ! -s "${source}" ]]; then
        echo "  WARNING: source not found or empty: $(basename "${source}") - skipping ${type_name}" >&2
        printf "0\t0"
        return 0
    fi

    local n_source
    n_source=$(seqkit stats -T "${source}" | tail -1 | cut -f4)

    if (( n_source < MIN_SOURCE_SEQS )); then
        echo "  WARNING: ${type_name} has only ${n_source} source sequences (< ${MIN_SOURCE_SEQS}) - skipping" >&2
        printf "%s\t0" "${n_source}"
        return 0
    fi

    # Resume: skip ISS if reads already exist
    if [[ -f "${out_fasta}" ]] && [[ "${FORCE}" == false ]]; then
        local n_reads
        n_reads=$(seqkit stats -T "${out_fasta}" | tail -1 | cut -f4)
        echo "  Already exists: ${type_name} (${n_reads} reads)" >&2
        printf "%s\t%s" "${n_source}" "${n_reads}"
        return 0
    fi

    mkdir -p "${type_dir}"

    echo "  Simulating ${N_READS_PER_TYPE} reads from ${n_source} ${type_name} members..." >&2

    local iss_prefix="${type_dir}/iss"
    iss generate \
        --genomes  "${source}" \
        --model    "${ISS_MODEL}" \
        --n_reads  "${N_READS_PER_TYPE}" \
        --cpus     "${THREADS}" \
        --seed     "${RAND_SEED}" \
        --output   "${iss_prefix}"

    cat "${iss_prefix}_R1.fastq" "${iss_prefix}_R2.fastq" \
        | seqkit fq2fa \
        | seqkit seq -w 0 \
        > "${out_fasta}"

    local n_reads
    n_reads=$(seqkit stats -T "${out_fasta}" | tail -1 | cut -f4)
    echo "    Saved: ${type_name}/reads.fasta (${n_reads} reads)" >&2
    printf "%s\t%s" "${n_source}" "${n_reads}"
}

# ---------------------------------------------------------------------------
# process_set: simulate all types for a Set and pool into one FASTA
# ---------------------------------------------------------------------------
process_set() {
    local set_num="$1"
    local set_dir="${OUTPUT_DIR}/set${set_num}"
    local pooled="${OUTPUT_DIR}/set${set_num}_rrna_reads.fasta"

    echo ""
    echo "============================================"
    echo "Set ${set_num}"
    echo "============================================"

    mkdir -p "${set_dir}"
    > "${pooled}"

    local html_rows=""
    local total_reads=0

    for type_name in "${RRNA_TYPES[@]}"; do
        local source
        source=$(set_source "${set_num}" "${type_name}")

        local counts
        counts=$(simulate_type "${set_dir}/${type_name}" "${type_name}" "${source}")
        local n_source n_reads
        n_source=$(echo "${counts}" | cut -f1)
        n_reads=$(echo "${counts}" | cut -f2)

        local type_fasta="${set_dir}/${type_name}/reads.fasta"
        if [[ -f "${type_fasta}" ]] && [[ -s "${type_fasta}" ]]; then
            cat "${type_fasta}" >> "${pooled}"
        fi

        total_reads=$(( total_reads + n_reads ))
        html_rows="${html_rows}      <tr><td>${type_name//_/ }</td><td>${n_source}</td><td>$(basename "${source}")</td><td>${n_reads}</td></tr>\n"
    done

    local n_pooled
    n_pooled=$(seqkit stats -T "${pooled}" | tail -1 | cut -f4)
    echo "  Pooled: set${set_num}_rrna_reads.fasta (${n_pooled} reads total)"

    # Save per-set data for HTML report
    printf "%s" "${html_rows}" > "${set_dir}/.html_rows.tmp"
    printf "%s" "${n_pooled}"  > "${set_dir}/.total.tmp"
}

process_set 1
process_set 2
process_set 3

################################################################################
# HTML SUMMARY REPORT
################################################################################

echo ""
echo "============================================"
echo "Writing HTML summary"
echo "============================================"

OUTPUT_HTML="${OUTPUT_DIR}/rrna_simulation_summary.html"

N_READS_PER_TYPE="${N_READS_PER_TYPE}" \
ISS_MODEL="${ISS_MODEL}" \
RAND_SEED="${RAND_SEED}" \
OUTPUT_HTML="${OUTPUT_HTML}" \
python3 - "${OUTPUT_DIR}" <<'PYEOF'
import sys, os, datetime
from pathlib import Path

output_dir = Path(sys.argv[1])
e = os.environ
n_per_type = e['N_READS_PER_TYPE']
iss_model  = e['ISS_MODEL']
rand_seed  = e['RAND_SEED']
date_str   = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

set_labels = {
    "1": "Set 1 - sensitive db (97% threshold for all types)",
    "2": "Set 2 - default db (bacteria SSU 90%, others 95%, Rfam 97%)",
    "3": "Set 3 - fast db (bacteria SSU 85%, others 90%, Rfam 97%)",
}

sections = ""
for sn in ["1", "2", "3"]:
    rows_file  = output_dir / f"set{sn}" / ".html_rows.tmp"
    total_file = output_dir / f"set{sn}" / ".total.tmp"
    if not rows_file.exists():
        continue
    rows  = rows_file.read_text().replace('\\n', '\n')
    total = total_file.read_text().strip() if total_file.exists() else "?"
    sections += f"""
<section>
<h2>{set_labels[sn]}</h2>
<div class="description">
<p>Source: non-seed cluster members at the matched threshold. Reads simulated with InSilicoSeq
({iss_model} model, 150 bp paired-end, seed {rand_seed}). Target: {n_per_type} reads per rRNA type (R1+R2 combined).
Types with fewer than 10 source sequences are skipped.</p>
</div>
<div class="table-wrap"><table>
  <thead>
    <tr><th>rRNA type</th><th>Source sequences</th><th>Source file</th><th>Reads simulated</th></tr>
  </thead>
  <tbody>
{rows}    <tr style="font-weight:bold; border-top: 2px solid #2c3e50;"><td>Total</td><td></td><td></td><td>{total}</td></tr>
  </tbody>
</table></div>
<p>Output: <code>set{sn}_rrna_reads.fasta</code></p>
</section>
"""

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>rRNA Read Simulation Summary</title>
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
<h1>rRNA Read Simulation Summary</h1>
<p class="meta">Generated: {date_str}</p>
{sections}
<p class="footer">Generated by simulate_rrna_reads.sh - SortMeRNA database pipeline</p>
</body>
</html>"""

Path(e['OUTPUT_HTML']).write_text(html)
print(f"  Done - {e['OUTPUT_HTML']}")
PYEOF

# Clean up temp files
rm -f "${OUTPUT_DIR}"/set*/.html_rows.tmp "${OUTPUT_DIR}"/set*/.total.tmp

echo ""
echo "============================================"
echo "Simulation complete"
echo "============================================"
echo ""
echo "Outputs:"
for sn in 1 2 3; do
    f="${OUTPUT_DIR}/set${sn}_rrna_reads.fasta"
    [[ -f "${f}" ]] && echo "  Set ${sn}: ${f}"
done
echo "  HTML summary: ${OUTPUT_DIR}/rrna_simulation_summary.html"
echo ""
echo "Next steps:"
echo "  Experiment 1 - scalability (rRNA sensitivity at scale):"
echo "    bash \$SMR_DB_ROOT_DIR/scripts/benchmarking/run_scalability.sh \\"
echo "        ${OUTPUT_DIR}/set2_rrna_reads.fasta \\"
echo "        ${OUTPUT_DIR}/scalability_rrna \\"
echo "        ${THREADS} \\"
echo "        --index-dir \$INDEX_DIR"
echo "  Experiment 2 - sensitivity across db configurations:"
echo "    Run set1 against smr_v\${SMR_VERSION}_sensitive_db"
echo "    Run set2 against smr_v\${SMR_VERSION}_default_db"
echo "    Run set3 against smr_v\${SMR_VERSION}_fast_db"
