#!/usr/bin/env bash
# run_smr_benchmark.sh  Run SortMeRNA (default db, e=1e-5) on the 5 paired-end
# benchmark datasets from the RiboDetector paper (Zhang et al. 2022, Nucleic
# Acids Research, doi:10.1093/nar/gkac112) + oral metatranscriptome dataset (MetaT).
#
# ============================================================================
# BENCHMARK DATASETS (simulated with ART_Illumina v2.3.7: -p -l 100 -ss HS25
# -m 150 -s 10, i.e. paired-end 100 bp HiSeq 2500 error model, except
# Amplicon_16S which is real sequencing data).
#
# Test column: FN = rRNA input, measures sensitivity (false negative rate);
#              FP = non-rRNA input, measures specificity (false positive rate);
#              FN+FP = mixed rRNA and non-rRNA.
#
# Dataset       Test   Pairs          Description
# -----------   -----  -------------  -----------------------------------------
# SILVA_rRNA    FN     20,000,000     SILVA SSU+LSU rRNA sequences
# OMA_CDS       FP     20,000,000     prokaryotic and eukaryotic mRNA
# oma_silva     FP     1,027,675      OMA mRNA CDSs with >=70% identity to rRNA genes;
#                                     estimates FPR on rRNA-similar mRNA
# homd_fp       FP     100,558        HOMD oral-microbe mRNA CDSs with >=70% identity to
#                                     high-FPR (>=0.5) rRNA hits in OMA_CDS
# ENA_virus     FP     27,206,792     Viral gene sequences from ENA
# Amplicon_16S  FN     7,917,920      Real 16S V1-V2 amplicon reads (oral microbiome study)
# Human_ncRNA   FP     6,330,381      Human non-coding RNA
# MetaT         FN+FP  9,165,829      Oral metatranscriptome: 4.7M prokaryotic mRNA,
#                                     2.5M human mRNA, 73K viral mRNA, 1.9M rRNA (21% rRNA fraction)
# ============================================================================
#
# Misclassifications:
#   rRNA datasets     - read pairs NOT classified as rRNA (false negatives)
#   non-rRNA datasets - read pairs incorrectly classified as rRNA (false positives)
#   mixed (MetaT)     - reported as NA
#
# IO (MB) is measured from /proc/PID/io (read_bytes + write_bytes). Requires
# Linux; on other systems IO is reported as NA.
#
# Usage:
#   bash run_smr_benchmark.sh <data_dir> <output_dir> [OPTIONS]
#
# Arguments:
#   data_dir     Directory containing paired read files named <dataset>.1.fq
#                and <dataset>.2.fq for each dataset listed above.
#   output_dir   Root directory for all outputs.
#
# Options:
#   --threads INT   Number of threads (default: 4).
#
# Required env vars:
#   SMR_BIN         Full path to the SortMeRNA binary.
#   INDEX_DIR       Directory containing SortMeRNA index subdirectories.
#   SMR_DB_ROOT_DIR Root of the sortmerna-database repository.
#
# Outputs under <output_dir>/<dataset>/smr/:
#   out/aligned.log    SortMeRNA summary log.
#   out/aligned.blast  BLAST-format alignments.
#   runtime_seconds.txt
#   peak_rss_mb.txt
#   io_mb.txt
# <output_dir>/summary.tsv

set -euo pipefail

DATA_DIR="$1"
OUTPUT_DIR="$2"
THREADS=4

shift 2 || shift $#
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

: "${SMR_BIN:?SMR_BIN env var not set}"
: "${INDEX_DIR:?INDEX_DIR env var not set}"
: "${SMR_DB_ROOT_DIR:?SMR_DB_ROOT_DIR env var not set}"

SMR_VERSION=$("${SMR_BIN}" --version 2>&1 | grep "^SortMeRNA version" | awk '{print $3}')
DB_CONFIG="smr_v${SMR_VERSION}_default_db"
DB_FASTA="${INDEX_DIR}/${DB_CONFIG}/${DB_CONFIG}.fasta"
DB_IDX="${INDEX_DIR}/${DB_CONFIG}/idx"

DATASET_NAMES=(oma_silva homd_fp OMA_CDS SILVA_rRNA Amplicon_16S Human_ncRNA MetaT ENA_virus)
declare -A DATASET_CLASS=(
    [SILVA_rRNA]=rrna     [OMA_CDS]=nonrrna    [ENA_virus]=nonrrna
    [Amplicon_16S]=rrna   [Human_ncRNA]=nonrrna [MetaT]=mixed
    [oma_silva]=nonrrna   [homd_fp]=nonrrna
)
declare -A DATASET_FILE=(
    [MetaT]=oral_metat.prokar_human_virus_mrna.rrna
)

mkdir -p "${OUTPUT_DIR}"

# _run_timed <runtime_file> <rss_file> <io_file> <cmd...>
# Runs cmd in background, tracks wall-clock time, peak RSS, and total IO (MB).
# IO tracking uses /proc/PID/io (Linux only); io_file gets NA on other systems.
_tree_rss_kb() {
    ps -e -o pid=,ppid=,rss= 2>/dev/null | awk -v root="$1" '
        { pid[NR]=$1; ppid[NR]=$2; rss[NR]=$3 }
        END {
            in_tree[root]=1; changed=1
            while (changed) {
                changed=0
                for (i=1;i<=NR;i++) {
                    if (!in_tree[pid[i]] && in_tree[ppid[i]]) {
                        in_tree[pid[i]]=1; changed=1
                    }
                }
            }
            total=0
            for (i=1;i<=NR;i++) if (in_tree[pid[i]]) total+=rss[i]
            print total
        }'
}

_run_timed() {
    local rt_f="$1" rss_f="$2" io_f="$3"
    shift 3
    local start end peak_rss_mb=0 last_io_mb=0 rss_kb rss_mb rb wb
    start=$(date +%s)
    "$@" &
    local pid=$!
    while kill -0 "${pid}" 2>/dev/null; do
        rss_kb=$(_tree_rss_kb "${pid}")
        rss_mb=$(( (rss_kb + 0) / 1024 ))
        (( rss_mb > peak_rss_mb )) && peak_rss_mb=${rss_mb}
        if [[ -f "/proc/${pid}/io" ]]; then
            rb=$(awk '/^read_bytes:/{print $2}'  "/proc/${pid}/io" 2>/dev/null || echo 0)
            wb=$(awk '/^write_bytes:/{print $2}' "/proc/${pid}/io" 2>/dev/null || echo 0)
            last_io_mb=$(( (rb + wb) / 1024 / 1024 ))
        fi
        sleep 5
    done
    wait "${pid}" || { echo "  ERROR: command failed"; exit 1; }
    end=$(date +%s)
    echo $(( end - start )) > "${rt_f}"
    echo "${peak_rss_mb}"   > "${rss_f}"
    if [[ -f "/proc/1/io" ]]; then
        echo "${last_io_mb}" > "${io_f}"
    else
        echo "NA" > "${io_f}"
    fi
}

echo "============================================"
echo "SortMeRNA benchmark"
echo "============================================"
echo "  Data dir:    ${DATA_DIR}"
echo "  Output:      ${OUTPUT_DIR}"
echo "  Threads:     ${THREADS}"
echo "  SMR DB:      ${DB_CONFIG}"
echo "  SMR E-value: 1e-5"
echo ""

summary_tsv="${OUTPUT_DIR}/summary.tsv"
printf 'Method\tDataset\tSequence_type\tTotal_sequences\tClassified_rRNA\tMisclassifications\tWall_time_s\tIO_MB\tMemory_MB\n' \
    > "${summary_tsv}"

for dataset in "${DATASET_NAMES[@]}"; do
    pfx="${DATASET_FILE[${dataset}]:-${dataset}}"
    R1="${DATA_DIR}/${pfx}.1.fq"
    R2="${DATA_DIR}/${pfx}.2.fq"

    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        echo "WARNING: skipping ${dataset} - ${pfx}.1.fq / .2.fq not found in ${DATA_DIR}" >&2
        continue
    fi

    out_dir="${OUTPUT_DIR}/${dataset}"
    mkdir -p "${out_dir}"
    class="${DATASET_CLASS[${dataset}]}"

    echo "--------------------------------------------"
    echo "Dataset: ${dataset} [${class}]"
    echo "--------------------------------------------"

    total_pairs=$(seqkit stats -T "${R1}" | tail -1 | cut -f4)

    smr_dir="${out_dir}/smr"
    aligned_log="${smr_dir}/out/aligned.log"

    if [[ ! -f "${aligned_log}" ]]; then
        echo "  Running SortMeRNA..."
        echo "  CMD: conda run -n sortmerna-bench ${SMR_BIN} --ref ${DB_FASTA} --reads ${R1} --reads ${R2} --idx-dir ${DB_IDX} --workdir ${smr_dir} --fastx --blast 1 --paired_out --zip-out 1 --threads ${THREADS} -e 1e-5"
        mkdir -p "${smr_dir}"
        _run_timed \
            "${smr_dir}/runtime_seconds.txt" \
            "${smr_dir}/peak_rss_mb.txt" \
            "${smr_dir}/io_mb.txt" \
            conda run -n sortmerna-bench "${SMR_BIN}" \
                --ref     "${DB_FASTA}" \
                --reads   "${R1}" \
                --reads   "${R2}" \
                --idx-dir "${DB_IDX}" \
                --workdir "${smr_dir}" \
                --fastx --blast 1 --paired_out \
                --zip-out 1 \
                --threads "${THREADS}" \
                -e 1e-5
        echo "  Done: $(cat "${smr_dir}/runtime_seconds.txt")s - peak RSS: $(cat "${smr_dir}/peak_rss_mb.txt") MB"
    else
        echo "  SortMeRNA: already done"
    fi

    smr_total=$( grep "Total reads = "                        "${aligned_log}" | grep -oE '[0-9]+' | head -1)
    smr_rrna=$(  grep "Total reads passing E-value threshold" "${aligned_log}" | grep -oE '[0-9]+' | head -1)
    smr_pairs=$(( smr_rrna / 2 ))
    if [[ "${class}" == "rrna" ]]; then
        smr_misclass=$(( total_pairs - smr_pairs ))
    elif [[ "${class}" == "nonrrna" ]]; then
        smr_misclass=${smr_pairs}
    else
        smr_misclass="NA"
    fi
    smr_rt=$( cat "${smr_dir}/runtime_seconds.txt" 2>/dev/null || echo NA)
    smr_io=$( cat "${smr_dir}/io_mb.txt"           2>/dev/null || echo NA)
    smr_ram=$(cat "${smr_dir}/peak_rss_mb.txt"     2>/dev/null || echo NA)

    printf 'SortMeRNA\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${dataset}" "${class}" "${total_pairs}" \
        "${smr_pairs}" "${smr_misclass}" "${smr_rt}" "${smr_io}" "${smr_ram}" \
        >> "${summary_tsv}"

    echo ""
done

echo "Done. Summary: ${summary_tsv}"

# --- HTML summary ---
OUTPUT_HTML="${OUTPUT_DIR}/summary.html"
SMR_VERSION="${SMR_VERSION}" DB_CONFIG="${DB_CONFIG}" \
python3 - "${summary_tsv}" "${OUTPUT_HTML}" <<'PYEOF'
import sys, csv, datetime, os

tsv_path, html_path = sys.argv[1], sys.argv[2]
smr_version = os.environ['SMR_VERSION']
db_config   = os.environ['DB_CONFIG']
date_str    = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

rows = ""
with open(tsv_path) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for r in reader:
        cls   = r['Sequence_type']
        total = r['Total_sequences']
        mis   = r['Misclassifications']
        if cls == 'rrna' and mis != 'NA' and total != 'NA':
            rate = f"{(1 - int(mis)/int(total))*100:.2f}%"
            metric = 'Sensitivity'
        elif cls == 'nonrrna' and mis != 'NA' and total != 'NA':
            rate = f"{int(mis)/int(total)*100:.2f}%"
            metric = 'FPR'
        else:
            rate, metric = 'NA', 'NA'

        classified = r['Classified_rRNA']
        if cls == 'mixed' and classified != 'NA' and total != 'NA':
            pct = f"{int(classified)/int(total)*100:.1f}%"
            mis_cell = f"{classified} ({pct} of reads classified as rRNA; ~21% expected)"
        else:
            mis_cell = mis

        rows += (
            f"<tr>"
            f"<td>{r['Dataset']}</td>"
            f"<td>{cls}</td>"
            f"<td>{total}</td>"
            f"<td>{mis_cell}</td>"
            f"<td>{metric}</td>"
            f"<td>{rate}</td>"
            f"<td>{r['Wall_time_s']}</td>"
            f"<td>{r['Memory_MB']}</td>"
            f"</tr>\n"
        )

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>SortMeRNA Benchmark Summary</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 40px; color: #333; }}
  h1   {{ color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 8px; }}
  .meta   {{ font-size: 0.9em; color: #666; margin-bottom: 1.2em; }}
  .description {{ background: #f8f9fa; border-left: 4px solid #3498db; padding: 10px 16px; margin: 12px 0; }}
  .table-wrap {{ overflow-x: auto; margin: 16px 0; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.88em; }}
  th, td {{ border: 1px solid #ddd; padding: 8px 14px; text-align: left; }}
  th {{ background: #2c3e50; color: #fff; }}
  tr:nth-child(even) {{ background: #f2f2f2; }}
  .footer {{ font-size: 0.8em; color: #999; margin-top: 40px; border-top: 1px solid #ddd; padding-top: 8px; }}
</style>
</head>
<body>
<h1>SortMeRNA Benchmark Summary</h1>
<p class="meta">Generated: {date_str} | SortMeRNA v{smr_version} | Database: {db_config} | E-value: 1e-5 | AWS EC2 r6i.16xlarge (64 vCPUs, 512 GB RAM), 40 threads</p>
<div class="description">
<p><b>Datasets:</b> Benchmark datasets from Deng et al. 2022 (<i>Nucleic Acids Research</i>, <a href="https://doi.org/10.1093/nar/gkac112">doi:10.1093/nar/gkac112</a>).
FN = rRNA input (measures sensitivity); FP = non-rRNA input (measures specificity); FN+FP = mixed.</p>
<table style="margin-top:8px;font-size:0.85em;width:auto">
<thead><tr><th>Dataset</th><th>Test</th><th>Pairs</th><th>Description</th></tr></thead>
<tbody>
<tr><td>SILVA_rRNA</td><td>FN</td><td>20,000,000</td><td>SILVA SSU+LSU rRNA sequences</td></tr>
<tr><td>OMA_CDS</td><td>FP</td><td>20,000,000</td><td>Prokaryotic and eukaryotic mRNA</td></tr>
<tr><td>oma_silva</td><td>FP</td><td>1,027,675</td><td>OMA mRNA CDSs with &ge;70% identity to rRNA genes; FPR on rRNA-similar mRNA</td></tr>
<tr><td>homd_fp</td><td>FP</td><td>100,558</td><td>HOMD oral-microbe mRNA CDSs with &ge;70% identity to high-FPR rRNA hits in OMA_CDS</td></tr>
<tr><td>ENA_virus</td><td>FP</td><td>27,206,792</td><td>Viral gene sequences from ENA</td></tr>
<tr><td>Amplicon_16S</td><td>FN</td><td>7,917,920</td><td>Real 16S V1-V2 amplicon reads (oral microbiome)</td></tr>
<tr><td>Human_ncRNA</td><td>FP</td><td>6,330,381</td><td>Human non-coding RNA</td></tr>
<tr><td>MetaT</td><td>FN+FP</td><td>9,165,829</td><td>Oral metatranscriptome: 4.7M prokaryotic mRNA, 2.5M human mRNA, 73K viral mRNA, 1.9M rRNA (21% rRNA fraction)</td></tr>
</tbody>
</table>
<p style="margin-top:10px"><b>Metrics:</b> Sensitivity = (total - misclassifications) / total for FN datasets.
FPR = misclassifications / total for FP datasets.
MetaT reports reads classified as rRNA vs. the expected ~21% rRNA fraction.</p>
</div>
<div class="table-wrap"><table>
<thead><tr>
  <th>Dataset</th><th>Type</th><th>Total pairs</th>
  <th>Misclassifications / rRNA classified</th><th>Metric</th><th>Value</th>
  <th>Wall time (s)</th><th>Memory (MB)</th>
</tr></thead>
<tbody>
{rows}</tbody>
</table></div>
<p class="footer">Generated by run_smr_benchmark.sh</p>
</body>
</html>"""

with open(html_path, 'w') as f:
    f.write(html)
print(f"  HTML summary: {html_path}")
PYEOF
