#!/usr/bin/env bash
# run_pacbio_sweep.sh  SortMeRNA --passes / --num_seeds parameter sweep for PacBio HiFi reads.
#
# Tests 21 values of --num_seeds against real rRNA reads
# (Karst et al. 2021, AGP, Qiita study 10317, ~4,500 bp 16S+ITS+23S operons) and
# PBSIM3-simulated non-rRNA reads from the masked T2T genome. Each combination
# reports sensitivity and FPR. Run SortMeRNA with -e 1e-5 (default db).
#
# Usage:
#   bash run_pacbio_sweep.sh <rrna_reads> <nonrrna_reads> <output_dir> [threads]
#
# Arguments:
#   rrna_reads     Real PacBio rRNA reads - FASTA/FASTQ, optionally gzipped
#                  ($PACBIO_DIR/karst2021_253k.fna.gz)
#   nonrrna_reads  PBSIM3-simulated non-rRNA reads - FASTA/FASTQ, optionally gzipped
#                  ($NON_RRNA_DIR/non_rrna_pacbio_253089_T2T.fastq.gz)
#   output_dir     Directory for all sweep outputs ($PACBIO_DIR/sweep)
#   threads        Number of threads (default: 4)
#
# Required env vars:  SMR_BIN, INDEX_DIR, SMR_VERSION
#
# Outputs under <output_dir>:
#   sweep_results.tsv               sensitivity, FPR, wall time per combination
#   <label>/rrna/out/aligned.log    SortMeRNA log per combination
#   <label>/nonrrna/out/aligned.log

set -euo pipefail

RRNA_READS="$1"
NONRRNA_READS="$2"
SWEEP_DIR="$3"
THREADS="${4:-4}"

SMR_BIN="${SMR_BIN:?Please set SMR_BIN (see README Set paths section)}"
INDEX_DIR="${INDEX_DIR:?Please set INDEX_DIR (see README Set paths section)}"
SMR_VERSION="${SMR_VERSION:?Please set SMR_VERSION (see README Set paths section)}"

DB_NAME="smr_v${SMR_VERSION}_default_db"
REF_DB="${INDEX_DIR}/${DB_NAME}/${DB_NAME}.fasta"
IDX_DIR="${INDEX_DIR}/${DB_NAME}/idx"

# Total reads in each dataset - used to compute sensitivity and FPR
TOTAL_RRNA=253089
TOTAL_NONRRNA=253089

# Quick-check mode: subsample to 10K reads for fast parameter iteration.
# Disable by setting QUICK=0.
QUICK="${QUICK:-1}"
QUICK_N=10000
if [[ "${QUICK}" == "1" ]]; then
    echo "Quick mode: subsampling to ${QUICK_N} reads..."
    RRNA_10K="${SWEEP_DIR}/rrna_10k.fasta.gz"
    NONRRNA_10K="${SWEEP_DIR}/nonrrna_10k.fastq.gz"
    [[ ! -f "${RRNA_10K}" ]]    && seqkit sample -n "${QUICK_N}" --rand-seed 42 "${RRNA_READS}"    | gzip -c > "${RRNA_10K}"
    [[ ! -f "${NONRRNA_10K}" ]] && seqkit sample -n "${QUICK_N}" --rand-seed 42 "${NONRRNA_READS}" | gzip -c > "${NONRRNA_10K}"
    RRNA_READS="${RRNA_10K}"
    NONRRNA_READS="${NONRRNA_10K}"
    TOTAL_RRNA=${QUICK_N}
    TOTAL_NONRRNA=${QUICK_N}
fi

echo "============================================"
echo "PacBio parameter sweep"
echo "============================================"
echo "  rRNA reads:     ${RRNA_READS}"
echo "  Non-rRNA reads: ${NONRRNA_READS}"
echo "  Output:         ${SWEEP_DIR}"
echo "  DB:             ${REF_DB}"
echo "  Index:          ${IDX_DIR}"
echo "  Threads:        ${THREADS}"
echo ""

mkdir -p "${SWEEP_DIR}"
RESULTS="${SWEEP_DIR}/sweep_results.tsv"
printf 'passes\tnum_seeds\trrna_aligned\tsensitivity\tnonrrna_aligned\tfpr\twall_sec\tpeak_rss_mb\n' > "${RESULTS}"

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

# Runs SortMeRNA in background, tracks wall time and peak RSS.
# Writes wall_sec and peak_rss_mb to stdout as "WALL=<n> RSS=<n>".
run_smr_timed() {
    local reads="$1" workdir="$2" passes="$3" num_seeds="$4"
    local start peak_rss_mb=0 rss_kb rss_mb
    start=$(date +%s)
    "${SMR_BIN}" \
        --ref "${REF_DB}" \
        --reads "${reads}" \
        --workdir "${workdir}" \
        --idx-dir "${IDX_DIR}" \
        --passes "${passes}" \
        --num_seeds "${num_seeds}" \
        --threads "${THREADS}" \
        --fastx --blast 1 \
        -e 1e-5 &
    local pid=$!
    while kill -0 "${pid}" 2>/dev/null; do
        rss_kb=$(_tree_rss_kb "${pid}")
        rss_mb=$(( (rss_kb + 0) / 1024 ))
        (( rss_mb > peak_rss_mb )) && peak_rss_mb=${rss_mb}
        sleep 5
    done
    wait "${pid}" || { echo "  ERROR: sortmerna failed" >&2; exit 1; }
    echo "WALL=$(( $(date +%s) - start )) RSS=${peak_rss_mb}"
}

# NOTE: --passes is currently ignored by SortMeRNA v6.0.2 (always runs default 18,9,3).
# Sweep is over --num_seeds only until --passes is fixed.
PASSES_LIST=("18,9,3")
SEEDS_LIST=(1000 900 800 700 600 500 400 300 200 100 50 45 40 35 30 25 20 15 10 5 2)

for passes in "${PASSES_LIST[@]}"; do
    for num_seeds in "${SEEDS_LIST[@]}"; do
        label="p${passes//,/_}_s${num_seeds}"
        rrna_log="${SWEEP_DIR}/${label}/rrna/out/aligned.log"
        nonrrna_log="${SWEEP_DIR}/${label}/nonrrna/out/aligned.log"

        echo "--------------------------------------------"
        echo "passes=${passes}  num_seeds=${num_seeds}"
        echo "--------------------------------------------"

        wall_rrna=0 rss_rrna=0 wall_nonrrna=0 rss_nonrrna=0

        if [[ -f "${rrna_log}" ]]; then
            echo "  rRNA run already exists - skipping"
        else
            mkdir -p "${SWEEP_DIR}/${label}/rrna"
            echo "  Running rRNA..."
            result=$(run_smr_timed "${RRNA_READS}" "${SWEEP_DIR}/${label}/rrna" "${passes}" "${num_seeds}")
            wall_rrna=$(echo "${result}" | grep -oP '(?<=WALL=)\d+')
            rss_rrna=$(echo "${result}"  | grep -oP '(?<=RSS=)\d+')
            echo "  rRNA done: ${wall_rrna}s peak RSS ${rss_rrna} MB"
        fi

        if [[ -f "${nonrrna_log}" ]]; then
            echo "  Non-rRNA run already exists - skipping"
        else
            mkdir -p "${SWEEP_DIR}/${label}/nonrrna"
            echo "  Running non-rRNA..."
            result=$(run_smr_timed "${NONRRNA_READS}" "${SWEEP_DIR}/${label}/nonrrna" "${passes}" "${num_seeds}")
            wall_nonrrna=$(echo "${result}" | grep -oP '(?<=WALL=)\d+')
            rss_nonrrna=$(echo "${result}"  | grep -oP '(?<=RSS=)\d+')
            echo "  Non-rRNA done: ${wall_nonrrna}s peak RSS ${rss_nonrrna} MB"
        fi

        wall=$(( wall_rrna + wall_nonrrna ))
        peak_rss=$(( rss_rrna > rss_nonrrna ? rss_rrna : rss_nonrrna ))

        rrna_aligned=$(grep -oP '(?<=Total reads passing E-value threshold = )\d+' "${rrna_log}")
        nonrrna_aligned=$(grep -oP '(?<=Total reads passing E-value threshold = )\d+' "${nonrrna_log}")
        sens=$(awk "BEGIN { printf \"%.4f\", ${rrna_aligned} / ${TOTAL_RRNA} }")
        fpr=$(awk "BEGIN { printf \"%.4f\", ${nonrrna_aligned} / ${TOTAL_NONRRNA} }")

        printf '%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\n' \
            "${passes}" "${num_seeds}" \
            "${rrna_aligned}" "${sens}" \
            "${nonrrna_aligned}" "${fpr}" \
            "${wall}" "${peak_rss}" \
            >> "${RESULTS}"

        echo "  sensitivity=${sens} (${rrna_aligned}/${TOTAL_RRNA})  fpr=${fpr} (${nonrrna_aligned}/${TOTAL_NONRRNA})  time=${wall}s  peak_rss=${peak_rss}MB"
    done
done

echo ""
echo "=== Sweep complete ==="
column -t -s $'\t' "${RESULTS}"

# --- ROC plot ---
PLOT_PNG="${SWEEP_DIR}/roc_num_seeds.png"
python3 - "${RESULTS}" "${PLOT_PNG}" <<'PYEOF'
import sys, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

tsv_path, plot_path = sys.argv[1], sys.argv[2]

rows = []
with open(tsv_path) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for r in reader:
        rows.append({
            'num_seeds': int(r['num_seeds']),
            'sensitivity': float(r['sensitivity']),
            'selectivity': 1.0 - float(r['fpr']),
        })

rows.sort(key=lambda r: r['num_seeds'])
seeds    = [r['num_seeds']    for r in rows]
sens     = [r['sensitivity']  for r in rows]
sel      = [r['selectivity']  for r in rows]

fig, ax = plt.subplots(figsize=(8, 6))

cmap = cm.get_cmap('viridis', len(rows))
for i, r in enumerate(rows):
    ax.scatter(r['selectivity'], r['sensitivity'], color=cmap(i), s=60, zorder=3)
    ax.annotate(str(r['num_seeds']),
                xy=(r['selectivity'], r['sensitivity']),
                xytext=(4, 4), textcoords='offset points',
                fontsize=7, color=cmap(i))

ax.plot(sel, sens, color='grey', linewidth=0.8, linestyle='--', zorder=2)

ax.set_xlabel('Selectivity (1 - FPR)', fontsize=12)
ax.set_ylabel('Sensitivity (TPR)', fontsize=12)
ax.set_title('SortMeRNA PacBio: num_seeds sweep\n(labels = num_seeds value)', fontsize=11)
ax.set_xlim(0, 1.02)
ax.set_ylim(0, 1.02)
ax.axhline(1.0, color='green', linewidth=0.6, linestyle=':')
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(plot_path, dpi=150)
print(f"  ROC plot: {plot_path}")
PYEOF
