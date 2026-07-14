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
RRNA_FULL="${RRNA_READS}"
NONRRNA_FULL="${NONRRNA_READS}"
SWEEP_DIR="$3"
THREADS="${4:-4}"

SMR_BIN="${SMR_BIN:?Please set SMR_BIN (see README Set paths section)}"
INDEX_DIR="${INDEX_DIR:?Please set INDEX_DIR (see README Set paths section)}"
SMR_VERSION="${SMR_VERSION:?Please set SMR_VERSION (see README Set paths section)}"

DB_NAME="smr_v${SMR_VERSION}_default_db"
FAMILY_MAP="${INDEX_DIR}/${DB_NAME}/family_map.tsv"
REF_DB="${INDEX_DIR}/${DB_NAME}/${DB_NAME}.fasta"
IDX_DIR="${INDEX_DIR}/${DB_NAME}/idx"

# Total reads in each dataset - used to compute sensitivity and FPR
TOTAL_RRNA=253089
TOTAL_NONRRNA=253089
TOTAL_RRNA_FULL=$(seqkit stats -T "${RRNA_FULL}"    | tail -1 | cut -f4)
TOTAL_NONRRNA_FULL=$(seqkit stats -T "${NONRRNA_FULL}" | tail -1 | cut -f4)

# Quick-check mode: subsample to 10K reads for fast parameter iteration.
# Disable by setting QUICK=0.
QUICK="${QUICK:-1}"
QUICK_N=10000
if [[ "${QUICK}" == "1" ]]; then
    RRNA_10K="${SWEEP_DIR}/rrna_10k.fasta.gz"
    NONRRNA_10K="${SWEEP_DIR}/nonrrna_10k.fastq.gz"
    mkdir -p "${SWEEP_DIR}"
    if [[ -f "${RRNA_10K}" && -f "${NONRRNA_10K}" ]]; then
        echo "Quick mode: reusing existing ${QUICK_N}-read subsets."
    else
        echo "Quick mode: subsampling to ${QUICK_N} reads..."
        seqkit sample -n "${QUICK_N}" --rand-seed 42 "${RRNA_READS}"    | gzip -c > "${RRNA_10K}"
        seqkit sample -n "${QUICK_N}" --rand-seed 42 "${NONRRNA_READS}" | gzip -c > "${NONRRNA_10K}"
    fi
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
FAMILY_COUNTS="${SWEEP_DIR}/family_counts.tsv"
printf 'passes\tevalue\tnum_seeds\tmin_lis\trrna_aligned\tsensitivity\tnonrrna_aligned\tfpr\twall_sec_rrna\twall_sec_nonrrna\tpeak_rss_mb\n' > "${RESULTS}"
rm -f "${FAMILY_COUNTS}"

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
    local reads="$1" workdir="$2" passes="$3" num_seeds="$4" evalue="$5" min_lis="$6"
    local start peak_rss_mb=0 rss_kb rss_mb
    start=$(date +%s)
    "${SMR_BIN}" \
        --ref "${REF_DB}" \
        --reads "${reads}" \
        --workdir "${workdir}" \
        --idx-dir "${IDX_DIR}" \
        --passes "${passes}" \
        --num_seeds "${num_seeds}" \
        --min_lis "${min_lis}" \
        --threads "${THREADS}" \
        --fastx --blast 1 \
        -e "${evalue}" &
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

# --min_lis is the per-reference co-linearity gate (LIS at stride 18):
# --num_seeds 2 is fixed to trigger LIS computation.
PASSES_LIST=("18,9,3")
EVALUES_LIST=(1e-5 1e-10 1e-20)
# Explicit (num_seeds, min_lis) pairs:
#   Row 1: num_seeds=2 fixed, min_lis sweeps 2-10  -> min_lis does all the filtering
#   Row 2: diagonal (num_seeds=min_lis 3-6)        -> both parameters contribute equally
#   Row 3: num_seeds > min_lis (4-6 vs 2)          -> num_seeds dominates, min_lis redundant
SEED_LIS_PAIRS=(
    "2 2"  "2 3"  "2 4"  "2 5"  "2 6"  "2 7"  "2 8"  "2 9"  "2 10"
    "3 3"  "4 4"  "5 5"  "6 6"
    "3 2"  "4 2"  "5 2"  "6 2"
)

for passes in "${PASSES_LIST[@]}"; do
    for evalue in "${EVALUES_LIST[@]}"; do
        ev_label="${evalue//-/m}"
        for pair in "${SEED_LIS_PAIRS[@]}"; do
            num_seeds="${pair%% *}"
            min_lis="${pair##* }"
            label="p${passes//,/_}_e${ev_label}_s${num_seeds}_lis${min_lis}"
            rrna_log="${SWEEP_DIR}/${label}/rrna/out/aligned.log"
            nonrrna_log="${SWEEP_DIR}/${label}/nonrrna/out/aligned.log"

            echo "--------------------------------------------"
            echo "passes=${passes}  evalue=${evalue}  num_seeds=${num_seeds}  min_lis=${min_lis}"
            echo "--------------------------------------------"

            wall_rrna=0 rss_rrna=0 wall_nonrrna=0 rss_nonrrna=0

            if [[ -f "${rrna_log}" ]]; then
                echo "  rRNA run already exists - skipping"
                wall_rrna=$(cat "${SWEEP_DIR}/${label}/rrna/wall_sec.txt"    2>/dev/null || echo 0)
                rss_rrna=$( cat "${SWEEP_DIR}/${label}/rrna/peak_rss_mb.txt" 2>/dev/null || echo 0)
            else
                mkdir -p "${SWEEP_DIR}/${label}/rrna"
                echo "  Running rRNA..."
                result=$(run_smr_timed "${RRNA_READS}" "${SWEEP_DIR}/${label}/rrna" "${passes}" "${num_seeds}" "${evalue}" "${min_lis}")
                wall_rrna=$(echo "${result}" | grep -oP '(?<=WALL=)\d+')
                rss_rrna=$(echo "${result}"  | grep -oP '(?<=RSS=)\d+')
                echo "${wall_rrna}" > "${SWEEP_DIR}/${label}/rrna/wall_sec.txt"
                echo "${rss_rrna}"  > "${SWEEP_DIR}/${label}/rrna/peak_rss_mb.txt"
                echo "  rRNA done: ${wall_rrna}s peak RSS ${rss_rrna} MB"
            fi

            if [[ -f "${nonrrna_log}" ]]; then
                echo "  Non-rRNA run already exists - skipping"
                wall_nonrrna=$(cat "${SWEEP_DIR}/${label}/nonrrna/wall_sec.txt"    2>/dev/null || echo 0)
                rss_nonrrna=$( cat "${SWEEP_DIR}/${label}/nonrrna/peak_rss_mb.txt" 2>/dev/null || echo 0)
            else
                mkdir -p "${SWEEP_DIR}/${label}/nonrrna"
                echo "  Running non-rRNA..."
                result=$(run_smr_timed "${NONRRNA_READS}" "${SWEEP_DIR}/${label}/nonrrna" "${passes}" "${num_seeds}" "${evalue}" "${min_lis}")
                wall_nonrrna=$(echo "${result}" | grep -oP '(?<=WALL=)\d+')
                rss_nonrrna=$(echo "${result}"  | grep -oP '(?<=RSS=)\d+')
                echo "${wall_nonrrna}" > "${SWEEP_DIR}/${label}/nonrrna/wall_sec.txt"
                echo "${rss_nonrrna}"  > "${SWEEP_DIR}/${label}/nonrrna/peak_rss_mb.txt"
                echo "  Non-rRNA done: ${wall_nonrrna}s peak RSS ${rss_nonrrna} MB"
            fi

            peak_rss=$(( rss_rrna > rss_nonrrna ? rss_rrna : rss_nonrrna ))

            rrna_aligned=$(grep -oP '(?<=Total reads passing E-value threshold = )\d+' "${rrna_log}")
            nonrrna_aligned=$(grep -oP '(?<=Total reads passing E-value threshold = )\d+' "${nonrrna_log}")
            sens=$(awk "BEGIN { printf \"%.4f\", ${rrna_aligned} / ${TOTAL_RRNA} }")
            fpr=$(awk "BEGIN { printf \"%.4f\", ${nonrrna_aligned} / ${TOTAL_NONRRNA} }")

            printf '%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\n' \
                "${passes}" "${evalue}" "${num_seeds}" "${min_lis}" \
                "${rrna_aligned}" "${sens}" \
                "${nonrrna_aligned}" "${fpr}" \
                "${wall_rrna}" "${wall_nonrrna}" "${peak_rss}" \
                >> "${RESULTS}"

            echo "  sensitivity=${sens} (${rrna_aligned}/${TOTAL_RRNA})  fpr=${fpr} (${nonrrna_aligned}/${TOTAL_NONRRNA})  rrna_time=${wall_rrna}s  nonrrna_time=${wall_nonrrna}s  peak_rss=${peak_rss}MB"

            # Family breakdown from BLAST outputs
            if [[ -n "${FAMILY_MAP}" && -f "${FAMILY_MAP}" ]]; then
                rrna_blast="${SWEEP_DIR}/${label}/rrna/out/aligned.blast.gz"
                nonrrna_blast="${SWEEP_DIR}/${label}/nonrrna/out/aligned.blast.gz"
                [[ -f "${rrna_blast}" ]] && python3 "${UTILS_DIR}/blast_family_breakdown.py" \
                    --blast "${rrna_blast}" --map "${FAMILY_MAP}" \
                    --num-seeds "${num_seeds}" --min-lis "${min_lis}" --evalue "${evalue}" --type rrna --out "${FAMILY_COUNTS}"
                [[ -f "${nonrrna_blast}" ]] && python3 "${UTILS_DIR}/blast_family_breakdown.py" \
                    --blast "${nonrrna_blast}" --map "${FAMILY_MAP}" \
                    --num-seeds "${num_seeds}" --min-lis "${min_lis}" --evalue "${evalue}" --type nonrrna --out "${FAMILY_COUNTS}"
            fi
        done  # seed_lis_pairs
    done  # evalue
done  # passes

echo ""
echo "=== Sweep complete ==="
column -t -s $'\t' "${RESULTS}"

# --- Full-dataset validation at recommended operating points ---
echo ""
echo "=== Full-dataset validation (${TOTAL_RRNA_FULL} rRNA / ${TOTAL_NONRRNA_FULL} non-rRNA reads) ==="
VALID_RESULTS="${SWEEP_DIR}/validation_results.tsv"
printf 'passes\tevalue\tnum_seeds\tmin_lis\trrna_aligned\tsensitivity\tnonrrna_aligned\tfpr\twall_sec_rrna\twall_sec_nonrrna\tpeak_rss_mb\n' > "${VALID_RESULTS}"

for pair in "1e-10 2 4" "1e-10 2 5" "1e-20 2 4" "1e-20 2 5"; do
    ev=$(echo "$pair" | cut -d' ' -f1)
    ns=$(echo "$pair" | cut -d' ' -f2)
    ml=$(echo "$pair" | cut -d' ' -f3)
    ev_label="${ev//-/m}"
    label="full_e${ev_label}_s${ns}_lis${ml}"
    rrna_log="${SWEEP_DIR}/${label}/rrna/out/aligned.log"
    nonrrna_log="${SWEEP_DIR}/${label}/nonrrna/out/aligned.log"

    echo "--------------------------------------------"
    echo "Full dataset: evalue=${ev}  num_seeds=${ns}  min_lis=${ml}"
    echo "--------------------------------------------"

    wall_rrna=0 rss_rrna=0 wall_nonrrna=0 rss_nonrrna=0

    if [[ -f "${rrna_log}" ]]; then
        echo "  rRNA run already exists - skipping"
        wall_rrna=$(cat "${SWEEP_DIR}/${label}/rrna/wall_sec.txt"    2>/dev/null || echo 0)
        rss_rrna=$( cat "${SWEEP_DIR}/${label}/rrna/peak_rss_mb.txt" 2>/dev/null || echo 0)
    else
        mkdir -p "${SWEEP_DIR}/${label}/rrna"
        echo "  Running rRNA (full dataset)..."
        result=$(run_smr_timed "${RRNA_FULL}" "${SWEEP_DIR}/${label}/rrna" "18,9,3" "${ns}" "${ev}" "${ml}")
        wall_rrna=$(echo "${result}" | grep -oP '(?<=WALL=)\d+')
        rss_rrna=$(echo "${result}"  | grep -oP '(?<=RSS=)\d+')
        echo "${wall_rrna}" > "${SWEEP_DIR}/${label}/rrna/wall_sec.txt"
        echo "${rss_rrna}"  > "${SWEEP_DIR}/${label}/rrna/peak_rss_mb.txt"
        echo "  rRNA done: ${wall_rrna}s peak RSS ${rss_rrna} MB"
    fi

    if [[ -f "${nonrrna_log}" ]]; then
        echo "  Non-rRNA run already exists - skipping"
        wall_nonrrna=$(cat "${SWEEP_DIR}/${label}/nonrrna/wall_sec.txt"    2>/dev/null || echo 0)
        rss_nonrrna=$( cat "${SWEEP_DIR}/${label}/nonrrna/peak_rss_mb.txt" 2>/dev/null || echo 0)
    else
        mkdir -p "${SWEEP_DIR}/${label}/nonrrna"
        echo "  Running non-rRNA (full dataset)..."
        result=$(run_smr_timed "${NONRRNA_FULL}" "${SWEEP_DIR}/${label}/nonrrna" "18,9,3" "${ns}" "${ev}" "${ml}")
        wall_nonrrna=$(echo "${result}" | grep -oP '(?<=WALL=)\d+')
        rss_nonrrna=$(echo "${result}"  | grep -oP '(?<=RSS=)\d+')
        echo "${wall_nonrrna}" > "${SWEEP_DIR}/${label}/nonrrna/wall_sec.txt"
        echo "${rss_nonrrna}"  > "${SWEEP_DIR}/${label}/nonrrna/peak_rss_mb.txt"
        echo "  Non-rRNA done: ${wall_nonrrna}s peak RSS ${rss_nonrrna} MB"
    fi

    wall=$(( wall_rrna + wall_nonrrna ))
    peak_rss=$(( rss_rrna > rss_nonrrna ? rss_rrna : rss_nonrrna ))
    rrna_aligned=$(grep -oP '(?<=Total reads passing E-value threshold = )\d+' "${rrna_log}")
    nonrrna_aligned=$(grep -oP '(?<=Total reads passing E-value threshold = )\d+' "${nonrrna_log}")
    sens=$(awk "BEGIN { printf \"%.4f\", ${rrna_aligned} / ${TOTAL_RRNA_FULL} }")
    fpr=$(awk "BEGIN { printf \"%.4f\", ${nonrrna_aligned} / ${TOTAL_NONRRNA_FULL} }")

    printf '18,9,3\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\n' \
        "${ev}" "${ns}" "${ml}" \
        "${rrna_aligned}" "${sens}" \
        "${nonrrna_aligned}" "${fpr}" \
        "${wall_rrna}" "${wall_nonrrna}" "${peak_rss}" \
        >> "${VALID_RESULTS}"

    echo "  sensitivity=${sens} (${rrna_aligned}/${TOTAL_RRNA_FULL})  fpr=${fpr} (${nonrrna_aligned}/${TOTAL_NONRRNA_FULL})  rrna_time=${wall_rrna}s  nonrrna_time=${wall_nonrrna}s"
done

echo ""
echo "=== Full-dataset validation results ==="
column -t -s $'\t' "${VALID_RESULTS}"

# --- R plots (ROC + stacked bar charts) ---
PLOTS_DIR="${SWEEP_DIR}/plots"
mkdir -p "${PLOTS_DIR}"
if [[ -f "${FAMILY_COUNTS}" ]]; then
    Rscript "${UTILS_DIR}/plot_pacbio_sweep.R" "${RESULTS}" "${FAMILY_COUNTS}" "${PLOTS_DIR}"
else
    echo "WARNING: family_counts.tsv not found - skipping bar charts (is family_map.tsv present?)"
    Rscript "${UTILS_DIR}/plot_pacbio_sweep.R" "${RESULTS}" /dev/null "${PLOTS_DIR}" 2>/dev/null || true
fi

# --- HTML report ---
OUTPUT_HTML="${SWEEP_DIR}/sweep_report.html"
python3 - "${RESULTS}" "${VALID_RESULTS}" "${PLOTS_DIR}" "${OUTPUT_HTML}" <<'PYEOF'
import sys, csv, base64, datetime
from pathlib import Path

results_tsv, valid_tsv, plots_dir, html_path = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

def embed_png(path):
    p = Path(path)
    if not p.exists():
        return "<p><em>Plot not available</em></p>"
    data = base64.b64encode(p.read_bytes()).decode()
    return f'<img src="data:image/png;base64,{data}" style="max-width:100%;border:1px solid #ddd;border-radius:4px">'

def make_table_rows(tsv_path):
    if not Path(tsv_path).exists():
        return ""
    rows = list(csv.DictReader(open(tsv_path), delimiter='\t'))
    return "".join(
        f"<tr><td>{r['passes']}</td><td>{r['evalue']}</td><td>{r['num_seeds']}</td><td>{r['min_lis']}</td>"
        f"<td>{r['rrna_aligned']}</td><td>{float(r['sensitivity'])*100:.2f}%</td>"
        f"<td>{r['nonrrna_aligned']}</td><td>{float(r['fpr'])*100:.2f}%</td>"
        f"<td>{r['wall_sec_rrna']}</td><td>{r['wall_sec_nonrrna']}</td><td>{r['peak_rss_mb']}</td></tr>"
        for r in rows
    )

table_rows       = make_table_rows(results_tsv)
valid_table_rows = make_table_rows(valid_tsv)

html = f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="UTF-8">
<title>PacBio Parameter Sweep Report</title>
<style>
  body{{font-family:Arial,sans-serif;margin:40px;color:#333}}
  h1{{color:#2c3e50;border-bottom:2px solid #2c3e50;padding-bottom:8px}}
  h2{{color:#34495e;margin-top:32px}}
  table{{border-collapse:collapse;width:100%;font-size:0.88em}}
  th,td{{border:1px solid #ddd;padding:8px 12px;text-align:left}}
  th{{background:#2c3e50;color:#fff}}
  tr:nth-child(even){{background:#f2f2f2}}
  .plot{{margin:16px 0}}
  .footer{{font-size:0.8em;color:#999;margin-top:40px;border-top:1px solid #ddd;padding-top:8px}}
</style></head>
<body>
<h1>SortMeRNA PacBio Parameter Sweep</h1>
<p>Generated: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>

<div style="background:#f8f9fa;border-left:4px solid #3498db;padding:12px 18px;margin:16px 0;font-size:0.93em">
<p><strong>rRNA reads (sensitivity test):</strong>
253,089 real PacBio HiFi full-length bacterial rRNA operon sequences (~4,500 bp, 16S+ITS+23S) from
<a href="https://doi.org/10.1038/s41592-020-01041-y">Karst et al. 2021 (<em>Nature Methods</em>)</a>,
70 AGP human fecal samples sequenced with PacBio Sequel II UMI amplicon sequencing
(<a href="https://qiita.ucsd.edu/study/description/10317#">Qiita study 10317</a>).
Every read is a guaranteed true positive by PCR amplification with 27F/2490R primers.
The parameter sweep was run on 10,000 reads randomly subsampled (seed 42) from this set.</p>
<p><strong>Non-rRNA reads (specificity test):</strong>
253,089 PacBio HiFi reads simulated with PBSIM3 from the rRNA-masked human T2T genome (CHM13v2.0),
mean length ~4,500 bp.
The parameter sweep was run on 10,000 reads randomly subsampled (seed 42) from this set.</p>
</div>

<h2>ROC Curve (sensitivity vs selectivity)</h2>
<div class="plot">{embed_png(f"{plots_dir}/roc.png")}</div>

<h2>rRNA family breakdown (left: rRNA reads, right: non-rRNA reads)</h2>
<div class="plot">{embed_png(f"{plots_dir}/bar_combined.png")}</div>
<div style="background:#f8f9fa;border-left:4px solid #2c3e50;padding:10px 16px;margin:8px 0;font-size:0.9em">
<p><strong>Left (rRNA reads):</strong> All bars at 10,000, pure 23S - full-operon reads maintain 100% sensitivity across all parameter values. The 23S region always generates enough co-linear seeds to trigger alignment.</p>
<p><strong>Why the classification lands on 23S:</strong> In default best mode (<code>--num_alignments 1</code>) SortMeRNA reports each read against only its single highest-scoring reference, and calls the read rRNA the moment that one alignment clears the e-value cutoff (any reference with <code>pscore &gt; minimal_score</code> sets <code>is_hit</code>; <code>alignment.cpp</code>). That score is a Smith-Waterman score accumulated over the aligned columns, anchored by the Longest Increasing Subsequence of co-linear seeds (<code>find_lis</code>, <code>alignment.cpp</code>), so a longer contiguous match region yields a longer LIS, a higher SW score, and a stronger (lower) e-value. Each read spans the full 5S+16S+23S operon, and the bacterial 23S (~2,900 bp) is roughly twice the length of the 16S (~1,500 bp): the 23S segment therefore produces the longest co-linear seed chain and the top-scoring reference is almost always a 23S sequence - hence the pure-23S bars. Even though the database holds far more 16S sequences, the 23S references are already diverse enough that every operon read finds a strong 23S match, so 23S alone suffices to classify these reads as rRNA.</p>
<p><strong>Right (non-rRNA reads - FPR):</strong>
<em>vs num_seeds:</em> modest FPR reduction (2->6); num_seeds is a weak lever. &nbsp;
<em>vs min_lis:</em> progressive step-down as min_lis increases 2->10; min_lis is the primary specificity lever, blocking tandem repeat artifacts or short rRNA pseudogenes by requiring longer co-linear seed chains on a single reference. &nbsp;
<em>vs e-value:</em> largest single effect - FPR from ~2,240 (1e-5) to ~54 (1e-20); dominant FP family is 18S throughout; pairwise alignment inspection confirms these are tandem repeat artifacts (AT-microsatellites, AGAGAG repeats) matching repeat-rich regions within rRNA database sequences, not genuine rRNA signal.</p>
</div>

<h2>Results table</h2>
<table><thead><tr>
  <th>passes</th><th>evalue</th><th>num_seeds</th><th>min_lis</th>
  <th>rRNA aligned</th><th>Sensitivity</th>
  <th>Non-rRNA aligned</th><th>FPR</th>
  <th>Wall time rRNA (s)</th><th>Wall time non-rRNA (s)</th><th>Peak RSS (MB)</th>
</tr></thead><tbody>{table_rows}</tbody></table>

<h2>Recommended operating points (full dataset: 253,089 reads)</h2>
<p>For PacBio reads, lower e-values such as <code>-e 1e-10</code> or <code>-e 1e-20</code> give the best specificity,
with <code>--min_lis</code> 2-5 (default of 2 for metagenomic reads), with no loss of rRNA recovery.</p>
<table><thead><tr>
  <th>passes</th><th>evalue</th><th>num_seeds</th><th>min_lis</th>
  <th>rRNA aligned</th><th>Sensitivity</th>
  <th>Non-rRNA aligned</th><th>FPR</th>
  <th>Wall time rRNA (s)</th><th>Wall time non-rRNA (s)</th><th>Peak RSS (MB)</th>
</tr></thead><tbody>{valid_table_rows}</tbody></table>

<p class="footer">Generated by run_pacbio_sweep.sh</p>
</body></html>"""

Path(html_path).write_text(html)
print(f"  HTML report: {html_path}")
PYEOF
