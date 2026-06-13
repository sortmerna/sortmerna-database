#!/usr/bin/env bash
# run_scalability.sh  Run SortMeRNA at increasing read volumes and generate plots.
#
# Usage:
#   bash run_scalability.sh <reads.fasta> <output_dir> <threads> [OPTIONS]
#
# Arguments:
#   reads.fasta    Input reads in FASTA format (must have >= max scale-point reads)
#   output_dir     Directory for all outputs
#   threads        Number of SortMeRNA threads (default: 4)
#
# Options:
#   --scale INT,...   Comma-separated read counts to test (default: 10000,100000,1000000,10000000)
#   --label STR       Label used in plot titles and output filenames
#                     (default: basename of reads.fasta without extension)
#   --config STR      SortMeRNA database config name under INDEX_DIR (default: smr_v<version>_default_db,
#                     where <version> is read from SMR_BIN at runtime)
#   --index-dir DIR   Directory containing SortMeRNA index subdirectories
#                     (overrides INDEX_DIR env var if set)
#   --seed INT        Random seed for seqkit subsampling (default: 42)
#   --score-split     Pass --score_split to SortMeRNA: compute S_min from per-thread chunk
#                     size rather than total dataset size, making the E-value threshold
#                     less dependent on total read count. Off by default (standard behavior).
#
# Required env vars (can be overridden with the matching option above):
#   SMR_BIN          Full path to SortMeRNA binary
#   INDEX_DIR        Directory containing SortMeRNA index subdirectories (overridable with --index-dir)
#   SMR_DB_ROOT_DIR  Root of the sortmerna-database repository
#
# Outputs (all under <output_dir>):
#   scale_<N>/reads_<N>.fasta           Subsampled reads
#   scale_<N>/smr_out/out/aligned.log   SortMeRNA summary log
#   scale_<N>/smr_out/out/aligned.blast BLAST-format alignments (--blast 1)
#   scale_<N>/runtime_seconds.txt       Wall-clock runtime in seconds
#   scale_<N>/peak_rss_mb.txt          Peak resident set size in MB
#   plots/                              Generated figures (PNG)

set -euo pipefail

READS="$1"
OUTPUT_DIR="$2"
THREADS="${3:-4}"

SCALE_POINTS_CSV="10000,100000,1000000,10000000"
RAND_SEED=42
DB_CONFIG="${DB_CONFIG:-}"
LABEL=""
INDEX_DIR_OPT=""
SCORE_SPLIT=false
EVALUE=""
READS_DIR=""

shift 3 || true
while [[ $# -gt 0 ]]; do
    case "$1" in
        --scale)       SCALE_POINTS_CSV="$2"; shift 2 ;;
        --label)       LABEL="$2";            shift 2 ;;
        --config)      DB_CONFIG="$2";        shift 2 ;;
        --index-dir)   INDEX_DIR_OPT="$2";    shift 2 ;;
        --seed)        RAND_SEED="$2";        shift 2 ;;
        --score-split) SCORE_SPLIT=true;      shift   ;;
        --evalue)      EVALUE="$2";           shift 2 ;;
        --reads-dir)   READS_DIR="$2";        shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

[[ -n "${INDEX_DIR_OPT}" ]] && INDEX_DIR="${INDEX_DIR_OPT}"
if [[ ! -f "${READS}" ]] && [[ -f "${READS}.gz" ]]; then READS="${READS}.gz"; fi
[[ -z "${LABEL}" ]] && LABEL="$(basename "${READS}" | sed 's/\.fasta\.gz$//; s/\.fasta$//')"

: "${SMR_BIN:?SMR_BIN env var not set}"
: "${INDEX_DIR:?INDEX_DIR not set - pass --index-dir or export INDEX_DIR}"
: "${SMR_DB_ROOT_DIR:?SMR_DB_ROOT_DIR env var not set}"

if [[ -z "${DB_CONFIG}" ]]; then
    SMR_VERSION=$("${SMR_BIN}" --version 2>&1 | grep "^SortMeRNA version" | awk '{print $3}')
    DB_CONFIG="smr_v${SMR_VERSION}_default_db"
fi

UTILS_DIR="${SMR_DB_ROOT_DIR}/scripts/utils"
DB_FASTA="${INDEX_DIR}/${DB_CONFIG}/${DB_CONFIG}.fasta"
DB_IDX="${INDEX_DIR}/${DB_CONFIG}/idx"

IFS=',' read -ra SCALE_POINTS <<< "${SCALE_POINTS_CSV}"

echo "============================================"
echo "Scalability benchmark: ${LABEL}"
echo "============================================"
echo "  Reads:       ${READS}"
echo "  Output:      ${OUTPUT_DIR}"
echo "  Scale points: ${SCALE_POINTS_CSV}"
echo "  DB config:   ${DB_CONFIG}"
echo "  Threads:     ${THREADS}"
echo "  Seed:        ${RAND_SEED}"
echo "  Score split: ${SCORE_SPLIT}"
echo "  E-value:     ${EVALUE:-default}"
echo "  Reads dir:   ${READS_DIR:-none (subsample from reads)}"
echo ""

mkdir -p "${OUTPUT_DIR}"
scale_dirs=()

for n in "${SCALE_POINTS[@]}"; do
    scale_dir="${OUTPUT_DIR}/scale_${n}"
    mkdir -p "${scale_dir}"
    scale_dirs+=("${scale_dir}")

    echo "--------------------------------------------"
    echo "Scale point: ${n} reads"
    echo "--------------------------------------------"

    if [[ -n "${READS_DIR}" ]]; then
        subset_fa="${READS_DIR}/scale_${n}/reads_${n}.fasta"
        if [[ ! -f "${subset_fa}" ]]; then
            echo "  ERROR: --reads-dir set but ${subset_fa} not found" >&2; exit 1
        fi
        echo "  Reusing reads: ${subset_fa}"
    else
        subset_fa="${scale_dir}/reads_${n}.fasta"
        if [[ ! -f "${subset_fa}" ]]; then
            echo "  Subsampling ${n} reads..."
            seqkit sample -2 -n "${n}" --rand-seed "${RAND_SEED}" "${READS}" \
                | seqkit seq -w 0 > "${subset_fa}"
            echo "  Saved: $(basename "${subset_fa}")"
        else
            echo "  Already exists: $(basename "${subset_fa}")"
        fi
    fi

    smr_workdir="${scale_dir}/smr_out"
    aligned_log="${smr_workdir}/out/aligned.log"

    if [[ ! -f "${aligned_log}" ]]; then
        echo "  Running SortMeRNA..."
        start=$(date +%s)
        smr_args=(
            --ref     "${DB_FASTA}"
            --reads   "${subset_fa}"
            --idx-dir "${DB_IDX}"
            --workdir "${smr_workdir}"
            --fastx --blast 1
            --threads "${THREADS}"
        )
        [[ "${SCORE_SPLIT}" == true ]] && smr_args+=(--score_split)
        [[ -n "${EVALUE}" ]] && smr_args+=(-e "${EVALUE}")
        "${SMR_BIN}" "${smr_args[@]}" &
        smr_pid=$!
        peak_rss_mb=0
        while kill -0 "${smr_pid}" 2>/dev/null; do
            rss=$(ps -p "${smr_pid}" -o rss --no-headers 2>/dev/null || echo 0)
            rss_mb=$(( (rss + 0) / 1024 ))
            if (( rss_mb > peak_rss_mb )); then peak_rss_mb=${rss_mb}; fi
            sleep 5
        done
        wait "${smr_pid}" || { echo "  ERROR: sortmerna failed"; exit 1; }
        end=$(date +%s)
        runtime=$(( end - start ))
        echo "${runtime}" > "${scale_dir}/runtime_seconds.txt"
        echo "${peak_rss_mb}" > "${scale_dir}/peak_rss_mb.txt"
        echo "  Done: ${runtime}s - peak RSS: ${peak_rss_mb} MB"
    else
        echo "  Already exists: aligned.log"
    fi
done

echo ""
echo "============================================"
echo "Generating plots..."
echo "============================================"

python3 "${UTILS_DIR}/plot_scalability.py" \
    --output-dir "${OUTPUT_DIR}/plots" \
    --scale-dirs "${scale_dirs[@]}" \
    --label "${LABEL}"

echo ""
echo "Done. Results in: ${OUTPUT_DIR}"
