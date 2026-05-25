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
#   --scale INT,...  Comma-separated read counts to test (default: 10000,100000,1000000,10000000)
#   --label STR      Label used in plot titles and output filenames
#                    (default: basename of reads.fasta without extension)
#   --config STR     SortMeRNA database config name (default: smr_v5.0.0_default_db)
#   --seed INT       Random seed for seqkit subsampling (default: 42)
#
# Required env vars:
#   SMR_BIN          Full path to SortMeRNA binary
#   INDEX_DIR        Directory containing SortMeRNA index subdirectories
#   SMR_DB_ROOT_DIR  Root of the sortmerna-database repository
#
# Outputs (all under <output_dir>):
#   scale_<N>/reads_<N>.fasta           Subsampled reads
#   scale_<N>/smr_out/out/aligned.log   SortMeRNA summary log
#   scale_<N>/smr_out/out/aligned.blast BLAST-format alignments (--blast 1)
#   scale_<N>/runtime_seconds.txt       Wall-clock runtime in seconds
#   plots/                              Generated figures (PNG)

set -euo pipefail

READS="$1"
OUTPUT_DIR="$2"
THREADS="${3:-4}"

SCALE_POINTS_CSV="10000,100000,1000000,10000000"
RAND_SEED=42
DB_CONFIG="${DB_CONFIG:-smr_v5.0.0_default_db}"
LABEL=""

shift 3 || true
while [[ $# -gt 0 ]]; do
    case "$1" in
        --scale)  SCALE_POINTS_CSV="$2"; shift 2 ;;
        --label)  LABEL="$2";            shift 2 ;;
        --config) DB_CONFIG="$2";        shift 2 ;;
        --seed)   RAND_SEED="$2";        shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

[[ -z "${LABEL}" ]] && LABEL="$(basename "${READS}" .fasta)"

: "${SMR_BIN:?SMR_BIN env var not set}"
: "${INDEX_DIR:?INDEX_DIR env var not set}"
: "${SMR_DB_ROOT_DIR:?SMR_DB_ROOT_DIR env var not set}"

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

    subset_fa="${scale_dir}/reads_${n}.fasta"
    if [[ ! -f "${subset_fa}" ]]; then
        echo "  Subsampling ${n} reads..."
        seqkit sample -n "${n}" --rand-seed "${RAND_SEED}" "${READS}" \
            | seqkit seq -w 0 > "${subset_fa}"
        echo "  Saved: $(basename "${subset_fa}")"
    else
        echo "  Already exists: $(basename "${subset_fa}")"
    fi

    smr_workdir="${scale_dir}/smr_out"
    aligned_log="${smr_workdir}/out/aligned.log"

    if [[ ! -f "${aligned_log}" ]]; then
        echo "  Running SortMeRNA..."
        start=$(date +%s)
        "${SMR_BIN}" \
            --ref "${DB_FASTA}" \
            --reads "${subset_fa}" \
            --idx-dir "${DB_IDX}" \
            --workdir "${smr_workdir}" \
            --fastx --blast 1 \
            --threads "${THREADS}"
        end=$(date +%s)
        runtime=$(( end - start ))
        echo "${runtime}" > "${scale_dir}/runtime_seconds.txt"
        echo "  Done: ${runtime}s"
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
