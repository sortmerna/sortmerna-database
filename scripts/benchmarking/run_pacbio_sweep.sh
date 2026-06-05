#!/usr/bin/env bash
# run_pacbio_sweep.sh  SortMeRNA --passes / --num_seeds parameter sweep for PacBio HiFi reads.
#
# Tests 4x4=16 combinations of --passes and --num_seeds against real rRNA reads
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

# Quick-check mode: iterate on a 10k subset first, then run the full dataset
# seqtk sample -s 42 "$RRNA_READS"    10000 > "${RRNA_READS%.fastq}_10k.fastq"
# seqtk sample -s 42 "$NONRRNA_READS" 10000 > "${NONRRNA_READS%.fastq}_10k.fastq"
# RRNA_READS="${RRNA_READS%.fastq}_10k.fastq"
# NONRRNA_READS="${NONRRNA_READS%.fastq}_10k.fastq"
# TOTAL_RRNA=10000; TOTAL_NONRRNA=10000

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
printf 'passes\tnum_seeds\trrna_aligned\tsensitivity\tnonrrna_aligned\tfpr\twall_sec\n' > "${RESULTS}"

run_smr() {
    local reads="$1"
    local workdir="$2"
    local passes="$3"
    local num_seeds="$4"
    "${SMR_BIN}" \
        --ref "${REF_DB}" \
        --reads "${reads}" \
        --workdir "${workdir}" \
        --idx-dir "${IDX_DIR}" \
        --passes "${passes}" \
        --num_seeds "${num_seeds}" \
        --threads "${THREADS}" \
        --fastx --blast 1 \
        -e 1e-5
}

PASSES_LIST=("18,9,3" "100,50,10" "200,100,20" "500,200,50")
SEEDS_LIST=(2 5 10 25)

for passes in "${PASSES_LIST[@]}"; do
    for num_seeds in "${SEEDS_LIST[@]}"; do
        label="p${passes//,/_}_s${num_seeds}"
        rrna_log="${SWEEP_DIR}/${label}/rrna/out/aligned.log"
        nonrrna_log="${SWEEP_DIR}/${label}/nonrrna/out/aligned.log"

        echo "--------------------------------------------"
        echo "passes=${passes}  num_seeds=${num_seeds}"
        echo "--------------------------------------------"

        start=$SECONDS

        if [[ -f "${rrna_log}" ]]; then
            echo "  rRNA run already exists - skipping"
        else
            mkdir -p "${SWEEP_DIR}/${label}/rrna"
            run_smr "${RRNA_READS}" "${SWEEP_DIR}/${label}/rrna" "${passes}" "${num_seeds}"
        fi

        if [[ -f "${nonrrna_log}" ]]; then
            echo "  Non-rRNA run already exists - skipping"
        else
            mkdir -p "${SWEEP_DIR}/${label}/nonrrna"
            run_smr "${NONRRNA_READS}" "${SWEEP_DIR}/${label}/nonrrna" "${passes}" "${num_seeds}"
        fi

        wall=$(( SECONDS - start ))

        rrna_aligned=$(grep -oP '(?<=Total reads passing E-value threshold = )\d+' "${rrna_log}")
        nonrrna_aligned=$(grep -oP '(?<=Total reads passing E-value threshold = )\d+' "${nonrrna_log}")
        sens=$(awk "BEGIN { printf \"%.4f\", ${rrna_aligned} / ${TOTAL_RRNA} }")
        fpr=$(awk "BEGIN { printf \"%.4f\", ${nonrrna_aligned} / ${TOTAL_NONRRNA} }")

        printf '%s\t%d\t%s\t%s\t%s\t%s\t%d\n' \
            "${passes}" "${num_seeds}" \
            "${rrna_aligned}" "${sens}" \
            "${nonrrna_aligned}" "${fpr}" \
            "${wall}" \
            >> "${RESULTS}"

        echo "  sensitivity=${sens} (${rrna_aligned}/${TOTAL_RRNA})  fpr=${fpr} (${nonrrna_aligned}/${TOTAL_NONRRNA})  time=${wall}s"
    done
done

echo ""
echo "=== Sweep complete ==="
column -t -s $'\t' "${RESULTS}"
