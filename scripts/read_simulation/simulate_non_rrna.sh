#!/usr/bin/env bash

################################################################################
# simulate_non_rrna.sh
#
# Mask rRNA loci in the T2T genome, simulate Illumina reads with InSilicoSeq,
# and filter Rfam non-rRNA sequences by length to produce two test sets:
#
#   non_rRNA_test_1M_T2T.fasta   - 1M simulated T2T genome reads (rRNA loci masked)
#   non_rRNA_test_Rfam.fasta     - Rfam non-rRNA sequences in read-length range (100-200 bp)
#
# Inputs (from download_non_rrna.sh):
#   t2t/${T2T_VERSION}.fa.gz         - T2T genome (compressed)
#   t2t/${T2T_VERSION}_rrna_loci.bed - rRNA loci BED (from extract_rrna_loci.py)
#   rfam_non_rrna_all.fasta          - All Rfam non-rRNA sequences
#
# Steps:
#   1. Decompress T2T genome
#   2. Mask rRNA loci with bedtools maskfasta
#   3. Simulate 1M Illumina PE reads with InSilicoSeq
#   4. Convert to FASTA -> non_rRNA_test_1M_T2T.fasta
#   5. Filter Rfam sequences 100-200 bp -> non_rRNA_test_Rfam.fasta
#
# Usage: bash simulate_non_rrna.sh [output_dir [threads]] [OPTIONS]
#
# Positional:
#   output_dir      Directory with download_non_rrna.sh outputs (default: $NON_RRNA_DIR or data/non_rrna)
#   threads         Number of threads for InSilicoSeq (default: 4)
#
# Options:
#   --t2t-reads INT   Number of simulated T2T reads to include (default: 1000000)
#                     InSilicoSeq is run with N_T2T/2 read pairs so that R1+R2
#                     combined equals exactly N_T2T individual reads.
#   --model STR       InSilicoSeq error model: HiSeq, NovaSeq, MiSeq (default: HiSeq)
#   --min-len INT     Minimum Rfam sequence length to include (default: 100)
#   --max-len INT     Maximum Rfam sequence length to include (default: 200)
#   --seed INT        Random seed for InSilicoSeq simulation (default: 42)
#   --skip-mask       Skip masking step, use existing masked genome
#   --skip-sim        Skip InSilicoSeq simulation, use existing output
#   -h, --help        Show help
#
# Environment variables:
#   T2T_VERSION   T2T genome version string (default: chm13v2.0)
#
# Requires: bedtools, iss (InSilicoSeq), seqkit
#
################################################################################

set -euo pipefail

POSITIONAL=()
N_T2T=1000000
ISS_MODEL=HiSeq
RFAM_MIN_LEN=100
RFAM_MAX_LEN=200
RAND_SEED=42
SKIP_MASK=false
SKIP_SIM=false

T2T_VERSION="${T2T_VERSION:-chm13v2.0}"

show_help() {
    grep '^#' "$0" | grep -v '#!/usr/bin/env bash' | sed 's/^# \?//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --t2t-reads) N_T2T="$2"; shift 2 ;;
        --model) ISS_MODEL="$2"; shift 2 ;;
        --min-len) RFAM_MIN_LEN="$2"; shift 2 ;;
        --max-len) RFAM_MAX_LEN="$2"; shift 2 ;;
        --seed) RAND_SEED="$2"; shift 2 ;;
        --skip-mask) SKIP_MASK=true; shift ;;
        --skip-sim) SKIP_SIM=true; shift ;;
        -h|--help) show_help ;;
        *) POSITIONAL+=("$1"); shift ;;
    esac
done

OUTPUT_DIR="${POSITIONAL[0]:-${NON_RRNA_DIR:-data/non_rrna}}"
THREADS="${POSITIONAL[1]:-4}"

OUTPUT_DIR="$(mkdir -p "${OUTPUT_DIR}" && cd "${OUTPUT_DIR}" && pwd)"

T2T_DIR="${OUTPUT_DIR}/t2t"
ISS_DIR="${OUTPUT_DIR}/iss"

T2T_GENOME_GZ="${T2T_DIR}/${T2T_VERSION}.fa.gz"
T2T_FA="${T2T_DIR}/${T2T_VERSION}.fa"
T2T_RRNA_BED="${T2T_DIR}/${T2T_VERSION}_rrna_loci.bed"
T2T_MASKED="${T2T_DIR}/${T2T_VERSION}_masked.fa"

RFAM_ALL="${OUTPUT_DIR}/rfam_non_rrna_all.fasta"
T2T_OUTPUT="${OUTPUT_DIR}/non_rRNA_test_1M_T2T.fasta"
RFAM_OUTPUT="${OUTPUT_DIR}/non_rRNA_test_Rfam.fasta"

echo "============================================"
echo "Non-rRNA Read Simulation Script"
echo "Output directory: ${OUTPUT_DIR}"
echo "T2T version:      ${T2T_VERSION}"
echo "T2T reads:        ${N_T2T}"
echo "ISS model:        ${ISS_MODEL}"
echo "Rfam length:      ${RFAM_MIN_LEN}-${RFAM_MAX_LEN} bp"
echo "Random seed:      ${RAND_SEED}"
echo "Threads:          ${THREADS}"
echo "============================================"
echo ""

for tool in bedtools iss seqkit; do
    if ! command -v "${tool}" &> /dev/null; then
        echo "Error: ${tool} not found. Install via: conda env create -f environment.yml"
        exit 1
    fi
done

for f in "${T2T_GENOME_GZ}" "${T2T_RRNA_BED}" "${RFAM_ALL}"; do
    if [[ ! -f "${f}" ]]; then
        echo "Error: required input not found: ${f}"
        echo "Run download_non_rrna.sh first."
        exit 1
    fi
done

################################################################################
# 1. DECOMPRESS T2T GENOME
################################################################################

echo "============================================"
echo "Step 1: Decompress T2T genome"
echo "============================================"

if [[ ! -f "${T2T_FA}" ]]; then
    echo "Decompressing ${T2T_VERSION}.fa.gz..."
    zcat "${T2T_GENOME_GZ}" > "${T2T_FA}"
    echo "  Saved: ${T2T_VERSION}.fa"
else
    echo "Already exists: ${T2T_VERSION}.fa"
fi

################################################################################
# 2. MASK rRNA LOCI
################################################################################

echo ""
echo "============================================"
echo "Step 2: Mask rRNA loci"
echo "============================================"

if [[ "${SKIP_MASK}" == false ]]; then
    if [[ ! -f "${T2T_MASKED}" ]]; then
        n_loci=$(wc -l < "${T2T_RRNA_BED}")
        echo "Masking ${n_loci} rRNA loci with bedtools maskfasta..."
        bedtools maskfasta \
            -fi "${T2T_FA}" \
            -bed "${T2T_RRNA_BED}" \
            -fo "${T2T_MASKED}"
        echo "  Saved: ${T2T_VERSION}_masked.fa"
    else
        echo "Already exists: ${T2T_VERSION}_masked.fa"
    fi
else
    echo "Skipping mask (--skip-mask)."
    if [[ ! -f "${T2T_MASKED}" ]]; then
        echo "Error: --skip-mask set but ${T2T_MASKED} does not exist"
        exit 1
    fi
fi

################################################################################
# 3. SIMULATE ILLUMINA READS WITH InSilicoSeq
################################################################################

echo ""
echo "============================================"
echo "Step 3: Simulate Illumina PE reads (InSilicoSeq)"
echo "============================================"

mkdir -p "${ISS_DIR}"
ISS_PREFIX="${ISS_DIR}/${T2T_VERSION}_masked"
ISS_R1="${ISS_PREFIX}_R1.fastq"
ISS_R2="${ISS_PREFIX}_R2.fastq"

if [[ "${SKIP_SIM}" == false ]]; then
    if [[ ! -f "${ISS_R1}" ]] || [[ ! -f "${ISS_R2}" ]]; then
        # iss generates n_reads PE pairs, so R1+R2 combined = 2 * n_reads.
        # Use N_T2T/2 pairs to get exactly N_T2T individual reads after combining.
        N_PAIRS=$(( N_T2T / 2 ))
        echo "Generating ${N_PAIRS} PE pairs (${N_T2T} individual reads) with InSilicoSeq..."
        iss generate \
            --genomes "${T2T_MASKED}" \
            --model "${ISS_MODEL}" \
            --n_reads "${N_PAIRS}" \
            --cpus "${THREADS}" \
            --seed "${RAND_SEED}" \
            --output "${ISS_PREFIX}"
        echo "  Saved: $(basename "${ISS_R1}"), $(basename "${ISS_R2}")"
    else
        echo "Already exists: InSilicoSeq output files"
    fi
else
    echo "Skipping simulation (--skip-sim)."
    if [[ ! -f "${ISS_R1}" ]] || [[ ! -f "${ISS_R2}" ]]; then
        echo "Error: --skip-sim set but ISS output not found: ${ISS_R1}"
        exit 1
    fi
fi

################################################################################
# 4. CONVERT T2T READS TO FASTA
################################################################################

echo ""
echo "============================================"
echo "Step 4: Convert T2T reads to FASTA"
echo "============================================"

echo "Converting ISS FASTQ to FASTA..."
cat "${ISS_R1}" "${ISS_R2}" | seqkit fq2fa > "${T2T_OUTPUT}"

n_t2t=$(seqkit stats -T "${T2T_OUTPUT}" | tail -1 | cut -f4)
echo "  Saved: non_rRNA_test_1M_T2T.fasta (${n_t2t} reads)"

################################################################################
# 5. FILTER Rfam SEQUENCES BY LENGTH
################################################################################

echo ""
echo "============================================"
echo "Step 5: Filter Rfam sequences (${RFAM_MIN_LEN}-${RFAM_MAX_LEN} bp)"
echo "============================================"

seqkit seq --min-len "${RFAM_MIN_LEN}" --max-len "${RFAM_MAX_LEN}" \
    "${RFAM_ALL}" > "${RFAM_OUTPUT}"

n_rfam=$(seqkit stats -T "${RFAM_OUTPUT}" | tail -1 | cut -f4)
n_rfam_total=$(seqkit stats -T "${RFAM_ALL}" | tail -1 | cut -f4)
echo "  Saved: non_rRNA_test_Rfam.fasta (${n_rfam} of ${n_rfam_total} sequences in ${RFAM_MIN_LEN}-${RFAM_MAX_LEN} bp range)"

echo ""
echo "============================================"
echo "Simulation complete"
echo "============================================"
echo ""
echo "Outputs:"
echo "  Masked genome:  ${T2T_MASKED}"
echo "  T2T test set:   ${T2T_OUTPUT} (${n_t2t} reads)"
echo "  Rfam test set:  ${RFAM_OUTPUT} (${n_rfam} sequences, ${RFAM_MIN_LEN}-${RFAM_MAX_LEN} bp)"
echo ""
echo "Next step: Run simulate_rrna_reads.sh to simulate rRNA reads for sensitivity testing"
