#!/usr/bin/env bash

################################################################################
# simulate_non_rrna.sh
#
# Mask rRNA loci in the T2T genome, simulate Illumina reads with ART, and
# combine with Rfam non-rRNA sequences into a ~1M-read non-rRNA test set.
#
# Inputs (from download_non_rrna.sh):
#   t2t/${T2T_VERSION}.fa.gz         - T2T genome (compressed)
#   t2t/${T2T_VERSION}_rrna_loci.bed - rRNA loci BED (from extract_rrna_loci.py)
#   rfam_non_rrna_sampled.fasta      - Rfam non-rRNA sequences (150,000)
#
# Steps:
#   1. Decompress T2T genome
#   2. Mask rRNA loci with bedtools maskfasta
#   3. Simulate Illumina PE reads with ART (art_illumina)
#   4. Subsample + combine T2T reads and Rfam sequences
#
# Output:
#   non_rRNA_test_1M.fasta           - ~1M non-rRNA sequences for specificity testing
#
# Usage: bash simulate_non_rrna.sh [output_dir [threads]] [OPTIONS]
#
# Positional:
#   output_dir    Directory with download_non_rrna.sh outputs (default: $NON_RRNA_DIR or data/non_rrna)
#   threads       Number of threads for ART (default: 4)
#
# Options:
#   --t2t-reads INT   Number of simulated T2T reads to include (default: 850000)
#   --seed INT        Random seed for subsampling (default: 42)
#   --skip-mask       Skip masking step, use existing masked genome
#   --skip-art        Skip ART simulation, use existing ART output
#   -h, --help        Show help
#
# Environment variables:
#   T2T_VERSION   T2T genome version string (default: chm13v2.0)
#
# Requires: bedtools, art_illumina, seqkit
#
################################################################################

set -euo pipefail

POSITIONAL=()
N_T2T=850000
RAND_SEED=42
SKIP_MASK=false
SKIP_ART=false

T2T_VERSION="${T2T_VERSION:-chm13v2.0}"

show_help() {
    grep '^#' "$0" | grep -v '#!/usr/bin/env bash' | sed 's/^# \?//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --t2t-reads) N_T2T="$2"; shift 2 ;;
        --seed) RAND_SEED="$2"; shift 2 ;;
        --skip-mask) SKIP_MASK=true; shift ;;
        --skip-art) SKIP_ART=true; shift ;;
        -h|--help) show_help ;;
        *) POSITIONAL+=("$1"); shift ;;
    esac
done

OUTPUT_DIR="${POSITIONAL[0]:-${NON_RRNA_DIR:-data/non_rrna}}"
THREADS="${POSITIONAL[1]:-4}"

OUTPUT_DIR="$(mkdir -p "${OUTPUT_DIR}" && cd "${OUTPUT_DIR}" && pwd)"

T2T_DIR="${OUTPUT_DIR}/t2t"
ART_DIR="${OUTPUT_DIR}/art"

T2T_GENOME_GZ="${T2T_DIR}/${T2T_VERSION}.fa.gz"
T2T_FA="${T2T_DIR}/${T2T_VERSION}.fa"
T2T_RRNA_BED="${T2T_DIR}/${T2T_VERSION}_rrna_loci.bed"
T2T_MASKED="${T2T_DIR}/${T2T_VERSION}_masked.fa"

RFAM_SAMPLED="${OUTPUT_DIR}/rfam_non_rrna_sampled.fasta"
OUTPUT_FASTA="${OUTPUT_DIR}/non_rRNA_test_1M.fasta"

echo "============================================"
echo "Non-rRNA Read Simulation Script"
echo "Output directory: ${OUTPUT_DIR}"
echo "T2T version:      ${T2T_VERSION}"
echo "T2T reads:        ${N_T2T}"
echo "Random seed:      ${RAND_SEED}"
echo "Threads:          ${THREADS}"
echo "============================================"
echo ""

for tool in bedtools art_illumina seqkit; do
    if ! command -v "${tool}" &> /dev/null; then
        echo "Error: ${tool} not found. Install via: conda env create -f environment.yml"
        exit 1
    fi
done

# Verify required inputs exist
for f in "${T2T_GENOME_GZ}" "${T2T_RRNA_BED}" "${RFAM_SAMPLED}"; do
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
# 3. SIMULATE ILLUMINA READS WITH ART
################################################################################

echo ""
echo "============================================"
echo "Step 3: Simulate Illumina PE reads (ART)"
echo "============================================"

mkdir -p "${ART_DIR}"
ART_PREFIX="${ART_DIR}/${T2T_VERSION}_masked"
ART_R1="${ART_PREFIX}1.fq"
ART_R2="${ART_PREFIX}2.fq"

if [[ "${SKIP_ART}" == false ]]; then
    if [[ ! -f "${ART_R1}" ]] || [[ ! -f "${ART_R2}" ]]; then
        # Calculate fold coverage to produce roughly N_T2T * 1.2 reads per end so that
        # after subsampling we have exactly N_T2T total (R1+R2 combined).
        # total_bases = N_T2T * 1.2 * 150  (150 bp reads, 20% oversample for subsampling)
        # fold = total_bases / genome_size
        genome_size=$(seqkit stats -T "${T2T_MASKED}" | tail -1 | cut -f5)
        # awk: compute fold, round up to 2 decimal places
        FOLD=$(awk -v n="${N_T2T}" -v gs="${genome_size}" \
            'BEGIN { printf "%.2f\n", (n * 1.2 * 150) / gs }')
        echo "Masked genome size: ${genome_size} bp"
        echo "ART fold coverage:  ${FOLD}x  (targeting ~$(( N_T2T * 12 / 10 )) reads x 2 ends)"
        echo "Running art_illumina (this may take several minutes)..."
        art_illumina \
            -ss HS25 \
            -i "${T2T_MASKED}" \
            -p \
            -l 150 \
            -f "${FOLD}" \
            -rs "${RAND_SEED}" \
            -o "${ART_PREFIX}"
        echo "  Saved: $(basename "${ART_R1}"), $(basename "${ART_R2}")"
    else
        echo "Already exists: ART output files"
    fi
else
    echo "Skipping ART simulation (--skip-art)."
    if [[ ! -f "${ART_R1}" ]] || [[ ! -f "${ART_R2}" ]]; then
        echo "Error: --skip-art set but ART output not found: ${ART_R1}"
        exit 1
    fi
fi

################################################################################
# 4. SUBSAMPLE + COMBINE
################################################################################

echo ""
echo "============================================"
echo "Step 4: Subsample T2T reads and combine with Rfam"
echo "============================================"

T2T_READS_FASTA="${ART_DIR}/${T2T_VERSION}_reads.fasta"

echo "Converting ART FASTQ to FASTA and subsampling to ${N_T2T} reads..."
cat "${ART_R1}" "${ART_R2}" \
    | seqkit sample --rand-seed "${RAND_SEED}" -n "${N_T2T}" \
    | seqkit fq2fa \
    > "${T2T_READS_FASTA}"

n_t2t_actual=$(seqkit stats -T "${T2T_READS_FASTA}" | tail -1 | cut -f4)
echo "  T2T reads: ${n_t2t_actual}"

n_rfam=$(seqkit stats -T "${RFAM_SAMPLED}" | tail -1 | cut -f4)
echo "  Rfam sequences: ${n_rfam}"

echo "Combining into ${OUTPUT_FASTA}..."
cat "${T2T_READS_FASTA}" "${RFAM_SAMPLED}" > "${OUTPUT_FASTA}"

n_total=$(seqkit stats -T "${OUTPUT_FASTA}" | tail -1 | cut -f4)
echo "  Total: ${n_total} sequences"

echo ""
echo "============================================"
echo "Simulation complete"
echo "============================================"
echo ""
echo "Outputs:"
echo "  Masked genome:     ${T2T_MASKED}"
echo "  T2T reads FASTA:   ${T2T_READS_FASTA}"
echo "  Non-rRNA test set: ${OUTPUT_FASTA}"
echo "    T2T reads:       ${n_t2t_actual}"
echo "    Rfam sequences:  ${n_rfam}"
echo "    Total:           ${n_total}"
echo ""
echo "Next step: Run simulate_rrna_reads.sh to simulate rRNA reads for sensitivity testing"
