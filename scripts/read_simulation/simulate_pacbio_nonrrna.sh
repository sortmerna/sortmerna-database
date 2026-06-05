#!/usr/bin/env bash

################################################################################
# simulate_pacbio_nonrrna.sh
#
# Simulate PacBio HiFi non-rRNA reads from the masked T2T genome using PBSIM3.
# The masked genome (rRNA loci replaced with Ns) is produced by
# download_non_rrna.sh. Rfam non-rRNA families are not used here because most
# sequences are shorter than typical PacBio HiFi read lengths.
#
# Output:
#   non_rrna_pacbio_<N>_T2T.fastq  - N simulated PacBio HiFi reads
#
# Input (from download_non_rrna.sh / simulate_non_rrna.sh):
#   t2t/${T2T_VERSION}_masked.fa   - T2T genome with rRNA loci masked
#
# Steps:
#   1. Compute --depth to produce approximately N_READS reads at the target
#      mean length from the masked genome (genome size estimated with seqkit)
#   2. Run PBSIM3 (--strategy wgs) across all chromosomes
#   3. Merge per-chromosome FASTQ outputs into a single file
#   4. Subsample to exactly N_READS with seqkit if PBSIM3 overshoots
#   5. Clean up PBSIM3 intermediate files
#
# Usage:
#   bash simulate_pacbio_nonrrna.sh [output_dir] [OPTIONS]
#
# Positional:
#   output_dir   Directory containing download_non_rrna.sh outputs and where
#                simulated non-rRNA PacBio reads will be written
#                (default: $NON_RRNA_DIR or data/non_rrna)
#
# Options:
#   --reads INT          Target number of simulated reads (default: 253089)
#                        Default matches the Karst et al. 2021 PacBio dataset
#                        used for the rRNA sensitivity test.
#   --length-mean INT    Mean read length in bp (default: 4500)
#                        Matches the ~4,500 bp 16S+ITS+23S amplicons in Karst
#                        et al. 2021, enabling a fair length-matched comparison.
#   --length-sd INT      Read length standard deviation in bp (default: 500)
#   --accuracy FLOAT     Mean read accuracy (default: 0.9999, matching Karst et al. 2021
#                        CCS error rate of ~0.0007-0.008%)
#   --model STR          PBSIM3 error model file path (default: auto-detect
#                        ERRHMM-SEQUEL.model from PBSIM3 installation)
#   --seed INT           Random seed (default: 42)
#   --keep-intermediates Keep per-chromosome FASTQ/MAF/ref files (default: remove)
#   -h, --help           Show help
#
# Environment variables:
#   T2T_VERSION   T2T genome version string (default: chm13v2.0)
#   PBSIM3_BIN    Full path to pbsim3 binary (default: pbsim3 on PATH)
#
# Requires: pbsim3, seqkit
#
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Defaults ---
OUTPUT_DIR="${NON_RRNA_DIR:-data/non_rrna}"
N_READS=253089
LENGTH_MEAN=4500
LENGTH_SD=500
ACCURACY=0.9999  # Karst et al. 2021 CCS error rate ~0.0007-0.008% -> accuracy 99.992-99.9993%
MODEL_PATH=""
SEED=42
KEEP_INTERMEDIATES=false
T2T_VERSION="${T2T_VERSION:-chm13v2.0}"
PBSIM3_BIN="${PBSIM3_BIN:-pbsim}"   # PBSIM3 required (yukiteruono/pbsim3); binary may be named pbsim or pbsim3

# --- Arg parsing ---
if [[ $# -gt 0 && "$1" != --* ]]; then
    OUTPUT_DIR="$1"; shift
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --reads)               N_READS="$2";        shift 2 ;;
        --length-mean)         LENGTH_MEAN="$2";    shift 2 ;;
        --length-sd)           LENGTH_SD="$2";      shift 2 ;;
        --accuracy)            ACCURACY="$2";       shift 2 ;;
        --model)               MODEL_PATH="$2";     shift 2 ;;
        --seed)                SEED="$2";           shift 2 ;;
        --keep-intermediates)  KEEP_INTERMEDIATES=true; shift ;;
        -h|--help)
            sed -n '3,60p' "$0" | grep '^#' | sed 's/^# \?//'
            exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

T2T_DIR="${OUTPUT_DIR}/t2t"
MASKED_FA="${T2T_DIR}/${T2T_VERSION}_masked.fa"
MASKED_FA_GZ="${MASKED_FA}.gz"
MASKED_FA_TMP=""

if [[ -f "${MASKED_FA}" ]]; then
    : # use as-is
elif [[ -f "${MASKED_FA_GZ}" ]]; then
    echo "Decompressing ${MASKED_FA_GZ}..."
    MASKED_FA_TMP="${T2T_DIR}/${T2T_VERSION}_masked_tmp.fa"
    pigz -dc "${MASKED_FA_GZ}" > "${MASKED_FA_TMP}" 2>/dev/null \
        || gzip -dc "${MASKED_FA_GZ}" > "${MASKED_FA_TMP}"
    MASKED_FA="${MASKED_FA_TMP}"
else
    echo "ERROR: masked genome not found (tried .fa and .fa.gz): ${MASKED_FA}" >&2
    echo "  Run download_non_rrna.sh and simulate_non_rrna.sh first." >&2
    exit 1
fi

# --- Auto-detect PBSIM3 model if not specified ---
# Note: ERRHMM-SEQUEL.model is from PacBio Sequel (not Sequel II CCS).
# No Sequel II CCS model exists in PBSIM3; ERRHMM-SEQUEL.model is the closest
# available. Combined with --accuracy-mean 0.9999 it approximates CCS error rates.
if [[ -z "${MODEL_PATH}" ]]; then
    _bin_dir="$(dirname "$(command -v "${PBSIM3_BIN}" 2>/dev/null || true)")"
    _env_root="${_bin_dir}/.."
    for _candidate in \
        "${_env_root}/data/ERRHMM-SEQUEL.model" \
        "${_env_root}/share/pbsim3/data/ERRHMM-SEQUEL.model" \
        "${_env_root}/share/pbsim/data/ERRHMM-SEQUEL.model"; do
        if [[ -f "${_candidate}" ]]; then
            MODEL_PATH="${_candidate}"
            break
        fi
    done
    if [[ -z "${MODEL_PATH}" ]]; then
        echo "ERROR: could not auto-detect PBSIM3 error model. Set --model explicitly." >&2
        echo "  Searched under: ${_env_root}" >&2
        exit 1
    fi
fi

N_LABEL="$(python3 -c "
n=${N_READS}
if n % 1000000 == 0: print(f'{n//1000000}M')
elif n % 1000 == 0:  print(f'{n//1000}K')
else:                print(str(n))
")"

OUTPUT_FASTQ="${OUTPUT_DIR}/non_rrna_pacbio_${N_LABEL}_T2T.fastq.gz"
PBSIM3_DIR="${OUTPUT_DIR}/pbsim3"
PBSIM3_PREFIX="${PBSIM3_DIR}/pbsim3_tmp"
mkdir -p "${PBSIM3_DIR}"

echo "============================================"
echo "PacBio non-rRNA simulation (PBSIM3)"
echo "============================================"
echo "  Masked genome: ${MASKED_FA}"
echo "  Target reads:  ${N_READS}"
echo "  Length mean:   ${LENGTH_MEAN} bp"
echo "  Length SD:     ${LENGTH_SD} bp"
echo "  Accuracy:      ${ACCURACY}"
echo "  Model:         ${MODEL_PATH}"
echo "  Seed:          ${SEED}"
echo "  Output:        ${OUTPUT_FASTQ}"
echo ""

if [[ -f "${OUTPUT_FASTQ}" ]]; then
    echo "Already exists: ${OUTPUT_FASTQ}"
    exit 0
fi

# Estimate genome size (excluding Ns) for depth calculation
echo "Estimating effective genome size..."
genome_bp=$(seqkit seq --seq-type dna "${MASKED_FA}" | tr -d 'nN\n' | wc -c)
echo "  Effective (non-N) bases: ${genome_bp}"

# depth = (N_READS * LENGTH_MEAN) / genome_bp, add 5% buffer
depth=$(python3 -c "print(f'{${N_READS} * ${LENGTH_MEAN} * 1.05 / ${genome_bp}:.4f}')")
echo "  PBSIM3 depth:  ${depth}x"
echo ""

echo "Running PBSIM3..."
"${PBSIM3_BIN}" \
    --strategy wgs \
    --method errhmm \
    --errhmm "${MODEL_PATH}" \
    --depth "${depth}" \
    --length-mean "${LENGTH_MEAN}" \
    --length-sd "${LENGTH_SD}" \
    --accuracy-mean "${ACCURACY}" \
    --seed "${SEED}" \
    --prefix "${PBSIM3_PREFIX}" \
    --genome "${MASKED_FA}"

echo ""
echo "Merging per-chromosome FASTQ files..."
cat "${PBSIM3_PREFIX}"_*.fq.gz > "${OUTPUT_FASTQ}.tmp.gz"

# Subsample to exactly N_READS if needed
actual=$(seqkit stats -T "${OUTPUT_FASTQ}.tmp.gz" | tail -1 | cut -f4)
echo "  Simulated reads: ${actual}, target: ${N_READS}"
if (( actual > N_READS )); then
    echo "  Subsampling to ${N_READS}..."
    seqkit sample -n "${N_READS}" --rand-seed "${SEED}" "${OUTPUT_FASTQ}.tmp.gz" \
        | gzip -c > "${OUTPUT_FASTQ}"
    rm "${OUTPUT_FASTQ}.tmp.gz"
else
    mv "${OUTPUT_FASTQ}.tmp.gz" "${OUTPUT_FASTQ}"
fi

if [[ "${KEEP_INTERMEDIATES}" == false ]]; then
    echo "Removing PBSIM3 intermediate files..."
    rm -rf "${PBSIM3_DIR}"
fi

[[ -n "${MASKED_FA_TMP}" && -f "${MASKED_FA_TMP}" ]] && rm -f "${MASKED_FA_TMP}"

final=$(seqkit stats -T "${OUTPUT_FASTQ}" | tail -1 | cut -f4)
echo ""
echo "Done."
echo "  Output: ${OUTPUT_FASTQ}"
echo "  Reads:  ${final}"
