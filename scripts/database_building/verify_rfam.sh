#!/bin/bash
#
# verify_rfam.sh - Verify Rfam rRNA sequences with Infernal cmsearch
#
# Verifies 5S (RF00001) and 5.8S (RF00002) sequences (full and seed sets)
# against their corresponding Rfam covariance models using --hmmonly --cut_ga,
# and outputs:
#   verified_<gene>_<type>.fasta  - sequences with an above-threshold hit
#   flagged_<gene>_<type>.fasta   - sequences with no qualifying hit
#   cmsearch_log_<gene>_<type>.tsv - hit coordinates and scores for kept sequences
#   <gene>_<type>_cmsearch.tblout - raw cmsearch output (kept for auditing)
#
# Requires: Infernal (cmsearch), seqkit, Python 3
# Usage: bash verify_rfam.sh [input_dir [output_dir [threads]]]

set -euo pipefail

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
  *) POSITIONAL+=("$1"); shift ;;
  esac
done

INPUT_DIR="${POSITIONAL[0]:-data}"
OUTPUT_DIR="${POSITIONAL[1]:-data/verified_rfam}"
THREADS="${POSITIONAL[2]:-4}"
RFAM_DIR="${RFAM_DIR:-${INPUT_DIR}/rfam}"
CMS_DIR="${CMS_DIR:-${INPUT_DIR}/cms}"
UTILS_DIR="${UTILS_DIR:-${SMR_DB_ROOT_DIR}/scripts/utils}"

RFAM_VERSION="${RFAM_VERSION:?Please set RFAM_VERSION}"

echo "============================================"
echo "Rfam Verification Script"
echo "Version: ${RFAM_VERSION}"
echo "Rfam sequences: ${RFAM_DIR}"
echo "CMs: ${CMS_DIR}"
echo "Output: ${OUTPUT_DIR}"
echo "Threads: ${THREADS}"
echo "============================================"

if ! command -v cmsearch &>/dev/null; then
  echo "Error: Infernal (cmsearch) not found. Please install Infernal first."
  exit 1
fi
if ! command -v seqkit &>/dev/null; then
  echo "Error: seqkit not found. Please install seqkit first."
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

STATS_TSV="${OUTPUT_DIR}/verification_stats.tsv"
rm -f "${STATS_TSV}"

# ── helpers ───────────────────────────────────────────────────────────────────

# Return path to base or base.gz, whichever exists; empty string if neither
find_rfam_fasta() {
  local base="$1"
  if [[ -f "${base}" ]]; then
    echo "${base}"
  elif [[ -f "${base}.gz" ]]; then
    echo "${base}.gz"
  fi
}

convert_u_to_t() {
  local input="$1" output="$2"
  if [[ ! -f "${output}" ]]; then
    echo "  Converting U -> T: $(basename "${output}")"
    seqkit seq --rna2dna "${input}" -o "${output}"
  else
    echo "  Already converted: $(basename "${output}")"
  fi
}

run_cmsearch() {
  local cm="$1" fasta="$2" tblout="$3"
  if [[ ! -f "${tblout}" ]]; then
    echo "  cmsearch: $(basename "${fasta}") vs $(basename "${cm}")"
    echo "  cmd: cmsearch --cpu ${THREADS} --tblout $(basename "${tblout}") --cut_ga --noali $(basename "${cm}") $(basename "${fasta}")"
    cmsearch --cpu "${THREADS}" --tblout "${tblout}" --cut_ga --noali \
      "${cm}" "${fasta}" > /dev/null
  else
    echo "  tblout exists, skipping: $(basename "${tblout}")"
  fi
}

filter_hits() {
  local fasta="$1" tblout="$2" clean="$3" flagged="$4" log="$5" gene="$6" type="$7"
  python3 "${UTILS_DIR}/parse_cmsearch.py" \
    --fasta   "${fasta}"    \
    --tblout  "${tblout}"   \
    --clean   "${clean}"    \
    --flagged "${flagged}"  \
    --log     "${log}"      \
    --gene    "${gene}"     \
    --domain  "${type}"     \
    --stats   "${STATS_TSV}"
}

verify_family() {
  local gene="$1" type="$2" fasta="$3" cm="$4"
  if [[ ! -f "${fasta}" ]] || [[ ! -s "${fasta}" ]]; then
    echo "  Skipping ${gene} ${type}: file not found or empty"
    return
  fi
  echo ""
  echo "  Verifying ${gene} ${type}..."
  local tblout="${OUTPUT_DIR}/${gene}_${type}_cmsearch.tblout"
  run_cmsearch "${cm}" "${fasta}" "${tblout}"
  filter_hits \
    "${fasta}" \
    "${tblout}" \
    "${OUTPUT_DIR}/verified_${gene}_${type}.fasta" \
    "${OUTPUT_DIR}/flagged_${gene}_${type}.fasta" \
    "${OUTPUT_DIR}/cmsearch_log_${gene}_${type}.tsv" \
    "${gene}" "${type}"
}

# ── 5S rRNA (RF00001) ─────────────────────────────────────────────────────────

echo ""
echo "============================================"
echo "Processing 5S rRNA (RF00001)"
echo "============================================"

RF5S_FULL_DNA="${OUTPUT_DIR}/rfam_5s_full_dna.fasta"
RF5S_SEED_DNA="${OUTPUT_DIR}/rfam_5s_seed_dna.fasta"

RF5S_FULL_SRC=$(find_rfam_fasta "${RFAM_DIR}/RF00001_5S_rRNA_full.fa")
RF5S_SEED_SRC=$(find_rfam_fasta "${RFAM_DIR}/RF00001_5S_rRNA_seed.fa")
[[ -n "${RF5S_FULL_SRC}" ]] && convert_u_to_t "${RF5S_FULL_SRC}" "${RF5S_FULL_DNA}" || echo "  Skipping 5S full: file not found in ${RFAM_DIR}"
[[ -n "${RF5S_SEED_SRC}" ]] && convert_u_to_t "${RF5S_SEED_SRC}" "${RF5S_SEED_DNA}" || echo "  Skipping 5S seed: file not found in ${RFAM_DIR}"

verify_family "5s" "full" "${RF5S_FULL_DNA}" "${CMS_DIR}/RF00001.cm"
verify_family "5s" "seed" "${RF5S_SEED_DNA}" "${CMS_DIR}/RF00001.cm"

# ── 5.8S rRNA (RF00002) ───────────────────────────────────────────────────────

echo ""
echo "============================================"
echo "Processing 5.8S rRNA (RF00002)"
echo "============================================"

RF58S_FULL_DNA="${OUTPUT_DIR}/rfam_58s_full_dna.fasta"
RF58S_SEED_DNA="${OUTPUT_DIR}/rfam_58s_seed_dna.fasta"

RF58S_FULL_SRC=$(find_rfam_fasta "${RFAM_DIR}/RF00002_5_8S_rRNA_full.fa")
RF58S_SEED_SRC=$(find_rfam_fasta "${RFAM_DIR}/RF00002_5_8S_rRNA_seed.fa")
[[ -n "${RF58S_FULL_SRC}" ]] && convert_u_to_t "${RF58S_FULL_SRC}" "${RF58S_FULL_DNA}" || echo "  Skipping 5.8S full: file not found in ${RFAM_DIR}"
[[ -n "${RF58S_SEED_SRC}" ]] && convert_u_to_t "${RF58S_SEED_SRC}" "${RF58S_SEED_DNA}" || echo "  Skipping 5.8S seed: file not found in ${RFAM_DIR}"

verify_family "5.8s" "full" "${RF58S_FULL_DNA}" "${CMS_DIR}/RF00002.cm"
verify_family "5.8s" "seed" "${RF58S_SEED_DNA}" "${CMS_DIR}/RF00002.cm"

# ── summary ───────────────────────────────────────────────────────────────────

echo ""
echo "============================================"
echo "Verification complete!"
echo "============================================"

HTML_SUMMARY="${OUTPUT_DIR}/verification_summary.html"
python3 "${UTILS_DIR}/generate_verification_summary.py" \
  "${STATS_TSV}" \
  --output "${HTML_SUMMARY}" \
  --title "Rfam Verification Summary" \
  --version "Rfam ${RFAM_VERSION}" \
  --tool "cmsearch --cut_ga"

echo ""
echo "Verified sequences: ${OUTPUT_DIR}/verified_*.fasta"
echo "Flagged sequences:  ${OUTPUT_DIR}/flagged_*.fasta"
echo "Hit logs:           ${OUTPUT_DIR}/cmsearch_log_*.tsv"
echo "HTML summary:       ${HTML_SUMMARY}"
echo ""
echo "Next step: Run cluster_sequences.sh to cluster verified sequences"
