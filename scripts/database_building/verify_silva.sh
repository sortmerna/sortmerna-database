#!/bin/bash
#
# verify_silva.sh - Verify SILVA rRNA sequences with Infernal cmsearch
#
# Splits SILVA SSU and LSU by domain, verifies each domain against the
# corresponding Rfam covariance model using --cut_ga, and outputs:
#   verified_<gene>_<domain>.fasta  — sequences with an above-threshold hit
#   flagged_<gene>_<domain>.fasta   — sequences with no qualifying hit (dropped by default)
#   cmsearch_log_<gene>_<domain>.tsv — hit coordinates and scores for kept sequences
#   <gene>_<domain>_cmsearch.tblout — raw cmsearch output (kept for auditing)
#
# CMs used:
#   SSU bacteria:           RF00177
#   SSU archaea:            RF01959
#   SSU eukaryota nuclear:  RF01960
#   SSU eukaryota mito:     RF02545
#   SSU eukaryota chloro:   RF00177  (plastids are bacterial-derived)
#   LSU bacteria:           RF02541
#   LSU archaea:            RF02540
#   LSU eukaryota:          RF02543
#
# Requires: Infernal (cmsearch), seqkit, Python 3
# Usage: bash verify_silva.sh [input_dir [output_dir [threads]]]

set -euo pipefail

# store all input arguments as an array for easier passing to Python scripts
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
  *) POSITIONAL+=("$1"); shift ;;
  esac
done

INPUT_DIR="${POSITIONAL[0]:-data}"
OUTPUT_DIR="${POSITIONAL[1]:-data/verified}"
THREADS="${POSITIONAL[2]:-4}"
SILVA_DIR="${SILVA_DIR:-${INPUT_DIR}/silva}"
CMS_DIR="${CMS_DIR:-${INPUT_DIR}/cms}"
UTILS_DIR="${UTILS_DIR:-${SMR_DB_ROOT_DIR}/scripts/utils}"

SILVA_SSU_VERSION="${SILVA_SSU_VERSION:?Please set SILVA_SSU_VERSION}"
SILVA_SSU_PATH="${SILVA_SSU_PATH:?Please set SILVA_SSU_PATH}"
SILVA_LSU_VERSION="${SILVA_LSU_VERSION:?Please set SILVA_LSU_VERSION}"
SILVA_LSU_PATH="${SILVA_LSU_PATH:?Please set SILVA_LSU_PATH}"

SILVA_SSU="${SILVA_DIR}/$(basename "${SILVA_SSU_PATH%.gz}")"
SILVA_LSU="${SILVA_DIR}/$(basename "${SILVA_LSU_PATH%.gz}")"

echo "============================================"
echo "SILVA Verification Script"
echo "SSU: ${SILVA_SSU}"
echo "LSU: ${SILVA_LSU}"
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

# ── helpers ──────────────────────────────────────────────────────────────────

convert_u_to_t() {
  local input="$1" output="$2"
  if [[ ! -f "${output}" ]]; then
  echo "  Converting U -> T: $(basename "${output}")"
  seqkit seq --rna2dna "${input}" -o "${output}"
  else
  echo "  Already converted: $(basename "${output}")"
  fi
}

split_by_domain() {
  local input="$1" prefix="$2"
  local total
  total=$(seqkit stats -T "${input}" | tail -1 | cut -f4)
  echo "  Splitting ${total} sequences by domain..."
  seqkit grep -r --by-name -p " Bacteria;"  "${input}" -o "${prefix}_bacteria.fasta"  2>/dev/null || true
  seqkit grep -r --by-name -p " Archaea;"   "${input}" -o "${prefix}_archaea.fasta"   2>/dev/null || true
  seqkit grep -r --by-name -p " Eukaryota;" "${input}" -o "${prefix}_eukaryota.fasta" 2>/dev/null || true
}

split_euk_organellar() {
  local euk_fasta="$1" prefix="$2"
  seqkit grep -r --by-name -p "Mitochondria" "${euk_fasta}" -o "${prefix}_mito.fasta"  2>/dev/null || true
  seqkit grep -r --by-name -p "Chloroplast"  "${euk_fasta}" -o "${prefix}_chloro.fasta" 2>/dev/null || true

  local organellar_ids="${prefix}_organellar_ids.tmp"
  cat "${prefix}_mito.fasta" "${prefix}_chloro.fasta" \
  | seqkit seq --name --only-id > "${organellar_ids}" 2>/dev/null || true

  if [[ -s "${organellar_ids}" ]]; then
  seqkit grep -v -f "${organellar_ids}" "${euk_fasta}" -o "${prefix}_nuclear.fasta"
  else
  cp "${euk_fasta}" "${prefix}_nuclear.fasta"
  fi
  rm -f "${organellar_ids}"
}

run_cmsearch() {
  local cm="$1" fasta="$2" tblout="$3"
  if [[ ! -f "${tblout}" ]]; then
  echo "  cmsearch: $(basename "${fasta}") vs $(basename "${cm}")"
  echo "  cmd: cmsearch --hmmonly --cpu ${THREADS} --tblout $(basename "${tblout}") --cut_ga --noali $(basename "${cm}") $(basename "${fasta}")"
  cmsearch --hmmonly --cpu "${THREADS}" --tblout "${tblout}" --cut_ga --noali \
    "${cm}" "${fasta}" > /dev/null
  else
  echo "  tblout exists, skipping: $(basename "${tblout}")"
  fi
}

filter_hits() {
  local fasta="$1" tblout="$2" clean="$3" flagged="$4" log="$5" gene="$6" domain="$7"
  python3 "${UTILS_DIR}/parse_cmsearch.py" \
  --fasta   "${fasta}"    \
  --tblout  "${tblout}"   \
  --clean   "${clean}"    \
  --flagged "${flagged}"  \
  --log     "${log}"      \
  --gene    "${gene}"     \
  --domain  "${domain}"   \
  --stats   "${STATS_TSV}"
}

verify_domain() {
  local gene="$1" domain="$2" fasta="$3" cm="$4"
  if [[ ! -f "${fasta}" ]] || [[ ! -s "${fasta}" ]]; then return; fi
  echo ""
  echo "  Verifying ${gene} ${domain}..."
  local tblout="${OUTPUT_DIR}/${gene}_${domain}_cmsearch.tblout"
  run_cmsearch "${cm}" "${fasta}" "${tblout}"
  filter_hits \
  "${fasta}" \
  "${tblout}" \
  "${OUTPUT_DIR}/verified_${gene}_${domain}.fasta" \
  "${OUTPUT_DIR}/flagged_${gene}_${domain}.fasta" \
  "${OUTPUT_DIR}/cmsearch_log_${gene}_${domain}.tsv" \
  "${gene}" "${domain}"
}

# ── SSU ──────────────────────────────────────────────────────────────────────

echo ""
echo "============================================"
echo "Processing SILVA SSU (${SILVA_SSU_VERSION})"
echo "============================================"

SSU_DNA="${OUTPUT_DIR}/silva_ssu_dna.fasta"
convert_u_to_t "${SILVA_SSU}" "${SSU_DNA}"

SSU_DOM="${OUTPUT_DIR}/silva_ssu_dom"
split_by_domain "${SSU_DNA}" "${SSU_DOM}"

# Split eukaryota into nuclear / mito / chloro
EUK_SSU_PREFIX="${OUTPUT_DIR}/silva_ssu_euk"
if [[ -f "${SSU_DOM}_eukaryota.fasta" ]]; then
  split_euk_organellar "${SSU_DOM}_eukaryota.fasta" "${EUK_SSU_PREFIX}"
fi

verify_domain "ssu" "bacteria" "${SSU_DOM}_bacteria.fasta" "${CMS_DIR}/RF00177.cm"
verify_domain "ssu" "archaea"  "${SSU_DOM}_archaea.fasta" "${CMS_DIR}/RF01959.cm"
verify_domain "ssu" "eukaryota_nuclear" "${EUK_SSU_PREFIX}_nuclear.fasta" "${CMS_DIR}/RF01960.cm"
verify_domain "ssu" "eukaryota_mito" "${EUK_SSU_PREFIX}_mito.fasta" "${CMS_DIR}/RF02545.cm"
verify_domain "ssu" "eukaryota_chloro" "${EUK_SSU_PREFIX}_chloro.fasta" "${CMS_DIR}/RF00177.cm"

# ── LSU ──────────────────────────────────────────────────────────────────────

echo ""
echo "============================================"
echo "Processing SILVA LSU (${SILVA_LSU_VERSION})"
echo "============================================"

LSU_DNA="${OUTPUT_DIR}/silva_lsu_dna.fasta"
convert_u_to_t "${SILVA_LSU}" "${LSU_DNA}"

LSU_DOM="${OUTPUT_DIR}/silva_lsu_dom"
split_by_domain "${LSU_DNA}" "${LSU_DOM}"

verify_domain "lsu" "bacteria"  "${LSU_DOM}_bacteria.fasta"  "${CMS_DIR}/RF02541.cm"
verify_domain "lsu" "archaea"   "${LSU_DOM}_archaea.fasta"   "${CMS_DIR}/RF02540.cm"
verify_domain "lsu" "eukaryota" "${LSU_DOM}_eukaryota.fasta" "${CMS_DIR}/RF02543.cm"

# ── summary ──────────────────────────────────────────────────────────────────

echo ""
echo "============================================"
echo "Verification complete!"
echo "============================================"

HTML_SUMMARY="${OUTPUT_DIR}/verification_summary.html"
python3 "${UTILS_DIR}/generate_verification_summary.py" \
  "${STATS_TSV}" \
  --output "${HTML_SUMMARY}" \
  --silva-ssu-version "${SILVA_SSU_VERSION}" \
  --silva-lsu-version "${SILVA_LSU_VERSION}"

echo ""
echo "Verified sequences: ${OUTPUT_DIR}/verified_*.fasta"
echo "Flagged sequences:  ${OUTPUT_DIR}/flagged_*.fasta"
echo "Hit logs:           ${OUTPUT_DIR}/cmsearch_log_*.tsv"
echo "HTML summary:       ${HTML_SUMMARY}"
echo ""
echo "Flagged files contain sequences SILVA kept that Rfam disagrees with."
echo "They are excluded from downstream clustering by default."
echo ""
echo "Next step: Run cluster_sequences.sh to cluster verified sequences"
