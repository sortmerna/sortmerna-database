#!/bin/bash
#
# cluster_sequences.sh - Cluster rRNA sequences using VSEARCH
#
# Clusters SILVA and RFAM sequences at multiple identity thresholds
# by kingdom, generating a summary table in markdown format.
#
# For each clustering threshold, outputs:
#   clustered_XX.fasta              - centroid/seed sequences (the database)
#   clustered_XX.uc                 - VSEARCH cluster membership
#   clustered_XX_test_members.fasta - non-seed members (for simulating test reads)
#   clustered_XX_cluster_mapping.txt - member-to-seed mapping

set -euo pipefail

# Argument parsing — named flags first, then positional
FORCE=false
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
  --force) FORCE=true; shift ;;
  *) POSITIONAL+=("$1"); shift ;;
  esac
done

# Configuration
INPUT_DIR="${POSITIONAL[0]:-data}"
OUTPUT_DIR="${POSITIONAL[1]:-data/clustered}"
THREADS="${POSITIONAL[2]:-4}"
VERIFIED_DIR="${VERIFIED_DIR:-${INPUT_DIR}/verified}"
CLUSTERED_DIR="${CLUSTERED_DIR:-${OUTPUT_DIR}}"
UTILS_DIR="${UTILS_DIR:-${SMR_DB_ROOT_DIR}/scripts/utils}"

SILVA_SSU_VERSION="${SILVA_SSU_VERSION:?Please set SILVA_SSU_VERSION}"
SILVA_LSU_VERSION="${SILVA_LSU_VERSION:?Please set SILVA_LSU_VERSION}"
RFAM_VERSION="${RFAM_VERSION:-15.1}"

# Clustering thresholds to test (as percentages)
THRESHOLDS=(97 95 90 85)

echo "============================================"
echo "Sequence Clustering Script"
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Threads: ${THREADS}"
echo "Thresholds: ${THRESHOLDS[*]}"
echo "Force re-run vsearch: ${FORCE}"
echo "============================================"

# Check for required tools
if ! command -v vsearch &> /dev/null; then
  echo "Error: VSEARCH not found. Please install VSEARCH first."
  exit 1
fi

if ! command -v seqkit &> /dev/null; then
  echo "Error: seqkit not found. Please install seqkit first."
  exit 1
fi

echo "Using VSEARCH version: $(vsearch --version 2>&1 | head -1)"

# Create output directories
mkdir -p "${CLUSTERED_DIR}"

# TSV file accumulating results for the summary table (written by add_result, read by generate_summary.py)
RESULTS_TSV="${CLUSTERED_DIR}/clustering_results_2.tsv"
rm -f "${RESULTS_TSV}"

# Function to get sequence stats
get_seq_stats() {
  local fasta="$1"
  if [[ -f "${fasta}" ]]; then
  seqkit stats -T "${fasta}" | tail -1 | cut -f4,5
  else
  echo "0	0"
  fi
}

# Function to add result row to TSV
add_result() {
  local ref_db="$1"
  local db="$2"
  local kingdom="$3"
  local clustering="$4"
  local num_seqs="$5"
  local total_nt="$6"
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
  "${ref_db}" "${db}" "${kingdom}" "${clustering}" "${num_seqs}" "${total_nt}" \
  >> "${RESULTS_TSV}"
}

# Function to cluster sequences
# Outputs:
#   <output>.fasta          - centroid/seed sequences (the clustered database)
#   <output>.uc             - VSEARCH cluster membership file
#   <output>_test_members.fasta   - non-seed cluster members (for simulating test reads)
#   <output>_cluster_mapping.txt  - mapping of member sequence IDs to their seed IDs
cluster_sequences() {
  local input="$1"
  local threshold="$2"
  local output="$3"

  if [[ ! -f "${input}" ]] || [[ ! -s "${input}" ]]; then
  return 1
  fi

  local identity=$(echo "scale=2; ${threshold}/100" | bc)
  local uc_file="${output%.fasta}.uc"
  local test_members="${output%.fasta}_test_members.fasta"
  local cluster_mapping="${output%.fasta}_cluster_mapping.txt"

  if [[ -f "${uc_file}" ]] && [[ -s "${output}" ]] && [[ "${FORCE}" != "true" ]]; then
  echo "  .uc and .fasta exist, skipping vsearch: $(basename "${uc_file}")"
  else
  local log_file="${output%.fasta}.log"
  echo "  Clustering at ${threshold}%: $(basename "${output}")"
  echo "  vsearch --cluster_fast ${input} --id ${identity} --centroids ${output} --uc ${uc_file} --log ${log_file} --threads ${THREADS} --strand both --notrunclabels --sizeorder --sizeout --quiet"

  vsearch \
    --cluster_fast "${input}" \
    --id "${identity}" \
    --centroids "${output}" \
    --uc "${uc_file}" \
    --log "${log_file}" \
    --threads "${THREADS}" \
    --strand both \
    --notrunclabels \
    --sizeorder \
    --sizeout \
    --quiet
  fi

  # Parse .uc file: extract member IDs and cluster mapping in one pass
  local member_ids="${output%.fasta}_member_ids.tmp"
  python3 "${UTILS_DIR}/parse_uc.py" "${uc_file}" \
  --member-ids "${member_ids}" \
  --mapping "${cluster_mapping}"

  # Extract non-seed member sequences from original input
  if [[ -s "${member_ids}" ]]; then
  seqkit grep -f "${member_ids}" "${input}" -o "${test_members}"
  local member_count=$(seqkit stats -T "${test_members}" | tail -1 | cut -f4)
  echo "  Non-seed test members: ${member_count} sequences"
  else
  echo "  No non-seed members at this threshold"
  touch "${test_members}"
  fi

  rm -f "${member_ids}"

  return 0
}

# Process SILVA SSU (pre-verified, pre-split by verify_silva.sh)
echo ""
echo "============================================"
echo "Processing SILVA SSU sequences"
echo "============================================"

for domain in bacteria archaea eukaryota_nuclear eukaryota_mito eukaryota_chloro; do
  kingdom_file="${VERIFIED_DIR}/verified_ssu_${domain}.fasta"
  if [[ ! -f "${kingdom_file}" ]] || [[ ! -s "${kingdom_file}" ]]; then continue; fi

  stats=$(get_seq_stats "${kingdom_file}")
  num_seqs=$(echo "${stats}" | cut -f1)
  total_nt=$(echo "${stats}" | cut -f2)
  add_result "SILVA ${SILVA_SSU_VERSION}" "SSU Ref NR 99" "${domain}" "verified" "${num_seqs}" "${total_nt}"

  for threshold in "${THRESHOLDS[@]}"; do
  output="${CLUSTERED_DIR}/silva_ssu_${domain}_${threshold}.fasta"
  if cluster_sequences "${kingdom_file}" "${threshold}" "${output}"; then
    stats=$(get_seq_stats "${output}")
    num_seqs=$(echo "${stats}" | cut -f1)
    total_nt=$(echo "${stats}" | cut -f2)
    add_result "SILVA ${SILVA_SSU_VERSION}" "SSU Ref NR 99" "${domain}" "${threshold}%" "${num_seqs}" "${total_nt}"
  fi
  python3 "${UTILS_DIR}/check_leakage.py" "${output}" "${output%.fasta}_test_members.fasta"
  done
done

# Process SILVA LSU (pre-verified, pre-split by verify_silva.sh)
echo ""
echo "============================================"
echo "Processing SILVA LSU sequences"
echo "============================================"

for domain in bacteria archaea eukaryota; do
  kingdom_file="${VERIFIED_DIR}/verified_lsu_${domain}.fasta"
  if [[ ! -f "${kingdom_file}" ]] || [[ ! -s "${kingdom_file}" ]]; then continue; fi

  stats=$(get_seq_stats "${kingdom_file}")
  num_seqs=$(echo "${stats}" | cut -f1)
  total_nt=$(echo "${stats}" | cut -f2)
  add_result "SILVA ${SILVA_LSU_VERSION}" "LSU Ref NR 99" "${domain}" "verified" "${num_seqs}" "${total_nt}"

  for threshold in "${THRESHOLDS[@]}"; do
  output="${CLUSTERED_DIR}/silva_lsu_${domain}_${threshold}.fasta"
  if cluster_sequences "${kingdom_file}" "${threshold}" "${output}"; then
    stats=$(get_seq_stats "${output}")
    num_seqs=$(echo "${stats}" | cut -f1)
    total_nt=$(echo "${stats}" | cut -f2)
    add_result "SILVA ${SILVA_LSU_VERSION}" "LSU Ref NR 99" "${domain}" "${threshold}%" "${num_seqs}" "${total_nt}"
  fi
  python3 "${UTILS_DIR}/check_leakage.py" "${output}" "${output%.fasta}_test_members.fasta"
  done
done

# Process RFAM 5S
echo ""
echo "============================================"
echo "Processing RFAM 5S sequences"
echo "============================================"

RFAM_5S_FULL=$(find "${INPUT_DIR}" -name "*RF00001*full*.fa" -type f 2>/dev/null | head -1)
RFAM_5S_SEED=$(find "${INPUT_DIR}" -name "*RF00001*seed*.fa" -type f 2>/dev/null | head -1)

# 5S seed sequences
if [[ -n "${RFAM_5S_SEED}" && -f "${RFAM_5S_SEED}" ]]; then
  stats=$(get_seq_stats "${RFAM_5S_SEED}")
  num_seqs=$(echo "${stats}" | cut -f1)
  total_nt=$(echo "${stats}" | cut -f2)
  add_result "RFAM ${RFAM_VERSION}" "5S seed" "root" "100%" "${num_seqs}" "${total_nt}"
fi

# 5S full sequences
if [[ -n "${RFAM_5S_FULL}" && -f "${RFAM_5S_FULL}" ]]; then
  stats=$(get_seq_stats "${RFAM_5S_FULL}")
  num_seqs=$(echo "${stats}" | cut -f1)
  total_nt=$(echo "${stats}" | cut -f2)
  add_result "RFAM ${RFAM_VERSION}" "5S" "root" "100%" "${num_seqs}" "${total_nt}"

  # Cluster at each threshold
  for threshold in "${THRESHOLDS[@]}"; do
  output="${CLUSTERED_DIR}/rfam_5s_${threshold}.fasta"
  if cluster_sequences "${RFAM_5S_FULL}" "${threshold}" "${output}"; then
    stats=$(get_seq_stats "${output}")
    num_seqs=$(echo "${stats}" | cut -f1)
    total_nt=$(echo "${stats}" | cut -f2)
    add_result "RFAM ${RFAM_VERSION}" "5S" "root" "${threshold}%" "${num_seqs}" "${total_nt}"
  fi
  python3 "${UTILS_DIR}/check_leakage.py" "${output}" "${output%.fasta}_test_members.fasta"
  done
fi

# Process RFAM 5.8S
echo ""
echo "============================================"
echo "Processing RFAM 5.8S sequences"
echo "============================================"

RFAM_5_8S_FULL=$(find "${INPUT_DIR}" -name "*RF00002*full*.fa" -type f 2>/dev/null | head -1)
RFAM_5_8S_SEED=$(find "${INPUT_DIR}" -name "*RF00002*seed*.fa" -type f 2>/dev/null | head -1)

# 5.8S seed sequences
if [[ -n "${RFAM_5_8S_SEED}" && -f "${RFAM_5_8S_SEED}" ]]; then
  stats=$(get_seq_stats "${RFAM_5_8S_SEED}")
  num_seqs=$(echo "${stats}" | cut -f1)
  total_nt=$(echo "${stats}" | cut -f2)
  add_result "RFAM ${RFAM_VERSION}" "5.8S seed" "eukaryota" "100%" "${num_seqs}" "${total_nt}"
fi

# 5.8S full sequences
if [[ -n "${RFAM_5_8S_FULL}" && -f "${RFAM_5_8S_FULL}" ]]; then
  stats=$(get_seq_stats "${RFAM_5_8S_FULL}")
  num_seqs=$(echo "${stats}" | cut -f1)
  total_nt=$(echo "${stats}" | cut -f2)
  add_result "RFAM ${RFAM_VERSION}" "5.8S" "eukaryota" "100%" "${num_seqs}" "${total_nt}"

  # Cluster at each threshold
  for threshold in "${THRESHOLDS[@]}"; do
  output="${CLUSTERED_DIR}/rfam_5_8s_${threshold}.fasta"
  if cluster_sequences "${RFAM_5_8S_FULL}" "${threshold}" "${output}"; then
    stats=$(get_seq_stats "${output}")
    num_seqs=$(echo "${stats}" | cut -f1)
    total_nt=$(echo "${stats}" | cut -f2)
    add_result "RFAM ${RFAM_VERSION}" "5.8S" "eukaryota" "${threshold}%" "${num_seqs}" "${total_nt}"
  fi
  python3 "${UTILS_DIR}/check_leakage.py" "${output}" "${output%.fasta}_test_members.fasta"
  done
fi

# Generate markdown table
echo ""
echo "============================================"
echo "Clustering complete!"
echo "============================================"

TABLE_FILE="${CLUSTERED_DIR}/clustering_summary.html"
python3 "${UTILS_DIR}/generate_summary.py" "${RESULTS_TSV}" --output "${TABLE_FILE}" --thresholds "${THRESHOLDS[@]}" --silva-ssu-version "${SILVA_SSU_VERSION}" --silva-lsu-version "${SILVA_LSU_VERSION}" --rfam-version "${RFAM_VERSION}"

echo ""
echo "Summary table written to: ${TABLE_FILE}"
echo ""
cat "${TABLE_FILE}"

echo ""
echo "Next step: Run build_sortmerna_index.sh to build SortMeRNA indices"
