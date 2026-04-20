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

# Configuration
INPUT_DIR="${1:-data}"
OUTPUT_DIR="${2:-data/clustered}"
THREADS="${3:-4}"
SILVA_DIR="${SILVA_DIR:-${INPUT_DIR}/silva}"
CLUSTERED_DIR="${CLUSTERED_DIR:-${OUTPUT_DIR}}"
UTILS_DIR="${UTILS_DIR:-${SMR_DB_ROOT_DIR}/scripts/utils}"

SILVA_VERSION="${SILVA_VERSION:-138.2}"
RFAM_VERSION="${RFAM_VERSION:-15.1}"

# Clustering thresholds to test (as percentages)
THRESHOLDS=(97 95 90 85)

echo "============================================"
echo "Sequence Clustering Script"
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Threads: ${THREADS}"
echo "Thresholds: ${THRESHOLDS[*]}"
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
mkdir -p "${SILVA_DIR}/by_kingdom"

# TSV file accumulating results for the summary table (written by add_result, read by generate_summary.py)
RESULTS_TSV="${CLUSTERED_DIR}/clustering_results.tsv"
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

# Function to split SILVA file by kingdom
split_by_kingdom() {
  local input="$1"
  local output_prefix="$2"

  local total
  total=$(seqkit stats -T "${input}" | tail -1 | cut -f4)
  echo "Splitting by kingdom: ${input} (${total} sequences)"

  # Extract Archaea
  seqkit grep -r --by-name -p " Archaea;" "${input}" -o "${output_prefix}_archaea.fasta" 2>/dev/null || true

  # Extract Bacteria
  seqkit grep -r --by-name -p " Bacteria;" "${input}" -o "${output_prefix}_bacteria.fasta" 2>/dev/null || true

  # Extract Eukaryota
  seqkit grep -r --by-name -p " Eukaryota;" "${input}" -o "${output_prefix}_eukaryota.fasta" 2>/dev/null || true

  # Verify split is complete
  local n_archaea n_bacteria n_eukaryota split_sum
  n_archaea=$(seqkit stats -T "${output_prefix}_archaea.fasta" 2>/dev/null | tail -1 | cut -f4)
  n_bacteria=$(seqkit stats -T "${output_prefix}_bacteria.fasta" 2>/dev/null | tail -1 | cut -f4)
  n_eukaryota=$(seqkit stats -T "${output_prefix}_eukaryota.fasta" 2>/dev/null | tail -1 | cut -f4)
  split_sum=$(( n_archaea + n_bacteria + n_eukaryota ))

  if [[ "${split_sum}" -eq "${total}" ]]; then
    echo "  Kingdom split OK: archaea=${n_archaea}, bacteria=${n_bacteria}, eukaryota=${n_eukaryota}"
  else
    echo "  Warning: kingdom split mismatch — input=${total}, split sum=${split_sum} (archaea=${n_archaea}, bacteria=${n_bacteria}, eukaryota=${n_eukaryota})"
    echo "  $(( total - split_sum )) sequences not assigned to any kingdom"
  fi
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

  if [[ -f "${uc_file}" ]]; then
    echo "  .uc file exists, skipping vsearch: $(basename "${uc_file}")"
  else
    echo "  Clustering at ${threshold}%: $(basename "${output}")"
    echo "  vsearch --cluster_fast ${input} --id ${identity} --centroids ${output} --uc ${uc_file} --threads ${THREADS} --strand both --notrunclabels --quiet"

    vsearch \
      --cluster_fast "${input}" \
      --id "${identity}" \
      --centroids "${output}" \
      --uc "${uc_file}" \
      --threads "${THREADS}" \
      --strand both \
      --notrunclabels \
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

# Process SILVA SSU
echo ""
echo "============================================"
echo "Processing SILVA SSU sequences"
echo "============================================"

SILVA_SSU="${SILVA_DIR}/SILVA_${SILVA_VERSION}_SSURef_NR99_tax_silva.fasta"
if [[ -f "${SILVA_SSU}" ]]; then
  SSU_PREFIX="${SILVA_DIR}/by_kingdom/silva_ssu"
  split_by_kingdom "${SILVA_SSU}" "${SSU_PREFIX}"

  for kingdom in archaea bacteria eukaryota; do
    kingdom_file="${SSU_PREFIX}_${kingdom}.fasta"
    if [[ -f "${kingdom_file}" ]] && [[ -s "${kingdom_file}" ]]; then
      # Add 99% (SILVA NR99 unclustered baseline) stats
      stats=$(get_seq_stats "${kingdom_file}")
      num_seqs=$(echo "${stats}" | cut -f1)
      total_nt=$(echo "${stats}" | cut -f2)
      add_result "SILVA ${SILVA_VERSION}" "SSU Ref NR 99" "${kingdom}" "99%" "${num_seqs}" "${total_nt}"

      # Cluster at each threshold
      for threshold in "${THRESHOLDS[@]}"; do
        output="${CLUSTERED_DIR}/silva_ssu_${kingdom}_${threshold}.fasta"
        if cluster_sequences "${kingdom_file}" "${threshold}" "${output}"; then
          stats=$(get_seq_stats "${output}")
          num_seqs=$(echo "${stats}" | cut -f1)
          total_nt=$(echo "${stats}" | cut -f2)
          add_result "SILVA ${SILVA_VERSION}" "SSU Ref NR 99" "${kingdom}" "${threshold}%" "${num_seqs}" "${total_nt}"
        fi
        python3 "${UTILS_DIR}/check_leakage.py" "${output}" "${output%.fasta}_test_members.fasta"
      done
    fi
  done
else
  echo "Warning: SILVA SSU file not found"
fi

# Process SILVA LSU
echo ""
echo "============================================"
echo "Processing SILVA LSU sequences"
echo "============================================"

SILVA_LSU="${SILVA_DIR}/SILVA_${SILVA_VERSION}_LSURef_NR99_tax_silva.fasta"
if [[ -f "${SILVA_LSU}" ]]; then
  LSU_PREFIX="${SILVA_DIR}/by_kingdom/silva_lsu"
  split_by_kingdom "${SILVA_LSU}" "${LSU_PREFIX}"

  for kingdom in archaea bacteria eukaryota; do
    kingdom_file="${LSU_PREFIX}_${kingdom}.fasta"
    if [[ -f "${kingdom_file}" ]] && [[ -s "${kingdom_file}" ]]; then
      # Add 99% (SILVA NR99 unclustered baseline) stats
      stats=$(get_seq_stats "${kingdom_file}")
      num_seqs=$(echo "${stats}" | cut -f1)
      total_nt=$(echo "${stats}" | cut -f2)
      add_result "SILVA ${SILVA_VERSION}" "LSU Ref NR 99" "${kingdom}" "99%" "${num_seqs}" "${total_nt}"

      # Cluster at each threshold
      for threshold in "${THRESHOLDS[@]}"; do
        output="${CLUSTERED_DIR}/silva_lsu_${kingdom}_${threshold}.fasta"
        if cluster_sequences "${kingdom_file}" "${threshold}" "${output}"; then
          stats=$(get_seq_stats "${output}")
          num_seqs=$(echo "${stats}" | cut -f1)
          total_nt=$(echo "${stats}" | cut -f2)
          add_result "SILVA ${SILVA_VERSION}" "LSU Ref NR 99" "${kingdom}" "${threshold}%" "${num_seqs}" "${total_nt}"
        fi
        python3 "${UTILS_DIR}/check_leakage.py" "${output}" "${output%.fasta}_test_members.fasta"
      done
    fi
  done
else
  echo "Warning: SILVA LSU file not found"
fi

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
python3 "${UTILS_DIR}/generate_summary.py" "${RESULTS_TSV}" --output "${TABLE_FILE}" --thresholds "${THRESHOLDS[@]}"

echo ""
echo "Summary table written to: ${TABLE_FILE}"
echo ""
cat "${TABLE_FILE}"

echo ""
echo "Next step: Run build_sortmerna_index.sh to build SortMeRNA indices"
