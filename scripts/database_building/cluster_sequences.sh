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

# SILVA version (update as needed)
SILVA_VERSION="138.2"

# RFAM version
RFAM_VERSION="15.1"

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
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/by_kingdom"

# Array to store results for markdown table
declare -a RESULTS=()
ROW_NUM=0

# Function to get sequence stats
get_seq_stats() {
    local fasta="$1"
    if [[ -f "${fasta}" ]]; then
        seqkit stats -T "${fasta}" | tail -1 | cut -f4,5
    else
        echo "0	0"
    fi
}

# Function to add result to table
add_result() {
    local ref_db="$1"
    local db="$2"
    local kingdom="$3"
    local clustering="$4"
    local num_seqs="$5"
    local total_nt="$6"

    ROW_NUM=$((ROW_NUM + 1))

    # Calculate log10(#seq)
    local log10_seq=""
    if [[ "${num_seqs}" -gt 0 ]]; then
        log10_seq=$(echo "scale=2; l(${num_seqs})/l(10)" | bc -l)
    fi

    RESULTS+=("${ROW_NUM}|${ref_db}|${db}|${kingdom}|${clustering}|${num_seqs}|${total_nt}|${log10_seq}")
}

# Function to split SILVA file by kingdom
split_by_kingdom() {
    local input="$1"
    local output_prefix="$2"

    echo "Splitting by kingdom: ${input}"

    # Extract Archaea
    seqkit grep -r -p "^Archaea" "${input}" -o "${output_prefix}_archaea.fasta" 2>/dev/null || true

    # Extract Bacteria
    seqkit grep -r -p "^Bacteria" "${input}" -o "${output_prefix}_bacteria.fasta" 2>/dev/null || true

    # Extract Eukaryota
    seqkit grep -r -p "^Eukary" "${input}" -o "${output_prefix}_eukaryota.fasta" 2>/dev/null || true
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

    if [[ -f "${output}" ]]; then
        echo "  Output exists, skipping: ${output}"
        return 0
    fi

    echo "  Clustering at ${threshold}%: $(basename "${output}")"

    vsearch \
        --cluster_fast "${input}" \
        --id "${identity}" \
        --centroids "${output}" \
        --uc "${uc_file}" \
        --threads "${THREADS}" \
        --strand both \
        --notrunclabels \
        --quiet

    # Extract non-seed member IDs from .uc file (H = hit/member records)
    # .uc columns: type, cluster#, length, %id, strand, _, _, cigar, query_label, target_label
    local member_ids="${output%.fasta}_member_ids.tmp"
    awk -F'\t' '$1 == "H" { print $9 }' "${uc_file}" > "${member_ids}"

    # Create cluster mapping: member_id <tab> seed_id
    {
        echo -e "member_id\tseed_id"
        awk -F'\t' '$1 == "H" { print $9 "\t" $10 }' "${uc_file}"
    } > "${cluster_mapping}"

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

SILVA_SSU=$(find "${INPUT_DIR}" -name "*SSU*.fasta" -type f 2>/dev/null | head -1)
if [[ -n "${SILVA_SSU}" && -f "${SILVA_SSU}" ]]; then
    SSU_PREFIX="${OUTPUT_DIR}/by_kingdom/silva_ssu"
    split_by_kingdom "${SILVA_SSU}" "${SSU_PREFIX}"

    for kingdom in archaea bacteria eukaryota; do
        kingdom_file="${SSU_PREFIX}_${kingdom}.fasta"
        if [[ -f "${kingdom_file}" ]] && [[ -s "${kingdom_file}" ]]; then
            # Add 100% (unclustered) stats
            stats=$(get_seq_stats "${kingdom_file}")
            num_seqs=$(echo "${stats}" | cut -f1)
            total_nt=$(echo "${stats}" | cut -f2)
            add_result "SILVA ${SILVA_VERSION}" "SSU Ref NR 99" "${kingdom}" "100%" "${num_seqs}" "${total_nt}"

            # Cluster at each threshold
            for threshold in "${THRESHOLDS[@]}"; do
                output="${OUTPUT_DIR}/silva_ssu_${kingdom}_${threshold}.fasta"
                if cluster_sequences "${kingdom_file}" "${threshold}" "${output}"; then
                    stats=$(get_seq_stats "${output}")
                    num_seqs=$(echo "${stats}" | cut -f1)
                    total_nt=$(echo "${stats}" | cut -f2)
                    add_result "SILVA ${SILVA_VERSION}" "SSU Ref NR 99" "${kingdom}" "${threshold}%" "${num_seqs}" "${total_nt}"
                fi
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

SILVA_LSU=$(find "${INPUT_DIR}" -name "*LSU*.fasta" -type f 2>/dev/null | head -1)
if [[ -n "${SILVA_LSU}" && -f "${SILVA_LSU}" ]]; then
    LSU_PREFIX="${OUTPUT_DIR}/by_kingdom/silva_lsu"
    split_by_kingdom "${SILVA_LSU}" "${LSU_PREFIX}"

    for kingdom in archaea bacteria eukaryota; do
        kingdom_file="${LSU_PREFIX}_${kingdom}.fasta"
        if [[ -f "${kingdom_file}" ]] && [[ -s "${kingdom_file}" ]]; then
            # Add 100% (unclustered) stats
            stats=$(get_seq_stats "${kingdom_file}")
            num_seqs=$(echo "${stats}" | cut -f1)
            total_nt=$(echo "${stats}" | cut -f2)
            add_result "SILVA ${SILVA_VERSION}" "LSU Ref NR 99" "${kingdom}" "100%" "${num_seqs}" "${total_nt}"

            # Cluster at each threshold
            for threshold in "${THRESHOLDS[@]}"; do
                output="${OUTPUT_DIR}/silva_lsu_${kingdom}_${threshold}.fasta"
                if cluster_sequences "${kingdom_file}" "${threshold}" "${output}"; then
                    stats=$(get_seq_stats "${output}")
                    num_seqs=$(echo "${stats}" | cut -f1)
                    total_nt=$(echo "${stats}" | cut -f2)
                    add_result "SILVA ${SILVA_VERSION}" "LSU Ref NR 99" "${kingdom}" "${threshold}%" "${num_seqs}" "${total_nt}"
                fi
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
        output="${OUTPUT_DIR}/rfam_5s_${threshold}.fasta"
        if cluster_sequences "${RFAM_5S_FULL}" "${threshold}" "${output}"; then
            stats=$(get_seq_stats "${output}")
            num_seqs=$(echo "${stats}" | cut -f1)
            total_nt=$(echo "${stats}" | cut -f2)
            add_result "RFAM ${RFAM_VERSION}" "5S" "root" "${threshold}%" "${num_seqs}" "${total_nt}"
        fi
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
        output="${OUTPUT_DIR}/rfam_5_8s_${threshold}.fasta"
        if cluster_sequences "${RFAM_5_8S_FULL}" "${threshold}" "${output}"; then
            stats=$(get_seq_stats "${output}")
            num_seqs=$(echo "${stats}" | cut -f1)
            total_nt=$(echo "${stats}" | cut -f2)
            add_result "RFAM ${RFAM_VERSION}" "5.8S" "eukaryota" "${threshold}%" "${num_seqs}" "${total_nt}"
        fi
    done
fi

# Generate markdown table
echo ""
echo "============================================"
echo "Clustering complete!"
echo "============================================"

TABLE_FILE="${OUTPUT_DIR}/clustering_summary.md"

{
    echo "# Clustering Summary"
    echo ""
    echo "| # | Ref DB | DB | kingdom | clustering | #sequences | total size (nt) | log10(#seq) |"
    echo "|---|--------|-----|---------|------------|------------|-----------------|-------------|"

    for row in "${RESULTS[@]}"; do
        IFS='|' read -r num ref_db db kingdom clustering num_seqs total_nt log10_seq <<< "${row}"
        echo "| ${num} | ${ref_db} | ${db} | ${kingdom} | ${clustering} | ${num_seqs} | ${total_nt} | ${log10_seq} |"
    done
} > "${TABLE_FILE}"

echo ""
echo "Summary table written to: ${TABLE_FILE}"
echo ""
cat "${TABLE_FILE}"

echo ""
echo "Next step: Run build_sortmerna_index.sh to build SortMeRNA indices"
