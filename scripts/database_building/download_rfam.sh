#!/bin/bash
#
# download_rfam.sh - Download RFAM rRNA families
#
# Downloads both full and seed FASTA sequences for 5S (RF00001) and 5.8S (RF00002) rRNA
# Source: https://rfam.org/

set -euo pipefail

# Configuration
RFAM_VERSION="15.1"
RFAM_FTP_URL="https://ftp.ebi.ac.uk/pub/databases/Rfam/${RFAM_VERSION}"
OUTPUT_DIR="${1:-data/rfam}"

# RFAM families to download
declare -A RFAM_FAMILIES=(
    ["RF00001"]="5S_rRNA"
    ["RF00002"]="5_8S_rRNA"
)

echo "============================================"
echo "RFAM Database Download Script"
echo "Version: ${RFAM_VERSION}"
echo "Output directory: ${OUTPUT_DIR}"
echo "============================================"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Function to download file
download_file() {
    local url="$1"
    local output="$2"
    local description="$3"

    echo ""
    echo "Downloading ${description}..."
    echo "URL: ${url}"

    if [[ -f "${output}" ]]; then
        echo "File already exists: ${output}"
        echo "Skipping download. Delete file to re-download."
        return 0
    fi

    if command -v wget &> /dev/null; then
        wget --progress=bar:force -O "${output}" "${url}"
    elif command -v curl &> /dev/null; then
        curl -L --progress-bar -o "${output}" "${url}"
    else
        echo "Error: Neither wget nor curl found. Please install one of them."
        exit 1
    fi

    echo "Downloaded: ${output}"
}

# Download full and seed sequences for each family
for family_id in "${!RFAM_FAMILIES[@]}"; do
    family_name="${RFAM_FAMILIES[$family_id]}"

    # Download full sequences (all members of the family) from FTP
    full_url="${RFAM_FTP_URL}/fasta_files/${family_id}.fa.gz"
    full_file="${OUTPUT_DIR}/${family_id}_${family_name}_full.fa.gz"
    download_file "${full_url}" "${full_file}" "${family_name} (${family_id}) full sequences"

    # Download seed sequences (curated representatives) from RFAM API
    seed_url="https://rfam.org/family/${family_id}/alignment/fastau"
    seed_file="${OUTPUT_DIR}/${family_id}_${family_name}_seed.fa"
    download_file "${seed_url}" "${seed_file}" "${family_name} (${family_id}) seed sequences"
done

# Decompress full sequence files
echo ""
echo "Decompressing files..."

for gz_file in "${OUTPUT_DIR}"/*_full.fa.gz; do
    if [[ -f "${gz_file}" ]]; then
        fa_file="${gz_file%.gz}"
        if [[ ! -f "${fa_file}" ]]; then
            echo "Decompressing: ${gz_file}"
            gunzip -k "${gz_file}"
        else
            echo "Already decompressed: ${fa_file}"
        fi
    fi
done

# Print summary
echo ""
echo "============================================"
echo "Download complete!"
echo "============================================"
echo "Files downloaded to: ${OUTPUT_DIR}"
echo ""

# Count sequences
echo "Sequence counts:"
echo ""
echo "Full sequences:"
for fa in "${OUTPUT_DIR}"/*_full.fa; do
    if [[ -f "${fa}" ]]; then
        count=$(seqkit stats -T "${fa}" | tail -1 | cut -f4)
        echo "  $(basename "${fa}"): ${count} sequences"
    fi
done

echo ""
echo "Seed sequences:"
for fa in "${OUTPUT_DIR}"/*_seed.fa; do
    if [[ -f "${fa}" ]]; then
        count=$(seqkit stats -T "${fa}" | tail -1 | cut -f4)
        echo "  $(basename "${fa}"): ${count} sequences"
    fi
done

echo ""
echo "Next step: Run cluster_sequences.sh to cluster the sequences"
