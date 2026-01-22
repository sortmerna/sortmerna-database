#!/bin/bash
#
# download_silva.sh - Download SILVA rRNA database
#
# Downloads SILVA Release 138.2 Ref NR 99 SSU and LSU sequences
# Source: https://www.arb-silva.de/

set -euo pipefail

# Configuration
SILVA_VERSION="138.2"
SILVA_BASE_URL="https://www.arb-silva.de/fileadmin/silva_databases/release_${SILVA_VERSION}/Exports"
OUTPUT_DIR="${1:-data/silva}"

# SILVA file names
SSU_FILE="SILVA_${SILVA_VERSION}_SSURef_NR99_tax_silva.fasta.gz"
LSU_FILE="SILVA_${SILVA_VERSION}_LSURef_NR99_tax_silva.fasta.gz"

echo "============================================"
echo "SILVA Database Download Script"
echo "Version: ${SILVA_VERSION}"
echo "Output directory: ${OUTPUT_DIR}"
echo "============================================"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Function to download and verify file
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

# Download SSU (16S/18S) sequences
download_file \
    "${SILVA_BASE_URL}/${SSU_FILE}" \
    "${OUTPUT_DIR}/${SSU_FILE}" \
    "SSU (16S/18S) rRNA sequences"

# Download LSU (23S/28S) sequences
download_file \
    "${SILVA_BASE_URL}/${LSU_FILE}" \
    "${OUTPUT_DIR}/${LSU_FILE}" \
    "LSU (23S/28S) rRNA sequences"

# Decompress files
echo ""
echo "Decompressing files..."

for gz_file in "${OUTPUT_DIR}"/*.fasta.gz; do
    if [[ -f "${gz_file}" ]]; then
        fasta_file="${gz_file%.gz}"
        if [[ ! -f "${fasta_file}" ]]; then
            echo "Decompressing: ${gz_file}"
            gunzip -k "${gz_file}"
        else
            echo "Already decompressed: ${fasta_file}"
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
echo "SSU sequences: ${OUTPUT_DIR}/${SSU_FILE%.gz}"
echo "LSU sequences: ${OUTPUT_DIR}/${LSU_FILE%.gz}"
echo ""

# Count sequences
echo "Sequence counts:"
for fasta in "${OUTPUT_DIR}"/*.fasta; do
    if [[ -f "${fasta}" ]]; then
        count=$(seqkit stats -T "${fasta}" | tail -1 | cut -f4)
        echo "  $(basename "${fasta}"): ${count} sequences"
    fi
done

echo ""
echo "Next step: Run download_rfam.sh to download RFAM sequences"
