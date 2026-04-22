#!/bin/bash
#
# download_cms.sh - Download and press Rfam covariance models for rRNA verification
#
# Requires: Infernal (cmpress), wget or curl
# Usage: bash download_cms.sh [output_dir]

set -euo pipefail

OUTPUT_DIR="${1:-data/cms}"

# Rfam CMs needed for SILVA verification
declare -A CMS=(
  [RF00177]="Bacterial SSU (16S)"
  [RF01959]="Archaeal SSU (16S)"
  [RF01960]="Eukaryotic nuclear SSU (18S)"
  [RF02541]="Bacterial LSU (23S)"
  [RF02540]="Archaeal LSU (23S)"
  [RF02543]="Eukaryotic LSU (28S)"
  [RF00001]="5S rRNA"
  [RF00002]="5.8S rRNA"
  [RF02545]="Mitochondrial SSU"
  [RF02546]="Mitochondrial LSU"
)

echo "============================================"
echo "Rfam CM Download Script"
echo "Output directory: ${OUTPUT_DIR}"
echo "============================================"

if ! command -v cmpress &>/dev/null; then
  echo "Error: Infernal (cmpress) not found. Please install Infernal first."
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

for rfam_id in "${!CMS[@]}"; do
  description="${CMS[$rfam_id]}"
  cm_file="${OUTPUT_DIR}/${rfam_id}.cm"

  if [[ -f "${cm_file}" ]]; then
  echo "Already downloaded: ${cm_file} (${description})"
  else
  echo "Downloading ${rfam_id} — ${description}..."
  url="https://rfam.org/family/${rfam_id}/cm"
  if command -v wget &>/dev/null; then
    wget -q -O "${cm_file}" "${url}"
  elif command -v curl &>/dev/null; then
    curl -sL -o "${cm_file}" "${url}"
  else
    echo "Error: Neither wget nor curl found."
    exit 1
  fi
  echo "  Saved: ${cm_file}"
  fi

  if [[ ! -f "${cm_file}.i1m" ]]; then
  echo "  Pressing: ${rfam_id}.cm"
  cmpress "${cm_file}"
  else
  echo "  Already pressed: ${rfam_id}.cm"
  fi
done

echo ""
echo "============================================"
echo "All CMs ready in: ${OUTPUT_DIR}"
echo "============================================"
echo ""
echo "Next step: Run verify_silva.sh to verify downloaded SILVA sequences"
