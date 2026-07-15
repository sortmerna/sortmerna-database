#!/usr/bin/env bash
################################################################################
# download_pacbio_metagenomics.sh
#
# Download PacBio HiFi metagenomic runs from the human gut study
# (BioProject PRJNA1139951, Table S28). Selects the "pb.concat.no_hsap"
# libraries (human reads removed, per-sample concatenated PacBio runs),
# excluding subsampled variants (Library_name containing ".sub_").
#
# SRA accessions are read from the bundled metadata file
# (assets/Table_S28_15947_PB.ONT.ILMN_metadata_SRA_v2.txt) for reproducibility:
#   Library_name = column 5, Run (SRR accession) = column 6.
#
# Usage: bash download_pacbio_metagenomics.sh [output_dir [threads [max_runs]]]
#
# Positional:
#   output_dir   Output directory (default: data/pacbio_metagenomics)
#   threads      Threads for fasterq-dump (default: 4)
#   max_runs     Max number of runs to download (default: 20, to conserve disk;
#                use 0 for all available runs)
#
# Outputs:
#   <output_dir>/<SRR>.fastq.gz   one gzipped FASTQ per run
#
# Requires: sra-tools (prefetch, fasterq-dump), pigz or gzip.
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
METADATA="${SCRIPT_DIR}/assets/Table_S28_15947_PB.ONT.ILMN_metadata_SRA_v2.txt"

OUTPUT_DIR="${1:-data/pacbio_metagenomics}"
THREADS="${2:-4}"
MAX_RUNS="${3:-20}"

# Cap fasterq-dump scratch space (overrides its premature estimate-based abort).
# Set below the free space on the temp volume, leaving headroom for the .sra and
# final .fastq. Override via env, e.g. DISK_LIMIT_TMP=100G.
DISK_LIMIT_TMP="${DISK_LIMIT_TMP:-60G}"

OUTPUT_DIR="$(mkdir -p "${OUTPUT_DIR}" && cd "${OUTPUT_DIR}" && pwd)"

for tool in prefetch fasterq-dump; do
    if ! command -v "${tool}" &> /dev/null; then
        echo "Error: ${tool} not found. Install sra-tools (e.g. conda install -c bioconda sra-tools)"
        exit 1
    fi
done

GZIP=$(command -v pigz || command -v gzip)

# Full-run pb.concat.no_hsap accessions (exclude .sub_ subsamples)
mapfile -t RUNS < <(awk -F'\t' 'NR>1 && $5 ~ /pb\.concat\.no_hsap/ && $5 !~ /\.sub_/ {print $6}' "${METADATA}")

# Cap the number of runs to conserve disk (max_runs=0 means all)
if [[ "${MAX_RUNS}" -gt 0 && "${#RUNS[@]}" -gt "${MAX_RUNS}" ]]; then
    RUNS=("${RUNS[@]:0:${MAX_RUNS}}")
fi

echo "============================================"
echo "PacBio metagenomics SRA download"
echo "Output directory: ${OUTPUT_DIR}"
echo "Runs to download:  ${#RUNS[@]}"
echo "============================================"

for run in "${RUNS[@]}"; do
    fastq_gz="${OUTPUT_DIR}/${run}.fastq.gz"
    if [[ -f "${fastq_gz}" ]]; then
        echo "Already exists: ${run}.fastq.gz"
        continue
    fi
    echo "Downloading ${run} ..."
    prefetch --output-directory "${OUTPUT_DIR}" "${run}"
    tmp_dir="${OUTPUT_DIR}/${run}.tmp"
    mkdir -p "${tmp_dir}"
    fasterq-dump --threads "${THREADS}" --outdir "${OUTPUT_DIR}" \
        --temp "${tmp_dir}" --disk-limit-tmp "${DISK_LIMIT_TMP}" \
        "${OUTPUT_DIR}/${run}/${run}.sra"
    rm -rf "${tmp_dir}"
    "${GZIP}" "${OUTPUT_DIR}/${run}.fastq"
    rm -rf "${OUTPUT_DIR}/${run}"
    echo "  Saved: ${run}.fastq.gz"
done

echo ""
echo "Download complete: ${#RUNS[@]} runs in ${OUTPUT_DIR}"
