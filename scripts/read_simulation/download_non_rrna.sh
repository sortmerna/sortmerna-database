#!/usr/bin/env bash

################################################################################
# download_non_rrna.sh
#
# Download source data for SortMeRNA specificity testing (false positive rate).
# Provides ~1M non-rRNA sequences from two sources:
#
#   1. Human T2T genome (version set via T2T_VERSION / T2T_ACCESSION / T2T_NAME env vars)
#      rRNA loci are masked using the T2T GFF annotations before read simulation.
#      850,000 150bp PE reads are simulated with ART (see simulate_non_rrna.sh).
#
#   2. Rfam non-rRNA families (150,000 sequences)
#      tRNA, SRP RNA, tmRNA, RNase P (bac+euk), U1-U6 spliceosomal
#      RF00005/RF00017/RF00023/RF00010/RF00009/RF00003/RF00004/RF00015/RF00020/RF00026
#
# Rfam sequences share structural features with rRNA and are the most likely
# source of false positives - making them the most challenging specificity test.
#
# Usage: bash download_non_rrna.sh [OPTIONS]
#
# Environment variables (set in README "Set paths" section):
#   T2T_ACCESSION       NCBI accession for CHM13v2.0 (default: GCA_009914755.4)
#   T2T_BASE            NCBI FTP base URL for the T2T assembly
#   RFAM_NON_RRNA_FTP   Rfam FASTA FTP base URL (default: CURRENT release)
#
# Usage: bash download_non_rrna.sh [output_dir [threads]] [OPTIONS]
#
# Positional:
#   output_dir          Output directory (default: $NON_RRNA_DIR or data/non_rrna)
#   threads             Number of threads (default: 4)
#
# Options:
#   --rfam INT          Number of Rfam non-rRNA sequences to sample (default: 150000)
#   --seed INT          Random seed for reproducibility (default: 42)
#   --skip-download     Skip download, use existing files
#   -h, --help          Show help
#
# Outputs:
#   t2t/${T2T_VERSION}.fa.gz           - T2T genome (input for ART simulation)
#   t2t/${T2T_VERSION}.gff.gz          - T2T GFF annotations (rRNA loci for masking)
#   t2t/${T2T_VERSION}_rrna_loci.bed   - rRNA loci in BED format (for bedtools maskfasta)
#   rfam/RF*.fa.gz                     - Rfam non-rRNA family FASTA files
#   rfam_non_rrna_sampled.fasta        - Sampled Rfam sequences ready for use
#
# Next step: Run simulate_non_rrna.sh to mask rRNA loci, simulate T2T reads
#            with ART, and combine with Rfam sequences into non_rRNA_test_1M.fasta
#
################################################################################

set -euo pipefail

POSITIONAL=()
N_Rfam=150000
RAND_SEED=42
SKIP_DOWNLOAD=false

T2T_ACCESSION="${T2T_ACCESSION:-GCA_009914755.4}"
T2T_NAME="${T2T_NAME:-CHM13_T2T_v2.0}"
T2T_VERSION="${T2T_VERSION:-chm13v2.0}"
T2T_BASE="${T2T_BASE:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/${T2T_ACCESSION}_${T2T_NAME}}"
RFAM_NON_RRNA_FTP="${RFAM_NON_RRNA_FTP:-https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files}"

show_help() {
    grep '^#' "$0" | grep -v '#!/usr/bin/env bash' | sed 's/^# \?//'
    exit 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --rfam) N_Rfam="$2"; shift 2 ;;
        --seed) RAND_SEED="$2"; shift 2 ;;
        --skip-download) SKIP_DOWNLOAD=true; shift ;;
        -h|--help) show_help ;;
        *) POSITIONAL+=("$1"); shift ;;
    esac
done

OUTPUT_DIR="${POSITIONAL[0]:-${NON_RRNA_DIR:-data/non_rrna}}"
THREADS="${POSITIONAL[1]:-4}"

# Use absolute path
OUTPUT_DIR="$(mkdir -p "${OUTPUT_DIR}" && cd "${OUTPUT_DIR}" && pwd)"

echo "============================================"
echo "Non-rRNA Source Download Script"
echo "Output directory: ${OUTPUT_DIR}"
echo "T2T accession:    ${T2T_ACCESSION}"
echo "T2T name:         ${T2T_NAME}"
echo "T2T version:      ${T2T_VERSION}"
echo "T2T base URL:     ${T2T_BASE}"
echo "Rfam FTP:         ${RFAM_NON_RRNA_FTP}"
echo "Rfam sample size: ${N_Rfam}"
echo "Random seed:      ${RAND_SEED}"
echo "============================================"
echo ""

# Check for required tools
for tool in wget seqkit; do
    if ! command -v "${tool}" &> /dev/null; then
        echo "Error: ${tool} not found. Install via: conda env create -f environment.yml"
        exit 1
    fi
done

################################################################################
# 1. HUMAN T2T GENOME
################################################################################

T2T_DIR="${OUTPUT_DIR}/t2t"

echo "============================================"
echo "Downloading T2T genome (${T2T_VERSION})"
echo "============================================"

mkdir -p "${T2T_DIR}"

T2T_GENOME_GZ="${T2T_DIR}/${T2T_VERSION}.fa.gz"
T2T_GFF_GZ="${T2T_DIR}/${T2T_VERSION}.gff.gz"

if [[ "${SKIP_DOWNLOAD}" == false ]]; then
    if [[ ! -f "${T2T_GENOME_GZ}" ]]; then
        echo "Downloading ${T2T_VERSION} genome (~889 MB compressed)..."
        wget -c --progress=bar "${T2T_BASE}/${T2T_ACCESSION}_${T2T_NAME}_genomic.fna.gz" \
            -O "${T2T_GENOME_GZ}"
        echo "  Saved: ${T2T_VERSION}.fa.gz"
    else
        echo "Already exists: ${T2T_VERSION}.fa.gz"
    fi

    if [[ ! -f "${T2T_GFF_GZ}" ]]; then
        echo "Downloading ${T2T_VERSION} GFF annotations (for rRNA locus masking)..."
        wget -c "${T2T_BASE}/${T2T_ACCESSION}_${T2T_NAME}_genomic.gff.gz" \
            -O "${T2T_GFF_GZ}"
        echo "  Saved: ${T2T_VERSION}.gff.gz"
    else
        echo "Already exists: ${T2T_VERSION}.gff.gz"
    fi
fi

echo "Extracting rRNA loci from GFF to BED (for masking in simulation step)..."
T2T_RRNA_BED="${T2T_DIR}/${T2T_VERSION}_rrna_loci.bed"
if [[ ! -f "${T2T_RRNA_BED}" ]]; then
    zcat "${T2T_GFF_GZ}" \
        | awk '$3 == "rRNA" {print $1"\t"$4-1"\t"$5}' \
        | sort -k1,1 -k2,2n \
        > "${T2T_RRNA_BED}"
    echo "  rRNA loci: $(wc -l < "${T2T_RRNA_BED}") regions -> ${T2T_RRNA_BED}"
else
    echo "Already exists: ${T2T_VERSION}_rrna_loci.bed ($(wc -l < "${T2T_RRNA_BED}") regions)"
fi

################################################################################
# 2. Rfam NON-rRNA FAMILIES
################################################################################

RFAM_DIR="${OUTPUT_DIR}/rfam"

# Non-rRNA families to download (family_id -> description)
declare -A RFAM_NON_RRNA=(
    ["RF00005"]="tRNA"
    ["RF00017"]="SRP_RNA"
    ["RF00023"]="tmRNA"
    ["RF00010"]="RNaseP_bacterial"
    ["RF00009"]="RNaseP_eukaryotic"
    ["RF00003"]="U1_spliceosomal"
    ["RF00004"]="U2_spliceosomal"
    ["RF00015"]="U4_spliceosomal"
    ["RF00020"]="U5_spliceosomal"
    ["RF00026"]="U6_spliceosomal"
)

if [[ "${SKIP_DOWNLOAD}" == false ]]; then
    echo ""
    echo "============================================"
    echo "Downloading Rfam non-rRNA families"
    echo "============================================"

    mkdir -p "${RFAM_DIR}"

    for family_id in "${!RFAM_NON_RRNA[@]}"; do
        family_name="${RFAM_NON_RRNA[$family_id]}"
        outfile="${RFAM_DIR}/${family_id}_${family_name}.fa.gz"
        if [[ -f "${outfile}" ]]; then
            echo "Already exists: $(basename "${outfile}")"
            continue
        fi
        echo "Downloading ${family_id} (${family_name})..."
        wget -c "${RFAM_NON_RRNA_FTP}/${family_id}.fa.gz" -O "${outfile}" || {
            echo "Warning: ${family_id} download failed, continuing..."
            rm -f "${outfile}"
        }
    done

    echo "Decompressing Rfam files..."
    for gz in "${RFAM_DIR}"/*.gz; do
        [[ -f "${gz}" ]] && gunzip -kf "${gz}"
    done
fi

echo ""
echo "Sampling Rfam non-rRNA sequences..."

cat "${RFAM_DIR}"/*.fa > "${OUTPUT_DIR}/rfam_non_rrna_all.fasta" 2>/dev/null || true

seqkit sample --rand-seed "${RAND_SEED}" -n "${N_Rfam}" \
    "${OUTPUT_DIR}/rfam_non_rrna_all.fasta" > "${OUTPUT_DIR}/rfam_non_rrna_sampled.fasta" 2>/dev/null

echo "Rfam non-rRNA: $(seqkit stats -T "${OUTPUT_DIR}/rfam_non_rrna_sampled.fasta" | tail -1 | cut -f4) sequences"

echo ""
echo "============================================"
echo "Download complete"
echo "============================================"
echo ""
echo "Outputs:"
echo "  T2T genome:   ${T2T_DIR}/chm13v2.0.fa.gz"
echo "  T2T GFF:      ${T2T_DIR}/chm13v2.0.gff.gz"
echo "  rRNA BED:     ${T2T_DIR}/chm13v2.0_rrna_loci.bed"
echo "  Rfam sampled: ${OUTPUT_DIR}/rfam_non_rrna_sampled.fasta"
echo ""
echo "Next step: Run simulate_non_rrna.sh to mask rRNA loci, simulate T2T reads"
echo "           with ART, and combine with Rfam sequences into non_rRNA_test_1M.fasta"
echo "  3. Target: <0.1% false positive rate"
