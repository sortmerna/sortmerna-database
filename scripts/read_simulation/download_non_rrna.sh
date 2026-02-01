#!/usr/bin/env bash

################################################################################
# download_non_rrna.sh
#
# Download and prepare non-rRNA sequences for SortMeRNA specificity testing
# (false positive rate). Creates a 1M sequence test set from diverse sources.
#
# Sources:
#   1. Bacterial mRNA   - RefSeq bacteria RNA, filtered for mRNA only
#   2. Eukaryotic mRNA  - Ensembl cDNA (human, mouse, zebrafish, C. elegans, Arabidopsis)
#   3. Rfam non-rRNA    - Individual non-rRNA Rfam families (tRNA, SRP, tmRNA, etc.)
#   4. Genomic fragments - Random fragments from bacterial/eukaryotic genomes
#
# Safety filters:
#   - Header keyword filter (seqkit grep) to remove labeled rRNA
#   - Barrnap (bac + euk) to detect rRNA genes in genomic fragments
#
# Usage: bash download_non_rrna.sh [OPTIONS]
#
# Options:
#   -o, --output DIR        Output directory (default: data/non_rrna)
#   --bacterial INT         Number of bacterial mRNA sequences (default: 500000)
#   --eukaryotic INT        Number of eukaryotic mRNA sequences (default: 300000)
#   --rfam INT              Number of Rfam non-rRNA sequences (default: 150000)
#   --genomic INT           Number of random genomic sequences (default: 50000)
#   --seed INT              Random seed for reproducibility (default: 42)
#   --skip-download         Skip download, use existing files
#   --threads INT           Number of threads (default: 8)
#   -h, --help              Show help
#
# Output:
#   non_rRNA_test_1M.fasta  - Combined 1M non-rRNA test sequences
#   non_rRNA_metadata.txt   - Composition and source information
#   non_rRNA_stats.txt      - Sequence statistics
#
################################################################################

set -euo pipefail

# Defaults
OUTPUT_DIR="data/non_rrna"
N_BACTERIAL=500000
N_EUKARYOTIC=300000
N_RFAM=150000
N_GENOMIC=50000
RAND_SEED=42
SKIP_DOWNLOAD=false
THREADS=8

show_help() {
    grep '^#' "$0" | grep -v '#!/usr/bin/env bash' | sed 's/^# \?//'
    exit 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
        --bacterial) N_BACTERIAL="$2"; shift 2 ;;
        --eukaryotic) N_EUKARYOTIC="$2"; shift 2 ;;
        --rfam) N_RFAM="$2"; shift 2 ;;
        --genomic) N_GENOMIC="$2"; shift 2 ;;
        --seed) RAND_SEED="$2"; shift 2 ;;
        --skip-download) SKIP_DOWNLOAD=true; shift ;;
        --threads) THREADS="$2"; shift 2 ;;
        -h|--help) show_help ;;
        *) echo "Error: Unknown option: $1"; show_help ;;
    esac
done

# Use absolute path
OUTPUT_DIR="$(mkdir -p "${OUTPUT_DIR}" && cd "${OUTPUT_DIR}" && pwd)"

TOTAL_SEQS=$((N_BACTERIAL + N_EUKARYOTIC + N_RFAM + N_GENOMIC))

echo "============================================"
echo "Non-rRNA Sequence Download Script"
echo "Output directory: ${OUTPUT_DIR}"
echo "Random seed: ${RAND_SEED}"
echo "============================================"
echo ""
echo "Target composition:"
echo "  Bacterial mRNA:    ${N_BACTERIAL}"
echo "  Eukaryotic mRNA:   ${N_EUKARYOTIC}"
echo "  Rfam non-rRNA:     ${N_RFAM}"
echo "  Genomic DNA:       ${N_GENOMIC}"
echo "  --------------------------------"
echo "  Total:             ${TOTAL_SEQS}"
echo ""

# Check for required tools
for tool in wget seqkit barrnap python3; do
    if ! command -v "${tool}" &> /dev/null; then
        echo "Error: ${tool} not found. Install via: conda env create -f environment.yml"
        exit 1
    fi
done

if ! python3 -c "from Bio import SeqIO" 2>/dev/null; then
    echo "Error: BioPython required: pip install biopython"
    exit 1
fi

################################################################################
# 1. BACTERIAL mRNA (RefSeq)
################################################################################

REFSEQ_DIR="${OUTPUT_DIR}/refseq_bacteria"

if [[ "${SKIP_DOWNLOAD}" == false ]]; then
    echo ""
    echo "============================================"
    echo "Downloading bacterial mRNA from RefSeq"
    echo "============================================"
    echo ""
    echo "NOTE: Each file is ~100-500MB. Downloading 10 files."

    mkdir -p "${REFSEQ_DIR}"

    REFSEQ_BASE="https://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria"

    for i in {1..10}; do
        file_num=$(printf "%03d" $i)
        outfile="${REFSEQ_DIR}/bacteria.${file_num}.rna.fna.gz"
        if [[ -f "${outfile}" ]]; then
            echo "Already exists: $(basename "${outfile}")"
            continue
        fi
        echo "Downloading bacteria.${file_num}.rna.fna.gz..."
        wget -q -c "${REFSEQ_BASE}/bacteria.${file_num}.rna.fna.gz" -O "${outfile}" || {
            echo "Warning: bacteria.${file_num}.rna.fna.gz not found, continuing..."
            rm -f "${outfile}"
            continue
        }
    done

    echo "Decompressing bacterial RNA files..."
    for gz in "${REFSEQ_DIR}"/*.gz; do
        [[ -f "${gz}" ]] && gunzip -kf "${gz}"
    done
fi

echo ""
echo "Processing bacterial mRNA sequences..."

python3 << PYTHON_BACTERIAL
from Bio import SeqIO
import random

random.seed(${RAND_SEED})

target_count = ${N_BACTERIAL}

# rRNA/tRNA keywords to exclude
exclude_terms = [
    'ribosomal rna', 'rrna', 'trna', 'transfer rna',
    '16s', '23s', '5s', '18s', '28s', '5.8s',
    'ribosomal', 'transfer'
]

import glob
input_files = sorted(glob.glob("${REFSEQ_DIR}/*.fna"))

mRNA_seqs = []
excluded = 0

print(f"Scanning {len(input_files)} RefSeq files...")
for fpath in input_files:
    for record in SeqIO.parse(fpath, 'fasta'):
        description = record.description.lower()

        is_excluded = any(term in description for term in exclude_terms)
        if is_excluded:
            excluded += 1
            continue

        mRNA_seqs.append(record)

        if len(mRNA_seqs) % 50000 == 0:
            print(f"  Collected {len(mRNA_seqs)} mRNA sequences...")

        if len(mRNA_seqs) >= target_count * 2:
            break
    if len(mRNA_seqs) >= target_count * 2:
        break

print(f"Collected {len(mRNA_seqs)} mRNA sequences (excluded {excluded} rRNA/tRNA)")

random.shuffle(mRNA_seqs)
sampled = mRNA_seqs[:target_count]

output_file = "${OUTPUT_DIR}/bacterial_mRNA_filtered.fasta"
with open(output_file, 'w') as out:
    SeqIO.write(sampled, out, 'fasta')

print(f"Saved {len(sampled)} sequences to {output_file}")
PYTHON_BACTERIAL

echo "Bacterial mRNA: $(seqkit stats -T "${OUTPUT_DIR}/bacterial_mRNA_filtered.fasta" | tail -1 | cut -f4) sequences"

################################################################################
# 2. EUKARYOTIC mRNA (Ensembl)
################################################################################

ENSEMBL_DIR="${OUTPUT_DIR}/ensembl"
ENSEMBL_RELEASE="112"

if [[ "${SKIP_DOWNLOAD}" == false ]]; then
    echo ""
    echo "============================================"
    echo "Downloading eukaryotic cDNA from Ensembl"
    echo "============================================"

    mkdir -p "${ENSEMBL_DIR}"

    ENSEMBL_BASE="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta"

    declare -A ENSEMBL_FILES=(
        ["human"]="${ENSEMBL_BASE}/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
        ["mouse"]="${ENSEMBL_BASE}/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
        ["zebrafish"]="${ENSEMBL_BASE}/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz"
        ["worm"]="${ENSEMBL_BASE}/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz"
    )

    # Arabidopsis is on Ensembl Plants (separate FTP)
    ENSEMBL_FILES["arabidopsis"]="https://ftp.ensemblgenomes.org/pub/plants/release-59/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"

    for organism in "${!ENSEMBL_FILES[@]}"; do
        url="${ENSEMBL_FILES[$organism]}"
        outfile="${ENSEMBL_DIR}/${organism}_cdna.fa.gz"
        if [[ -f "${outfile}" ]]; then
            echo "Already exists: $(basename "${outfile}")"
            continue
        fi
        echo "Downloading ${organism} cDNA..."
        wget -q -c "${url}" -O "${outfile}" || {
            echo "Warning: ${organism} download failed, continuing..."
            rm -f "${outfile}"
        }
    done

    echo "Decompressing Ensembl files..."
    for gz in "${ENSEMBL_DIR}"/*.gz; do
        [[ -f "${gz}" ]] && gunzip -kf "${gz}"
    done
fi

echo ""
echo "Sampling eukaryotic mRNA sequences..."

seqkit sample --rand-seed "${RAND_SEED}" -n "${N_EUKARYOTIC}" \
    "${ENSEMBL_DIR}"/*.fa > "${OUTPUT_DIR}/eukaryotic_mRNA_sampled.fasta" 2>/dev/null

echo "Eukaryotic mRNA: $(seqkit stats -T "${OUTPUT_DIR}/eukaryotic_mRNA_sampled.fasta" | tail -1 | cut -f4) sequences"

################################################################################
# 3. RFAM NON-rRNA FAMILIES
################################################################################

RFAM_DIR="${OUTPUT_DIR}/rfam"
RFAM_VERSION="15.1"
RFAM_FTP="https://ftp.ebi.ac.uk/pub/databases/Rfam/${RFAM_VERSION}/fasta_files"

# Non-rRNA families to download (family_id -> description)
declare -A RFAM_NON_RRNA=(
    ["RF00005"]="tRNA"
    ["RF00017"]="SRP_RNA"
    ["RF00023"]="tmRNA"
    ["RF00010"]="RNaseP"
    ["RF00373"]="RNaseMRP"
    ["RF00004"]="U2_spliceosomal_RNA"
    ["RF00026"]="U6_spliceosomal_RNA"
    ["RF00030"]="RNaseP_bact_a"
    ["RF01854"]="bacterial_RNaseP_class_B"
    ["RF00169"]="bacterial_SRP"
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
        wget -q -c "${RFAM_FTP}/${family_id}.fa.gz" -O "${outfile}" || {
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

seqkit sample --rand-seed "${RAND_SEED}" -n "${N_RFAM}" \
    "${OUTPUT_DIR}/rfam_non_rrna_all.fasta" > "${OUTPUT_DIR}/rfam_non_rrna_sampled.fasta" 2>/dev/null

echo "Rfam non-rRNA: $(seqkit stats -T "${OUTPUT_DIR}/rfam_non_rrna_sampled.fasta" | tail -1 | cut -f4) sequences"

################################################################################
# 4. RANDOM GENOMIC DNA
################################################################################

GENOME_DIR="${OUTPUT_DIR}/genomes"

if [[ "${SKIP_DOWNLOAD}" == false ]]; then
    echo ""
    echo "============================================"
    echo "Downloading genomes for random fragments"
    echo "============================================"

    mkdir -p "${GENOME_DIR}"

    declare -A GENOMES=(
        ["ecoli"]="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
        ["bsubtilis"]="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Bacillus_subtilis/reference/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz"
        ["scerevisiae"]="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/reference/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz"
    )

    for name in "${!GENOMES[@]}"; do
        url="${GENOMES[$name]}"
        outfile="${GENOME_DIR}/${name}.fna.gz"
        if [[ -f "${outfile}" ]]; then
            echo "Already exists: $(basename "${outfile}")"
            continue
        fi
        echo "Downloading ${name} genome..."
        wget -q -c "${url}" -O "${outfile}" || {
            echo "Warning: ${name} download failed, continuing..."
            rm -f "${outfile}"
        }
    done

    echo "Decompressing genome files..."
    for gz in "${GENOME_DIR}"/*.gz; do
        [[ -f "${gz}" ]] && gunzip -kf "${gz}"
    done
fi

echo ""
echo "Generating random genomic DNA fragments..."

python3 << PYTHON_GENOMIC
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob
import random

random.seed(${RAND_SEED})

target_count = ${N_GENOMIC}
genome_dir = "${GENOME_DIR}"

# Load genomes
genomes = []
for gfile in sorted(glob.glob(f"{genome_dir}/*.fna")):
    for record in SeqIO.parse(gfile, 'fasta'):
        seq_str = str(record.seq).upper()
        if len(seq_str) > 1000:
            genomes.append((record.id, seq_str))
            print(f"Loaded: {record.id} ({len(seq_str):,} bp)")

if not genomes:
    print("Error: No genomes loaded")
    exit(1)

# Extract random fragments (500-3000 bp)
fragments = []
for i in range(target_count):
    genome_id, genome_seq = random.choice(genomes)
    length = random.randint(500, 3000)

    if len(genome_seq) < length:
        continue

    start = random.randint(0, len(genome_seq) - length)
    fragment = genome_seq[start:start + length]

    # Skip fragments with too many Ns
    if fragment.count('N') / len(fragment) > 0.1:
        continue

    record = SeqRecord(
        Seq(fragment),
        id=f"genomic_{genome_id}_{start}_{start+length}",
        description=f"random fragment from {genome_id}"
    )
    fragments.append(record)

    if (i + 1) % 10000 == 0:
        print(f"  Generated {i + 1} fragments...")

print(f"Generated {len(fragments)} genomic fragments")

output_file = "${OUTPUT_DIR}/genomic_random.fasta"
with open(output_file, 'w') as out:
    SeqIO.write(fragments, out, 'fasta')

print(f"Saved to {output_file}")
PYTHON_GENOMIC

echo "Genomic fragments: $(seqkit stats -T "${OUTPUT_DIR}/genomic_random.fasta" | tail -1 | cut -f4) sequences"

################################################################################
# 5. COMBINE AND LABEL
################################################################################

echo ""
echo "============================================"
echo "Combining and labeling sequences"
echo "============================================"

python3 << PYTHON_LABEL
from Bio import SeqIO

sources = {
    "${OUTPUT_DIR}/bacterial_mRNA_filtered.fasta": "bacterial_mRNA",
    "${OUTPUT_DIR}/eukaryotic_mRNA_sampled.fasta": "eukaryotic_mRNA",
    "${OUTPUT_DIR}/rfam_non_rrna_sampled.fasta": "rfam_ncRNA",
    "${OUTPUT_DIR}/genomic_random.fasta": "genomic_DNA",
}

all_seqs = []
counts = {}

for source_file, label in sources.items():
    count = 0
    for record in SeqIO.parse(source_file, 'fasta'):
        count += 1
        record.id = f"{label}_{count}"
        record.description = f"{label} | {record.description}"
        all_seqs.append(record)
    counts[label] = count
    print(f"  {label}: {count}")

output_file = "${OUTPUT_DIR}/non_rRNA_test_combined.fasta"
with open(output_file, 'w') as out:
    SeqIO.write(all_seqs, out, 'fasta')

print(f"\nCombined total: {len(all_seqs)}")
PYTHON_LABEL

################################################################################
# 6. SAFETY FILTERS - remove any rRNA contamination
################################################################################

echo ""
echo "============================================"
echo "Running rRNA safety filters"
echo "============================================"

COMBINED="${OUTPUT_DIR}/non_rRNA_test_combined.fasta"
FILTERED="${OUTPUT_DIR}/non_rRNA_test_1M.fasta"

# --- Filter 1: Header keyword filter ---
echo ""
echo "Filter 1: Header keyword scan..."

RRNA_KEYWORDS="ribosomal|rRNA|rrna|rDNA| 16S | 23S | 5S | 18S | 28S | 5\.8S "

seqkit grep -v -r -n -p "${RRNA_KEYWORDS}" "${COMBINED}" \
    > "${OUTPUT_DIR}/non_rRNA_keyword_filtered.fasta" 2>/dev/null

BEFORE=$(seqkit stats -T "${COMBINED}" | tail -1 | cut -f4)
AFTER_KW=$(seqkit stats -T "${OUTPUT_DIR}/non_rRNA_keyword_filtered.fasta" | tail -1 | cut -f4)
KW_REMOVED=$((BEFORE - AFTER_KW))
echo "  Keyword filter removed: ${KW_REMOVED} sequences"

# --- Filter 2: Barrnap (rRNA gene prediction) ---
echo ""
echo "Filter 2: Barrnap rRNA gene prediction..."

# Run barrnap for bacteria and eukaryotes
BARRNAP_HITS="${OUTPUT_DIR}/barrnap_hits.txt"

echo "  Running barrnap --kingdom bac..."
barrnap --kingdom bac --threads "${THREADS}" --quiet \
    "${OUTPUT_DIR}/non_rRNA_keyword_filtered.fasta" 2>/dev/null \
    | grep -v "^#" | cut -f1 | sort -u > "${OUTPUT_DIR}/barrnap_bac_ids.txt" || true

echo "  Running barrnap --kingdom euk..."
barrnap --kingdom euk --threads "${THREADS}" --quiet \
    "${OUTPUT_DIR}/non_rRNA_keyword_filtered.fasta" 2>/dev/null \
    | grep -v "^#" | cut -f1 | sort -u > "${OUTPUT_DIR}/barrnap_euk_ids.txt" || true

# Combine barrnap hits
cat "${OUTPUT_DIR}/barrnap_bac_ids.txt" "${OUTPUT_DIR}/barrnap_euk_ids.txt" \
    | sort -u > "${BARRNAP_HITS}"

BARRNAP_COUNT=$(wc -l < "${BARRNAP_HITS}" | tr -d ' ')
echo "  Barrnap detected rRNA in: ${BARRNAP_COUNT} sequences"

# Remove barrnap hits
if [[ -s "${BARRNAP_HITS}" ]]; then
    seqkit grep -v -f "${BARRNAP_HITS}" \
        "${OUTPUT_DIR}/non_rRNA_keyword_filtered.fasta" > "${FILTERED}" 2>/dev/null
else
    cp "${OUTPUT_DIR}/non_rRNA_keyword_filtered.fasta" "${FILTERED}"
fi

FINAL_COUNT=$(seqkit stats -T "${FILTERED}" | tail -1 | cut -f4)
TOTAL_REMOVED=$((BEFORE - FINAL_COUNT))

echo ""
echo "Safety filter summary:"
echo "  Before filtering:     ${BEFORE}"
echo "  Keyword filter:       -${KW_REMOVED}"
echo "  Barrnap filter:       -${BARRNAP_COUNT}"
echo "  Final clean set:      ${FINAL_COUNT}"

# Clean up intermediate files
rm -f "${OUTPUT_DIR}/non_rRNA_keyword_filtered.fasta"
rm -f "${OUTPUT_DIR}/barrnap_bac_ids.txt" "${OUTPUT_DIR}/barrnap_euk_ids.txt"

################################################################################
# 7. GENERATE METADATA AND STATISTICS
################################################################################

echo ""
echo "============================================"
echo "Generating metadata and statistics"
echo "============================================"

# Metadata
python3 << PYTHON_META
from Bio import SeqIO
from collections import Counter

fasta_file = "${FILTERED}"
meta_file = "${OUTPUT_DIR}/non_rRNA_metadata.txt"

label_counts = Counter()
total = 0

for record in SeqIO.parse(fasta_file, 'fasta'):
    total += 1
    # Label is the prefix before the first underscore+digit
    label = record.id.rsplit('_', 1)[0] if '_' in record.id else record.id
    label_counts[label] += 1

with open(meta_file, 'w') as f:
    f.write("Non-rRNA Test Set Metadata\n")
    f.write("==========================\n\n")
    f.write(f"Total sequences: {total}\n")
    f.write(f"Random seed: ${RAND_SEED}\n\n")
    f.write("Composition:\n")
    for label, count in sorted(label_counts.items()):
        pct = (count / total) * 100
        f.write(f"  {label}: {count} ({pct:.1f}%)\n")

print(f"Metadata written to: {meta_file}")
PYTHON_META

# Statistics
seqkit stats "${FILTERED}" > "${OUTPUT_DIR}/non_rRNA_stats.txt"

echo ""
echo "============================================"
echo "Non-rRNA test set complete!"
echo "============================================"
echo ""
echo "Final sequences: ${FINAL_COUNT}"
echo "Output: ${FILTERED}"
echo "Metadata: ${OUTPUT_DIR}/non_rRNA_metadata.txt"
echo "Stats: ${OUTPUT_DIR}/non_rRNA_stats.txt"
echo ""
seqkit stats "${FILTERED}"
echo ""
echo "Next steps:"
echo "  1. Simulate reads from these sequences"
echo "  2. Run SortMeRNA to measure specificity (false positive rate)"
echo "  3. Target: <0.1% false positive rate"
