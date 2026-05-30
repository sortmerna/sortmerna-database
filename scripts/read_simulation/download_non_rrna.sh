#!/usr/bin/env bash

################################################################################
# download_non_rrna.sh
#
# Download source data for SortMeRNA specificity testing (false positive rate).
# Provides ~1M non-rRNA sequences from two sources:
#
#   1. Human T2T genome (version set via T2T_VERSION / T2T_ACCESSION / T2T_NAME env vars)
#      rRNA loci are masked using the RefSeq GFF3 annotation before read simulation.
#      850,000 150bp PE reads are simulated with ART (see simulate_non_rrna.sh).
#
#   2. Rfam non-rRNA families (150,000 sequences)
#      tRNA, SRP RNA, tmRNA, RNase P (bac+euk), U1-U6 spliceosomal
#      RF00005/RF00017/RF00023/RF00010/RF00009/RF00003/RF00004/RF00015/RF00020/RF00026
#
# Rfam sequences share structural features with rRNA and are the most likely
# source of false positives - making them the most challenging specificity test.
#
# Usage: bash download_non_rrna.sh [output_dir [threads]] [OPTIONS]
#
# Positional:
#   output_dir          Output directory (default: $NON_RRNA_DIR or data/non_rrna)
#   threads             Number of threads (default: 4)
#
# Options:
#   --margin INT        Bp to extend each rRNA locus on each side (default: 100)
#   -h, --help          Show help
#
# Environment variables (set in README "Set paths" section):
#   T2T_ACCESSION       GCA accession for CHM13v2.0 genome FASTA (default: GCA_009914755.4)
#   T2T_GCF_ACCESSION   GCF accession for RefSeq rRNA annotation (default: GCF_009914755.1)
#   T2T_NAME            Assembly name (default: T2T-CHM13v2.0)
#   T2T_BASE            NCBI FTP base URL for the GCA assembly (genome FASTA)
#   T2T_GCF_BASE        NCBI FTP base URL for the GCF assembly (rRNA annotation)
#   RFAM_NON_RRNA_FTP   Rfam FASTA FTP base URL (default: CURRENT release)
#   CMS_DIR             Directory containing pressed Rfam CMs (from download_cms.sh)
#                       If set, cmsearch is run to supplement GFF3 rRNA annotation.
#
# Note: genome FASTA is from GCA (GenBank accession names, e.g. CP068277.2);
#       rRNA annotation GFF3 is from GCF (RefSeq, NC_ accession names).
#       The assembly report maps NC_ names to GenBank accession names to match the FASTA.
#       cmsearch uses RF01960 (18S), RF02543 (28S), RF00001 (5S), RF00002 (5.8S);
#       RF00177 (mt-12S) and RF02541 (mt-16S) are searched against chrM only.
#
# Outputs:
#   t2t/${T2T_VERSION}.fa.gz                    - T2T genome (input for ISS simulation)
#   t2t/${T2T_VERSION}_annotation.gff.gz        - RefSeq GFF3 annotation (~76 MB)
#   t2t/${T2T_VERSION}_assembly_report.txt      - Chromosome name mapping (NC_ -> GenBank)
#   t2t/${T2T_VERSION}_gff3_loci.bed            - rRNA loci from GFF3 annotation
#   t2t/RF*_hits.tbl                            - cmsearch tblout files (if CMS_DIR set)
#   t2t/${T2T_VERSION}_cmsearch_loci.bed        - rRNA loci from cmsearch (if CMS_DIR set)
#   t2t/RF*_extra_vs_gff3.bed                  - per-family regions found by cmsearch but not GFF3
#   t2t/RF*_extra_vs_gff3.fasta                - sequences for those regions (for BLAST verification)
#   t2t/${T2T_VERSION}_rrna_loci.bed            - final merged BED (for bedtools maskfasta)
#   rfam/RF*.fa.gz                              - Rfam non-rRNA family FASTA files
#   rfam_non_rrna_all.fasta                     - All Rfam non-rRNA sequences combined
#
# Next step: Run simulate_non_rrna.sh to mask rRNA loci, simulate T2T reads
#            with InSilicoSeq, and filter Rfam sequences by length
#
################################################################################

set -euo pipefail

POSITIONAL=()
RNA_LOCI_MARGIN=100

T2T_ACCESSION="${T2T_ACCESSION:-GCA_009914755.4}"
T2T_GCF_ACCESSION="${T2T_GCF_ACCESSION:-GCF_009914755.1}"
T2T_NAME="${T2T_NAME:-T2T-CHM13v2.0}"
T2T_VERSION="${T2T_VERSION:-chm13v2.0}"
T2T_BASE="${T2T_BASE:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/${T2T_ACCESSION}_${T2T_NAME}}"
T2T_GCF_BASE="${T2T_GCF_BASE:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/${T2T_GCF_ACCESSION}_${T2T_NAME}}"
RFAM_NON_RRNA_FTP="${RFAM_NON_RRNA_FTP:-https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files}"
CMS_DIR="${CMS_DIR:-}"

show_help() {
    grep '^#' "$0" | grep -v '#!/usr/bin/env bash' | sed 's/^# \?//'
    exit 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --margin) RNA_LOCI_MARGIN="$2"; shift 2 ;;
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
echo "T2T GCA accession: ${T2T_ACCESSION}"
echo "T2T GCF accession: ${T2T_GCF_ACCESSION}"
echo "T2T name:          ${T2T_NAME}"
echo "T2T version:       ${T2T_VERSION}"
echo "Rfam FTP:          ${RFAM_NON_RRNA_FTP}"
echo "rRNA loci margin:  ${RNA_LOCI_MARGIN} bp"
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
# 1. HUMAN T2T GENOME + ANNOTATION
################################################################################

T2T_DIR="${OUTPUT_DIR}/t2t"

echo "============================================"
echo "Downloading T2T genome and annotation (${T2T_VERSION})"
echo "============================================"

mkdir -p "${T2T_DIR}"

T2T_GENOME_GZ="${T2T_DIR}/${T2T_VERSION}.fa.gz"
T2T_GFF_GZ="${T2T_DIR}/${T2T_VERSION}_annotation.gff.gz"
T2T_ASSEMBLY_REPORT="${T2T_DIR}/${T2T_VERSION}_assembly_report.txt"

if [[ ! -f "${T2T_GENOME_GZ}" ]]; then
    echo "Downloading ${T2T_VERSION} genome (~889 MB compressed)..."
    wget -c --progress=bar "${T2T_BASE}/${T2T_ACCESSION}_${T2T_NAME}_genomic.fna.gz" \
        -O "${T2T_GENOME_GZ}"
    echo "  Saved: ${T2T_VERSION}.fa.gz"
else
    echo "Already exists: ${T2T_VERSION}.fa.gz"
fi

if [[ ! -f "${T2T_GFF_GZ}" ]]; then
    echo "Downloading ${T2T_VERSION} RefSeq GFF3 annotation (~76 MB, for rRNA locus masking)..."
    wget -c --progress=bar "${T2T_GCF_BASE}/${T2T_GCF_ACCESSION}_${T2T_NAME}_genomic.gff.gz" \
        -O "${T2T_GFF_GZ}"
    echo "  Saved: ${T2T_VERSION}_annotation.gff.gz"
else
    echo "Already exists: ${T2T_VERSION}_annotation.gff.gz"
fi

if [[ ! -f "${T2T_ASSEMBLY_REPORT}" ]]; then
    echo "Downloading assembly report (chromosome name mapping)..."
    wget -c "${T2T_GCF_BASE}/${T2T_GCF_ACCESSION}_${T2T_NAME}_assembly_report.txt" \
        -O "${T2T_ASSEMBLY_REPORT}"
    echo "  Saved: ${T2T_VERSION}_assembly_report.txt"
else
    echo "Already exists: ${T2T_VERSION}_assembly_report.txt"
fi

T2T_FA="${T2T_DIR}/${T2T_VERSION}.fa"
T2T_GFF3_BED="${T2T_DIR}/${T2T_VERSION}_gff3_loci.bed"
T2T_CMSEARCH_BED="${T2T_DIR}/${T2T_VERSION}_cmsearch_loci.bed"
T2T_RRNA_BED="${T2T_DIR}/${T2T_VERSION}_rrna_loci.bed"

if [[ ! -f "${T2T_RRNA_BED}" ]] || [[ ! -s "${T2T_RRNA_BED}" ]]; then

    echo "Extracting rRNA loci from GFF3..."
    python3 "${UTILS_DIR}/extract_rrna_loci.py" "${T2T_GFF_GZ}" "${T2T_GFF3_BED}" \
        --margin "${RNA_LOCI_MARGIN}" --name-map "${T2T_ASSEMBLY_REPORT}"
    sort -k1,1 -k2,2n "${T2T_GFF3_BED}" \
        | bedtools merge > "${T2T_GFF3_BED}.tmp" \
        && mv "${T2T_GFF3_BED}.tmp" "${T2T_GFF3_BED}"
    gff3_regions=$(wc -l < "${T2T_GFF3_BED}")
    gff3_bp=$(awk '{sum+=$3-$2} END{print sum}' "${T2T_GFF3_BED}")
    echo "  GFF3 loci: ${gff3_regions} regions, ${gff3_bp} bp"

    if [[ -n "${CMS_DIR}" ]] && command -v cmsearch &>/dev/null; then
        if [[ ! -f "${T2T_FA}" ]]; then
            echo "Decompressing ${T2T_VERSION}.fa.gz for cmsearch..."
            zcat "${T2T_GENOME_GZ}" > "${T2T_FA}"
        fi

        echo ""
        echo "============================================"
        echo "Running cmsearch for rRNA loci (RF01960, RF02543, RF00001, RF00002, RF00177, RF02541)"
        echo "============================================"

        for cm in RF01960 RF02543 RF00001 RF00002; do
            tblout="${T2T_DIR}/${cm}_hits.tbl"
            if [[ ! -f "${tblout}" ]]; then
                case "${cm}" in
                    RF01960|RF02543) hmmonly="--hmmonly" ;;
                    *)               hmmonly="" ;;
                esac
                echo "Running cmsearch ${cm}${hmmonly:+ (--hmmonly)} ..."
                cmsearch ${hmmonly} --cut_ga --cpu "${THREADS}" \
                    --tblout "${tblout}" \
                    "${CMS_DIR}/${cm}.cm" "${T2T_FA}" > /dev/null
                echo "  Saved: $(basename "${tblout}")"
            else
                echo "Already exists: ${cm}_hits.tbl"
            fi
        done

        # Mitochondrial rRNA (12S = RF00177, 16S = RF02541) - search chrM only
        chrM_acc=$(awk -F'\t' '$3 == "MT" {print $5}' "${T2T_ASSEMBLY_REPORT}")
        CHRM_FA="${T2T_DIR}/chrM.fa"
        if [[ -n "${chrM_acc}" ]]; then
            if [[ ! -f "${CHRM_FA}" ]]; then
                echo "Extracting chrM (${chrM_acc}) for mt-rRNA cmsearch..."
                seqkit grep -p "${chrM_acc}" "${T2T_FA}" -w 0 > "${CHRM_FA}"
            fi
            for cm in RF00177 RF02541; do
                tblout="${T2T_DIR}/${cm}_hits.tbl"
                if [[ ! -f "${tblout}" ]]; then
                    echo "Running cmsearch ${cm} (mt-rRNA, chrM only)..."
                    cmsearch --cut_ga --cpu "${THREADS}" \
                        --tblout "${tblout}" \
                        "${CMS_DIR}/${cm}.cm" "${CHRM_FA}" > /dev/null
                    echo "  Saved: $(basename "${tblout}")"
                else
                    echo "Already exists: ${cm}_hits.tbl"
                fi
            done
        else
            echo "Warning: chrM accession not found in assembly report - skipping mt-rRNA masking"
        fi

        cms_tblouts=(
            "${T2T_DIR}/RF01960_hits.tbl" "${T2T_DIR}/RF02543_hits.tbl"
            "${T2T_DIR}/RF00001_hits.tbl" "${T2T_DIR}/RF00002_hits.tbl"
        )
        for cm in RF00177 RF02541; do
            [[ -f "${T2T_DIR}/${cm}_hits.tbl" ]] && cms_tblouts+=("${T2T_DIR}/${cm}_hits.tbl")
        done

        awk -v margin="${RNA_LOCI_MARGIN}" '
            /^#/ { next }
            $17 != "!" { next }
            {
                s = ($8 < $9) ? $8 : $9
                e = ($8 < $9) ? $9 : $8
                start = (s - 1 - margin < 0) ? 0 : s - 1 - margin
                print $1 "\t" start "\t" (e + margin)
            }
        ' "${cms_tblouts[@]}" \
        | sort -k1,1 -k2,2n \
        | bedtools merge \
        > "${T2T_CMSEARCH_BED}"

        cms_regions=$(wc -l < "${T2T_CMSEARCH_BED}")
        cms_bp=$(awk '{sum+=$3-$2} END{print sum}' "${T2T_CMSEARCH_BED}")
        extra_bp=$(bedtools subtract -a "${T2T_CMSEARCH_BED}" -b "${T2T_GFF3_BED}" \
            | awk '{sum+=$3-$2} END{print sum+0}')
        echo "  cmsearch loci: ${cms_regions} regions, ${cms_bp} bp"
        echo "  Extra vs GFF3: ${extra_bp} bp not covered by GFF3 annotation"

        echo ""
        echo "Extracting per-family sequences found by cmsearch but not in GFF3..."
        for cm in RF01960 RF02543 RF00001 RF00002 RF00177 RF02541; do
            tblout="${T2T_DIR}/${cm}_hits.tbl"
            [[ ! -f "${tblout}" ]] && continue
            extra_bed="${T2T_DIR}/${cm}_extra_vs_gff3.bed"
            extra_fa="${T2T_DIR}/${cm}_extra_vs_gff3.fasta"
            awk -v margin="${RNA_LOCI_MARGIN}" '
                /^#/ { next }
                $17 != "!" { next }
                {
                    s = ($8 < $9) ? $8 : $9
                    e = ($8 < $9) ? $9 : $8
                    start = (s - 1 - margin < 0) ? 0 : s - 1 - margin
                    print $1 "\t" start "\t" (e + margin)
                }
            ' "${tblout}" \
            | sort -k1,1 -k2,2n \
            | bedtools merge \
            | bedtools subtract -a - -b "${T2T_GFF3_BED}" \
            > "${extra_bed}"
            n_regions=$(wc -l < "${extra_bed}")
            if [[ "${n_regions}" -gt 0 ]]; then
                bedtools getfasta -fi "${T2T_FA}" -bed "${extra_bed}" -fo "${extra_fa}"
                echo "  ${cm}: ${n_regions} regions not in GFF3 -> $(basename "${extra_fa}")"
            else
                echo "  ${cm}: no regions outside GFF3"
                rm -f "${extra_bed}"
            fi
        done

        # Combine GFF3 and cmsearch loci, sort by chrom/start (required by bedtools merge),
        # then collapse overlapping intervals into a non-redundant union BED.
        cat "${T2T_GFF3_BED}" "${T2T_CMSEARCH_BED}" \
            | sort -k1,1 -k2,2n \
            | bedtools merge \
            > "${T2T_RRNA_BED}"
    else
        [[ -z "${CMS_DIR}" ]] && echo "CMS_DIR not set - using GFF3 annotation only."
        cp "${T2T_GFF3_BED}" "${T2T_RRNA_BED}"
    fi

    total_regions=$(wc -l < "${T2T_RRNA_BED}")
    total_bp=$(awk '{sum+=$3-$2} END{print sum}' "${T2T_RRNA_BED}")
    echo "  Final BED: ${total_regions} regions, ${total_bp} bp -> ${T2T_RRNA_BED}"
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

echo ""
echo "Combining Rfam non-rRNA sequences..."

find "${RFAM_DIR}" -name "*.fa" | sort | xargs cat > "${OUTPUT_DIR}/rfam_non_rrna_all.fasta"

echo "Rfam non-rRNA: $(seqkit stats -T "${OUTPUT_DIR}/rfam_non_rrna_all.fasta" | tail -1 | cut -f4) sequences total"

echo ""
echo "============================================"
echo "Download complete"
echo "============================================"
echo ""
echo "Outputs:"
echo "  T2T genome:        ${T2T_DIR}/${T2T_VERSION}.fa.gz"
echo "  T2T annotation:    ${T2T_DIR}/${T2T_VERSION}_annotation.gff.gz"
echo "  Assembly report:   ${T2T_DIR}/${T2T_VERSION}_assembly_report.txt"
echo "  GFF3 loci BED:     ${T2T_GFF3_BED}"
[[ -f "${T2T_CMSEARCH_BED}" ]] && echo "  cmsearch loci BED: ${T2T_CMSEARCH_BED}"
echo "  Final rRNA BED:    ${T2T_RRNA_BED}"
echo "  Rfam all:          ${OUTPUT_DIR}/rfam_non_rrna_all.fasta"
echo ""
echo "Next step: Run simulate_non_rrna.sh to mask rRNA loci, simulate T2T reads"
echo "           with InSilicoSeq, and filter Rfam sequences by length"
