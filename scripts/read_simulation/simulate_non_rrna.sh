#!/usr/bin/env bash

################################################################################
# simulate_non_rrna.sh
#
# Mask rRNA loci in the T2T genome, simulate Illumina reads with InSilicoSeq,
# and sample Rfam non-rRNA sequences to produce two test sets:
#
#   non_rRNA_test_<N>_T2T.fasta  - N simulated T2T genome reads (rRNA loci masked, default 10M)
#   non_rRNA_test_Rfam.fasta     - Rfam non-rRNA sequences sampled evenly across families
#
# Inputs (from download_non_rrna.sh):
#   t2t/${T2T_VERSION}.fa.gz         - T2T genome (compressed)
#   t2t/${T2T_VERSION}_rrna_loci.bed - rRNA loci BED (from extract_rrna_loci.py)
#   rfam/RF*.fa                      - Rfam non-rRNA family FASTA files (decompressed)
#
# Steps:
#   1. Decompress T2T genome
#   2. Mask rRNA loci with bedtools maskfasta
#   3. Simulate N Illumina PE reads with InSilicoSeq (default: 10M)
#   4. Convert to FASTA -> non_rRNA_test_<N>_T2T.fasta
#   5. Sample N_RFAM/n_families sequences from each Rfam family -> non_rRNA_test_Rfam.fasta
#   6. Write HTML summary report
#
# Usage: bash simulate_non_rrna.sh [output_dir [threads]] [OPTIONS]
#
# Positional:
#   output_dir      Directory with download_non_rrna.sh outputs (default: $NON_RRNA_DIR or data/non_rrna)
#   threads         Number of threads for InSilicoSeq (default: 4)
#
# Options:
#   --t2t-reads INT   Number of simulated T2T reads to include (default: 10000000)
#                     ISS --n_reads splits evenly between R1 and R2, so R1+R2 combined = N_T2T.
#                     10M is the default so that subsampling covers all scalability scale points
#                     (10K, 100K, 1M, 10M) without re-running ISS.
#   --model STR       InSilicoSeq error model: HiSeq, NovaSeq, MiSeq (default: NovaSeq)
#   --rfam-reads INT  Total Rfam sequences to sample across all families (default: 500000)
#                     Fair-share allocation: quota redistributed from small families to large
#                     ones so the target is met exactly (see fair_share_rfam.py).
#   --seed INT        Random seed for InSilicoSeq simulation and Rfam sampling (default: 42)
#   -h, --help        Show help
#
# Environment variables:
#   T2T_VERSION         T2T genome version string (default: chm13v2.0)
#   T2T_ACCESSION       GCA accession (default: GCA_009914755.4)
#   T2T_GCF_ACCESSION   GCF accession (default: GCF_009914755.1)
#   T2T_NAME            Assembly name (default: T2T-CHM13v2.0)
#
# Requires: bedtools, iss (InSilicoSeq), seqkit
#
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UTILS_DIR="${SCRIPT_DIR}/../utils"

POSITIONAL=()
N_T2T=10000000
ISS_MODEL=NovaSeq
N_RFAM_READS=500000
RAND_SEED=42

T2T_VERSION="${T2T_VERSION:-chm13v2.0}"
T2T_ACCESSION="${T2T_ACCESSION:-GCA_009914755.4}"
T2T_GCF_ACCESSION="${T2T_GCF_ACCESSION:-GCF_009914755.1}"
T2T_NAME="${T2T_NAME:-T2T-CHM13v2.0}"

show_help() {
    grep '^#' "$0" | grep -v '#!/usr/bin/env bash' | sed 's/^# \?//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --t2t-reads) N_T2T="$2"; shift 2 ;;
        --model) ISS_MODEL="$2"; shift 2 ;;
        --rfam-reads) N_RFAM_READS="$2"; shift 2 ;;
        --seed) RAND_SEED="$2"; shift 2 ;;
        -h|--help) show_help ;;
        *) POSITIONAL+=("$1"); shift ;;
    esac
done

OUTPUT_DIR="${POSITIONAL[0]:-${NON_RRNA_DIR:-data/non_rrna}}"
THREADS="${POSITIONAL[1]:-4}"

OUTPUT_DIR="$(mkdir -p "${OUTPUT_DIR}" && cd "${OUTPUT_DIR}" && pwd)"

T2T_DIR="${OUTPUT_DIR}/t2t"
ISS_DIR="${OUTPUT_DIR}/iss"
RFAM_DIR="${OUTPUT_DIR}/rfam"

T2T_GENOME_GZ="${T2T_DIR}/${T2T_VERSION}.fa.gz"
T2T_FA="${T2T_DIR}/${T2T_VERSION}.fa"
T2T_RRNA_BED="${T2T_DIR}/${T2T_VERSION}_rrna_loci.bed"
T2T_MASKED="${T2T_DIR}/${T2T_VERSION}_masked.fa"

if   (( N_T2T % 1000000 == 0 )); then N_T2T_LABEL="$(( N_T2T / 1000000 ))M"
elif (( N_T2T % 1000    == 0 )); then N_T2T_LABEL="$(( N_T2T / 1000 ))K"
else                                   N_T2T_LABEL="${N_T2T}"
fi
T2T_OUTPUT="${OUTPUT_DIR}/non_rRNA_test_${N_T2T_LABEL}_T2T.fasta"

# Simulate 5% more reads than needed to absorb N-heavy reads removed by the ambiguity filter.
# Round up to the nearest even number (ISS splits --n_reads evenly between R1 and R2).
N_T2T_SIM=$(( N_T2T + N_T2T / 20 ))
(( N_T2T_SIM % 2 != 0 )) && (( N_T2T_SIM++ ))
RFAM_OUTPUT="${OUTPUT_DIR}/non_rRNA_test_Rfam.fasta"
OUTPUT_HTML="${OUTPUT_DIR}/non_rrna_test_set_summary.html"

echo "============================================"
echo "Non-rRNA Read Simulation Script"
echo "Output directory: ${OUTPUT_DIR}"
echo "T2T version:      ${T2T_VERSION}"
echo "T2T reads:        ${N_T2T}"
echo "ISS model:        ${ISS_MODEL}"
echo "Rfam reads:       ${N_RFAM_READS} total (evenly split across families)"
echo "Random seed:      ${RAND_SEED}"
echo "Threads:          ${THREADS}"
echo "============================================"
echo ""

for tool in bedtools iss seqkit fastp; do
    if ! command -v "${tool}" &> /dev/null; then
        echo "Error: ${tool} not found. Install via: conda env create -f environment.yml"
        exit 1
    fi
done

for f in "${T2T_GENOME_GZ}" "${T2T_RRNA_BED}"; do
    if [[ ! -f "${f}" ]]; then
        echo "Error: required input not found: ${f}"
        echo "Run download_non_rrna.sh first."
        exit 1
    fi
done

if ! ls "${RFAM_DIR}"/RF*.fa >/dev/null 2>&1; then
    echo "Error: Rfam family files not found in ${RFAM_DIR}"
    echo "Run download_non_rrna.sh first."
    exit 1
fi

################################################################################
# 1. DECOMPRESS T2T GENOME
################################################################################

echo "============================================"
echo "Step 1: Decompress T2T genome"
echo "============================================"

if [[ ! -f "${T2T_FA}" ]]; then
    echo "Decompressing ${T2T_VERSION}.fa.gz..."
    zcat "${T2T_GENOME_GZ}" > "${T2T_FA}"
    echo "  Saved: ${T2T_VERSION}.fa"
else
    echo "Already exists: ${T2T_VERSION}.fa"
fi

################################################################################
# 2. MASK rRNA LOCI
################################################################################

echo ""
echo "============================================"
echo "Step 2: Mask rRNA loci"
echo "============================================"

missing=0
while IFS= read -r chrom; do
    if ! grep -qPm1 "^>${chrom}(\s|$)" "${T2T_FA}"; then
        echo "Error: BED chromosome '${chrom}' not found in ${T2T_FA}."
        missing=1
    fi
done < <(cut -f1 "${T2T_RRNA_BED}" | sort -u)
if [[ "${missing}" -eq 1 ]]; then
    echo "Chromosome names in the BED must match FASTA headers."
    echo "Delete ${T2T_RRNA_BED} and re-run download_non_rrna.sh to regenerate."
    exit 1
fi

n_loci=$(wc -l < "${T2T_RRNA_BED}")
masked_bp=$(awk '{sum += $3-$2} END{print sum}' "${T2T_RRNA_BED}")

if [[ ! -f "${T2T_MASKED}" ]]; then
    echo "Masking ${n_loci} rRNA loci with bedtools maskfasta..."
    n_before=$(grep -v "^>" "${T2T_FA}" | tr -cd 'Nn' | wc -c)
    bedtools maskfasta \
        -fi "${T2T_FA}" \
        -bed "${T2T_RRNA_BED}" \
        -fo "${T2T_MASKED}"
    echo "  Saved: ${T2T_VERSION}_masked.fa"
    n_after=$(grep -v "^>" "${T2T_MASKED}" | tr -cd 'Nn' | wc -c)
    n_added=$(( n_after - n_before ))
    if [[ "${n_added}" -ne "${masked_bp}" ]]; then
        echo "Warning: newly masked bases (${n_added}) != BED total (${masked_bp})."
        echo "  This may indicate overlapping BED regions or intervals clamped at chromosome ends."
    fi
else
    echo "Already exists: ${T2T_VERSION}_masked.fa"
fi

################################################################################
# 3. SIMULATE ILLUMINA READS WITH InSilicoSeq
################################################################################

echo ""
echo "============================================"
echo "Step 3: Simulate Illumina PE reads (InSilicoSeq)"
echo "============================================"

mkdir -p "${ISS_DIR}"
ISS_PREFIX="${ISS_DIR}/${T2T_VERSION}_masked"
ISS_R1="${ISS_PREFIX}_R1.fastq"
ISS_R2="${ISS_PREFIX}_R2.fastq"

if [[ ! -f "${ISS_R1}" ]] || [[ ! -f "${ISS_R2}" ]]; then
    echo "Generating ${N_T2T_SIM} reads (${N_T2T} target + 5% buffer) with InSilicoSeq..."
    iss generate \
        --genomes "${T2T_MASKED}" \
        --model "${ISS_MODEL}" \
        --n_reads "${N_T2T_SIM}" \
        --cpus "${THREADS}" \
        --seed "${RAND_SEED}" \
        --gc_bias \
        --output "${ISS_PREFIX}"
    echo "  Saved: $(basename "${ISS_R1}"), $(basename "${ISS_R2}")"
else
    echo "Already exists: InSilicoSeq output files"
fi

################################################################################
# 3b. FILTER AND QC WITH FASTP
################################################################################

echo ""
echo "============================================"
echo "Step 3b: Filter reads with fastp"
echo "============================================"

FASTP_R1="${ISS_DIR}/t2t_filtered_R1.fastq"
FASTP_R2="${ISS_DIR}/t2t_filtered_R2.fastq"
FASTP_JSON="${ISS_DIR}/fastp.json"

fastp \
    --in1 "${ISS_R1}" --in2 "${ISS_R2}" \
    --out1 "${FASTP_R1}" --out2 "${FASTP_R2}" \
    --n_base_limit 10 \
    --disable_adapter_trimming \
    --json "${FASTP_JSON}" \
    --thread "${THREADS}"

read -r FASTP_READS_BEFORE FASTP_READS_AFTER FASTP_REMOVED_N FASTP_REMOVED_Q \
         FASTP_AVG_LEN FASTP_Q20 FASTP_Q30 FASTP_GC \
    < <(python3 - "${FASTP_JSON}" <<'PYEOF'
import json, sys
d = json.load(open(sys.argv[1]))
bf = d['summary']['before_filtering']
fr = d['filtering_result']
print(
    bf['total_reads'],
    fr['passed_filter_reads'],
    fr['too_many_N_reads'],
    fr['low_quality_reads'],
    round(bf['total_bases'] / bf['total_reads'], 1),
    round(bf['q20_rate'] * 100, 2),
    round(bf['q30_rate'] * 100, 2),
    round(bf['gc_content'] * 100, 2),
)
PYEOF
)
echo "  Before: ${FASTP_READS_BEFORE}  After: ${FASTP_READS_AFTER}"
echo "  Removed - N content: ${FASTP_REMOVED_N}  Low quality: ${FASTP_REMOVED_Q}"
echo "  Q20: ${FASTP_Q20}%  Q30: ${FASTP_Q30}%  GC: ${FASTP_GC}%"

################################################################################
# 4. CONVERT T2T READS TO FASTA
################################################################################

echo ""
echo "============================================"
echo "Step 4: Convert T2T reads to FASTA"
echo "============================================"

ISS_FILTERED_FA="${ISS_DIR}/t2t_filtered.fasta"
echo "Converting filtered FASTQ to FASTA..."
cat "${FASTP_R1}" "${FASTP_R2}" | seqkit fq2fa > "${ISS_FILTERED_FA}"

if (( FASTP_READS_AFTER < N_T2T )); then
    echo "Error: only ${FASTP_READS_AFTER} reads after filtering, need ${N_T2T}. Increase buffer."
    exit 1
fi

echo "Sampling ${N_T2T} reads..."
seqkit sample -n "${N_T2T}" --rand-seed "${RAND_SEED}" "${ISS_FILTERED_FA}" -o "${T2T_OUTPUT}"
n_t2t=$(seqkit stats -T "${T2T_OUTPUT}" | tail -1 | cut -f4)
echo "  Saved: $(basename "${T2T_OUTPUT}") (${n_t2t} reads)"

################################################################################
# 5. SAMPLE Rfam FAMILIES
################################################################################

echo ""
echo "============================================"
echo "Step 5: Sample ${N_RFAM_READS} total sequences from Rfam families (fair-share)"
echo "============================================"

# Compute per-family allocations: evenly distributed, with remainder
# redistributed from small families to large ones (see fair_share_rfam.py).
declare -A alloc_map
while IFS=$'\t' read -r stem n; do
    alloc_map["${stem}"]="${n}"
done < <(python3 "${UTILS_DIR}/fair_share_rfam.py" "${N_RFAM_READS}" "${RFAM_DIR}")

> "${RFAM_OUTPUT}"
rfam_html_rows=""
rfam_total_seqs=0

for fa in "${RFAM_DIR}"/RF*.fa; do
    [[ -f "${fa}" ]] || continue
    family_file=$(basename "${fa}" .fa)
    rfam_id="${family_file%%_*}"
    family_name="${family_file#*_}"
    n_alloc="${alloc_map[${family_file}]:-0}"

    stats=$(seqkit stats -T "${fa}" | tail -1)
    n_total=$(echo "${stats}" | cut -f4)
    min_len=$(echo "${stats}" | cut -f6)
    avg_len=$(echo "${stats}" | cut -f7)
    max_len=$(echo "${stats}" | cut -f8)

    rfam_tmp=$(mktemp)
    seqkit sample -n "${n_alloc}" --rand-seed "${RAND_SEED}" "${fa}" | seqkit seq -w 0 > "${rfam_tmp}"
    n_sampled=$(seqkit stats -T "${rfam_tmp}" | tail -1 | cut -f4)
    cat "${rfam_tmp}" >> "${RFAM_OUTPUT}"
    rm -f "${rfam_tmp}"

    rfam_total_seqs=$(( rfam_total_seqs + n_total ))
    echo "  ${rfam_id} (${family_name}): ${n_total} total, ${min_len}-${max_len} bp (avg ${avg_len}) -> ${n_sampled} sampled"
    rfam_html_rows="${rfam_html_rows}      <tr><td>${family_name//_/ }</td><td>${rfam_id}</td><td>${n_total}</td><td>${min_len}-${max_len}</td><td>${avg_len}</td><td>${n_sampled}</td></tr>\n"
done

n_rfam=$(seqkit stats -T "${RFAM_OUTPUT}" | tail -1 | cut -f4)
echo "  Total: ${n_rfam} sequences -> non_rRNA_test_Rfam.fasta"

################################################################################
# 6. HTML SUMMARY REPORT
################################################################################

echo ""
echo "============================================"
echo "Step 6: Write HTML summary"
echo "============================================"

# Gather cmsearch stats for HTML (if cmsearch was run in download_non_rrna.sh)
T2T_GFF3_BED="${T2T_DIR}/${T2T_VERSION}_gff3_loci.bed"
T2T_CMSEARCH_BED="${T2T_DIR}/${T2T_VERSION}_cmsearch_loci.bed"

CMS_USED=false
GFF3_REGIONS=0; GFF3_BP=0; CMS_REGIONS=0; CMS_BP=0; EXTRA_BP=0
cms_family_rows=""

if [[ -f "${T2T_CMSEARCH_BED}" ]] && [[ -f "${T2T_GFF3_BED}" ]]; then
    CMS_USED=true
    GFF3_REGIONS=$(wc -l < "${T2T_GFF3_BED}")
    GFF3_BP=$(awk '{sum+=$3-$2} END{print sum}' "${T2T_GFF3_BED}")
    CMS_REGIONS=$(wc -l < "${T2T_CMSEARCH_BED}")
    CMS_BP=$(awk '{sum+=$3-$2} END{print sum}' "${T2T_CMSEARCH_BED}")
    EXTRA_BP=$(bedtools subtract -a "${T2T_CMSEARCH_BED}" -b "${T2T_GFF3_BED}" \
        | awk '{sum+=$3-$2} END{print sum+0}')

    for cm_id in RF01960 RF02543 RF00001 RF00002 RF00177 RF02541; do
        extra_bed="${T2T_DIR}/${cm_id}_extra_vs_gff3.bed"
        [[ ! -f "${extra_bed}" ]] && continue
        n_extra=$(wc -l < "${extra_bed}")
        case "${cm_id}" in
            RF01960) cm_name="18S SSU rRNA"  ;;
            RF02543) cm_name="28S LSU rRNA"  ;;
            RF00001) cm_name="5S rRNA"       ;;
            RF00002) cm_name="5.8S rRNA"     ;;
            RF00177) cm_name="12S mt-rRNA"   ;;
            RF02541) cm_name="16S mt-rRNA"   ;;
        esac
        cms_family_rows="${cms_family_rows}      <tr><td>${cm_name}</td><td>${cm_id}</td><td>${n_extra}</td></tr>\n"
    done
fi

T2T_ACCESSION="${T2T_ACCESSION}" \
T2T_GCF_ACCESSION="${T2T_GCF_ACCESSION}" \
T2T_NAME="${T2T_NAME}" \
T2T_VERSION="${T2T_VERSION}" \
ISS_MODEL="${ISS_MODEL}" \
RAND_SEED="${RAND_SEED}" \
N_T2T_SIM="${FASTP_READS_BEFORE}" \
FASTP_REMOVED_N="${FASTP_REMOVED_N}" \
FASTP_REMOVED_Q="${FASTP_REMOVED_Q}" \
FASTP_AVG_LEN="${FASTP_AVG_LEN}" \
FASTP_Q20="${FASTP_Q20}" \
FASTP_Q30="${FASTP_Q30}" \
FASTP_GC="${FASTP_GC}" \
N_T2T="${n_t2t}" \
N_T2T_LABEL="${N_T2T_LABEL}" \
N_LOCI="${n_loci}" \
MASKED_BP="${masked_bp}" \
N_RFAM_READS="${N_RFAM_READS}" \
N_RFAM="${n_rfam}" \
RFAM_TOTAL_SEQS="${rfam_total_seqs}" \
CMS_USED="${CMS_USED}" \
GFF3_REGIONS="${GFF3_REGIONS}" \
GFF3_BP="${GFF3_BP}" \
CMS_REGIONS="${CMS_REGIONS}" \
CMS_BP="${CMS_BP}" \
EXTRA_BP="${EXTRA_BP}" \
OUTPUT_HTML="${OUTPUT_HTML}" \
python3 - "${rfam_html_rows}" "${cms_family_rows}" <<'PYEOF'
import sys, os, datetime

rfam_rows_raw = sys.argv[1]
cms_rows_raw  = sys.argv[2]

e = os.environ

# Build optional cmsearch section (inserted as a pre-built string to avoid
# brace-escaping issues with Python's .format() in the main template)
cms_used = e.get('CMS_USED', 'false') == 'true'
if cms_used:
    cms_rows = cms_rows_raw.replace('\\n', '\n')
    cms_section = (
        '<section>\n'
        '<h2>cmsearch rRNA Annotation Supplement</h2>\n'
        '<div class="description">\n'
        '<p>The RefSeq GFF3 annotation was supplemented with Infernal cmsearch to identify '
        'unannotated rRNA copies and rRNA pseudogenes not individually annotated in RefSeq. '
        'Covariance models RF01960 (18S), RF02543 (28S), RF00001 (5S), and RF00002 (5.8S) '
        'were run with --cut_ga; --hmmonly was applied for the two large subunit models.</p>\n'
        '</div>\n'
        '<div class="table-wrap"><table>\n'
        '  <thead><tr><th>Source</th><th>Regions</th><th>Bases (bp)</th></tr></thead>\n'
        '  <tbody>\n'
        f'    <tr><td>GFF3 annotation only</td><td>{e["GFF3_REGIONS"]}</td><td>{e["GFF3_BP"]}</td></tr>\n'
        f'    <tr><td>cmsearch (all families combined)</td><td>{e["CMS_REGIONS"]}</td><td>{e["CMS_BP"]}</td></tr>\n'
        f'    <tr><td>Extra bp found by cmsearch vs GFF3</td><td>-</td><td>{e["EXTRA_BP"]}</td></tr>\n'
        f'    <tr style="font-weight:bold; border-top: 2px solid #2c3e50;"><td>Final merged BED (GFF3 + cmsearch)</td><td>{e["N_LOCI"]}</td><td>{e["MASKED_BP"]}</td></tr>\n'
        '  </tbody>\n'
        '</table></div>\n'
        '<h3>Per-family regions found by cmsearch but absent from GFF3</h3>\n'
        '<div class="table-wrap"><table>\n'
        '  <thead><tr><th>Family</th><th>Rfam ID</th><th>Regions not in GFF3</th></tr></thead>\n'
        '  <tbody>\n'
        + cms_rows +
        '  </tbody>\n'
        '</table></div>\n'
        '</section>\n'
    )
else:
    cms_section = ''

html = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Non-rRNA Reference Sources Summary for Reads Simulation</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 40px; color: #333; }}
  h1   {{ color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 8px; }}
  h2   {{ color: #34495e; margin-top: 32px; }}
  h3   {{ color: #34495e; margin-top: 20px; }}
  .description {{ background: #f8f9fa; border-left: 4px solid #3498db; padding: 10px 16px; margin: 12px 0; }}
  .table-wrap  {{ overflow-x: auto; margin: 16px 0; }}
  table  {{ border-collapse: collapse; min-width: 600px; }}
  th, td {{ border: 1px solid #ddd; padding: 8px 14px; text-align: left; }}
  th     {{ background: #2c3e50; color: #fff; }}
  tr:nth-child(even) {{ background: #f2f2f2; }}
  .meta  {{ font-size: 0.9em; color: #666; margin-top: 4px; }}
  .footer {{ font-size: 0.8em; color: #999; margin-top: 40px; border-top: 1px solid #ddd; padding-top: 8px; }}
</style>
</head>
<body>
<h1>Non-rRNA Reference Sources Summary for Reads Simulation</h1>
<p class="meta">Generated: {date}</p>

<section>
<h2>T2T Genome (Human CHM13v2.0)</h2>
<div class="description">
<p>Illumina 150bp PE reads simulated from the T2T CHM13v2.0 human genome with rRNA loci
masked prior to simulation. Masking prevents false negatives where a simulated read
happens to span an rRNA region.</p>
</div>
<div class="table-wrap"><table>
  <thead><tr><th>Parameter</th><th>Value</th></tr></thead>
  <tbody>
    <tr><td>Assembly name</td><td>{t2t_name}</td></tr>
    <tr><td>GCA accession</td><td>{t2t_accession}</td></tr>
    <tr><td>GCF accession (annotation source)</td><td>{t2t_gcf_accession}</td></tr>
    <tr><td>Version</td><td>{t2t_version}</td></tr>
    <tr><td>rRNA loci masked</td><td>{n_loci}</td></tr>
    <tr><td>Total bases masked (including margins)</td><td>{masked_bp}</td></tr>
    <tr><td>Simulation tool</td><td>InSilicoSeq</td></tr>
    <tr><td>Error model</td><td>{iss_model}</td></tr>
    <tr><td>Read length</td><td>150 bp (paired-end)</td></tr>
    <tr><td>Random seed</td><td>{rand_seed}</td></tr>
    <tr><td>Reads simulated (with 5% buffer)</td><td>{n_t2t_sim}</td></tr>
    <tr><td>Reads removed (>10 N bases)</td><td>{fastp_removed_n}</td></tr>
    <tr><td>Reads removed (low quality)</td><td>{fastp_removed_q}</td></tr>
    <tr><td>Reads in test set</td><td>{n_t2t}</td></tr>
    <tr><td>Output file</td><td>non_rRNA_test_{n_t2t_label}_T2T.fasta</td></tr>
  </tbody>
</table></div>
</section>

<section>
<h2>Read Quality (fastp, before filtering)</h2>
<div class="description">
<p>Quality statistics from fastp computed on R1+R2 before filtering. fastp filters read pairs
where either read has &gt;10 ambiguous bases or fails quality thresholds.</p>
</div>
<div class="table-wrap"><table>
  <thead><tr><th>Metric</th><th>Value</th></tr></thead>
  <tbody>
    <tr><td>Average read length (bp)</td><td>{fastp_avg_len}</td></tr>
    <tr><td>Q20 bases (%)</td><td>{fastp_q20}</td></tr>
    <tr><td>Q30 bases (%)</td><td>{fastp_q30}</td></tr>
    <tr><td>GC content (%)</td><td>{fastp_gc}</td></tr>
  </tbody>
</table></div>
</section>

{cms_section}
<section>
<h2>Rfam Non-rRNA Families</h2>
<div class="description">
<p>Target total: {n_rfam_reads} sequences sampled using fair-share allocation across families.
Families smaller than the per-family quota contribute all their sequences; the remainder
is redistributed to larger families so the target is met exactly (see fair_share_rfam.py).
No read simulation applied - sequences are used as-is to test whether SortMeRNA
correctly rejects structurally complex non-rRNA sequences.</p>
</div>
<div class="table-wrap"><table>
  <thead>
    <tr>
      <th>Family</th><th>Rfam ID</th><th>Total sequences</th>
      <th>Length range (bp)</th><th>Average length (bp)</th><th>Sequences sampled</th>
    </tr>
  </thead>
  <tbody>
{rfam_rows}    <tr style="font-weight:bold; border-top: 2px solid #2c3e50;"><td>Total</td><td></td><td>{rfam_total_seqs}</td><td></td><td></td><td>{n_rfam}</td></tr>
  </tbody>
</table></div>
</section>

<p class="footer">Generated by simulate_non_rrna.sh - SortMeRNA database pipeline</p>
</body>
</html>
""".format(
    date=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    t2t_name=e['T2T_NAME'],
    t2t_accession=e['T2T_ACCESSION'],
    t2t_gcf_accession=e['T2T_GCF_ACCESSION'],
    t2t_version=e['T2T_VERSION'],
    n_loci=e['N_LOCI'],
    masked_bp=e['MASKED_BP'],
    iss_model=e['ISS_MODEL'],
    rand_seed=e['RAND_SEED'],
    n_t2t_sim=e['N_T2T_SIM'],
    fastp_removed_n=e['FASTP_REMOVED_N'],
    fastp_removed_q=e['FASTP_REMOVED_Q'],
    fastp_avg_len=e['FASTP_AVG_LEN'],
    fastp_q20=e['FASTP_Q20'],
    fastp_q30=e['FASTP_Q30'],
    fastp_gc=e['FASTP_GC'],
    n_t2t=e['N_T2T'],
    n_t2t_label=e['N_T2T_LABEL'],
    n_rfam_reads=e['N_RFAM_READS'],
    rfam_rows=rfam_rows_raw.replace('\\n', '\n'),
    rfam_total_seqs=e['RFAM_TOTAL_SEQS'],
    n_rfam=e['N_RFAM'],
    cms_section=cms_section,
)

open(e['OUTPUT_HTML'], 'w').write(html)
PYEOF

echo "  Done - ${OUTPUT_HTML}"

echo ""
echo "============================================"
echo "Simulation complete"
echo "============================================"
echo ""
echo "Outputs:"
echo "  Masked genome:  ${T2T_MASKED}"
echo "  T2T test set:   ${T2T_OUTPUT} (${n_t2t} reads)"
echo "  Rfam test set:  ${RFAM_OUTPUT} (${n_rfam} sequences)"
echo "  HTML summary:   ${OUTPUT_HTML}"
echo ""
echo "Next step: Run simulate_rrna_reads.sh to simulate rRNA reads for sensitivity testing"
