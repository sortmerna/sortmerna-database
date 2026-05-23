#!/usr/bin/env bash

################################################################################
# simulate_non_rrna.sh
#
# Mask rRNA loci in the T2T genome, simulate Illumina reads with InSilicoSeq,
# and sample Rfam non-rRNA sequences to produce two test sets:
#
#   non_rRNA_test_1M_T2T.fasta   - 1M simulated T2T genome reads (rRNA loci masked)
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
#   3. Simulate 1M Illumina PE reads with InSilicoSeq
#   4. Convert to FASTA -> non_rRNA_test_1M_T2T.fasta
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
#   --t2t-reads INT   Number of simulated T2T reads to include (default: 1000000)
#                     InSilicoSeq is run with N_T2T/2 read pairs so that R1+R2
#                     combined equals exactly N_T2T individual reads.
#   --model STR       InSilicoSeq error model: HiSeq, NovaSeq, MiSeq (default: HiSeq)
#   --rfam-reads INT  Total Rfam sequences to sample across all families (default: 500000)
#                     Fair-share allocation: quota redistributed from small families to large
#                     ones so the target is met exactly (see fair_share_rfam.py).
#   --seed INT        Random seed for InSilicoSeq simulation and Rfam sampling (default: 42)
#   --skip-mask       Skip masking step, use existing masked genome
#   --skip-sim        Skip InSilicoSeq simulation, use existing output
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
N_T2T=1000000
ISS_MODEL=HiSeq
N_RFAM_READS=500000
RAND_SEED=42
SKIP_MASK=false
SKIP_SIM=false

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
        --skip-mask) SKIP_MASK=true; shift ;;
        --skip-sim) SKIP_SIM=true; shift ;;
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

T2T_OUTPUT="${OUTPUT_DIR}/non_rRNA_test_1M_T2T.fasta"
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

for tool in bedtools iss seqkit; do
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

if [[ "${SKIP_MASK}" == false ]]; then
    if [[ ! -f "${T2T_MASKED}" ]]; then
        n_loci=$(wc -l < "${T2T_RRNA_BED}")
        echo "Masking ${n_loci} rRNA loci with bedtools maskfasta..."
        bedtools maskfasta \
            -fi "${T2T_FA}" \
            -bed "${T2T_RRNA_BED}" \
            -fo "${T2T_MASKED}"
        echo "  Saved: ${T2T_VERSION}_masked.fa"
    else
        echo "Already exists: ${T2T_VERSION}_masked.fa"
    fi
else
    echo "Skipping mask (--skip-mask)."
    if [[ ! -f "${T2T_MASKED}" ]]; then
        echo "Error: --skip-mask set but ${T2T_MASKED} does not exist"
        exit 1
    fi
fi

n_loci=$(wc -l < "${T2T_RRNA_BED}")
masked_bp=$(awk '{sum += $3-$2} END{print sum}' "${T2T_RRNA_BED}")

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

if [[ "${SKIP_SIM}" == false ]]; then
    if [[ ! -f "${ISS_R1}" ]] || [[ ! -f "${ISS_R2}" ]]; then
        # iss generates n_reads PE pairs, so R1+R2 combined = 2 * n_reads.
        # Use N_T2T/2 pairs to get exactly N_T2T individual reads after combining.
        N_PAIRS=$(( N_T2T / 2 ))
        echo "Generating ${N_PAIRS} PE pairs (${N_T2T} individual reads) with InSilicoSeq..."
        iss generate \
            --genomes "${T2T_MASKED}" \
            --model "${ISS_MODEL}" \
            --n_reads "${N_PAIRS}" \
            --cpus "${THREADS}" \
            --seed "${RAND_SEED}" \
            --output "${ISS_PREFIX}"
        echo "  Saved: $(basename "${ISS_R1}"), $(basename "${ISS_R2}")"
    else
        echo "Already exists: InSilicoSeq output files"
    fi
else
    echo "Skipping simulation (--skip-sim)."
    if [[ ! -f "${ISS_R1}" ]] || [[ ! -f "${ISS_R2}" ]]; then
        echo "Error: --skip-sim set but ISS output not found: ${ISS_R1}"
        exit 1
    fi
fi

################################################################################
# 4. CONVERT T2T READS TO FASTA
################################################################################

echo ""
echo "============================================"
echo "Step 4: Convert T2T reads to FASTA"
echo "============================================"

echo "Converting ISS FASTQ to FASTA..."
cat "${ISS_R1}" "${ISS_R2}" | seqkit fq2fa > "${T2T_OUTPUT}"

n_t2t=$(seqkit stats -T "${T2T_OUTPUT}" | tail -1 | cut -f4)
echo "  Saved: non_rRNA_test_1M_T2T.fasta (${n_t2t} reads)"

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
    max_len=$(echo "${stats}" | cut -f8)

    rfam_tmp=$(mktemp)
    seqkit sample -n "${n_alloc}" --rand-seed "${RAND_SEED}" "${fa}" > "${rfam_tmp}"
    n_sampled=$(seqkit stats -T "${rfam_tmp}" | tail -1 | cut -f4)
    cat "${rfam_tmp}" >> "${RFAM_OUTPUT}"
    rm -f "${rfam_tmp}"

    rfam_total_seqs=$(( rfam_total_seqs + n_total ))
    echo "  ${rfam_id} (${family_name}): ${n_total} total, ${min_len}-${max_len} bp -> ${n_sampled} sampled"
    rfam_html_rows="${rfam_html_rows}      <tr><td>${family_name//_/ }</td><td>${rfam_id}</td><td>${n_total}</td><td>${min_len}-${max_len}</td><td>${n_sampled}</td></tr>\n"
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

T2T_ACCESSION="${T2T_ACCESSION}" \
T2T_GCF_ACCESSION="${T2T_GCF_ACCESSION}" \
T2T_NAME="${T2T_NAME}" \
T2T_VERSION="${T2T_VERSION}" \
ISS_MODEL="${ISS_MODEL}" \
RAND_SEED="${RAND_SEED}" \
N_T2T="${n_t2t}" \
N_LOCI="${n_loci}" \
MASKED_BP="${masked_bp}" \
N_RFAM_READS="${N_RFAM_READS}" \
N_RFAM="${n_rfam}" \
RFAM_TOTAL_SEQS="${rfam_total_seqs}" \
OUTPUT_HTML="${OUTPUT_HTML}" \
python3 - "${rfam_html_rows}" <<'PYEOF'
import sys, os, datetime

rfam_rows_raw = sys.argv[1]

e = os.environ

html = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Non-rRNA Test Set Summary</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 40px; color: #333; }}
  h1   {{ color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 8px; }}
  h2   {{ color: #34495e; margin-top: 32px; }}
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
<h1>Non-rRNA Test Set Summary</h1>
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
    <tr><td>Reads generated</td><td>{n_t2t}</td></tr>
    <tr><td>Output file</td><td>non_rRNA_test_1M_T2T.fasta</td></tr>
  </tbody>
</table></div>
</section>

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
      <th>Length range (bp)</th><th>Sequences sampled</th>
    </tr>
  </thead>
  <tbody>
{rfam_rows}    <tr style="font-weight:bold; border-top: 2px solid #2c3e50;"><td>Total</td><td></td><td>{rfam_total_seqs}</td><td></td><td>{n_rfam}</td></tr>
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
    n_t2t=e['N_T2T'],
    n_rfam_reads=e['N_RFAM_READS'],
    rfam_rows=rfam_rows_raw.replace('\\n', '\n'),
    rfam_total_seqs=e['RFAM_TOTAL_SEQS'],
    n_rfam=e['N_RFAM'],
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
