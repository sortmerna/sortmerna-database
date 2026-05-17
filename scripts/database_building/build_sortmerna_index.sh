#!/bin/bash
#
# build_sortmerna_index.sh - Assemble and index SortMeRNA database configurations
#
# Concatenates per-domain clustered FASTA files into three database configurations
# and builds a SortMeRNA index for each:
#
#   smr_v4.7_sensitive_db - all SILVA at 97%, RFAM full at 97%
#   smr_v4.7_default_db   - SILVA at 95% (bacteria SSU at 90%), RFAM seed
#   smr_v4.7_fast_db      - SILVA at 90% (bacteria SSU at 85%), RFAM seed

set -euo pipefail

FORCE=false
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
  --force) FORCE=true; shift ;;
  *) POSITIONAL+=("$1"); shift ;;
  esac
done

INPUT_DIR="${POSITIONAL[0]:-data}"
OUTPUT_DIR="${POSITIONAL[1]:-data/index}"
THREADS="${POSITIONAL[2]:-4}"

CLUSTERED_DIR="${CLUSTERED_DIR:-${INPUT_DIR}/clustered}"
VERIFIED_RFAM_DIR="${VERIFIED_RFAM_DIR:-${INPUT_DIR}/verified_rfam}"

echo "============================================"
echo "SortMeRNA Index Building Script"
echo "Clustered dir:    ${CLUSTERED_DIR}"
echo "Verified RFAM dir: ${VERIFIED_RFAM_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Threads:          ${THREADS}"
echo "Force rebuild:    ${FORCE}"
echo "============================================"

if ! command -v sortmerna &>/dev/null; then
  echo "Error: SortMeRNA not found. Please install SortMeRNA first."
  exit 1
fi

if ! command -v seqkit &>/dev/null; then
  echo "Error: seqkit not found. Please install seqkit first."
  exit 1
fi

echo "Using SortMeRNA version: $(sortmerna --version 2>&1 | head -1)"

mkdir -p "${OUTPUT_DIR}"

# Concatenate a list of FASTA files into a combined reference and build the index.
# Args: config_name file1 file2 ...
build_config() {
  local name="$1"
  shift
  local files=("$@")

  local db_dir="${OUTPUT_DIR}/${name}"
  local combined="${db_dir}/${name}.fasta"
  local stats_file="${db_dir}/index.stats"

  mkdir -p "${db_dir}"

  echo ""
  echo "============================================"
  echo "Configuration: ${name}"
  echo "============================================"

  if [[ "${FORCE}" == false ]] && [[ -f "${stats_file}" ]]; then
    echo "  Index already exists - skipping (use --force to rebuild)"
    return 0
  fi

  # Assemble combined reference
  > "${combined}"
  local total=0
  for f in "${files[@]}"; do
    if [[ ! -f "${f}" ]] || [[ ! -s "${f}" ]]; then
      echo "  WARNING: missing or empty file: ${f}"
      continue
    fi
    local n
    n=$(seqkit stats -T "${f}" | tail -1 | cut -f4)
    echo "  + $(basename "${f}"): ${n} sequences"
    cat "${f}" >> "${combined}"
    total=$((total + n))
  done
  echo "  Total: ${total} sequences"

  # Build index
  echo "  Building SortMeRNA index..."
  local start
  start=$(date +%s)

  sortmerna \
    --ref "${combined}" \
    --idx-dir "${db_dir}/idx" \
    --index 1 \
    --threads "${THREADS}"

  local duration=$(( $(date +%s) - start ))
  local index_size
  index_size=$(du -sh "${db_dir}/idx" | cut -f1)

  cat > "${stats_file}" <<EOF
config:          ${name}
combined_fasta:  ${combined}
total_sequences: ${total}
build_time_sec:  ${duration}
index_size:      ${index_size}
threads:         ${THREADS}
build_date:      $(date -Iseconds)
sortmerna:       $(sortmerna --version 2>&1 | head -1)
EOF

  echo "  Done in ${duration}s - index size: ${index_size}"
  echo "  Index: ${db_dir}/idx"
}

# smr_v4.7_sensitive_db: all SILVA at 97%, RFAM full at 97%
build_config "smr_v4.7_sensitive_db" \
  "${CLUSTERED_DIR}/silva_ssu_bacteria_97.fasta" \
  "${CLUSTERED_DIR}/silva_ssu_archaea_97.fasta" \
  "${CLUSTERED_DIR}/silva_ssu_eukaryota_97.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_bacteria_97.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_archaea_97.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_eukaryota_97.fasta" \
  "${CLUSTERED_DIR}/rfam_5s_97.fasta" \
  "${CLUSTERED_DIR}/rfam_5_8s_97.fasta"

# smr_v4.7_default_db: SILVA at 95% (bacteria SSU at 90%), RFAM seed
build_config "smr_v4.7_default_db" \
  "${CLUSTERED_DIR}/silva_ssu_bacteria_90.fasta" \
  "${CLUSTERED_DIR}/silva_ssu_archaea_95.fasta" \
  "${CLUSTERED_DIR}/silva_ssu_eukaryota_95.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_bacteria_95.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_archaea_95.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_eukaryota_95.fasta" \
  "${VERIFIED_RFAM_DIR}/verified_5s_seed.fasta" \
  "${VERIFIED_RFAM_DIR}/verified_5.8s_seed.fasta"

# smr_v4.7_fast_db: SILVA at 90% (bacteria SSU at 85%), RFAM seed
build_config "smr_v4.7_fast_db" \
  "${CLUSTERED_DIR}/silva_ssu_bacteria_85.fasta" \
  "${CLUSTERED_DIR}/silva_ssu_archaea_90.fasta" \
  "${CLUSTERED_DIR}/silva_ssu_eukaryota_90.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_bacteria_90.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_archaea_90.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_eukaryota_90.fasta" \
  "${VERIFIED_RFAM_DIR}/verified_5s_seed.fasta" \
  "${VERIFIED_RFAM_DIR}/verified_5.8s_seed.fasta"

echo ""
echo "============================================"
echo "All indices built in: ${OUTPUT_DIR}"
echo "============================================"
echo ""
echo "Summary:"
printf "  %-30s  %10s  %10s  %s\n" "Configuration" "Sequences" "Build time" "Index size"
for name in smr_v4.7_sensitive_db smr_v4.7_default_db smr_v4.7_fast_db; do
  stats="${OUTPUT_DIR}/${name}/index.stats"
  if [[ -f "${stats}" ]]; then
    seqs=$(grep "total_sequences" "${stats}" | awk '{print $2}')
    secs=$(grep "build_time_sec"  "${stats}" | awk '{print $2}')
    size=$(grep "index_size"      "${stats}" | awk '{print $2}')
    printf "  %-30s  %10s  %9ss  %s\n" "${name}" "${seqs}" "${secs}" "${size}"
  fi
done
echo ""
echo "Next step: Run bash $SMR_DB_ROOT_DIR/scripts/read_simulation/simulate_rrna_reads.sh"
