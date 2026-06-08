#!/bin/bash
#
# build_sortmerna_index.sh - Assemble and index SortMeRNA database configurations
#
# Concatenates per-domain clustered FASTA files into three database configurations
# and builds a SortMeRNA index for each:
#
#   smr_v<SMR_VERSION>_sensitive_db - all SILVA at 97%, Rfam full at 97%
#   smr_v<SMR_VERSION>_default_db   - SILVA at 95% (bacteria SSU at 90%), Rfam full at 90%
#   smr_v<SMR_VERSION>_fast_db      - SILVA at 90% (bacteria SSU at 85%), Rfam full at 85%

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

CLUSTERED_DIR="${CLUSTERED_DIR:-${INPUT_DIR}/clustered}"
VERIFIED_RFAM_DIR="${VERIFIED_RFAM_DIR:-${INPUT_DIR}/verified_rfam}"

echo "============================================"
echo "SortMeRNA Index Building Script"
echo "Clustered dir:    ${CLUSTERED_DIR}"
echo "Verified Rfam dir: ${VERIFIED_RFAM_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Force rebuild:    ${FORCE}"
echo "============================================"

SMR_BIN="${SMR_BIN:?Please set SMR_BIN in your environment (see README Set paths section)}"

if [[ ! -x "${SMR_BIN}" ]]; then
  echo "Error: SortMeRNA not found or not executable: ${SMR_BIN}"
  exit 1
fi

if ! command -v seqkit &>/dev/null; then
  echo "Error: seqkit not found. Please install seqkit first."
  exit 1
fi

SMR_VERSION="${SMR_VERSION:?Please set SMR_VERSION in your environment (see README Set paths section)}"
ACTUAL_VERSION=$("${SMR_BIN}" --version 2>&1 | grep "^SortMeRNA version" | awk '{print $3}')
if [[ "${SMR_VERSION}" != "${ACTUAL_VERSION}" ]]; then
  echo "Error: SMR_VERSION=${SMR_VERSION} does not match ${SMR_BIN} version ${ACTUAL_VERSION}"
  exit 1
fi
echo "Using SortMeRNA version: ${SMR_VERSION}"
SMR_PREFIX="smr_v${SMR_VERSION}"

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

  # Always recompute masking stats (fast, no index rebuild needed)
  local masking_tsv="${db_dir}/masking_stats.tsv"
  rm -f "${masking_tsv}"
  local dust_flag=()
  command -v dustmasker &>/dev/null && dust_flag=(--run-dustmasker)
  for f in "${files[@]}"; do
    [[ -f "${f}" ]] || continue
    local subunit domain base
    base=$(basename "${f}" | sed 's/_[0-9]*\.fasta$//')
    subunit=$(echo "${base}" | sed 's/silva_\(ssu\|lsu\)_.*/\1/' | tr '[:lower:]' '[:upper:]')
    domain=$(echo "${base}"  | sed 's/silva_\(ssu\|lsu\)_//' | sed 's/rfam_.*/rfam/')
    [[ "${base}" == rfam_5_8s* ]] && subunit="Rfam 5.8S" && domain="all"
    [[ "${base}" == rfam_5s*    ]] && subunit="Rfam 5S"   && domain="all"
    python3 "${UTILS_DIR}/compute_masking_stats.py" \
      --fasta "${f}" --subunit "${subunit}" --domain "${domain}" \
      --out "${masking_tsv}" "${dust_flag[@]}"
  done

  if [[ "${FORCE}" == false ]] && [[ -f "${stats_file}" ]]; then
    echo "  Index already exists - skipping (use --force to rebuild)"
    return 0
  fi

  [[ "${FORCE}" == true ]] && rm -rf "${db_dir}/idx"

  # Assemble combined reference and family map
  local family_map="${db_dir}/family_map.tsv"
  > "${combined}"
  local total=0
  local fam_args=()
  for f in "${files[@]}"; do
    if [[ ! -f "${f}" ]] || [[ ! -s "${f}" ]]; then
      echo "  WARNING: missing or empty file: ${f}"
      continue
    fi
    local n family_label
    n=$(seqkit stats -T "${f}" | tail -1 | cut -f4)
    family_label=$(basename "${f}" | sed 's/_[0-9]*\.fasta$//' | sed 's/_/ /g')
    echo "  + $(basename "${f}"): ${n} sequences"
    cat "${f}" >> "${combined}"
    fam_args+=("${f}:${family_label}")
    total=$((total + n))
  done

  python3 "${UTILS_DIR}/build_family_map.py" "${family_map}" "${fam_args[@]}"
  echo "  Total: ${total} sequences"

  # Build index
  echo "  Building SortMeRNA index..."
  local start
  start=$(date +%s)

  "${SMR_BIN}" \
    --ref "${combined}" \
    --workdir "${db_dir}" \
    --idx-dir "${db_dir}/idx" \
    --task 5 &
  local smr_pid=$!

  # Poll the sortmerna process every 5 seconds to track peak CPU% and peak RSS (resident
  # set size in MB - the portion of RAM the process actually holds). The loop exits once
  # the process terminates; wait then captures its exit code so set -e still triggers on failure.
  local peak_cpu=0
  local peak_rss_mb=0
  while kill -0 "${smr_pid}" 2>/dev/null; do
    local cpu rss
    read -r cpu rss < <(ps -p "${smr_pid}" -o %cpu,rss --no-headers 2>/dev/null || echo "0 0")
    peak_cpu=$(awk -v a="${cpu}" -v b="${peak_cpu}" 'BEGIN{print (a+0>b+0)?a:b}')
    local rss_mb=$(( (rss + 0) / 1024 ))
    (( rss_mb > peak_rss_mb )) && peak_rss_mb=${rss_mb}
    sleep 5
  done
  wait "${smr_pid}" || { echo "  ERROR: sortmerna failed"; return 1; }

  local duration=$(( $(date +%s) - start ))
  local index_size
  index_size=$(du -sh "${db_dir}/idx" | cut -f1)

  cat > "${stats_file}" <<EOF
config:          ${name}
combined_fasta:  ${combined}
total_sequences: ${total}
build_time_sec:  ${duration}
index_size:      ${index_size}
peak_cpu_pct:    ${peak_cpu}
peak_rss_mb:     ${peak_rss_mb}
build_date:      $(date -Iseconds)
sortmerna:       ${SMR_VERSION}
EOF

  echo "  Done in ${duration}s - index size: ${index_size} - peak CPU: ${peak_cpu}% - peak RSS: ${peak_rss_mb} MB"
  echo "  Index: ${db_dir}/idx"
}

# sensitive: all SILVA at 97%, Rfam full at 97%
build_config "${SMR_PREFIX}_sensitive_db" \
  "${CLUSTERED_DIR}/silva_ssu_bacteria_97.fasta" \
  "${CLUSTERED_DIR}/silva_ssu_archaea_97.fasta" \
  "${CLUSTERED_DIR}/silva_ssu_eukaryota_97.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_bacteria_97.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_archaea_97.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_eukaryota_97.fasta" \
  "${CLUSTERED_DIR}/rfam_5s_97.fasta" \
  "${CLUSTERED_DIR}/rfam_5_8s_97.fasta"

# default: SILVA at 95% (bacteria SSU at 90%), Rfam full at 90%
build_config "${SMR_PREFIX}_default_db" \
  "${CLUSTERED_DIR}/silva_ssu_bacteria_90.fasta" \
  "${CLUSTERED_DIR}/silva_ssu_archaea_95.fasta" \
  "${CLUSTERED_DIR}/silva_ssu_eukaryota_95.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_bacteria_95.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_archaea_95.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_eukaryota_95.fasta" \
  "${CLUSTERED_DIR}/rfam_5s_90.fasta" \
  "${CLUSTERED_DIR}/rfam_5_8s_90.fasta"

# fast: SILVA at 90% (bacteria SSU at 85%), Rfam full at 85%
build_config "${SMR_PREFIX}_fast_db" \
  "${CLUSTERED_DIR}/silva_ssu_bacteria_85.fasta" \
  "${CLUSTERED_DIR}/silva_ssu_archaea_90.fasta" \
  "${CLUSTERED_DIR}/silva_ssu_eukaryota_90.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_bacteria_90.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_archaea_90.fasta" \
  "${CLUSTERED_DIR}/silva_lsu_eukaryota_90.fasta" \
  "${CLUSTERED_DIR}/rfam_5s_85.fasta" \
  "${CLUSTERED_DIR}/rfam_5_8s_85.fasta"

echo ""
echo "============================================"
echo "All indices built in: ${OUTPUT_DIR}"
echo "============================================"
echo ""
echo "Summary:"
printf "  %-30s  %10s  %10s  %10s  %11s  %s\n" "Configuration" "Sequences" "Build time" "Peak CPU%" "Peak RSS MB" "Index size"
for name in "${SMR_PREFIX}_sensitive_db" "${SMR_PREFIX}_default_db" "${SMR_PREFIX}_fast_db"; do
  stats="${OUTPUT_DIR}/${name}/index.stats"
  if [[ -f "${stats}" ]]; then
    seqs=$(grep "total_sequences" "${stats}" | awk '{print $2}')
    secs=$(grep "build_time_sec"  "${stats}" | awk '{print $2}')
    size=$(grep "index_size"      "${stats}" | awk '{print $2}')
    cpu=$(grep  "peak_cpu_pct"    "${stats}" | awk '{print $2}')
    rss=$(grep  "peak_rss_mb"     "${stats}" | awk '{print $2}')
    printf "  %-30s  %10s  %9ss  %10s  %11s  %s\n" "${name}" "${seqs}" "${secs}" "${cpu}" "${rss}" "${size}"
  fi
done

# Create index_build_summary.html as a copy of clustering_summary.html with a new
# "SortMeRNA Index Build Summary" section injected after "Recommended Database Configurations".
CLUSTER_HTML="${CLUSTERED_DIR}/clustering_summary.html"
INDEX_HTML="${OUTPUT_DIR}/index_build_summary.html"
echo ""
echo "Writing index build summary: ${INDEX_HTML}"

rows=""
for name in "${SMR_PREFIX}_sensitive_db" "${SMR_PREFIX}_default_db" "${SMR_PREFIX}_fast_db"; do
  stats="${OUTPUT_DIR}/${name}/index.stats"
  [[ -f "${stats}" ]] || continue
  s_seqs=$(grep "total_sequences" "${stats}" | awk '{print $2}')
  s_secs=$(grep "build_time_sec"  "${stats}" | awk '{print $2}')
  s_size=$(grep "index_size"      "${stats}" | awk '{print $2}')
  s_cpu=$(grep  "peak_cpu_pct"    "${stats}" | awk '{print $2}')
  s_rss=$(grep  "peak_rss_mb"     "${stats}" | awk '{print $2}')
  s_date=$(grep "build_date"      "${stats}" | awk '{print $2}')
  rows="${rows}      <tr><td>${name}</td><td>${s_seqs}</td><td>${s_secs}</td><td>${s_size}</td><td>${s_cpu}</td><td>${s_rss}</td><td>${s_date}</td></tr>\n"
done

CLUSTER_HTML="${CLUSTER_HTML}" INDEX_HTML="${INDEX_HTML}" \
OUTPUT_DIR="${OUTPUT_DIR}" SMR_PREFIX="${SMR_PREFIX}" \
python3 - "${rows}" <<'PYEOF'
import sys, os, csv
from pathlib import Path

rows_raw = sys.argv[1]
output_dir = os.environ['OUTPUT_DIR']
smr_prefix = os.environ['SMR_PREFIX']

def load_masking(config_name):
    p = Path(output_dir) / config_name / "masking_stats.tsv"
    if not p.exists():
        return []
    with open(p) as f:
        return list(csv.DictReader(f, delimiter='\t'))

def masking_table(configs, col_key, col_label):
    """Build an HTML table for one masking source (silva or dust)."""
    # collect all subunit/domain combos
    rows_by_sd = {}
    for cfg in configs:
        for row in load_masking(cfg):
            k = (row['subunit'], row['domain'])
            rows_by_sd.setdefault(k, {})[cfg] = row

    if not rows_by_sd:
        return "<p><em>Masking stats not available.</em></p>"

    masked_key = f"{col_key}_masked"
    pct_key    = f"{col_key}_pct"
    cfg_headers = "".join(f"<th>{c}</th>" for c in configs)
    thead = f"<tr><th>Subunit</th><th>Domain</th>{cfg_headers}</tr>"
    tbody = ""
    for (subunit, domain), by_cfg in sorted(rows_by_sd.items()):
        cells = ""
        for cfg in configs:
            r = by_cfg.get(cfg, {})
            if r:
                m, p = r.get(masked_key, 'NA'), r.get(pct_key, 'NA')
                t = r.get('total_seqs', '')
                cells += f"<td>{m}/{t} ({p}%)</td>"
            else:
                cells += "<td>-</td>"
        tbody += f"<tr><td>{subunit}</td><td>{domain}</td>{cells}</tr>\n"

    return (
        f"<table><thead>{thead}</thead><tbody>{tbody}</tbody></table>"
    )

configs = [f"{smr_prefix}_sensitive_db", f"{smr_prefix}_default_db", f"{smr_prefix}_fast_db"]
silva_table_html  = masking_table(configs, "silva", "SILVA soft masking")
dust_table_html   = masking_table(configs, "dust",  "dustmasker")
has_dust = any(
    load_masking(c) and any(r.get('dust_masked', 'NA') != 'NA' for r in load_masking(c))
    for c in configs
)

masking_section = (
    "<h2>Soft-masking in SortMeRNA Reference Sequences</h2>\n"
    "<div class=\"description\">\n"
    "<p>SILVA sequences carry soft-masked (lowercase) bases marking low-complexity regions "
    "such as tandem repeats. SortMeRNA's indexer does not distinguish case, so these regions "
    "are indexed identically to unmasked bases. To prevent tandem-repeat artifacts from "
    "generating spurious seed matches, soft-masked bases are converted to hard-masked N before "
    "indexing. N bases cannot participate in 18-mer seeds, so masked regions produce no hits.</p>\n"
    "<p><strong>Practical consequence:</strong> sequences with N's in the reference generate "
    "fewer seed hits than they should (only matching reads that happen to have the same N at "
    "those positions), but any alignment that does form is scored correctly against the "
    "unmasked positions. For rRNA databases this is generally fine since reference sequences "
    "are high quality with very few masked bases. For reads with many N's (e.g. low-quality "
    "regions of PacBio reads pre-CCS), the seed stage may miss some true matches that "
    "Smith-Waterman would have accepted.</p>\n"
    "</div>\n"
    "<h3>Table 1: SILVA soft masking (sequences with at least one lowercase base)</h3>\n"
    f"<div class=\"table-wrap\">{silva_table_html}</div>\n"
    "<h3>Table 2: Independent dustmasker low-complexity masking</h3>\n"
    "<div class=\"description\"><p>dustmasker is run independently of SILVA to provide a "
    "second estimate and guard against SILVA changing or disabling its masking in future "
    "releases.</p></div>\n"
    + (f"<div class=\"table-wrap\">{dust_table_html}</div>\n"
       if has_dust else
       "<p><em>dustmasker not available during this build.</em></p>\n")
)

section = (
    "<h2>SortMeRNA Index Build Summary</h2>\n"
    "<div class=\"description\"><p>"
    "Resource usage measured by polling <code>ps</code> every 5 seconds during index construction. "
    "Peak CPU% is the highest single-sample value; Peak RSS MB is the highest resident set size "
    "(RAM actually held by the process, excluding swap) observed during the build."
    "</p></div>\n"
    "<div class=\"table-wrap\"><table>\n"
    "  <thead>\n"
    "    <tr>\n"
    "      <th>Configuration</th><th>Sequences</th><th>Build time (s)</th><th>Index size</th>\n"
    "      <th>Peak CPU%</th><th>Peak RSS MB</th><th>Build date</th>\n"
    "    </tr>\n"
    "  </thead>\n"
    "  <tbody>\n"
    + rows_raw.replace('\\n', '\n')
    + "  </tbody>\n</table></div>\n"
    "<p class=\"footer\">Generated by build_sortmerna_index.sh - SortMeRNA database pipeline</p>"
)

cluster_html = os.environ['CLUSTER_HTML']
index_html   = os.environ['INDEX_HTML']

html = open(cluster_html).read() if os.path.isfile(cluster_html) else '<html><body></body></html>'

marker = 'Recommended Database Configurations'
pos = html.find(marker)
if pos != -1:
    close = html.find('</section>', pos)
    insert_at = close + len('</section>') if close != -1 else len(html)
else:
    insert_at = html.find('</body>')
    if insert_at == -1:
        insert_at = len(html)

html = html[:insert_at] + '\n<section>\n' + masking_section + '\n</section>\n<section>\n' + section + '\n</section>\n' + html[insert_at:]
open(index_html, 'w').write(html)
PYEOF

echo "  Done - ${INDEX_HTML}"

echo ""
echo "Next step: Run bash $SMR_DB_ROOT_DIR/scripts/read_simulation/simulate_rrna_reads.sh"
