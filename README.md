# SortMeRNA Database Benchmark

[![CI](https://github.com/sortmerna/sortmerna-database/actions/workflows/tests.yml/badge.svg)](https://github.com/sortmerna/sortmerna-database/actions/workflows/tests.yml)
[![Coverage](https://codecov.io/gh/sortmerna/sortmerna-database/branch/main/graph/badge.svg)](https://codecov.io/gh/sortmerna/sortmerna-database)

A comprehensive benchmarking framework for building and evaluating optimized rRNA databases for SortMeRNA across modern sequencing platforms.

**Status**: 🚧 Active Development - this repository is undergoing daily changes and is not yet ready for use. Scripts, workflows, and outputs may change without notice. Do not use this repository in production or depend on its outputs until this status is updated to indicate a stable release.

## Overview

This repository contains code and workflows to:
1. Download and independently verify rRNA sequences using Infernal covariance models
2. Build clustered rRNA databases from verified SILVA and Rfam sequences
3. Benchmark database performance (accuracy vs. size tradeoffs)
4. Optimize clustering parameters for different use cases
5. Validate databases against simulated and real datasets
   
## Goals

### Phase 1: Database Construction
- [x] Download latest SILVA database
- [x] Download Rfam rRNA families
- [x] Download Rfam covariance models (Infernal)
- [x] Independently verify rRNA sequences with Infernal cmsearch (--cut_ga)
- [x] Implement clustering pipeline
- [x] Test multiple clustering thresholds
- [x] Build SortMeRNA indices for each clustered database
- [x] Generate database statistics and metadata

### Phase 2: Validation and Benchmarking
- [x] Simulate Illumina 150bp rRNA using non-seed cluster members and non-rRNA reads using T2T genome + Rfam non-rRNA families (InSilicoSeq)
- [x] Experiment 1: Scalability - run SortMeRNA on T2T non-rRNA and rRNA reads at 10K, 100K, 1M, 10M scale points; generate runtime, FP rate, sensitivity, and E-value plots
- [x] Experiment 2: Simulate rRNA reads from non-seed members at each clustering threshold (simulate_rrna_reads.sh)
- [x] Experiment 2: Run SortMeRNA against matched database configuration and measure sensitivity
- [ ] Download real Illumina metatranscriptomics data
- [ ] Download real PacBio amplicon data (Karst et al. 2021)
- [ ] Download real PacBio metagenomics data
- [ ] Benchmark latest SortMeRNA vs. SortMeRNA v2.1b (original paper) and other tools: sensitivity, specificity, runtime, memory


## Installation

### Requirements

**Core tools:**
- [SortMeRNA](https://github.com/sortmerna/sortmerna) v6.0.1 - rRNA filtering (installed from source, see Quick Start)
- [VSEARCH](https://github.com/torognes/vsearch) >= 2.22 - sequence clustering
- [Infernal](http://eddylab.org/infernal/) >= 1.1.4 - covariance model search (cmsearch, cmpress)
- [SeqKit](https://bioinf.shenwei.me/seqkit/) >= 2.5 - sequence statistics

**Languages:**
- Python >= 3.8
- R >= 4.0

**Python libraries:**
- pandas, numpy, matplotlib, seaborn, biopython

**Simulation tools (for benchmarking):**
- [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) - Illumina read simulation
- [PBSIM3](https://github.com/yukiteruono/pbsim3) - PacBio and Oxford Nanopore read simulation

### Quick Start

**Install SortMeRNA from source**

```bash
export SMR_VERSION=6.0.2

# Download and extract
wget https://github.com/sortmerna/sortmerna/archive/refs/tags/v${SMR_VERSION}.tar.gz
tar xvf v${SMR_VERSION}.tar.gz && cd sortmerna-${SMR_VERSION}

# Pin GCC 11
sudo apt install -y gcc-11 g++-11
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 100
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 100

# Create a build environment and install Python build dependencies
conda create -n sortmerna_${SMR_VERSION} -c conda-forge -y python pyyaml jinja2 requests ninja cmake
conda activate sortmerna_${SMR_VERSION}

python setup.py all

conda deactivate
```

The binary is at `dist/bin/sortmerna` inside the source directory.

**Install benchmarking dependencies:**

```bash
# Clone repository
git clone https://github.com/sortmerna/sortmerna-database.git
cd sortmerna-database

# Create conda environment
conda env create -f environment.yml
conda activate sortmerna-bench
```

## Phase 1: Database Construction

Download and verify rRNA sequences from SILVA and Rfam, cluster them at multiple identity thresholds, and build SortMeRNA indices for three database configurations (sensitive, default, fast). This phase produces the reference databases used in Phase 2.

### Database Sources

#### SILVA rRNA Database
- **Version**: Release 138.2 Ref NR 99
- **URL**: https://www.arb-silva.de/
- **Content**: SSU (16S/18S) and LSU (23S/28S) rRNA sequences
- **Taxonomy**: Bacteria, Archaea, Eukarya
- **Size (raw)**: 510,495 sequences (SSU) and 95,279 sequences (LSU)
- **File type**: `_trunc` variants - sequences have been truncated so that all nucleotides not aligned to the SILVA reference alignment are removed. This produces clean, alignment-bounded rRNA sequences and avoids including flanking genomic context that would inflate database size and reduce SortMeRNA specificity.

#### Rfam
- **Version**: Rfam 15.1 (January 2026, 4227 families)
- **URL**: https://rfam.org/
- **Families**: RF00001 (5S), RF00002 (5.8S)
- **Size (raw)**: 712 seed alignment sequences (5S) and 61 seed alignment sequences (5.8S)

### Set paths & create working directory

```bash
# SortMeRNA version and binary - update when upgrading
export SMR_VERSION=6.0.2
export SMR_BIN=$HOME/sortmerna-${SMR_VERSION}/dist/bin/sortmerna

# Path to the cloned sortmerna-database repository
export SMR_DB_ROOT_DIR=$HOME/sortmerna-database
export UTILS_DIR=$SMR_DB_ROOT_DIR/scripts/utils

# Working directory - all data will be written here
export WORK_DIR=$HOME/working
export DATA_DIR=$WORK_DIR/data
mkdir -p $DATA_DIR && cd $WORK_DIR

# Data directories (adjust if needed)
export SILVA_DIR=$DATA_DIR/silva
export RFAM_DIR=$DATA_DIR/rfam
export CMS_DIR=$DATA_DIR/cms
export VERIFIED_DIR=$DATA_DIR/verified
export VERIFIED_RFAM_DIR=$DATA_DIR/verified_rfam
export CLUSTERED_DIR=$DATA_DIR/clustered
export INDEX_DIR=$DATA_DIR/index
export NON_RRNA_DIR=$DATA_DIR/non_rrna
export RRNA_SIM_DIR=$DATA_DIR/rrna_sim
export SCALABILITY_DIR=$DATA_DIR/scalability_test

# SILVA versions and full download URLs (update to use a different release or file type)
export SILVA_SSU_VERSION=138.2
export SILVA_SSU_PATH=https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz
export SILVA_LSU_VERSION=138.2
export SILVA_LSU_PATH=https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_LSURef_NR99_tax_silva_trunc.fasta.gz
export RFAM_VERSION=15.1

# Human T2T genome - used as non-rRNA specificity test source
export T2T_ACCESSION=GCA_009914755.4
# RefSeq (GCF) accession - used to download rRNA annotation GFF3
export T2T_GCF_ACCESSION=GCF_009914755.1
# NCBI assembly name used in FTP filenames (${T2T_ACCESSION}_${T2T_NAME}_genomic.fna.gz)
export T2T_NAME=T2T-CHM13v2.0
# human-readable version used for output filenames
export T2T_VERSION=chm13v2.0
# full FTP directory URL for GCA assembly (genome FASTA)
export T2T_BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/${T2T_ACCESSION}_${T2T_NAME}
# full FTP directory URL for GCF assembly (rRNA annotation GFF3)
export T2T_GCF_BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/${T2T_GCF_ACCESSION}_${T2T_NAME}

# Rfam non-rRNA families FTP (tied to RFAM_VERSION; change to CURRENT to always pull the latest)
export RFAM_NON_RRNA_FTP=https://ftp.ebi.ac.uk/pub/databases/Rfam/$RFAM_VERSION/fasta_files

```

### Download Source Databases

```bash
# Download SILVA
bash $SMR_DB_ROOT_DIR/scripts/database_building/download_silva.sh $SILVA_DIR

# Download Rfam
bash $SMR_DB_ROOT_DIR/scripts/database_building/download_rfam.sh $RFAM_DIR

# Download and press Rfam covariance models
bash $SMR_DB_ROOT_DIR/scripts/database_building/download_cms.sh $CMS_DIR
```

### Clustering Strategy

#### Reference Verification

Before clustering, all sequences were independently verified as rRNA using Infernal `cmsearch --cut_ga` against Rfam covariance models (see Verify Sequences for details). Models used per gene/domain:

- **SSU**: RF00177 (Bacteria), RF01959 (Archaea), RF01960 (Eukaryota)
- **LSU**: RF02541 (Bacteria), RF02540 (Archaea), RF02543 (Eukaryota)
- **5.8S**: RF00002 (Eukaryota)
- **5S**: RF00001 (all domains)

Trimming to the cmsearch hit coordinates was strategic for SortMeRNA: the tool builds a k-mer index over every reference sequence, so any non-rRNA nucleotides in a reference - flanking genomic DNA, assembly context, or sequence that SILVA's own truncation missed - would be indexed alongside the rRNA and could produce false-positive matches against non-rRNA reads.

#### Tools Considered

[VSEARCH](https://github.com/torognes/vsearch) was used for sequence clustering. VSEARCH is a high-performance open-source tool for sequence analysis that implements a centroid-based clustering algorithm. It is the same tool used by the SILVA project to generate their non-redundant (NR) databases ([SILVA 138.2 release notes](https://www.arb-silva.de/documentation/release-1382)), ensuring compatibility and reproducibility with established rRNA reference workflows.

Key features:
- Fast clustering with `--cluster_fast` using identity thresholds
- Supports both strands for rRNA sequences
- Memory-efficient for large databases
- Produces representative centroid sequences

#### Clustering Thresholds to Test

- **97%** - Species-level clustering (commonly used for 16S)
- **95%** - Genus-level approximation
- **90%** - Family-level approximation
- **85%** - Higher-order taxonomic grouping

### Verify Sequences

All sequences were independently verified as rRNA using Infernal's `cmsearch --cut_ga`. SILVA SSU/LSU sequences were screened with `--hmmonly` (profile-HMM only, no secondary-structure scoring), which is sufficient for long ribosomal RNAs and substantially faster; Rfam 5S and 5.8S sequences used a full covariance-model search because these shorter RNAs benefit more from secondary-structure information. Each verified sequence was trimmed to the cmsearch hit coordinates `[seq_from, seq_to]`; sequences with no above-threshold hit were written to `flagged_*.fasta` and excluded from downstream steps.

**Multi-hit merging.** Some sequences produce more than one `!`-threshold hit against the same covariance model, typically because cmsearch reports separate local alignments for different segments of the same rRNA gene. When this happens, the trimming window is expanded to `[min(seq_from), max(seq_to)]` across all hits - but only if that merged span covers at least **85%** of the original sequence length (`COVERAGE_THRESHOLD = 0.85`). If coverage is below this threshold the hits are non-contiguous enough to suggest a chimera, assembly artifact, or non-rRNA insertion between the matching regions; the sequence is instead trimmed to the single highest-scoring hit and its log entry is marked `_REVIEW`. The merged span intentionally includes any gap between hits: for rRNA, inter-hit gaps typically correspond to hypervariable regions that score poorly against the profile HMM but are still genuine rRNA sequence. Including a gap is therefore usually correct, while excluding it risks discarding valid sequence; the 85% threshold provides a safeguard against merging across large non-rRNA insertions.

**Minimum trimmed length.** After trimming, sequences shorter than **50%** of the original SILVA sequence length are routed to `flagged_*.fasta` rather than `clean`. This catches cases where only a small fragment of the rRNA was recovered after trimming - for example, a sequence where cmsearch aligned only to a short terminal region while the bulk of the sequence had no rRNA match. Such fragments are too partial to serve as reliable SortMeRNA reference sequences and could reduce specificity by matching non-rRNA reads that happen to resemble that short region. The note column in the TSV log records `too_short` for these cases.

**LSU Eukaryota is most affected.** Multi-hit and short-fragment cases were disproportionately observed in LSU eukaryotic sequences (RF02543, model length ~3400 bp). Among the ~16,000 LSU eukaryota sequences verified from SILVA 138.2, approximately 1.7% were flagged for low-coverage multi-hit patterns requiring manual review.

```bash
# Verify SILVA SSU and LSU sequences (bacteria, archaea, eukaryota) using --hmmonly --cut_ga
bash $SMR_DB_ROOT_DIR/scripts/database_building/verify_silva.sh \
    $WORK_DIR/data \
    $VERIFIED_DIR \
    4
```

```bash
# Verify Rfam 5S and 5.8S seed sequences using full CM search with --cut_ga
bash $SMR_DB_ROOT_DIR/scripts/database_building/verify_rfam.sh \
    $WORK_DIR/data \
    $VERIFIED_RFAM_DIR \
    4
```

Outputs per domain/type in each verified directory:

| File pattern | Description |
|---|---|
| `verified_<gene>_<domain>.fasta` | Sequences confirmed as rRNA - input to clustering |
| `flagged_<gene>_<domain>.fasta` | Sequences with no qualifying hit - excluded |
| `cmsearch_log_<gene>_<domain>.tsv` | Best hit coordinates and score for each kept sequence |
| `<gene>_<domain>_cmsearch.tblout` | Raw cmsearch output (kept for auditing) |
| `verification_summary.html` | Interactive summary with flagged counts and offset histogram |

**Verification Summary**

- SILVA 138.2: <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/verified/verification_summary.html" target="_blank">verification_summary.html</a>

- Rfam: <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/verified_rfam/verification_summary.html" target="_blank">verification_summary.html</a>

### Build Clustered Databases

```bash
# Cluster verified SILVA and Rfam sequences at multiple identity thresholds (97%, 95%, 90%, 85%)
bash $SMR_DB_ROOT_DIR/scripts/database_building/cluster_sequences.sh \
    $WORK_DIR/data \
    $CLUSTERED_DIR \
    4
```

By default, vsearch clustering is skipped for any threshold where the `.uc` output file already exists - only the downstream steps (member extraction, leakage check, summary table) are re-run. To force vsearch to re-cluster and overwrite existing `.uc` files, pass `--force`:

```bash
bash $SMR_DB_ROOT_DIR/scripts/database_building/cluster_sequences.sh \
    --force $WORK_DIR/data \
    $CLUSTERED_DIR \
    4
```

For each database and clustering threshold, the script outputs four files:

| File | Description |
|------|-------------|
| `*_XX.fasta` | Centroid/seed sequences - the clustered rRNA database used by SortMeRNA |
| `*_XX.uc` | VSEARCH cluster membership file |
| `*_XX_test_members.fasta` | Non-seed cluster members - used as source sequences for simulating test reads |
| `*_XX_cluster_mapping.txt` | Tab-delimited mapping of each member sequence ID to its seed sequence ID |

The centroid sequences (`*_XX.fasta`) become the SortMeRNA reference databases. The non-seed members (`*_XX_test_members.fasta`) are sequences that were clustered away at each threshold, providing a natural source of reads for benchmarking - since they are real rRNA sequences not present in the database, they test whether SortMeRNA can still identify similar but non-identical rRNA.

- Summary of total sequences per clustering threshold and SortMeRNA reference databases: <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/clustered/clustering_summary.html" target="_blank">clustering_summary.html</a>

### Build SortMeRNA Indices

Three database configurations are assembled from the clustered FASTA files and indexed with SortMeRNA:

| Configuration | SILVA SSU bacteria | SILVA SSU other | SILVA LSU | Rfam 5S / 5.8S |
|---|---|---|---|---|
| `smr_v<version>_sensitive_db` | 97% | 97% | 97% | 97% (full) |
| `smr_v<version>_default_db` | 90% | 95% | 95% | seed |
| `smr_v<version>_fast_db` | 85% | 90% | 90% | seed |

The sensitive database maximises recall; the default and fast databases trade a small number of sequences for a smaller index and faster runtime. Bacteria SSU uses a one-step lower threshold than the other SILVA domains in the default and fast configurations because bacterial 16S sequences are the largest and most diverse group - a slightly lower threshold helps limit index size without meaningfully reducing sensitivity.

```bash
bash $SMR_DB_ROOT_DIR/scripts/database_building/build_sortmerna_index.sh \
    $WORK_DIR/data \
    $INDEX_DIR
```

- Per-configuration index build report (sequence count, build time, index size, peak CPU%, peak RAM): <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/index/index_build_summary.html" target="_blank">index_build_summary.html</a>

Pass `--force` to rebuild an index that already exists:

```bash
bash $SMR_DB_ROOT_DIR/scripts/database_building/build_sortmerna_index.sh \
    --force $WORK_DIR/data \
    $INDEX_DIR
```

For each configuration the script produces:

| File / directory | Description |
|---|---|
| `<config>/<config>.fasta` | Combined reference FASTA (all domains concatenated) |
| `<config>/idx/` | SortMeRNA index files |
| `<config>/index.stats` | Sequence count, build time, index size, SortMeRNA version |

**Index building is done once.** The index only needs to be built once per database configuration. When you later run SortMeRNA alignment with the same `--ref` and `--idx-dir` paths, SortMeRNA finds the existing index and skips rebuilding - even across separate invocations or tmux sessions. Use `--force` to explicitly trigger a rebuild.

### Non-rRNA Test Sets

Two separate test sets from different sources, used independently to measure specificity:

| File | Source | Description |
|------|--------|-------------|
| `non_rRNA_test_10M_T2T.fasta` | Human T2T genome (CHM13v2.0) | 10M simulated 150bp PE reads, rRNA loci masked prior to simulation |
| `non_rRNA_test_Rfam.fasta` | Rfam non-rRNA families | 500K sequences sampled evenly across 10 families that share structural features with rRNA (the most challenging specificity test) |

```bash
bash $SMR_DB_ROOT_DIR/scripts/read_simulation/download_non_rrna.sh $NON_RRNA_DIR 4
```

**Preparation (T2T):**
1. Download CHM13v2.0 genome and RefSeq GFF3 annotation
2. Supplement the GFF3 with Infernal cmsearch (RF01960, RF02543, RF00001, RF00002 against the full genome; RF00177 and RF02541 against chrM only) to find rRNA copies not individually annotated in RefSeq (see below)
3. Mask rRNA loci with bedtools maskfasta
4. Simulate 10M 150bp PE reads with InSilicoSeq (NovaSeq error model)

**Supplementing GFF3 with cmsearch:**
The T2T-CHM13v2.0 assembly resolves the Nucleolar Organizer Regions (NORs) on the short arms of acrocentric chromosomes (chr13, chr14, chr15, chr21, chr22), assembling the large tandem arrays of rDNA repeat units encoding 18S/5.8S/28S rRNA that were largely absent or unsequenced in GRCh37/GRCh38. Comparison of cmsearch-identified rDNA loci against the RefSeq GFF3 annotation reveals rDNA sequence not covered by the annotation (exact figures are reported in <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/non_rrna/non_rrna_test_set_summary.html" target="_blank">non_rrna_test_set_summary.html</a> generated by `simulate_non_rrna.sh`), demonstrating that annotation-based masking alone is insufficient for this assembly. To run cmsearch, set `CMS_DIR` to a directory of pressed Rfam covariance models (see `download_non_rrna.sh` for details); without it the script falls back to GFF3-only masking.

BLAST analysis of the RF01960 (18S) regions not covered by GFF3 shows that the large majority are rRNA-derived: the most common top hits are 45S pre-ribosomal RNA entries, followed by 18S ribosomal pseudogene entries. All are included in the masking BED - if a sequence is annotated as a ribosomal pseudogene it is rRNA-derived and SortMeRNA, which operates purely on sequence content, would correctly flag reads from it.

Note: RF00001 (5S rRNA) hits include 5S pseudogenes distributed across many chromosomes in addition to functional copies - BLAST analysis confirms these are RNA5SP* pseudogene entries. These are included in the masking BED for the same reason.

Note: The mitochondrial 12S rRNA (homologue of the SSU rRNA: RF00177) and 16S rRNA (homologue of the LSU rRNA: RF02541) were identified by running cmsearch against chrM only and added to the masking BED. The two mt-rRNA genes are separated by ~70 bp and merge into a single BED region after adding the 100 bp margin. Without this step, reads simulated from mitochondrial rRNA loci would not be masked and could contribute artifactual alignments that inflate apparent false-positive signal.

**Preparation (Rfam):**
1. Download 10 non-rRNA families (tRNA, SRP RNA, tmRNA, RNase P, spliceosomal RNAs)
2. Sample evenly across families

Sequences are used as-is without read simulation. Most Rfam families consist of short RNAs (tRNA averages 73 bp, spliceosomal snRNAs 100-200 bp); simulating 150 bp reads from these sequences would discard the majority of them due to length. The natural length distribution is also the point: these structurally complex RNAs should be rejected by SortMeRNA regardless of length.

Sampling evenly across families (fair-share allocation, ~102K per large family at the 500K default) ensures all families are equally represented regardless of family size. A fixed random seed (`--seed 42`) makes the output reproducible.

```bash
bash $SMR_DB_ROOT_DIR/scripts/read_simulation/simulate_non_rrna.sh \
    $NON_RRNA_DIR \
    4 \
    --t2t-reads 10000000 \
    --rfam-reads 500000
```

- Non-rRNA reference sources summary for read simulation: <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/non_rrna/non_rrna_test_set_summary.html" target="_blank">non_rrna_test_set_summary.html</a>

### rRNA Test Sets

Three Sets of simulated Illumina rRNA reads, one per database configuration. Each Set uses non-seed cluster members at the matching clustering threshold as source sequences, so every simulated read is a real rRNA sequence absent from the database being tested:

| Set | Source members | Matched database | Reads |
|-----|---------------|-----------------|-------|
| Set 1 | 97% non-seeds (all types) | `smr_v${SMR_VERSION}_sensitive_db` | 12,500 per rRNA type x 8 types = 100,000 |
| Set 2 | bacteria SSU 90%, others 95%, Rfam 90% | `smr_v${SMR_VERSION}_default_db` | 12,500 per rRNA type x 8 types = 100,000 |
| Set 3 | bacteria SSU 85%, others 90%, Rfam 85% | `smr_v${SMR_VERSION}_fast_db` | 12,500 per rRNA type x 8 types = 100,000 |

Rfam 5S and 5.8S use 97% non-seed members for all Sets because the default and fast databases use seed-only Rfam (no threshold-based clustering), so 97% members are absent from all three databases.

**Preparation:**
1. Cluster sequences with `cluster_sequences.sh` (produces `*_test_members.fasta` files)
2. Run InSilicoSeq on each rRNA type's non-seed members (150bp PE, NovaSeq model)
3. Pool all types per Set into a single FASTA

```bash
bash $SMR_DB_ROOT_DIR/scripts/read_simulation/simulate_rrna_reads.sh \
    $RRNA_SIM_DIR \
    4 \
    --clustered-dir $CLUSTERED_DIR
```

- rRNA read simulation summary: <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/rrna_sim/rrna_simulation_summary.html" target="_blank">rrna_simulation_summary.html</a>

**Scalability pool (Experiment 1):** `simulate_rrna_reads.sh` also produces `rRNA_test_10M.fasta` alongside the three sensitivity sets. The pool is built from Set 2 non-seed sources (default database, 90-95% clustering thresholds) and combines two source types:
- **SILVA types (6 types, SSU + LSU):** IUPAC-cleaned and run through ISS (NovaSeq model, 150bp PE). ISS generates 10M reads.
- **Rfam 5S (avg 117bp) and 5.8S (avg 150bp):** included directly as-is (IUPAC-cleaned, no min-length filter). These sequences are shorter than the ISS read length and cannot be reliably simulated; including them directly ensures short rRNA - the hardest sequences to detect as `S_min` rises with read count - are represented in the scalability pool.

The SILVA ISS reads and Rfam direct sequences are shuffled together before output.

## Phase 2: Validation and Benchmarking

Evaluate the databases built in Phase 1 for sensitivity (rRNA detection) and specificity (false positive rate), and benchmark the latest SortMeRNA against v2.1b. Simulated Illumina reads from non-seed cluster members test sensitivity at each clustering threshold; non-rRNA reads from the T2T genome and Rfam test specificity. Real PacBio amplicon data provides an independent validation.

### Benchmarking Approach

#### Simulated Data (Illumina)
Generate synthetic Illumina reads with known rRNA/non-rRNA composition:
- **Tool**: InSilicoSeq (150bp paired-end, NovaSeq error model)
- **rRNA source**: Non-seed cluster members (`*_test_members.fasta`) - real rRNA sequences not present in the clustered database
- **Non-rRNA source**: `non_rRNA_test_10M_T2T.fasta` and `non_rRNA_test_Rfam.fasta` - tested separately

#### SortMeRNA is a filter, not a classifier

SortMeRNA is an rRNA **filter** - it answers "is this read rRNA?" not "which type of rRNA is this?". A read simulated from a bacterial SSU sequence may be detected via a k-mer match and Smith-Waterman extension against an archaeal SSU reference, and this is expected behavior given how conserved rRNA sequences are across domains. Per-type sensitivity is therefore not a meaningful metric. All rRNA types are pooled into a single reads file per Set, and sensitivity is measured as a single aggregate: detected / total.

#### E-value Filtering

SortMeRNA uses the Karlin-Altschul framework (`E = K · m · n · exp(-λ · S)`) with Gumbel parameters (λ, K) [6,7] computed via the ALP library, where `m` is the query length, `n` is the total length of the reference database, `S` is the alignment score, and `E` is the expected number of alignments with score >= S by chance. Rather than computing a floating-point E-value per alignment, SortMeRNA inverts the formula once at startup to derive a minimum Smith-Waterman score (`S_min = ln(E / (K·m·n)) / (-λ)`), then filters reads with a single integer comparison during alignment.

**Distinction from BLAST**: in BLAST, `m` is the length of the individual query sequence, giving a per-query E-value. In SortMeRNA, `m` is the total nucleotide count across all reads in the dataset. This means `S_min` sets a run-level threshold: the expected number of spurious alignments across all reads against the database is <= E, not <= E per read. The effective per-read E-value is therefore `E / total_reads` - with 1M reads and E=1.0, the per-read threshold is 0.000001, far stricter than BLAST's default. A consequence is that the threshold is dataset-size dependent: adding more reads makes filtering more stringent, and short reads face the same absolute score threshold as long reads regardless of their length.

This design is also motivated by the difference in database scale. SortMeRNA's rRNA reference databases are ~124M total bases (~69K sequences), roughly 10,000x smaller than BLAST's nt database (~1.3 trillion bases, ~96M sequences) (July 2023) [8]. If SortMeRNA used BLAST's per-query approach with `m` = individual read length (e.g. 150 bp), the search space `K·m·n` would be so small that many alignments would pass the filter. By setting `m` to the total reads length, SortMeRNA trades BLAST's per-query statistical framing for a run-level one in order to produce a meaningful threshold despite having a reference database that is ~10,000x smaller than BLAST's.

**Parallelism and `--score_split`**: SortMeRNA processes reads in parallel by dividing the input file into per-thread byte-range chunks - no physical splitting of the reads file occurs. By default, `m` is the total nucleotide count across all reads regardless of thread count, so `S_min` is identical across threads and the run-level threshold holds. The `--score_split` option changes this: it divides `m` by the number of threads, computing `S_min` as if each thread's chunk were an independent dataset. This lowers `S_min` (more lenient threshold), with an effect equivalent to increasing the E-value. It is off by default and should be used with care, as it makes the threshold dependent on the number of threads rather than the size of the dataset.

#### Experiment 1: Scalability

**Goal:** How do runtime and memory scale with read volume, and do sensitivity and false positive rate remain stable at large scale?

- **Design:** Fixed configuration (SMR default db) across increasing read volumes, run separately for rRNA reads and T2T non-rRNA reads
- **Read volumes:** 10,000 -> 100,000 -> 1,000,000 -> 10,000,000 reads
- **rRNA reads:** Subsampled from Set 2 non-seeds (90-95% identity members: default db) - tests sensitivity and the E-value threshold scaling effect (`S_min` increases with total read count, so sensitivity may shift across scale points)
- **Non-rRNA reads:** Subsampled from `non_rRNA_test_10M_T2T.fasta` (10M reads simulated from the masked T2T genome to cover all scale points without re-running InSilicoSeq) - tests false positive rate at scale. The Rfam non-rRNA sequences were not used here because the T2T genome provides a much larger pool to sample from, and InSilicoSeq can generate consistent 150bp paired-end reads at any volume. Rfam families have limited sequence counts and many members are shorter than 150bp, making it impractical to simulate uniform-length reads at the scales needed for this experiment.
- **Rationale for default db:** The configuration most users would deploy in practice; runtime numbers are directly interpretable for real-world use
- **Metrics:**
  - Wall-clock runtime at each scale point
  - Peak RAM usage at each scale point
  - Sensitivity (rRNA runs) and false positive rate (T2T runs) at each scale point
  - S_min (score threshold) vs. read count - illustrates the E-value scaling effect
  - E-value and % identity distributions of aligned reads at each scale point
- **`--score_split` comparison:** Re-run at each scale point with `--score_split` enabled. This option computes `S_min` from the per-thread chunk size rather than the total dataset size, making the threshold less sensitive to total read count. Comparing the two runs directly shows how much sensitivity and false positive rate shift when the E-value threshold is decoupled from dataset scale.

Run for T2T non-rRNA reads (false positive rate at scale), Rfam non-rRNA reads (hardest specificity test - structurally similar to rRNA), and rRNA reads (sensitivity at scale) across four E-value thresholds. E-value 1 is SortMeRNA's default. Requires `non_rRNA_test_10M_T2T.fasta` and `non_rRNA_test_Rfam.fasta` from `simulate_non_rrna.sh` and `rRNA_test_10M.fasta` from `simulate_rrna_reads.sh`. The Rfam non-rRNA set has 500K reads so its scale points are capped at 500K:

Subsample reads once at E-value 1 (the default), then reuse those reads for the remaining thresholds with `--reads-dir`:

```bash
for ev in 1 0.1 0.05 0.01 0.001 0.0001 0.00001; do
    reads_dir_t2t=$( [[ "${ev}" == "1" ]] && echo "" || echo "--reads-dir $SCALABILITY_DIR/scalability_t2t_ev1" )
    reads_dir_rfam=$( [[ "${ev}" == "1" ]] && echo "" || echo "--reads-dir $SCALABILITY_DIR/scalability_rfam_ev1" )
    reads_dir_rrna=$( [[ "${ev}" == "1" ]] && echo "" || echo "--reads-dir $SCALABILITY_DIR/scalability_rrna_ev1" )

    bash $SMR_DB_ROOT_DIR/scripts/benchmarking/run_scalability.sh \
        $NON_RRNA_DIR/non_rRNA_test_10M_T2T.fasta \
        $SCALABILITY_DIR/scalability_t2t_ev${ev} \
        4 \
        --index-dir $INDEX_DIR \
        --config smr_v${SMR_VERSION}_default_db \
        --evalue ${ev} ${reads_dir_t2t}

    bash $SMR_DB_ROOT_DIR/scripts/benchmarking/run_scalability.sh \
        $NON_RRNA_DIR/non_rRNA_test_Rfam.fasta \
        $SCALABILITY_DIR/scalability_rfam_ev${ev} \
        4 \
        --index-dir $INDEX_DIR \
        --config smr_v${SMR_VERSION}_default_db \
        --scale 10000,100000,500000 \
        --evalue ${ev} ${reads_dir_rfam}

    bash $SMR_DB_ROOT_DIR/scripts/benchmarking/run_scalability.sh \
        $RRNA_SIM_DIR/rRNA_test_10M.fasta \
        $SCALABILITY_DIR/scalability_rrna_ev${ev} \
        4 \
        --index-dir $INDEX_DIR \
        --config smr_v${SMR_VERSION}_default_db \
        --evalue ${ev} ${reads_dir_rrna}
done
```

Generate the ROC plot once all E-value runs are complete. Each point on the curve is one E-value threshold, showing the sensitivity vs. false positive rate tradeoff. T2T and Rfam non-rRNA appear as separate series - both use the same rRNA runs for TPR but differ in what non-rRNA reads are used for FPR:

```bash
python3 $SMR_DB_ROOT_DIR/scripts/utils/plot_roc_evalue.py \
    --output-dir $SCALABILITY_DIR/plots \
    --label scalability_benchmark \
    --evalues 1 0.1 0.05 0.01 0.001 0.0001 0.00001 \
    --rrna-dirs \
        $SCALABILITY_DIR/scalability_rrna_ev1 \
        $SCALABILITY_DIR/scalability_rrna_ev0.1 \
        $SCALABILITY_DIR/scalability_rrna_ev0.05 \
        $SCALABILITY_DIR/scalability_rrna_ev0.01 \
        $SCALABILITY_DIR/scalability_rrna_ev0.001 \
        $SCALABILITY_DIR/scalability_rrna_ev0.0001 \
        $SCALABILITY_DIR/scalability_rrna_ev0.00001 \
    --nonrrna-dirs \
        $SCALABILITY_DIR/scalability_t2t_ev1 \
        $SCALABILITY_DIR/scalability_t2t_ev0.1 \
        $SCALABILITY_DIR/scalability_t2t_ev0.05 \
        $SCALABILITY_DIR/scalability_t2t_ev0.01 \
        $SCALABILITY_DIR/scalability_t2t_ev0.001 \
        $SCALABILITY_DIR/scalability_t2t_ev0.0001 \
        $SCALABILITY_DIR/scalability_t2t_ev0.00001 \
    --series-labels "T2T non-rRNA" \
    --rrna-dirs \
        $SCALABILITY_DIR/scalability_rrna_ev1 \
        $SCALABILITY_DIR/scalability_rrna_ev0.1 \
        $SCALABILITY_DIR/scalability_rrna_ev0.05 \
        $SCALABILITY_DIR/scalability_rrna_ev0.01 \
        $SCALABILITY_DIR/scalability_rrna_ev0.001 \
        $SCALABILITY_DIR/scalability_rrna_ev0.0001 \
        $SCALABILITY_DIR/scalability_rrna_ev0.00001 \
    --nonrrna-dirs \
        $SCALABILITY_DIR/scalability_rfam_ev1 \
        $SCALABILITY_DIR/scalability_rfam_ev0.1 \
        $SCALABILITY_DIR/scalability_rfam_ev0.05 \
        $SCALABILITY_DIR/scalability_rfam_ev0.01 \
        $SCALABILITY_DIR/scalability_rfam_ev0.001 \
        $SCALABILITY_DIR/scalability_rfam_ev0.0001 \
        $SCALABILITY_DIR/scalability_rfam_ev0.00001 \
    --series-labels "Rfam non-rRNA" \
    --rrna-family-tsv $RRNA_SIM_DIR/rRNA_test_10M_family.tsv \
    --silva-version $SILVA_SSU_VERSION \
    --rfam-version $RFAM_VERSION \
    --smr-db-label "v${SMR_VERSION} default db"
```

A `scalability_benchmark_summary.html` file is written alongside the ROC plot containing: the embedded ROC curve and runtime/RAM plots; per-E-value and per-scale-point tables for T2T non-rRNA reads and Rfam non-rRNA reads (reads classified as rRNA); and a per-family sensitivity breakdown for SILVA rRNA reads showing how many reads from each rRNA family (silva ssu bacteria, rfam 5s, etc.) were assigned to rRNA by SortMeRNA.

- Scalability benchmark with plots: <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/scalability_test/plots/scalability_benchmark_summary.html" target="_blank">scalability_benchmark_summary.html</a>

#### Experiment 2: Sensitivity across database configurations

**Goal:** Does SortMeRNA correctly identify rRNA reads as divergence increases between the read set and the database?

- **Design:** Diagonal 3x3 - each read set tested against its matched database configuration only (same divergence level for reads and database)
- **E-value:** Use the optimal e-value identified from the Experiment 1 ROC curve (the point with the best sensitivity/specificity tradeoff). Run Experiment 1 first and inspect the ROC plot before running Experiment 2.
- **Read sets:**
  - Set 1: 12,500 reads x 8 rRNA types = 100,000 total reads simulated from 97% non-seed members -> tested against SMR sensitive db (97%)
  - Set 2: 12,500 reads x 8 rRNA types = 100,000 total reads simulated from 90-95% non-seed members -> tested against SMR default db (90-95%)
  - Set 3: 12,500 reads x 8 rRNA types = 100,000 total reads simulated from 85-90% non-seed members -> tested against SMR fast db (85-90%)
- **Rationale for 100,000 reads per Set:** Round number that is easy to interpret (1 missed read = 0.001% sensitivity drop), and fast to run. Kept at 100K rather than 1M to avoid conflating database-configuration differences with the E-value scaling effect: at higher read counts `S_min` rises (because `m` = total nucleotide count), which can suppress sensitivity independently of the database being tested. Experiment 1 characterizes that scaling effect separately; Experiment 2 holds read count fixed so results reflect only the database configuration.
- **Metric:** Sensitivity = detected / total (aggregate across all rRNA types)

Run SortMeRNA for each Set against its matched database configuration (rRNA reads produced in Phase 1 - see rRNA Test Sets). Replace `<optimal_evalue>` with the value identified from the Experiment 1 ROC curve:

```bash
export SENSITIVITY_DIR=$DATA_DIR/sensitivity_test

bash $SMR_DB_ROOT_DIR/scripts/benchmarking/run_sensitivity.sh \
    $SENSITIVITY_DIR \
    4 \
    --evalue 1e-5
```

Outputs a per-Set sensitivity table and HTML summary at `$SENSITIVITY_DIR/sensitivity_summary.html`.

#### Real Benchmark Datasets (PacBio)
Sensitivity test using real PacBio long-read amplicon data:
- **Source**: [Karst et al. (2021, *Nature Methods*)](https://doi.org/10.1038/s41592-020-01041-y) - 253,089 high-quality,
  full-length bacterial rRNA operon sequences (~4,500 bp, 16S+ITS+23S) from 70
  AGP human fecal samples, generated using PacBio Sequel II UMI amplicon
  sequencing
- **Rationale**: Every read is a guaranteed true positive by virtue of PCR
  amplification with 27F/2490R primers. Tests SortMeRNA's ability to handle ~4,500 bp long reads.
- **Non-rRNA source**: PBSIM3 simulated reads from the masked T2T genome and Rfam non-rRNA family sequences - tested separately
- **Experiments**:
  - **Sensitivity**: Run all 253,089 operon sequences through SortMeRNA;
    expected classification rate = 100% per database configuration
  - **Specificity**: Run PBSIM3-simulated non-rRNA reads through SortMeRNA;
    measure false positive rate per database configuration

#### PacBio Parameter Optimisation Sweep

SortMeRNA's default seeding parameters (`--passes 18,9,3`, `--num_seeds 2`) were
designed for Illumina and 454 reads in the 100-500 bp range. For ~4,500 bp
PacBio HiFi reads the seed window length (`--L 18`) and Levenshtein distance (1)
remain appropriate, but the pass strides can be made much sparser - reducing the
number of seed lookups per read from ~1,495 (default) down to ~90 - without
sacrificing sensitivity.

**Parameter grid** (4 x 4 = 16 combinations):

| `--passes` | Approx. stride ratio vs default | P1 windows / read | Total windows / read |
|---|---|---|---|
| `18,9,3`    | 1x (default)  | 250 | 1,495 |
| `100,50,10` | ~6x           |  45 |   449 |
| `200,100,20`| ~12x          |  23 |   225 |
| `500,200,50`| ~28x          |   9 |    90 |

`--num_seeds` tested: **2** (default), **5**, **10**, **25**

> **Degenerate combinations**: `--passes 500,200,50 --num_seeds 25` - P1 yields
> only ~9 windows so the seed threshold cannot be met in Pass 1; the code always
> cascades through all three passes. Such combinations appear in the grid to
> confirm they offer no speed benefit over the default.

**Metrics collected per run:**
- **Sensitivity** = aligned reads / 253,089 (expected ~100% for all combos - any drop is a regression)
- **Wall-clock time** - primary signal for selecting the winner
- **Mean SW alignment score** - proxy for alignment quality; should be stable across sparse-seeding combos

**Parameter sweep script:**

```bash
bash $SMR_DB_ROOT_DIR/scripts/benchmarking/run_pacbio_sweep.sh \
    /path/to/karst2021_253k.fastq \
    /path/to/pbsim3_nonrrna_253k.fastq \
    $WORK_DIR/results/pacbio_sweep \
    4
```

**Expected results table** (`sweep_results.tsv`):

```
passes        num_seeds  total   aligned  sensitivity  wall_sec
18,9,3        2          253089  ...      ...          ...
18,9,3        5          253089  ...      ...          ...
...
500,200,50    25         253089  ...      ...          ...
```

**Interpreting results**: sensitivity should remain flat at ~1.0000 across all
combinations (HiFi error rate ~0.1% means sparse seeds still land in error-free
stretches). The winning combination is the sparsest `--passes` / highest
`--num_seeds` pair that maintains full sensitivity - expect this to be around
`--passes 200,100,20 --num_seeds 5` based on the window-count analysis above.

#### Performance Metrics
- **Sensitivity**: True positives / (True positives + False negatives)
- **Specificity**: True negatives / (True negatives + False positives)
- **Precision**: True positives / (True positives + False positives)
- **F1 Score**: Harmonic mean of precision and recall
- **Runtime**: Wall-clock time for various read counts
- **Memory**: Peak RAM usage
- **Disk space**: Database + index size

### Benchmark SortMeRNA versions

Benchmarking compares the latest SortMeRNA against v2.1b (the version used in the original paper). Since both binaries are named `sortmerna`, install v2.1b separately and reference it by full path.

**Install SortMeRNA v2.1b:**

```bash
wget https://github.com/sortmerna/sortmerna/releases/download/2.1b/sortmerna-2.1b-linux.tar.gz
tar -xzf sortmerna-2.1b-linux.tar.gz -C /home/ubuntu/
sudo apt-get install -y libgomp1
# binary: /home/ubuntu/sortmerna-2.1b/sortmerna
export SMR_V2_BIN=/home/ubuntu/sortmerna-2.1b/sortmerna
```

*(benchmarking scripts coming soon)*

## Expected Outputs

### Database Statistics
- Sequence count reduction per clustering level
- Size reduction (GB → MB)
- Taxonomic coverage analysis
- Representative sequence distribution

### Performance Reports
- Sensitivity/Specificity curves
- Runtime vs. database size plots
- Memory usage profiles
- ROC curves for different thresholds

### Recommendations
- Optimal clustering threshold for general use
- Platform-specific database recommendations
- Trade-off analysis (accuracy vs. efficiency)

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with clear description

## Citation

If you use this benchmark in your research, please cite:

```
[Citation to be added upon publication]
```

## License

LGPL-3.0

## Contact

- **Issues**: https://github.com/sortmerna/sortmerna-database/issues

## References

1. Kopylova E, Noé L, Touzet H. SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. Bioinformatics. 2012 Dec 15;28(24):3211-7. doi: 10.1093/bioinformatics/bts611. Epub 2012 Oct 15
2. Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Res. 2013 Jan;41(Database issue):D590-6. doi: 10.1093/nar/gks1219. Epub 2012 Nov 28
3. Kalvari I, Nawrocki EP, Ontiveros-Palacios N, Argasinska J, Lamkiewicz K, Marz M, Griffiths-Jones S, Toffano-Nioche C, Gautheret D, Weinberg Z, Rivas E, Eddy SR, Finn RD, Bateman A, Petrov AI. Rfam 14: expanded coverage of metagenomic, viral and microRNA families. Nucleic Acids Res. 2021 Jan 8;49(D1):D192-D200. doi: 10.1093/nar/gkaa1047
4. Karst SM, Ziels RM, Kirkegaard RH et al. High-accuracy long-read amplicon sequences using unique molecular identifiers with Nanopore or PacBio sequencing. Nat Methods 18, 165-169 (2021). doi: [10.1038/s41592-020-01041-y](https://doi.org/10.1038/s41592-020-01041-y)
5. Gourlé H, Karlsson-Lindsjö O, Hayer J, Bongcam-Rudloff E. Simulating Illumina metagenomic data with InSilicoSeq. Bioinformatics. 2019 Feb 1;35(3):521-522. doi: [10.1093/bioinformatics/bty630](https://doi.org/10.1093/bioinformatics/bty630)
6. Karlin S, Altschul SF. Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes. Proc Natl Acad Sci U S A. 1990 Mar;87(6):2264-2268. doi: [10.1073/pnas.87.6.2264](https://doi.org/10.1073/pnas.87.6.2264)
7. Madden T. The BLAST Sequence Analysis Tool. 2013 Mar 15. In: The NCBI Handbook [Internet]. 2nd edition. Bethesda (MD): National Center for Biotechnology Information (US); 2013-. Available from: https://www.ncbi.nlm.nih.gov/books/NBK153387/
8. National Library of Medicine. BLAST Databases. July 2023. Available from: https://www.nlm.nih.gov/ncbi/workshops/2023-08_BLAST_evol/databases.html
