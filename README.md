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
- [ ] Simulate Illumina 150bp rRNA and non-rRNA reads using non-seed cluster members and non-rRNA sequences (ART)
- [ ] Download real Illumina metatranscriptomics data
- [ ] Download real PacBio amplicon data (Karst et al. 2021)
- [ ] Download real PacBio metagenomics data
- [ ] Benchmark latest SortMeRNA vs. SortMeRNA v2.1b (original paper) and other tools: sensitivity, specificity, runtime, memory


## Installation

### Requirements

**Core tools:**
- [SortMeRNA](https://github.com/sortmerna/sortmerna) v5.0.0 - rRNA filtering (installed separately, see below)
- [VSEARCH](https://github.com/torognes/vsearch) >= 2.22 - sequence clustering
- [Infernal](http://eddylab.org/infernal/) >= 1.1.4 - covariance model search (cmsearch, cmpress)
- [SeqKit](https://bioinf.shenwei.me/seqkit/) >= 2.5 - sequence statistics

**Languages:**
- Python >= 3.8
- R >= 4.0

**Python libraries:**
- pandas, numpy, matplotlib, seaborn, biopython

**Simulation tools (for benchmarking):**
- [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) - Illumina read simulation
- [PBSIM3](https://github.com/yukiteruono/pbsim3) - PacBio and Oxford Nanopore read simulation

### Quick Start

```bash
# Clone repository
git clone https://github.com/sortmerna/sortmerna-database.git
cd sortmerna-database

# Create conda environment
conda env create -f environment.yml
conda activate sortmerna-bench
```

**Install SortMeRNA separately:**

```bash
export SMR_VERSION=5.0.0
wget https://github.com/sortmerna/sortmerna/releases/download/v${SMR_VERSION}/sortmerna-${SMR_VERSION}-Linux.tar.gz
tar -xzf sortmerna-${SMR_VERSION}-Linux.tar.gz -C /home/ubuntu/
# binary: /home/ubuntu/sortmerna-${SMR_VERSION}-Linux/bin/sortmerna
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

# SILVA versions and full download URLs (update to use a different release or file type)
export SILVA_SSU_VERSION=138.2
export SILVA_SSU_PATH=https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz
export SILVA_LSU_VERSION=138.2
export SILVA_LSU_PATH=https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_LSURef_NR99_tax_silva_trunc.fasta.gz
export RFAM_VERSION=15.1

# Human T2T genome - used as non-rRNA specificity test source
export T2T_ACCESSION=GCA_009914755.4
# NCBI assembly name used in FTP filenames (${T2T_ACCESSION}_${T2T_NAME}_genomic.fna.gz)
export T2T_NAME=CHM13_T2T_v2.0
# human-readable version used for output filenames
export T2T_VERSION=chm13v2.0
# full FTP directory URL (derived from accession + name)
export T2T_BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/${T2T_ACCESSION}_${T2T_NAME}

# Rfam non-rRNA families FTP (tied to RFAM_VERSION; change to CURRENT to always pull the latest)
export RFAM_NON_RRNA_FTP=https://ftp.ebi.ac.uk/pub/databases/Rfam/$RFAM_VERSION/fasta_files

# SortMeRNA binary - full path
export SMR_VERSION=5.0.0
export SMR_BIN=/home/ubuntu/sortmerna-${SMR_VERSION}-Linux/bin/sortmerna

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
bash $SMR_DB_ROOT_DIR/scripts/database_building/verify_silva.sh $WORK_DIR/data $VERIFIED_DIR 4
```

```bash
# Verify Rfam 5S and 5.8S seed sequences using full CM search with --cut_ga
bash $SMR_DB_ROOT_DIR/scripts/database_building/verify_rfam.sh $WORK_DIR/data $VERIFIED_RFAM_DIR 4
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
bash $SMR_DB_ROOT_DIR/scripts/database_building/cluster_sequences.sh $WORK_DIR/data $CLUSTERED_DIR 4
```

By default, vsearch clustering is skipped for any threshold where the `.uc` output file already exists - only the downstream steps (member extraction, leakage check, summary table) are re-run. To force vsearch to re-cluster and overwrite existing `.uc` files, pass `--force`:

```bash
bash $SMR_DB_ROOT_DIR/scripts/database_building/cluster_sequences.sh --force $WORK_DIR/data $CLUSTERED_DIR 4
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
bash $SMR_DB_ROOT_DIR/scripts/database_building/build_sortmerna_index.sh $WORK_DIR/data $INDEX_DIR 4
```

- Per-configuration index build report (sequence count, build time, index size, peak CPU%, peak RAM): <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/index/index_build_summary.html" target="_blank">index_build_summary.html</a>

Pass `--force` to rebuild an index that already exists:

```bash
bash $SMR_DB_ROOT_DIR/scripts/database_building/build_sortmerna_index.sh --force $WORK_DIR/data $INDEX_DIR 4
```

For each configuration the script produces:

| File / directory | Description |
|---|---|
| `<config>/<config>.fasta` | Combined reference FASTA (all domains concatenated) |
| `<config>/idx/` | SortMeRNA index files |
| `<config>/index.stats` | Sequence count, build time, index size, SortMeRNA version |

**Index building is done once.** The index only needs to be built once per database configuration. When you later run SortMeRNA alignment with the same `--ref` and `--idx-dir` paths, SortMeRNA finds the existing index and skips rebuilding - even across separate invocations or tmux sessions. Use `--force` to explicitly trigger a rebuild.

### Non-rRNA Test Set (`non_rRNA_test_1M.fasta`)

~1M non-rRNA sequences from two sources that together test both straightforward
and challenging specificity cases:

| Source | Description | Count |
|--------|-------------|-------|
| **Human T2T genome (CHM13v2.0)** | Simulated 150bp PE reads, rRNA loci masked prior to simulation | 850,000 |
| **Rfam non-rRNA families** | 10 families chosen because they share structural features with rRNA and are the most likely source of false positives (see table below) | 150,000 |

```bash
bash $SMR_DB_ROOT_DIR/scripts/read_simulation/download_non_rrna.sh $NON_RRNA_DIR 4
```

**Preparation:**
1. Mask rRNA loci in CHM13v2.0 using known T2T annotations
2. Simulate 850,000 150bp PE reads with ART
3. Sample 150,000 Rfam non-rRNA sequences
4. Apply header keyword filter to remove any rRNA-related sequences

**Rfam non-rRNA families:**

| Family | Rfam ID |
|--------|---------|
| tRNA | RF00005 |
| SRP RNA | RF00017 |
| tmRNA | RF00023 |
| RNase P (bacterial) | RF00010 |
| RNase P (eukaryotic) | RF00009 |
| U1 spliceosomal | RF00003 |
| U2 spliceosomal | RF00004 |
| U4 spliceosomal | RF00015 |
| U5 spliceosomal | RF00020 |
| U6 spliceosomal | RF00026 |

All sampling uses a fixed random seed (`--seed 42`) for reproducibility.

**Output files:**

| File | Description |
|------|-------------|
| `non_rRNA_test_1M.fasta` | Final clean non-rRNA test sequences |
| `non_rRNA_metadata.txt` | Composition breakdown by source |
| `non_rRNA_stats.txt` | Sequence length and count statistics |

## Phase 2: Validation and Benchmarking

Evaluate the databases built in Phase 1 for sensitivity (rRNA detection) and specificity (false positive rate), and benchmark the latest SortMeRNA against v2.1b. Simulated Illumina reads from non-seed cluster members test sensitivity at each clustering threshold; non-rRNA reads from the T2T genome and Rfam test specificity. Real PacBio amplicon data provides an independent validation.

### Benchmarking Approach

#### Simulated Data (Illumina)
Generate synthetic Illumina reads with known rRNA/non-rRNA composition:
- **Tool**: ART (150bp paired-end, Illumina-specific error profiles)
- **rRNA source**: Non-seed cluster members (`*_test_members.fasta`) - real rRNA sequences not present in the clustered database
- **Non-rRNA source**: `non_rRNA_test_1M.fasta` - simulated T2T genome reads (rRNA loci masked) and Rfam non-rRNA families (tRNA, SRP RNA, tmRNA, RNase P, spliceosomal RNAs)

**Important note on SortMeRNA as a filter, not a classifier**

SortMeRNA is an rRNA **filter** - it answers "is this read rRNA?" not "which type of rRNA is this?". A read simulated from a bacterial SSU sequence may be detected via a k-mer match and Smith-Waterman extension against an archaeal SSU reference, and this is expected behavior given how conserved rRNA sequences are across domains. Per-type sensitivity is therefore not a meaningful metric. All rRNA types are pooled into a single reads file per Set, and sensitivity is measured as a single aggregate: detected / total.

**Experiment 1: Sensitivity across database configurations**

**Goal:** Does SortMeRNA correctly identify rRNA reads as divergence increases between the read set and the database?

- **Design:** Diagonal 3x3 - each read set tested against its matched database configuration only (same divergence level for reads and database)
- **Read sets:**
  - Set 1: 12,500 reads x 8 rRNA types = 100,000 total reads simulated from 97% non-seed members -> tested against SMR sensitive db (97%)
  - Set 2: 12,500 reads x 8 rRNA types = 100,000 total reads simulated from 90-95% non-seed members -> tested against SMR default db (90-95%)
  - Set 3: 12,500 reads x 8 rRNA types = 100,000 total reads simulated from 85-90% non-seed members -> tested against SMR fast db (85-90%)
- **Rationale for 100,000 reads per Set:** Round number that is easy to interpret (1 missed read = 0.001% sensitivity drop), aligns naturally with the first data point of the scalability experiment, and is fast to run
- **Metric:** Sensitivity = detected / total (aggregate across all rRNA types)

**Experiment 2: Scalability**

**Goal:** How do runtime and memory scale with read volume, and does sensitivity remain stable at large scale?

- **Design:** Fixed configuration (Set 2 non-seeds -> SMR default db) - the most representative real-world deployment scenario - across increasing read volumes
- **Read volumes:** 10,000 -> 100,000 -> 1,000,000 -> 10,000,000 reads
- **Rationale for default db:** The configuration most users would deploy in practice; runtime numbers are directly interpretable for real-world use
- **Metrics:**
  - Wall-clock runtime at each scale point
  - Peak RAM usage at each scale point
  - Sensitivity at each scale point - should remain flat since SortMeRNA processes reads independently; any degradation at 10M reads would indicate memory pressure or a bug worth flagging

#### Real Benchmark Datasets (PacBio)
Sensitivity test using real PacBio long-read amplicon data:
- **Source**: [Karst et al. (2021, *Nature Methods*)](https://doi.org/10.1038/s41592-020-01041-y) - 253,089 high-quality, 
  full-length bacterial rRNA operon sequences (~4,500 bp, 16S+ITS+23S) from 70 
  AGP human fecal samples, generated using PacBio Sequel II UMI amplicon 
  sequencing
- **Rationale**: Every read is a guaranteed true positive by virtue of PCR 
  amplification with 27F/2490R primers. Tests SortMeRNA's ability to handle ~4,500 bp long reads.
- **Non-rRNA source**: PBSIM3 simulated reads from `non_rRNA_test_1M.fasta` (simulated T2T genome reads and Rfam non-rRNA families)
- **Experiments**:
  - **Sensitivity**: Run all 253,089 operon sequences through SortMeRNA; 
    expected classification rate = 100% per database configuration
  - **Specificity**: Run PBSIM3-simulated non-rRNA reads through SortMeRNA; 
    measure false positive rate per database configuration
  - **Runtime/scaling**: Mix real rRNA + simulated non-rRNA reads at 0%, 10%, 
    25%, 50%, 75%, 90% rRNA composition

#### Performance Metrics
- **Sensitivity**: True positives / (True positives + False negatives)
- **Specificity**: True negatives / (True negatives + False positives)
- **Precision**: True positives / (True positives + False positives)
- **F1 Score**: Harmonic mean of precision and recall
- **Runtime**: Wall-clock time for various read counts
- **Memory**: Peak RAM usage
- **Disk space**: Database + index size

### Simulate Reads

*(coming soon)*

### Download PacBio Amplicon Data

*(coming soon)*

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

1. Kopylova E, Noé L, Touzet H. SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. Bioinformatics. 2012 Dec 15;28(24):3211-7. doi: 10.1093/bioinformatics/bts611. Epub 2012 Oct 15. PMID: 23071270.
2. Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Res. 2013 Jan;41(Database issue):D590-6. doi: 10.1093/nar/gks1219. Epub 2012 Nov 28. PMID: 23193283; PMCID: PMC3531112.
3. Kalvari I, Nawrocki EP, Ontiveros-Palacios N, Argasinska J, Lamkiewicz K, Marz M, Griffiths-Jones S, Toffano-Nioche C, Gautheret D, Weinberg Z, Rivas E, Eddy SR, Finn RD, Bateman A, Petrov AI. Rfam 14: expanded coverage of metagenomic, viral and microRNA families. Nucleic Acids Res. 2021 Jan 8;49(D1):D192-D200. doi: 10.1093/nar/gkaa1047. PMID: 33211869; PMCID: PMC7779021.
4. Karst SM, Ziels RM, Kirkegaard RH et al. High-accuracy long-read amplicon sequences using unique molecular identifiers with Nanopore or PacBio sequencing. Nat Methods 18, 165-169 (2021). https://doi.org/10.1038/s41592-020-01041-y

