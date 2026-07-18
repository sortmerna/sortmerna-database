# SortMeRNA Database Benchmark

[![CI](https://github.com/sortmerna/sortmerna-database/actions/workflows/tests.yml/badge.svg)](https://github.com/sortmerna/sortmerna-database/actions/workflows/tests.yml)
[![Coverage](https://codecov.io/gh/sortmerna/sortmerna-database/branch/main/graph/badge.svg)](https://codecov.io/gh/sortmerna/sortmerna-database)

A benchmarking framework for building and evaluating rRNA databases for SortMeRNA across modern sequencing platforms.

**Status**: 🚧 Active Development - this repository is undergoing daily changes and is not yet ready for use. Scripts, workflows, and outputs may change without notice. Do not use this repository in production or depend on its outputs until this status is updated to indicate a stable release.

## Table of Contents

- [Latest databases](#latest-databases)
- [Overview](#overview)
- [Installation](#installation)
- [Phase 1: Database Construction](#phase-1-database-construction-1)
  - [Download Source Databases](#download-source-databases)
  - [Verify Sequences](#verify-sequences)
  - [Build Clustered Databases](#build-clustered-databases)
  - [Build SortMeRNA Indices](#build-sortmerna-indices)
  - [Non-rRNA Test Sets](#non-rrna-test-sets)
  - [rRNA Test Sets](#rrna-test-sets)
- [Phase 2: Validation and Benchmarking](#phase-2-validation-and-benchmarking-1)
  - [Experiment 1: Scalability](#experiment-1-scalability)
  - [Experiment 2: Sensitivity across database configurations](#experiment-2-sensitivity-across-database-configurations)
  - [Experiment 3: Benchmark on Deng et al. 2022 datasets](#experiment-3-benchmark-on-deng-et-al-2022-datasets)
  - [Experiment 4: Short-read Illumina Metatranscriptomics](#experiment-4-short-read-illumina-metatranscriptomics)
  - [Experiment 5: PacBio rRNA Operon](#experiment-5-pacbio-rrna-operon-karst-et-al-2021)
  - [Experiment 6: PacBio Metagenomics](#experiment-6-pacbio-metagenomics-minich-et-al-2025)
  - [PacBio Parameter Optimisation Sweep](#pacbio-parameter-optimisation-sweep)
- [AWS Instances](#aws-instances)
- [References](#references)

## Latest databases

| Database | Sequences | Index size | Clustering | Recommended for | Link |
|---|---|---|---|---|---|
| `smr_v6.0.2_sensitive_db` | 515,371 | 3.8 GB | 97% all | Maximum sensitivity | coming soon
| `smr_v6.0.2_default_db` | 240,397 | 2.0 GB | 90-95% SILVA, 90% Rfam | General use (recommended) | coming soon
| `smr_v6.0.2_fast_db` | 137,179 | 1.1 GB | 85-90% SILVA, 85% Rfam | Speed-critical workflows | coming soon

> [!NOTE]
> Full build report: <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/index/index_build_summary.html" target="_blank">index_build_summary.html</a>

These databases differ from the raw SILVA/Rfam releases in three ways:

- **Trimmed to cmsearch coordinates** - flanking non-rRNA sequence is removed before indexing. SortMeRNA indexes every nucleotide in each reference, so untrimmed flanks would increase the false-positive rate.
- **Clustered with VSEARCH** - redundant sequences are merged at identity thresholds of 85-97%, reducing database size and index memory without meaningful sensitivity loss.
- **Hard-masked low-complexity regions** - VSEARCH DUST soft-masking (lowercase) is converted to Ns before indexing, preventing low-complexity kmers from contributing to spurious alignments.

These databases can be used with any SortMeRNA version. For SortMeRNA v6.0.0 and later we recommend `-e 1e-5` for Illumina/454 data and `-e 1e-20` for PacBio instead of the previous default of `-e 1`. Version `v6.0.0` switched the hash index from CMPH to BBHash and alignment from SSW to Parasail, and latest benchmarks show a minor improvement in selectivity at this threshold with no impact on sensitivity or runtime.

**Step 1 - Build the index (once per database):**
```bash
sortmerna --ref smr_v6.0.2_default_db.fasta --idx-dir idx/ --task 5
```

**Step 2 - Filter reads (reuse the index for every run):**
```bash
sortmerna --ref smr_v6.0.2_default_db.fasta --reads reads_R1.fq.gz --reads reads_R2.fq.gz \
    --idx-dir idx/ --workdir run/ --fastx --paired_in --threads 4 -e 1e-5
```

`--workdir` is the base directory for everything that doesn't have its own explicit path override. Under it, SortMeRNA creates:

| Subdirectory | Contents | Override option |
|---|---|---|
| `idx/` | Reference index | `--idx-dir` |
| `kvdb/` | Alignment key-value DB (intermediate) | `--kvdb` |
| `readb/` | Pre-processed/split reads | `--readb` |
| `out/` | Final output files (`aligned.*`, `other.*`) | `--aligned` / `--other` |

In the command above, `--workdir run/` combined with `--idx-dir idx/` means the index is read from `idx/` (explicit override) and everything else (kvdb, output files) lands under `run/`.

> **Note:** SortMeRNA will exit with an error if `run/kvdb/` already exists and is non-empty. Remove or empty it between runs.

**Paired-end output behavior (`--paired_in` / `--paired_out`):**

The default behavior (neither flag set) is per-read: each read is written to the aligned output only if it individually aligns. The pair is not kept together - one read can end up in the aligned output while its mate is absent (or in the `--other` output if specified).

| Scenario | `--paired_in` | `--paired_out` | Default (neither) |
|---|---|---|---|
| Both reads align | both -> aligned | both -> aligned | both -> aligned |
| Only R1 aligns | both -> aligned | both -> non-aligned | R1 -> aligned, R2 dropped |
| Only R2 aligns | both -> aligned | both -> non-aligned | R2 -> aligned, R1 dropped |

`--paired_in` is the conservative choice for most rRNA filtering workflows: it keeps pairs intact, which downstream tools (assemblers, quantifiers) require. Use `--paired_out` if you want the non-rRNA output to be pair-intact instead (e.g. for downstream assembly of the non-rRNA fraction).

## Overview

This repository contains code and workflows to:
1. Download and independently verify rRNA sequences using Infernal covariance models
2. Build clustered rRNA databases from verified SILVA and Rfam sequences
3. Benchmark database performance (accuracy vs. size tradeoffs)

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
- [x] Simulate rRNA reads from non-seed members at each clustering threshold
- [x] Experiment 1: Run SortMeRNA on T2T non-rRNA (10K-10M), Rfam ncRNA (10K-500K), and rRNA reads (10K-10M); generate runtime, FP rate, sensitivity, and E-value plots
- [x] Experiment 2: Run SortMeRNA against matched database configuration and measure sensitivity
- [x] Experiment 3: Run SortMeRNA on Deng et al. 2022 benchmark datasets (sensitivity and specificity)
- [ ] Experiment 4: Short-read Illumina metatranscriptomics data
- [x] Experiment 5: PacBio rRNA operon sensitivity/specificity (Karst et al. 2021, parameter sweep)
- [ ] Experiment 6: PacBio metagenomics classification (Minich et al. 2025)


## Installation

### AWS Instances

Benchmarking was run on AWS EC2:
- **c6i.xlarge** - 4 vCPUs, 8 GB RAM (compute-intensive runs)
- **r6i.xlarge** - 4 vCPUs, 32 GB RAM (memory-intensive runs)
- **r6i.16xlarge** - 64 vCPUs, 512 GB RAM (compute + memory-intensive runs)

### Requirements

**Core tools:**
- [SortMeRNA](https://github.com/sortmerna/sortmerna) v6.0.2 - rRNA filtering (installed from source, see Quick Start)
- [VSEARCH](https://github.com/torognes/vsearch) >= 2.22 - sequence clustering
- [Infernal](http://eddylab.org/infernal/) >= 1.1.4 - covariance model search (cmsearch, cmpress)
- [SeqKit](https://bioinf.shenwei.me/seqkit/) >= 2.5 - sequence statistics

**Languages:**
- Python >= 3.8
- R >= 4.0

**Simulation tools (for benchmarking):**
- [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) - Illumina read simulation
- [PBSIM3](https://github.com/yukiteruono/pbsim3) - PacBio and Oxford Nanopore read simulation

### Quick Start

**Install SortMeRNA**

Here SortMeRNA was installed from source, but it's quicker to use the binaries when available in the release.

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
export SENSITIVITY_DIR=$DATA_DIR/sensitivity_test
export PACBIO_DIR=$DATA_DIR/pacbio
export THIRD_PARTY_ILLUMINA_BENCHMARK_DIR=/home/ubuntu/results-r6i.16xlarge

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

Trimming to cmsearch coordinates reduces false positives: SortMeRNA indexes every nucleotide in each reference sequence, so non-rRNA flanking sequence would be indexed alongside rRNA and increase the false-positive rate.

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

**Multi-hit merging.** When cmsearch reports multiple above-threshold hits for the same sequence (typically separate local alignments over different segments of the same rRNA gene), the trimming window is expanded to `[min(seq_from), max(seq_to)]` - but only if the merged span covers at least **85%** of the original sequence length. Below that threshold the sequence is trimmed to the single highest-scoring hit and flagged `_REVIEW` (suggesting a chimera or non-rRNA insertion). The 85% threshold is permissive enough to include hypervariable inter-hit gaps, which are genuine rRNA, while guarding against merging across large non-rRNA insertions.

**Minimum trimmed length.** Sequences shorter than **50%** of the original length after trimming are routed to `flagged_*.fasta` (`too_short` in the log).

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

> [!NOTE]
> SILVA 138.2: <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/verified/verification_summary.html" target="_blank">verification_summary.html</a>
> Rfam 15.1: <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/verified_rfam/verification_summary.html" target="_blank">verification_summary.html</a>

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

For each database and clustering threshold, the script outputs five files:

| File | Description |
|------|-------------|
| `*_XX.fasta` | Centroid/seed sequences with VSEARCH DUST soft masking (lowercase) |
| `*_XX_masked.fasta.gz` | Hard-masked copy of centroids - lowercase bases converted to N |
| `*_XX.uc` | VSEARCH cluster membership file |
| `*_XX_test_members.fasta` | Non-seed cluster members - used as source sequences for simulating test reads |
| `*_XX_cluster_mapping.txt` | Tab-delimited mapping of each member sequence ID to its seed sequence ID |

During clustering, VSEARCH applies DUST low-complexity masking (`--qmask dust`, the default), which lowercases low-complexity regions (e.g. tandem repeats such as AT-microsatellites) in the centroid output. A hard-masked copy (`*_XX_masked.fasta.gz`) is immediately produced by converting all lowercase bases to N using `vsearch --fastx_mask --hardmask`, and verified against the soft-masked file. The hard-masked file is what gets indexed by SortMeRNA.

The centroid sequences (`*_XX_masked.fasta`) become the SortMeRNA reference databases. The non-seed members (`*_XX_test_members.fasta`) are sequences that were clustered away at each threshold, providing a natural source of reads for benchmarking - since they are real rRNA sequences not present in the database, they test whether SortMeRNA can still identify similar but non-identical rRNA.

> [!NOTE]
> Summary of total sequences per clustering threshold and SortMeRNA reference databases: <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/clustered/clustering_summary.html" target="_blank">clustering_summary.html</a>

### Build SortMeRNA Indices

Three database configurations are assembled from the clustered FASTA files and indexed with SortMeRNA:

| Configuration | SILVA SSU bacteria | SILVA SSU other | SILVA LSU | Rfam 5S / 5.8S |
|---|---|---|---|---|
| `smr_v<version>_sensitive_db` | 97% | 97% | 97% | 97% |
| `smr_v<version>_default_db` | 90% | 95% | 95% | 90% |
| `smr_v<version>_fast_db` | 85% | 90% | 90% | 85% |

The sensitive database maximises recall; the default and fast databases trade a small number of sequences for a smaller index and faster runtime.

```bash
bash $SMR_DB_ROOT_DIR/scripts/database_building/build_sortmerna_index.sh \
    $WORK_DIR/data \
    $INDEX_DIR
```

> [!NOTE]
> Per-configuration index build report (sequence count, build time, index size, peak CPU%, peak RAM): <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/index/index_build_summary.html" target="_blank">index_build_summary.html</a>

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
| `<config>/family_map.tsv` | Tab-separated mapping of sequence ID to rRNA family (used for BLAST hit breakdown) |

**Index building is done once.** The index only needs to be built once per database configuration. When you later run SortMeRNA alignment with the same `--ref` and `--idx-dir` paths, SortMeRNA finds the existing index and skips rebuilding - even across separate invocations or tmux sessions. Use `--force` to explicitly trigger a rebuild.

### Non-rRNA Test Sets

Two separate test sets from different sources, used independently to measure specificity:

| File | Source | Description |
|------|--------|-------------|
| `non_rRNA_test_10M_T2T.fasta` | Human T2T genome (CHM13v2.0) | 10M simulated 150bp PE reads, rRNA loci masked prior to simulation |
| `non_rRNA_test_Rfam.fasta` | Rfam non-rRNA families | 500K sequences sampled evenly across 10 families that share structural features with rRNA |

```bash
bash $SMR_DB_ROOT_DIR/scripts/read_simulation/download_non_rrna.sh $NON_RRNA_DIR 4
```

**Preparation (T2T):**
1. Download CHM13v2.0 genome and RefSeq GFF3 annotation
2. Supplement the GFF3 with Infernal cmsearch (RF01960, RF02543, RF00001, RF00002 against the full genome; RF00177 and RF02541 against chrM only) to find rRNA copies not individually annotated in RefSeq (see below)
3. Mask rRNA loci with bedtools maskfasta
4. Simulate 10M 150bp PE reads with InSilicoSeq (NovaSeq error model)

**Supplementing GFF3 with cmsearch:**
The T2T-CHM13v2.0 assembly resolves the Nucleolar Organizer Regions (NORs) on the short arms of acrocentric chromosomes (chr13, chr14, chr15, chr21, chr22), assembling the large tandem arrays of rDNA repeat units encoding 18S/5.8S/28S rRNA. Comparison of cmsearch-identified rDNA loci against the RefSeq GFF3 annotation reveals rDNA sequence not covered by the annotation, demonstrating that cmsearch supplements the GFF3 by recovering rDNA loci not individually annotated in RefSeq. To run cmsearch, set `CMS_DIR` to a directory of pressed Rfam covariance models (see `download_non_rrna.sh` for details); without it the script falls back to GFF3-only masking.

BLAST analysis of the RF01960 (18S) regions not covered by GFF3 shows that the large majority are rRNA-derived: the most common top hits are 45S pre-ribosomal RNA entries, followed by 18S ribosomal pseudogene entries. All are included in the masking BED - if a sequence is annotated as a ribosomal pseudogene it is rRNA-derived and SortMeRNA, which operates purely on sequence content, would flag reads from it.

Note: Mitochondrial 12S (RF00177) and 16S (RF02541) rRNA genes are identified by cmsearch against chrM only and added to the masking BED.

**Preparation (Rfam):**
1. Download 10 non-rRNA families (tRNA, SRP RNA, tmRNA, RNase P, spliceosomal RNAs)
2. Sample evenly across families

Sequences are used as-is without read simulation. Most Rfam families consist of short RNAs (tRNA averages 73 bp, spliceosomal snRNAs 100-200 bp); simulating 150 bp reads from these sequences would discard the majority of them due to length. The natural length distribution is also the point: these structurally complex RNAs should be rejected by SortMeRNA regardless of length.

Sampling uses fair-share allocation: the per-family quota starts at `target / n_families` and is recomputed after each round as small families are exhausted and removed, so larger families absorb the slack. If total available sequences across all families is less than the target, all sequences are used. A fixed random seed (`--seed 42`) makes the output reproducible.

```bash
bash $SMR_DB_ROOT_DIR/scripts/read_simulation/simulate_non_rrna.sh \
    $NON_RRNA_DIR \
    4 \
    --t2t-reads 10000000 \
    --rfam-reads 500000
```

> [!NOTE]
> Non-rRNA reference sources summary for read simulation: <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/non_rrna/non_rrna_test_set_summary.html" target="_blank">non_rrna_test_set_summary.html</a>

### rRNA Test Sets

Three Sets of simulated Illumina rRNA reads, one per database configuration. Each Set uses non-seed cluster members at the matching clustering threshold as source sequences, so every simulated read is a real rRNA sequence absent from the database being tested:

| Set | Source members | Matched database | Reads |
|-----|---------------|-----------------|-------|
| Set 1 | 97% non-seeds (all types) | `smr_v${SMR_VERSION}_sensitive_db` | 12,500 per rRNA type x 8 types = 100,000 |
| Set 2 | bacteria SSU 90%, others 95%, Rfam 90% | `smr_v${SMR_VERSION}_default_db` | 12,500 per rRNA type x 8 types = 100,000 |
| Set 3 | bacteria SSU 85%, others 90%, Rfam 85% | `smr_v${SMR_VERSION}_fast_db` | 12,500 per rRNA type x 8 types = 100,000 |

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

> [!NOTE]
> rRNA read simulation summary: <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/rrna_sim/rrna_simulation_summary.html" target="_blank">rrna_simulation_summary.html</a>

**Scalability pool (Experiment 1):** `simulate_rrna_reads.sh` also produces `rRNA_test_10M.fasta` alongside the three sensitivity sets. The pool is built from Set 2 non-seed sources (default database, 90-95% clustering thresholds) and combines two source types:
- **SILVA types (6 types, SSU + LSU):** IUPAC-cleaned and run through ISS (NovaSeq model, 150bp PE). ISS generates 10M reads.
- **Rfam 5S (avg 117bp) and 5.8S (avg 150bp):** included directly as-is (IUPAC-cleaned, no min-length filter). These sequences are shorter than the ISS read length and cannot be reliably simulated; including them directly ensures short rRNA are represented in the scalability pool.

The SILVA ISS reads and Rfam direct sequences are shuffled together before output.

## Phase 2: Validation and Benchmarking

Six experiments evaluate the databases built in Phase 1 across read types and experimental contexts:

- **Experiments 1-2** use simulated Illumina reads - scalability across read volumes and sensitivity across database configurations
- **Experiment 3** benchmarks on published, independently curated short-read datasets ([Deng et al. 2022](https://doi.org/10.1093/nar/gkac112))
- **Experiment 4** tests on real short-read Illumina metatranscriptomics (pending)
- **Experiment 5** evaluates sensitivity and specificity on real PacBio HiFi full-length rRNA operon sequences ([Karst et al. 2021](https://doi.org/10.1038/s41592-020-01041-y)) and includes a parameter sweep for long-read settings
- **Experiment 6** tests on real PacBio HiFi metagenomics where rRNA is an unknown fraction of the community ([Minich et al. 2025](https://www.cell.com/cell/fulltext/S0092-8674%2825%2900975-4)) (pending)

### Benchmarking Approach

#### Simulated Data (Illumina)
Generate synthetic Illumina reads with known rRNA/non-rRNA composition:
- **Tool**: InSilicoSeq (150bp paired-end, NovaSeq error model)
- **rRNA source**: Non-seed cluster members (`*_test_members.fasta`) - real rRNA sequences not present in the clustered database
- **Non-rRNA source**: `non_rRNA_test_10M_T2T.fasta` and `non_rRNA_test_Rfam.fasta` - tested separately

#### E-value Filtering

SortMeRNA uses the Karlin-Altschul framework (`E = K · m · n · exp(-λ · S)`) with Gumbel parameters (λ, K) [6,7] computed via the ALP library, where `m` and `n` are the effective query and database lengths (raw lengths minus an edge-effect correction), `S` is the alignment score, and `E` is the expected number of alignments with score >= S by chance. Rather than computing a floating-point E-value per alignment, SortMeRNA inverts the formula once at startup to derive a minimum Smith-Waterman score (`S_min = ln(E / (K·m·n)) / (-λ)`), then filters reads with a single integer comparison during alignment.

**Distinction from BLAST**: in BLAST, `m` is the length of the individual query sequence, giving a per-query E-value. In SortMeRNA, `m` is the total nucleotide count across all reads in the dataset (every read is filtered against the same score threshold regardless of its own length). This means `S_min` sets a run-level threshold: the expected number of spurious alignments across all reads against the database is <= E, not <= E per read.

This design is motivated by the difference in database scale. SortMeRNA's rRNA reference databases are ~143M total bases (~240K sequences), roughly 9,000x smaller than BLAST's nt database (~1.3 trillion bases, ~96M sequences) (National Library of Medicine, 2023). If SortMeRNA used BLAST's per-query approach with `m` = individual read length (e.g. 150 bp), the search space `K·m·n` would be so small that many alignments would pass the filter. By setting `m` to the total reads length, SortMeRNA trades BLAST's per-query statistical framing for a run-level one in order to produce a meaningful threshold despite having a reference database that is ~9,000x smaller than BLAST's.

#### Experiment 1: Scalability

**Goal:** How do runtime and memory scale with read volume, and do sensitivity and false positive rate remain stable at large scale?

- **Design:** Fixed configuration (SMR default db) across increasing read volumes, run separately for rRNA reads, T2T non-rRNA reads, and Rfam non-rRNA reads
- **Read volumes:** T2T and rRNA: 10,000 -> 100,000 -> 1,000,000 -> 10,000,000; Rfam: 10,000 -> 100,000 -> 500,000 (capped by family size)
- **rRNA reads:** Subsampled from Set 2 non-seeds (90-95% identity members: default db) - tests sensitivity and the E-value threshold scaling effect (`S_min` increases with total read count, so sensitivity may shift across scale points)
- **Non-rRNA reads (T2T):** Subsampled from `non_rRNA_test_10M_T2T.fasta` - tests false positive rate at scale up to 10M reads
- **Non-rRNA reads (Rfam):** Subsampled from `non_rRNA_test_Rfam.fasta` - tests specificity against structurally complex ncRNAs, capped at 500K due to family size limits
- **Rationale for default db:** The configuration most users would deploy in practice; runtime numbers are directly interpretable for real-world use
- **Metrics:**
  - Wall-clock runtime at each scale point
  - Peak RAM usage at each scale point
  - Sensitivity (rRNA runs) and false positive rate (T2T runs) at each scale point
  - S_min (score threshold) vs. read count - illustrates the E-value scaling effect
  - E-value and % identity distributions of aligned reads at each scale point

Run for T2T non-rRNA reads + Rfam non-rRNA reads (false positive rate at scale) and rRNA reads (sensitivity at scale) across four E-value thresholds. Requires `non_rRNA_test_10M_T2T.fasta` and `non_rRNA_test_Rfam.fasta` from `simulate_non_rrna.sh` and `rRNA_test_10M.fasta` from `simulate_rrna_reads.sh`. The Rfam non-rRNA set has 500K reads so its scale points are capped at 500K:

Subsample reads once, then reuse those reads for the remaining thresholds with `--reads-dir`:

```bash
for ev in 1 0.1 0.05 0.01 0.001 0.0001 1e-5 1e-10 1e-20; do
    reads_dir_t2t=$( [[ "${ev}" == "1" ]] && echo "" || echo "--reads-dir $SCALABILITY_DIR/scalability_t2t_ev1" )
    reads_dir_rfam=$( [[ "${ev}" == "1" ]] && echo "" || echo "--reads-dir $SCALABILITY_DIR/scalability_rfam_ev1" )
    reads_dir_rrna=$( [[ "${ev}" == "1" ]] && echo "" || echo "--reads-dir $SCALABILITY_DIR/scalability_rrna_ev1" )

    bash $SMR_DB_ROOT_DIR/scripts/benchmarking/run_scalability.sh \
        $NON_RRNA_DIR/non_rRNA_test_10M_T2T.fasta.gz \
        $SCALABILITY_DIR/scalability_t2t_ev${ev} \
        4 \
        --index-dir $INDEX_DIR \
        --config smr_v${SMR_VERSION}_default_db \
        --evalue ${ev} ${reads_dir_t2t}

    bash $SMR_DB_ROOT_DIR/scripts/benchmarking/run_scalability.sh \
        $NON_RRNA_DIR/non_rRNA_test_Rfam.fasta.gz \
        $SCALABILITY_DIR/scalability_rfam_ev${ev} \
        4 \
        --index-dir $INDEX_DIR \
        --config smr_v${SMR_VERSION}_default_db \
        --scale 10000,100000,500000 \
        --evalue ${ev} ${reads_dir_rfam}

    bash $SMR_DB_ROOT_DIR/scripts/benchmarking/run_scalability.sh \
        $RRNA_SIM_DIR/rRNA_test_10M.fasta.gz \
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
    --evalues 1 0.1 0.05 0.01 0.001 0.0001 1e-5 1e-10 1e-20 \
    --rrna-dirs \
        $SCALABILITY_DIR/scalability_rrna_ev1 \
        $SCALABILITY_DIR/scalability_rrna_ev0.1 \
        $SCALABILITY_DIR/scalability_rrna_ev0.05 \
        $SCALABILITY_DIR/scalability_rrna_ev0.01 \
        $SCALABILITY_DIR/scalability_rrna_ev0.001 \
        $SCALABILITY_DIR/scalability_rrna_ev0.0001 \
        $SCALABILITY_DIR/scalability_rrna_ev1e-5 \
        $SCALABILITY_DIR/scalability_rrna_ev1e-10 \
        $SCALABILITY_DIR/scalability_rrna_ev1e-20 \
    --nonrrna-dirs \
        $SCALABILITY_DIR/scalability_t2t_ev1 \
        $SCALABILITY_DIR/scalability_t2t_ev0.1 \
        $SCALABILITY_DIR/scalability_t2t_ev0.05 \
        $SCALABILITY_DIR/scalability_t2t_ev0.01 \
        $SCALABILITY_DIR/scalability_t2t_ev0.001 \
        $SCALABILITY_DIR/scalability_t2t_ev0.0001 \
        $SCALABILITY_DIR/scalability_t2t_ev1e-5 \
        $SCALABILITY_DIR/scalability_t2t_ev1e-10 \
        $SCALABILITY_DIR/scalability_t2t_ev1e-20 \
    --series-labels "T2T non-rRNA" \
    --rrna-dirs \
        $SCALABILITY_DIR/scalability_rrna_ev1 \
        $SCALABILITY_DIR/scalability_rrna_ev0.1 \
        $SCALABILITY_DIR/scalability_rrna_ev0.05 \
        $SCALABILITY_DIR/scalability_rrna_ev0.01 \
        $SCALABILITY_DIR/scalability_rrna_ev0.001 \
        $SCALABILITY_DIR/scalability_rrna_ev0.0001 \
        $SCALABILITY_DIR/scalability_rrna_ev1e-5 \
        $SCALABILITY_DIR/scalability_rrna_ev1e-10 \
        $SCALABILITY_DIR/scalability_rrna_ev1e-20 \
    --nonrrna-dirs \
        $SCALABILITY_DIR/scalability_rfam_ev1 \
        $SCALABILITY_DIR/scalability_rfam_ev0.1 \
        $SCALABILITY_DIR/scalability_rfam_ev0.05 \
        $SCALABILITY_DIR/scalability_rfam_ev0.01 \
        $SCALABILITY_DIR/scalability_rfam_ev0.001 \
        $SCALABILITY_DIR/scalability_rfam_ev0.0001 \
        $SCALABILITY_DIR/scalability_rfam_ev1e-5 \
        $SCALABILITY_DIR/scalability_rfam_ev1e-10 \
        $SCALABILITY_DIR/scalability_rfam_ev1e-20 \
    --series-labels "Rfam non-rRNA" \
    --rrna-family-tsv $RRNA_SIM_DIR/rRNA_test_10M_family.tsv \
    --silva-version $SILVA_SSU_VERSION \
    --rfam-version $RFAM_VERSION \
    --smr-db-label "v${SMR_VERSION} default db"
```

A `scalability_benchmark_summary.html` file is written alongside the ROC plot containing: the embedded ROC curve and runtime/RAM plots; per-E-value and per-scale-point tables for T2T non-rRNA reads and Rfam non-rRNA reads (reads classified as rRNA); and a per-family sensitivity breakdown for SILVA rRNA reads showing how many reads from each rRNA family (silva ssu bacteria, rfam 5s, etc.) were assigned to rRNA by SortMeRNA.

> [!NOTE]
> Scalability benchmark with plots: <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/scalability_test/plots/scalability_benchmark_summary.html" target="_blank">scalability_benchmark_summary.html</a>

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
bash $SMR_DB_ROOT_DIR/scripts/benchmarking/run_sensitivity.sh \
    $SENSITIVITY_DIR \
    4 \
    --evalue 1e-5
```

> [!NOTE]
> Outputs a per-Set sensitivity table and HTML summary: <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/sensitivity_test/sensitivity_summary.html" target="_blank">sensitivity_summary.html</a>.

#### Experiment 3: Benchmark on Deng et al. 2022 datasets

**Goal:** Evaluate SortMeRNA (default db, `-e 1e-5`) on the benchmark datasets published in [Deng et al. 2022 (*Nucleic Acids Research*)](https://doi.org/10.1093/nar/gkac112), using all 8 datasets (source from manuscript: [Zenodo](https://zenodo.org/records/5547691)). Note both `RiboDetector_benchmark_datasets.tar.gz` and `RiboDetector_metaT_dataset.tar.gz` were downloaded and `RiboDetector_metaT_dataset`/* reads moved into `RiboDetector_benchmark_datasets` for benchmarking; in addition `oma_silva` and `homd_fp` were provided directly by the authors.

Most datasets were simulated with ART_Illumina v2.3.7 (`-p -l 100 -ss HS25 -m 150 -s 10`, paired-end 100 bp HiSeq 2500 model). FN datasets (rRNA input) measure sensitivity; FP datasets (non-rRNA input) measure specificity; MetaT (mixed) is reported as NA since true labels are unknown at the read level.

| Dataset | Test | Pairs | Description |
|---|---|---|---|
| SILVA_rRNA | FN | 20,000,000 | SILVA SSU+LSU rRNA sequences |
| OMA_CDS | FP | 20,000,000 | prokaryotic and eukaryotic mRNA |
| oma_silva | FP | 1,027,675 | OMA mRNA CDSs with &ge;70% identity to rRNA genes; estimates FPR on rRNA-similar mRNA |
| homd_fp | FP | 100,558 | HOMD oral-microbe mRNA CDSs with &ge;70% identity to high-FPR (&ge;0.5) rRNA hits in OMA_CDS |
| ENA_virus | FP | 27,206,792 | Viral gene sequences from ENA |
| Amplicon_16S | FN | 7,917,920 | Real 16S V1-V2 amplicon reads (oral microbiome study) |
| Human_ncRNA | FP | 6,330,381 | Human non-coding RNA |
| MetaT | FN+FP | 9,165,829 | Oral metatranscriptome: 4.7M prokaryotic mRNA, 2.5M human mRNA, 73K viral mRNA, 1.9M rRNA (21% rRNA fraction) |

Run on AWS EC2 r6i.16xlarge (64 vCPUs, 512 GB RAM) using 40 threads to match the 40 CPU cores used in Deng et al. 2022, allowing a closer comparison of runtime results. Only SortMeRNA results are reported here; benchmarking against other tools is out of scope for this repository and will be covered in the upcoming publication. To run the benchmark:

```bash
bash $SMR_DB_ROOT_DIR/scripts/benchmarking/run_smr_benchmark.sh \
    /home/ubuntu/TESTS/RiboDetector_benchmark_datasets \
    $THIRD_PARTY_ILLUMINA_BENCHMARK_DIR \
    --threads 40
```

> [!NOTE]
> Outputs a per-Set summary HTML: <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/results-r6i.16xlarge/summary.html" target="_blank">summary.html</a>.

#### Experiment 4: Short-read Illumina Metatranscriptomics

- **Status**: [ ] Data download and analysis pending.

#### Experiment 5: PacBio rRNA Operon (Karst et al. 2021)

Sensitivity and specificity test using real PacBio HiFi full-length rRNA operon reads:
- **Source**: [Karst et al. (2021, *Nature Methods*)](https://doi.org/10.1038/s41592-020-01041-y) - 253,089 high-quality,
  full-length bacterial rRNA operon sequences (~4,500 bp, 16S+ITS+23S) from 70
  AGP human fecal samples, generated using PacBio Sequel II UMI amplicon
  sequencing. Raw data: [Qiita study 10317](https://qiita.ucsd.edu/study/description/10317#) (American Gut Project).
- **Rationale**: Every read is a guaranteed true positive by virtue of PCR
  amplification with 27F/2490R primers. Tests SortMeRNA's ability to handle ~4,500 bp long reads.
- **Non-rRNA source**: PBSIM3 simulated reads from the masked T2T genome (253,089 reads at ~4,500 bp mean to match the Karst dataset).

```bash
# Count total reads
seqkit stats -T $PACBIO_DIR/karst2021_253k.fna.gz
file	format	type	num_seqs	sum_len	min_len	avg_len	max_len
/home/ubuntu/working/data/pacbio/karst2021_253k.fna.gz	FASTA	DNA	253089	1117512051	3830	4415.5	5510


bash $SMR_DB_ROOT_DIR/scripts/read_simulation/simulate_pacbio_nonrrna.sh \
    $NON_RRNA_DIR \
    2>&1 | tee $NON_RRNA_DIR/pbsim3.log

seqkit stats -T $NON_RRNA_DIR/non_rrna_pacbio_253089_T2T.fastq.gz
file	format	type	num_seqs	sum_len	min_len	avg_len	max_len
/home/ubuntu/working/data/non_rrna/non_rrna_pacbio_253089_T2T.fastq.gz	FASTQ	DNA	253089	1141508536	98	4510.3	7252
```

- **Experiments**:
  - **Sensitivity**: Run all 253,089 operon sequences through SortMeRNA;
    expected classification rate = 100% per database configuration
  - **Specificity**: Run PBSIM3-simulated non-rRNA reads through SortMeRNA;
    measure false positive rate per database configuration

##### PacBio Parameter Optimisation Sweep

SortMeRNA's default parameters (`--passes 18,9,3`, `--num_seeds 2`, `--min_lis 2`, `-e 1e-3`) were designed for Illumina/454 short reads, where the challenge is sensitivity: a read may contain only a small fragment of one rRNA gene, so the seed and LIS thresholds are kept permissive (2 seeds, LIS >= 2) to avoid missing weakly-represented rRNA signal. The E-value threshold similarly reflects the short-read regime, where a high-scoring alignment over 150 bp is already strong evidence.

`--passes 18,9,3`: the strides are defined relative to the seed length (L, L/2, 3), not read length, so the seed-matching mechanism is identical across read lengths. On long reads, pass 1 alone saturates the seed threshold and early-exit fires immediately, making passes 2 and 3 largely redundant in practice.

**Sweep findings:** For PacBio long reads (5-15 kb), random k-mer matches and degraded pseudogene fragments accumulate enough seeds to pass default thresholds, driving FPR as high as 13.05% against human non-rRNA reads. The parameter sweep reveals:

- **`-e` is the dominant specificity lever** - tightening from `1e-5` to `1e-20` alone reduces FPR from 13.05% to 0.02% with no sensitivity loss
- **`--min_lis` provides targeted additional filtering** - requires co-linear seed chains on a single reference, directly blocking short pseudogene fragments (~50 bp, max LIS ~3) while preserving detection of 5S (~120 bp, max LIS ~7) and larger subunits
- **`--num_seeds` contributes modestly** and is largely subsumed by `--min_lis` at the same threshold value
- **Sensitivity stays 100% throughout** - full-operon reads generate LIS >> 10 against 23S references regardless of parameter choice

```bash
bash $SMR_DB_ROOT_DIR/scripts/benchmarking/run_pacbio_sweep.sh \
    $PACBIO_DIR/karst2021_253k.fna.gz \
    $NON_RRNA_DIR/non_rrna_pacbio_253089_T2T.fastq.gz \
    $PACBIO_DIR/sweep_lis \
    4
```

Sweeps 17 `(num_seeds, min_lis)` pairs x 3 e-values on 10K subsampled reads. Outputs `sweep_results.tsv`, ROC curve (3 panels by e-value), and family breakdown bar charts showing which rRNA subunit drives sensitivity and FPR at each parameter combination.

For PacBio reads, lower e-values such as `-e 1e-10` or better yet `-e 1e-20` give the best specificity, with `--min_lis 2-6`, with no loss of rRNA recovery. For PacBio reads 5-45 kb in length, `--min_lis 6` (`--num_seeds 2`) is recommended.

> [!NOTE]
> PacBio Parameter Sweep results: <a href="https://sortmerna.github.io/sortmerna-database/results/silva_138.2_Rfam_15.1/working/data/pacbio/sweep_report.html" target="_blank">sweep_report.html</a>.

#### Experiment 6: PacBio Metagenomics (Minich et al. 2025)

- **Source**: [Minich et al. (2025, *Cell*)](https://www.cell.com/cell/fulltext/S0092-8674%2825%2900975-4) - PacBio HiFi (gDNA fragmented (megaruptor3 speed31) SMARTBELL3.0 with 5kb size selection at end) metagenomics from 47 fecal samples, mean read length N50 ~9,663 bp (SD += 1,868). Unlike Karst et al., reads are shotgun metagenomic - no guaranteed rRNA content per read. Raw data: [PRJNA1139951](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1139951). Metadata: [Supplementary Table S28](https://data.mendeley.com/datasets/ks5tvfzbzr/1).
- **Rationale**: Tests SortMeRNA on long-read metagenomics where rRNA reads are a small unknown fraction of a mixed community. True labels are unknown at the read level; results reported as fraction of reads classified as rRNA and family breakdown.

Download the first 20 of the 47 available per-sample PacBio HiFi runs, to conserve disk (human reads removed, full-run `pb.concat.no_hsap` libraries; SRA accessions are read from the bundled `assets/Table_S28_15947_PB.ONT.ILMN_metadata_SRA_v2.txt`). Requires sra-tools (`prefetch`, `fasterq-dump`) and `pigz`/`gzip`:

```bash
bash $SMR_DB_ROOT_DIR/scripts/read_simulation/download_pacbio_metagenomics.sh \
    $PACBIO_DIR/pacbio_metagenomics \
    4
```

Writes one gzipped FASTQ per run (`<SRR>.fastq.gz`) to the output directory. Positional arguments are the output directory (default `data/pacbio_metagenomics`), thread count for `fasterq-dump` (default 4), and max runs to download (default `20`; use `0` to download all 47).

Run SortMeRNA and/or Infernal cmsearch on the downloaded samples. SortMeRNA is run once per sample, at the operating point recommended by the Experiment 5 sweep (`-e 1e-20 --num_seeds 2 --min_lis 6`). cmsearch (`--hmmonly --cut_ga`) is run once per sample against the rRNA covariance models (16S/18S/23S/28S/5S/5.8S) as an orthogonal reference - there is no read-level ground truth for shotgun metagenomes, so cmsearch serves as a second, structure-aware predictor to compare against. Select the method with `--tool sortmerna|cmsearch|both` (default `both`):

```bash
bash $SMR_DB_ROOT_DIR/scripts/benchmarking/run_pacbio_metagenomics.sh \
    --tool both \
    $PACBIO_DIR/pacbio_metagenomics \
    $PACBIO_DIR/metagenomics_results \
    4
```

Both methods record reads classified as rRNA, wall time, and peak RSS. Results are written to `metagenomics_results.tsv` (SortMeRNA, one row per sample) and `cmsearch_results.tsv` (one row per sample). The SortMeRNA blast output (read + reference coordinates) and cmsearch `rrna.tblout` (read + model coordinates) are kept per sample so rRNA read coordinates can be compared between the two methods. SortMeRNA requires `SMR_BIN`, `INDEX_DIR`, `SMR_VERSION`; cmsearch requires `CMS_DIR` (pressed Rfam rRNA CMs from `download_cms.sh`) and Infernal (`cmsearch`, `cmpress`).

## Expected Outputs

### Performance Metrics
- **Sensitivity**: True positives / (True positives + False negatives)
- **Specificity**: True negatives / (True negatives + False positives)
- **Precision**: True positives / (True positives + False positives)
- **F1 Score**: Harmonic mean of precision and recall
- **Runtime**: Wall-clock time for various read counts
- **Memory**: Peak RAM usage
- **Disk space**: Database + index size

### Database Statistics
- Sequence count reduction per clustering level
- Size reduction (GB -> MB)
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

1. Kopylova E, Noé L, Touzet H. SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. Bioinformatics. 2012 Dec 15;28(24):3211-7. doi: [10.1093/bioinformatics/bts611](https://doi.org/10.1093/bioinformatics/bts611). Epub 2012 Oct 15
2. Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Res. 2013 Jan;41(Database issue):D590-6. doi: [10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219). Epub 2012 Nov 28
3. Kalvari I, Nawrocki EP, Ontiveros-Palacios N, Argasinska J, Lamkiewicz K, Marz M, Griffiths-Jones S, Toffano-Nioche C, Gautheret D, Weinberg Z, Rivas E, Eddy SR, Finn RD, Bateman A, Petrov AI. Rfam 14: expanded coverage of metagenomic, viral and microRNA families. Nucleic Acids Res. 2021 Jan 8;49(D1):D192-D200. doi: [10.1093/nar/gkaa1047](https://doi.org/10.1093/nar/gkaa1047)
4. Gourlé H, Karlsson-Lindsjö O, Hayer J, Bongcam-Rudloff E. Simulating Illumina metagenomic data with InSilicoSeq. Bioinformatics. 2019 Feb 1;35(3):521-522. doi: [10.1093/bioinformatics/bty630](https://doi.org/10.1093/bioinformatics/bty630)
5. Karlin S, Altschul SF. Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes. Proc Natl Acad Sci U S A. 1990 Mar;87(6):2264-2268. doi: [10.1073/pnas.87.6.2264](https://doi.org/10.1073/pnas.87.6.2264)
6. Madden T. The BLAST Sequence Analysis Tool. 2013 Mar 15. In: The NCBI Handbook [Internet]. 2nd edition. Bethesda (MD): National Center for Biotechnology Information (US); 2013-. Available from: https://www.ncbi.nlm.nih.gov/books/NBK153387/
7. National Library of Medicine. BLAST Databases. July 2023. Available from: https://www.nlm.nih.gov/ncbi/workshops/2023-08_BLAST_evol/databases.html
8. Deng ZL, Münch PC, Mreches R, McHardy AC. Rapid and accurate identification of ribosomal RNA sequences via deep learning. Nucleic Acids Res. 2022 Jun 10;50(10):e60. doi: [10.1093/nar/gkac112](https://doi.org/10.1093/nar/gkac112)
9. Karst SM, Ziels RM, Kirkegaard RH et al. High-accuracy long-read amplicon sequences using unique molecular identifiers with Nanopore or PacBio sequencing. Nat Methods 18, 165-169 (2021). doi: [10.1038/s41592-020-01041-y](https://doi.org/10.1038/s41592-020-01041-y)
10. Minich JJ, Allsing N, Din MO, Tisza MJ, Maleta K, McDonald D, Hartwick N, Mamerto A, Brennan C, Hansen L, Shaffer J, Murray ER, Duong T, Knight R, Stephenson K, Manary MJ, Michael TP. Culture-independent meta-pangenomics enabled by long-read metagenomics reveals associations with pediatric undernutrition. Cell. 2025 Nov 13;188(23):6666-6686.e25. doi: [10.1016/j.cell.2025.08.020](https://doi.org/10.1016/j.cell.2025.08.020). PMID: 40930091.
