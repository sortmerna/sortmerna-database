# SortMeRNA Database Benchmark

A comprehensive benchmarking framework for building and evaluating optimized rRNA databases for SortMeRNA across modern sequencing platforms.

## Overview

This repository contains code and workflows to:
1. Build clustered rRNA databases from SILVA and RFAM
2. Benchmark database performance (accuracy vs. size tradeoffs)
3. Optimize clustering parameters for different use cases
4. Validate databases against simulated and real datasets
   
## Repository Structure

```
.
├── scripts/
│   ├── database_building/
│   │   ├── download_silva.sh      # Download SILVA SSU/LSU rRNA sequences
│   │   ├── download_rfam.sh       # Download RFAM 5S/5.8S rRNA sequences
│   │   ├── cluster_sequences.sh   # Cluster sequences by kingdom at multiple thresholds
```

## Goals

### Phase 1: Database Construction
- [x] Download latest SILVA database
- [x] Download RFAM rRNA families
- [x] Implement clustering pipeline
- [ ] Test multiple clustering thresholds
- [ ] Build SortMeRNA indices for each clustered database
- [ ] Generate database statistics and metadata

### Phase 2: Database Validation
- [ ] Create simulated datasets with known rRNA content
- [ ] Illumina short reads (75bp, 150bp, 250bp)
- [ ] PacBio HiFi reads (~15kb)
- [ ] Oxford Nanopore reads (10-50kb)
- [ ] Validate against real benchmark datasets
- [ ] Measure sensitivity and specificity per clustering level

### Phase 3: Performance Benchmarking
- [ ] Runtime analysis across database sizes
- [ ] Memory footprint measurements
- [ ] Index building time vs. search time
- [ ] Scalability testing (1M, 10M, 50M, 100M reads)

### Phase 4: Comparison & Recommendations
- [ ] Identify optimal clustering threshold
- [ ] Generate platform-specific recommendations
- [ ] Document database build protocols

## Database Sources

### SILVA rRNA Database
- **Version**: Release 138.2 Ref NR 99
- **URL**: https://www.arb-silva.de/
- **Content**: SSU (16S/18S) and LSU (23S/28S) rRNA sequences
- **Taxonomy**: Bacteria, Archaea, Eukarya
- **Size (raw)**: 510,495 sequences (SSU) and 95,279 sequences (LSU)

### RFAM
- **Version**: Rfam 15.1 (January 2026, 4227 families)
- **URL**: https://rfam.org/
- **Families**: RF00001 (5S), RF00002 (5.8S)
- **Size (raw)**: 712 seed alignment sequences (5S) and 61 seed alignment sequences (5.8S)

## Clustering Strategy

### Tools Considered

We use [VSEARCH](https://github.com/torognes/vsearch) for sequence clustering. VSEARCH is a high-performance open-source tool for sequence analysis that implements a centroid-based clustering algorithm. It is the same tool used by the SILVA project to generate their non-redundant (NR) databases, ensuring compatibility and reproducibility with established rRNA reference workflows.

Key features:
- Fast clustering with `--cluster_fast` using identity thresholds
- Supports both strands for rRNA sequences
- Memory-efficient for large databases
- Produces representative centroid sequences

### Clustering Thresholds to Test

- **97%** - Species-level clustering (commonly used for 16S)
- **95%** - Genus-level approximation
- **90%** - Family-level approximation
- **85%** - Higher-order taxonomic grouping

## Benchmarking Approach

### Simulated Data
Generate synthetic reads with known rRNA/non-rRNA composition:
- **Tools**: ART (Illumina), PBSIM3 (PacBio/Nanopore)
- **Composition**: 0%, 10%, 25%, 50%, 75%, 90% rRNA content
- **Organisms**: E. coli, S. cerevisiae, H. sapiens, mixed metatranscriptomic
- **Error profiles**: Platform-specific error rates

### Real Benchmark Datasets
- Public RNA-seq datasets with validated rRNA content
- Spike-in experiments
- Metatranscriptomic samples with known composition

### Performance Metrics
- **Sensitivity**: True positives / (True positives + False negatives)
- **Specificity**: True negatives / (True negatives + False positives)
- **Precision**: True positives / (True positives + False positives)
- **F1 Score**: Harmonic mean of precision and recall
- **Runtime**: Wall-clock time for various read counts
- **Memory**: Peak RAM usage
- **Disk space**: Database + index size

## Installation

### Requirements

**Core tools:**
- [SortMeRNA](https://github.com/sortmerna/sortmerna) v4.3.7 - rRNA filtering (installed separately from GitHub release)
- [VSEARCH](https://github.com/torognes/vsearch) >= 2.22 - sequence clustering
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
git clone https://github.com/yourusername/sortmerna-database-benchmark.git
cd sortmerna-database-benchmark

# Create conda environment
conda env create -f environment.yml
conda activate sortmerna-bench
```

**Install SortMeRNA separately:** Download v4.3.7 binaries from https://github.com/sortmerna/sortmerna/releases

## Usage

### 1. Download Source Databases

```bash
# Download SILVA
bash scripts/database_building/download_silva.sh

# Download RFAM
bash scripts/database_building/download_rfam.sh
```

### 2. Build Clustered Databases

```bash
# Args: input_dir output_dir threads
bash scripts/database_building/cluster_sequences.sh data data/clustered 4
```

For each database and clustering threshold, the script outputs four files:

| File | Description |
|------|-------------|
| `*_XX.fasta` | Centroid/seed sequences — the clustered rRNA database used by SortMeRNA |
| `*_XX.uc` | VSEARCH cluster membership file |
| `*_XX_test_members.fasta` | Non-seed cluster members — used as source sequences for simulating test reads |
| `*_XX_cluster_mapping.txt` | Tab-delimited mapping of each member sequence ID to its seed sequence ID |

The centroid sequences (`*_XX.fasta`) become the SortMeRNA reference databases. The non-seed members (`*_XX_test_members.fasta`) are sequences that were clustered away at each threshold, providing a natural source of reads for benchmarking — since they are real rRNA sequences not present in the database, they test whether SortMeRNA can still identify similar but non-identical rRNA.

### 3. Download Non-rRNA Test Sequences (Specificity Testing)

To measure the false positive rate (specificity), we need a large set of sequences that are definitively **not** rRNA. SortMeRNA should reject all of these; any that are classified as rRNA are false positives.

```bash
bash scripts/read_simulation/download_non_rrna.sh -o data/non_rrna --threads 8
```

The script downloads ~1M non-rRNA sequences from four sources:

| Source | Description | Default count |
|--------|-------------|---------------|
| **RefSeq bacterial mRNA** | Protein-coding transcripts from NCBI RefSeq bacteria RNA, filtered to exclude rRNA/tRNA | 500,000 |
| **Ensembl eukaryotic cDNA** | cDNA from human, mouse, zebrafish, *C. elegans*, and *Arabidopsis* | 300,000 |
| **Rfam non-rRNA families** | Non-coding RNA from Rfam families: tRNA, SRP RNA, tmRNA, RNase P, spliceosomal RNAs | 150,000 |
| **Random genomic fragments** | 500–3000 bp fragments from *E. coli*, *B. subtilis*, and *S. cerevisiae* genomes | 50,000 |

**Safety filters** are applied to the combined set to remove any rRNA contamination:

1. **Header keyword filter** — `seqkit grep` removes sequences with rRNA-related terms (ribosomal, rRNA, 16S, 23S, etc.)
2. **Barrnap rRNA prediction** — runs HMM-based rRNA gene detection for both bacterial (`--kingdom bac`) and eukaryotic (`--kingdom euk`) models, catching unlabeled rRNA especially in genomic fragments

Output files:

| File | Description |
|------|-------------|
| `non_rRNA_test_1M.fasta` | Final clean non-rRNA test sequences |
| `non_rRNA_metadata.txt` | Composition breakdown by source |
| `non_rRNA_stats.txt` | Sequence length and count statistics |

All sampling uses a fixed random seed (`--seed 42`) for reproducibility.

### 4. Run Benchmarks

### 5. Compare Databases

```bash
# Compare all clustering levels
python scripts/benchmarking/compare_databases.py \
  --results results/*_metrics.json \
  --output results/comparison_report.html
```

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

GPL-3.0

## Contact

- **Issues**: https://github.com/yourusername/sortmerna-database-benchmark/issues
- **Email**: jenya.kopylov@gmail.com

## References

1. Kopylova E, Noé L, Touzet H. SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. Bioinformatics. 2012 Dec 15;28(24):3211-7. doi: 10.1093/bioinformatics/bts611. Epub 2012 Oct 15. PMID: 23071270.
2. Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Res. 2013 Jan;41(Database issue):D590-6. doi: 10.1093/nar/gks1219. Epub 2012 Nov 28. PMID: 23193283; PMCID: PMC3531112.
3. Kalvari I, Nawrocki EP, Ontiveros-Palacios N, Argasinska J, Lamkiewicz K, Marz M, Griffiths-Jones S, Toffano-Nioche C, Gautheret D, Weinberg Z, Rivas E, Eddy SR, Finn RD, Bateman A, Petrov AI. Rfam 14: expanded coverage of metagenomic, viral and microRNA families. Nucleic Acids Res. 2021 Jan 8;49(D1):D192-D200. doi: 10.1093/nar/gkaa1047. PMID: 33211869; PMCID: PMC7779021.

**Status**: 🚧 Active Development
