import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from build_family_map import parse_header, write_family_map, resolve_db_version, map_taxonomy, N_TAX_LEVELS


class TestParseHeader:
    def test_full_silva_header(self):
        line = ">KU725478.45631.48554 Bacteria;Pseudomonadota;Alphaproteobacteria;Rickettsiales;Mitochondria;Incertae Sedis;Sphagnum angustifolium;size=5"
        seq_id, original, levels = parse_header(line)
        assert seq_id == "KU725478.45631.48554"
        assert original == "KU725478.45631.48554 Bacteria;Pseudomonadota;Alphaproteobacteria;Rickettsiales;Mitochondria;Incertae Sedis;Sphagnum angustifolium;size=5"
        assert levels[0] == "Bacteria"
        assert levels[1] == "Pseudomonadota"
        assert levels[2] == "Alphaproteobacteria"
        assert levels[3] == "Rickettsiales"
        assert levels[4] == "Mitochondria"
        assert levels[5] == "Incertae Sedis"
        assert levels[6] == "Sphagnum angustifolium"
        assert len(levels) == N_TAX_LEVELS

    def test_size_tag_stripped(self):
        # 2 tokens: domain=first, species=last, middle empty
        line = ">AB123.1.500 Bacteria;Firmicutes;size=10"
        seq_id, original, levels = parse_header(line)
        assert original == "AB123.1.500 Bacteria;Firmicutes;size=10"
        assert levels[0] == "Bacteria"   # domain
        assert levels[6] == "Firmicutes"  # species
        assert levels[1:6] == [""] * 5

    def test_bare_header_no_taxonomy(self):
        line = ">RF00001.1"
        seq_id, original, levels = parse_header(line)
        assert seq_id == "RF00001.1"
        assert original == "RF00001.1"
        assert levels == [""] * N_TAX_LEVELS

    def test_partial_taxonomy_first_and_last_only(self):
        # <7 tokens: only domain (first) and species (last) filled
        line = ">SEQ1 Archaea;Euryarchaeota"
        seq_id, _, levels = parse_header(line)
        assert levels[0] == "Archaea"        # domain
        assert levels[6] == "Euryarchaeota"  # species
        assert levels[1:6] == [""] * 5
        assert len(levels) == N_TAX_LEVELS

    def test_deep_lineage_first_and_last_only(self):
        # >7 tokens (8 here): real species preserved as last, middle empty
        line = (">AF106036.1.3725 Eukaryota;Discoba;Discicristata;Euglenozoa;"
                "Euglenida;Aphagea;Distigma;Distigma proteus;size=1")
        seq_id, _, levels = parse_header(line)
        assert levels[0] == "Eukaryota"          # domain
        assert levels[6] == "Distigma proteus"   # species (last token, not truncated)
        assert levels[1:6] == [""] * 5

    def test_single_token_domain_only(self):
        assert map_taxonomy(["Eukaryota"]) == ["Eukaryota"] + [""] * 6

    def test_exactly_seven_fills_all(self):
        toks = ["Bacteria", "Pseudomonadota", "Alphaproteobacteria",
                "Rickettsiales", "Mitochondria", "Incertae Sedis", "Sphagnum angustifolium"]
        assert map_taxonomy(toks) == toks

    def test_leading_gt_stripped(self):
        line = ">SEQ1 Bacteria;Firmicutes;size=1"
        seq_id, _, _ = parse_header(line)
        assert not seq_id.startswith(">")

    def test_parse_tax_false_leaves_levels_empty(self):
        line = ">QKKF02033054.1/420106-419872 Laodelphax striatellus isolate; whole genome shotgun sequence;size=1"
        seq_id, original, levels = parse_header(line, parse_tax=False)
        assert seq_id == "QKKF02033054.1/420106-419872"
        assert original == "QKKF02033054.1/420106-419872 Laodelphax striatellus isolate; whole genome shotgun sequence;size=1"
        assert levels == [""] * N_TAX_LEVELS


class TestResolveDbVersion:
    def test_silva_ssu(self):
        assert resolve_db_version("silva ssu bacteria", "138.2", "138.2", "15.1") == "SILVA 138.2"

    def test_silva_lsu_uses_lsu_version(self):
        assert resolve_db_version("silva lsu archaea", "138.2", "138.1", "15.1") == "SILVA 138.1"

    def test_rfam(self):
        assert resolve_db_version("rfam 5s", "138.2", "138.2", "15.1") == "Rfam 15.1"

    def test_missing_version_empty(self):
        assert resolve_db_version("silva ssu bacteria", "", "", "") == ""
        assert resolve_db_version("rfam 5s", "138.2", "138.2", "") == ""


class TestWriteFamilyMap:
    HEADER = ("seq_id\trrna_family\tdatabase_version\tseq_length\toriginal\t"
              "domain\tphylum\tclass\torder\ttax_family\tgenus\tspecies")

    def test_output_columns(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(
            ">SEQ1 Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;L. acidophilus;size=2\n"
            "ACGT\n"
            ">SEQ2 Archaea;Euryarchaeota\n"
            "TTTT\n"
        )
        out = tmp_path / "family_map.tsv"
        total = write_family_map(str(out), [(str(fasta), "silva ssu bacteria")],
                                 silva_ssu_version="138.2")
        assert total == 2
        lines = out.read_text().splitlines()
        assert lines[0] == self.HEADER
        row1 = lines[1].split("\t")
        assert row1[0] == "SEQ1"
        assert row1[1] == "silva ssu bacteria"
        assert row1[2] == "SILVA 138.2"
        assert row1[3] == "4"  # ACGT
        assert row1[4] == "SEQ1 Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;L. acidophilus;size=2"
        assert row1[5] == "Bacteria"          # domain
        assert row1[9] == "Lactobacillaceae"  # tax_family
        row2 = lines[2].split("\t")
        assert row2[0] == "SEQ2"
        assert row2[5] == "Archaea"  # domain
        assert row2[7] == ""         # class empty

    def test_seq_length_multiline(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">SEQ1 Bacteria;Firmicutes\nACGTACGT\nACG\n>SEQ2 Archaea\nAC\n")
        out = tmp_path / "map.tsv"
        write_family_map(str(out), [(str(fasta), "silva ssu bacteria")], silva_ssu_version="138.2")
        lines = out.read_text().splitlines()
        assert lines[1].split("\t")[3] == "11"  # 8 + 3
        assert lines[2].split("\t")[3] == "2"

    def test_multiple_fasta_files(self, tmp_path):
        f1 = tmp_path / "ssu.fasta"
        f1.write_text(">SSU1 Bacteria;Firmicutes\nACGT\n")
        f2 = tmp_path / "lsu.fasta"
        f2.write_text(">LSU1 Bacteria;Proteobacteria\nTTTT\n")
        out = tmp_path / "map.tsv"
        total = write_family_map(str(out), [
            (str(f1), "silva ssu bacteria"),
            (str(f2), "silva lsu bacteria"),
        ])
        assert total == 2
        lines = out.read_text().splitlines()
        assert lines[1].split("\t")[1] == "silva ssu bacteria"
        assert lines[2].split("\t")[1] == "silva lsu bacteria"

    def test_rfam_taxonomy_columns_empty(self, tmp_path):
        fasta = tmp_path / "rfam.fasta"
        fasta.write_text(
            ">QKKF02033054.1/1-100 Laodelphax striatellus, whole genome shotgun sequence;size=1\n"
            "ACGT\n"
        )
        out = tmp_path / "map.tsv"
        write_family_map(str(out), [(str(fasta), "rfam 5s")], rfam_version="15.1")
        row = out.read_text().splitlines()[1].split("\t")
        assert row[0] == "QKKF02033054.1/1-100"
        assert row[1] == "rfam 5s"
        assert row[2] == "Rfam 15.1"
        assert row[3] == "4"
        assert row[4] == "QKKF02033054.1/1-100 Laodelphax striatellus, whole genome shotgun sequence;size=1"
        assert row[5:] == [""] * N_TAX_LEVELS  # all 7 taxonomy columns empty

    def test_non_header_lines_skipped(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">SEQ1 Bacteria\nACGT\nACGT\n>SEQ2 Archaea\nTTTT\n")
        out = tmp_path / "map.tsv"
        total = write_family_map(str(out), [(str(fasta), "SSU")])
        assert total == 2
