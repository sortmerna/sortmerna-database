import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from build_family_map import parse_header, write_family_map, N_TAX_LEVELS


class TestParseHeader:
    def test_full_silva_header(self):
        line = ">KU725478.45631.48554 Bacteria;Pseudomonadota;Alphaproteobacteria;Rickettsiales;Mitochondria;Incertae Sedis;Sphagnum angustifolium;size=5"
        seq_id, levels = parse_header(line)
        assert seq_id == "KU725478.45631.48554"
        assert levels[0] == "Bacteria"
        assert levels[1] == "Pseudomonadota"
        assert levels[2] == "Alphaproteobacteria"
        assert levels[3] == "Rickettsiales"
        assert levels[4] == "Mitochondria"
        assert levels[5] == "Incertae Sedis"
        assert levels[6] == "Sphagnum angustifolium"
        assert len(levels) == N_TAX_LEVELS

    def test_size_tag_stripped(self):
        line = ">AB123.1.500 Bacteria;Firmicutes;size=10"
        seq_id, levels = parse_header(line)
        assert levels[0] == "Bacteria"
        assert levels[1] == "Firmicutes"
        assert levels[2] == ""

    def test_bare_header_no_taxonomy(self):
        line = ">RF00001.1"
        seq_id, levels = parse_header(line)
        assert seq_id == "RF00001.1"
        assert levels == [""] * N_TAX_LEVELS

    def test_partial_taxonomy_padded(self):
        line = ">SEQ1 Archaea;Euryarchaeota"
        seq_id, levels = parse_header(line)
        assert levels[0] == "Archaea"
        assert levels[1] == "Euryarchaeota"
        assert levels[2:] == [""] * (N_TAX_LEVELS - 2)
        assert len(levels) == N_TAX_LEVELS

    def test_leading_gt_stripped(self):
        line = ">SEQ1 Bacteria;Firmicutes;size=1"
        seq_id, _ = parse_header(line)
        assert not seq_id.startswith(">")


class TestWriteFamilyMap:
    def test_output_columns(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(
            ">SEQ1 Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;L. acidophilus;size=2\n"
            "ACGT\n"
            ">SEQ2 Archaea;Euryarchaeota\n"
            "TTTT\n"
        )
        out = tmp_path / "family_map.tsv"
        total = write_family_map(str(out), [(str(fasta), "SSU Bacteria (16S)")])
        assert total == 2
        lines = out.read_text().splitlines()
        assert lines[0] == "seq_id\trna_family\tdomain\tphylum\tclass\torder\ttax_family\tgenus\tspecies"
        row1 = lines[1].split("\t")
        assert row1[0] == "SEQ1"
        assert row1[1] == "SSU Bacteria (16S)"
        assert row1[2] == "Bacteria"
        assert row1[6] == "Lactobacillaceae"
        row2 = lines[2].split("\t")
        assert row2[0] == "SEQ2"
        assert row2[2] == "Archaea"
        assert row2[4] == ""  # class empty

    def test_multiple_fasta_files(self, tmp_path):
        f1 = tmp_path / "ssu.fasta"
        f1.write_text(">SSU1 Bacteria;Firmicutes\nACGT\n")
        f2 = tmp_path / "lsu.fasta"
        f2.write_text(">LSU1 Bacteria;Proteobacteria\nTTTT\n")
        out = tmp_path / "map.tsv"
        total = write_family_map(str(out), [
            (str(f1), "SSU Bacteria (16S)"),
            (str(f2), "LSU Bacteria (23S)"),
        ])
        assert total == 2
        lines = out.read_text().splitlines()
        assert lines[1].split("\t")[1] == "SSU Bacteria (16S)"
        assert lines[2].split("\t")[1] == "LSU Bacteria (23S)"

    def test_non_header_lines_skipped(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">SEQ1 Bacteria\nACGT\nACGT\n>SEQ2 Archaea\nTTTT\n")
        out = tmp_path / "map.tsv"
        total = write_family_map(str(out), [(str(fasta), "SSU")])
        assert total == 2
