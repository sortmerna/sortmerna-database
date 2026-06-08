import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from hard_mask_fasta import hard_mask_sequence, hard_mask_fasta, MaskingStats


class TestHardMaskSequence:
    def test_all_uppercase_unchanged(self):
        seq, n = hard_mask_sequence("ACGTACGT")
        assert seq == "ACGTACGT"
        assert n == 0

    def test_all_lowercase_to_N(self):
        seq, n = hard_mask_sequence("acgtacgt")
        assert seq == "NNNNNNNN"
        assert n == 8

    def test_mixed_case(self):
        seq, n = hard_mask_sequence("ACGTacgt")
        assert seq == "ACGTNNN N".replace(" ", "")
        assert seq == "ACGTNNNN"
        assert n == 4

    def test_existing_N_preserved(self):
        seq, n = hard_mask_sequence("ACGTNNacgt")
        assert seq == "ACGTNNNNNN"
        assert n == 4  # only the lowercase ones

    def test_empty_sequence(self):
        seq, n = hard_mask_sequence("")
        assert seq == ""
        assert n == 0

    def test_all_N_no_new_masking(self):
        seq, n = hard_mask_sequence("NNNNNN")
        assert seq == "NNNNNN"
        assert n == 0

    def test_iupac_uppercase_preserved(self):
        # Uppercase IUPAC ambiguity codes (R,Y,W,S,...) are NOT soft masking
        seq, n = hard_mask_sequence("ACGTRYWSMK")
        assert seq == "ACGTRYWSMK"
        assert n == 0

    def test_iupac_lowercase_masked(self):
        # Lowercase IUPAC codes (r,y,w,...) ARE soft masking → N
        seq, n = hard_mask_sequence("ACGTrywsmk")
        assert seq == "ACGTNNNNN N".replace(" ", "")
        assert seq == "ACGTNNNNNN"
        assert n == 6

    def test_tandem_repeat_pattern(self):
        # AT-microsatellite as seen in SILVA soft-masked sequences
        # "atatatatat at" = 12 lowercase bases
        seq, n = hard_mask_sequence("GCATatatatatatatGCAT")
        assert seq == "GCATNNNNNNNNNNNNGCAT"
        assert n == 12
        assert seq.startswith("GCAT")
        assert seq.endswith("GCAT")
        assert "N" * 12 in seq


class TestHardMaskFasta:
    def test_no_masking_needed(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text(">seq1\nACGTACGT\n>seq2\nGGGGCCCC\n")
        out = tmp_path / "out.fasta"
        stats = hard_mask_fasta(str(fasta), str(out))
        assert stats.total_sequences == 2
        assert stats.masked_sequences == 0
        assert stats.masked_bases == 0
        assert "ACGTACGT" in out.read_text()

    def test_full_masking(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text(">seq1\nacgtacgt\n")
        out = tmp_path / "out.fasta"
        stats = hard_mask_fasta(str(fasta), str(out))
        assert stats.total_sequences == 1
        assert stats.masked_sequences == 1
        assert stats.masked_bases == 8
        assert "NNNNNNNN" in out.read_text()

    def test_partial_masking(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text(">seq1 header preserved\nACGTatatatGCAT\n")
        out = tmp_path / "out.fasta"
        stats = hard_mask_fasta(str(fasta), str(out))
        assert stats.total_sequences == 1
        assert stats.masked_sequences == 1
        assert stats.masked_bases == 6
        content = out.read_text()
        assert ">seq1 header preserved" in content
        assert "ACGTNNNNNNGCAT" in content

    def test_header_preserved_unchanged(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        header = ">AB123.1.2654 Bacteria;Firmicutes;size=5"
        fasta.write_text(f"{header}\nACGTatGCAT\n")
        out = tmp_path / "out.fasta"
        hard_mask_fasta(str(fasta), str(out))
        assert header in out.read_text()

    def test_multiline_sequence(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text(">seq1\nACGTacgt\nGGGGcccc\n")
        out = tmp_path / "out.fasta"
        stats = hard_mask_fasta(str(fasta), str(out))
        assert stats.total_sequences == 1
        assert stats.masked_bases == 8
        content = out.read_text()
        assert "ACGTNNNNGGGGNNN" in content.replace("\n", "")

    def test_gzipped_input_output(self, tmp_path):
        fasta_gz = tmp_path / "in.fasta.gz"
        with gzip.open(fasta_gz, "wt") as f:
            f.write(">seq1\nACGTatGCAT\n")
        out_gz = tmp_path / "out.fasta.gz"
        stats = hard_mask_fasta(str(fasta_gz), str(out_gz))
        assert stats.masked_bases == 2
        with gzip.open(out_gz, "rt") as f:
            content = f.read()
        assert "ACGTNNNGCAT".replace("NNN", "NN") in content.replace("\n", "")
        assert "ACGTNNGCAT" in content.replace("\n", "")

    def test_multiple_sequences_stats(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text(
            ">seq1\nACGT\n"       # no masking
            ">seq2\nacgt\n"       # fully masked
            ">seq3\nACGTacgt\n"  # partial
        )
        out = tmp_path / "out.fasta"
        stats = hard_mask_fasta(str(fasta), str(out))
        assert stats.total_sequences == 3
        assert stats.masked_sequences == 2
        assert stats.total_bases == 16
        assert stats.masked_bases == 8

    def test_windows_line_endings(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        fasta.write_bytes(b">seq1\r\nACGTatGCAT\r\n>seq2\r\nGGGG\r\n")
        out = tmp_path / "out.fasta"
        stats = hard_mask_fasta(str(fasta), str(out))
        assert stats.total_sequences == 2
        assert stats.masked_bases == 2
        assert "ACGTNNGCAT" in out.read_text().replace("\n", "")

    def test_trailing_whitespace_in_sequence(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text(">seq1\nACGTatGCAT   \n")
        out = tmp_path / "out.fasta"
        stats = hard_mask_fasta(str(fasta), str(out))
        assert stats.total_bases == 10  # trailing spaces stripped
        assert stats.masked_bases == 2

    def test_empty_fasta(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text("")
        out = tmp_path / "out.fasta"
        stats = hard_mask_fasta(str(fasta), str(out))
        assert stats.total_sequences == 0
        assert stats.masked_bases == 0
