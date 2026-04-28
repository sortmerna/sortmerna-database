import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from parse_cmsearch import trim_to_hit, parse_tblout, read_fasta


# ── trim_to_hit ───────────────────────────────────────────────────────────────

class TestTrimToHit:
  def test_no_trim_full_sequence(self):
    assert trim_to_hit("ACGT", 1, 4) == "ACGT"

  def test_5prime_trim(self):
    assert trim_to_hit("NNNACGT", 4, 7) == "ACGT"

  def test_3prime_trim(self):
    assert trim_to_hit("ACGTNNN", 1, 4) == "ACGT"

  def test_both_ends_trimmed(self):
    assert trim_to_hit("NNNACGTNNN", 4, 7) == "ACGT"

  def test_single_base(self):
    assert trim_to_hit("ACGT", 2, 2) == "C"

  def test_seq_from_1_seq_to_1(self):
    assert trim_to_hit("ACGT", 1, 1) == "A"

  def test_coordinates_are_1based_inclusive(self):
    seq = "ABCDE"
    assert trim_to_hit(seq, 2, 4) == "BCD"


# ── parse_tblout ──────────────────────────────────────────────────────────────

TBLOUT_HEADER = "#target name accession query name accession mdl mdl from mdl to seq from seq to strand trunc pass gc bias score E-value inc description of target\n"

def make_tblout_line(seq_id, seq_from, seq_to, score, evalue, inc="!"):
  # fields[0]=seq_id [7]=seq_from [8]=seq_to [9]=strand [14]=score [15]=evalue [16]=inc
  return f"{seq_id} - query - hmm 1 100 {seq_from} {seq_to} + - 1 0.5 0.0 {score} {evalue} {inc} description\n"


class TestParseTblout:
  def test_single_hit(self, tmp_path):
    f = tmp_path / "hits.tblout"
    f.write_text(TBLOUT_HEADER + make_tblout_line("seq1", 1, 500, 900.0, "0"))
    hits = parse_tblout(str(f))
    assert "seq1" in hits
    assert hits["seq1"] == (900.0, "0", 1, 500, "+")

  def test_best_hit_kept_for_duplicate_seq_id(self, tmp_path):
    f = tmp_path / "hits.tblout"
    f.write_text(
      TBLOUT_HEADER +
      make_tblout_line("seq1", 347, 1576, 990.5, "6.1e-296") +
      make_tblout_line("seq1", 1,   180,  61.7, "2.4e-15")
    )
    hits = parse_tblout(str(f))
    assert len(hits) == 1
    assert hits["seq1"][0] == 990.5
    assert hits["seq1"][2] == 347

  def test_non_inclusion_hit_ignored(self, tmp_path):
    f = tmp_path / "hits.tblout"
    f.write_text(TBLOUT_HEADER + make_tblout_line("seq1", 1, 500, 30.0, "0.5", inc="?"))
    hits = parse_tblout(str(f))
    assert hits == {}

  def test_comment_lines_skipped(self, tmp_path):
    f = tmp_path / "hits.tblout"
    f.write_text("# this is a comment\n# another comment\n")
    assert parse_tblout(str(f)) == {}

  def test_short_lines_skipped(self, tmp_path):
    f = tmp_path / "hits.tblout"
    f.write_text("seq1 incomplete line\n")
    assert parse_tblout(str(f)) == {}

  def test_empty_file(self, tmp_path):
    f = tmp_path / "empty.tblout"
    f.write_text("")
    assert parse_tblout(str(f)) == {}

  def test_multiple_sequences(self, tmp_path):
    f = tmp_path / "hits.tblout"
    f.write_text(
      TBLOUT_HEADER +
      make_tblout_line("seq1", 1,   1500, 1100.0, "0") +
      make_tblout_line("seq2", 540, 1697,  931.6, "3.9e-278")
    )
    hits = parse_tblout(str(f))
    assert set(hits.keys()) == {"seq1", "seq2"}
    assert hits["seq2"][2] == 540


# ── read_fasta ────────────────────────────────────────────────────────────────

class TestReadFasta:
  def test_single_sequence(self, tmp_path):
    f = tmp_path / "seqs.fasta"
    f.write_text(">seq1 some description\nACGT\n")
    records = list(read_fasta(str(f)))
    assert records == [("seq1 some description", "ACGT")]

  def test_multiple_sequences(self, tmp_path):
    f = tmp_path / "seqs.fasta"
    f.write_text(">seq1\nACGT\n>seq2\nTTTT\n")
    records = list(read_fasta(str(f)))
    assert len(records) == 2
    assert records[0] == ("seq1", "ACGT")
    assert records[1] == ("seq2", "TTTT")

  def test_multiline_sequence(self, tmp_path):
    f = tmp_path / "seqs.fasta"
    f.write_text(">seq1\nACGT\nACGT\nACGT\n")
    records = list(read_fasta(str(f)))
    assert records == [("seq1", "ACGTACGTACGT")]

  def test_empty_file(self, tmp_path):
    f = tmp_path / "empty.fasta"
    f.write_text("")
    assert list(read_fasta(str(f))) == []

  def test_header_without_description(self, tmp_path):
    f = tmp_path / "seqs.fasta"
    f.write_text(">seq1\nACGT\n")
    records = list(read_fasta(str(f)))
    assert records[0][0] == "seq1"
