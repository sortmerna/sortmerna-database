import pytest
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from parse_cmsearch import trim_to_hit, parse_tblout, read_fasta

TEST_DIR = Path(__file__).parent / "data"
SCRIPT   = Path(__file__).parent.parent / "parse_cmsearch.py"


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
    assert hits["seq1"] == (900.0, "0", 1, 500, "+", 1, 1, 500)

  def test_multi_hit_merged_span_and_best_coords_stored(self, tmp_path):
    # Two hits: best has score 990.5 at [347,1576]; second at [1,180]
    # merged span = [1, 1576]; best coords = [347, 1576]
    f = tmp_path / "hits.tblout"
    f.write_text(
      TBLOUT_HEADER +
      make_tblout_line("seq1", 347, 1576, 990.5, "6.1e-296") +
      make_tblout_line("seq1", 1,   180,  61.7, "2.4e-15")
    )
    hits = parse_tblout(str(f))
    assert len(hits) == 1
    score, evalue, merged_from, merged_to, strand, n_hits, best_sf, best_st = hits["seq1"]
    assert score      == 990.5
    assert n_hits     == 2
    assert merged_from == 1      # min(347, 1)
    assert merged_to   == 1576   # max(1576, 180)
    assert best_sf     == 347    # best hit's own sf
    assert best_st     == 1576   # best hit's own st

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
    assert hits["seq2"][2] == 540   # single hit — merged_from == sf


# ── read_fasta ────────────────────────────────────────────────────────────────

class TestReadFasta:
  def test_single_sequence(self, tmp_path):
    f = tmp_path / "seqs.fasta"
    f.write_text(">seq1 some description\nACGT\n")
    records = list(read_fasta(str(f)))
    assert records[0][0] == "seq1 some description"
    assert str(records[0][1]) == "ACGT"

  def test_multiple_sequences(self, tmp_path):
    f = tmp_path / "seqs.fasta"
    f.write_text(">seq1\nACGT\n>seq2\nTTTT\n")
    records = list(read_fasta(str(f)))
    assert len(records) == 2
    assert records[0][0] == "seq1"
    assert str(records[0][1]) == "ACGT"
    assert records[1][0] == "seq2"
    assert str(records[1][1]) == "TTTT"

  def test_multiline_sequence(self, tmp_path):
    f = tmp_path / "seqs.fasta"
    f.write_text(">seq1\nACGT\nACGT\nACGT\n")
    records = list(read_fasta(str(f)))
    assert len(records) == 1
    assert records[0][0] == "seq1"
    assert str(records[0][1]) == "ACGTACGTACGT"

  @pytest.mark.filterwarnings("ignore::skbio.io.registry.FormatIdentificationWarning")
  def test_empty_file(self, tmp_path):
    f = tmp_path / "empty.fasta"
    f.write_text("")
    assert list(read_fasta(str(f))) == []

  def test_header_without_description(self, tmp_path):
    f = tmp_path / "seqs.fasta"
    f.write_text(">seq1\nACGT\n")
    records = list(read_fasta(str(f)))
    assert records[0][0] == "seq1"
    assert str(records[0][1]) == "ACGT"


# ── helpers for end-to-end tests ─────────────────────────────────────────────

def run_script(tmp_path, fasta_text, tblout_text):
  fasta   = tmp_path / "in.fasta"
  tblout  = tmp_path / "in.tblout"
  clean   = tmp_path / "clean.fasta"
  flagged = tmp_path / "flagged.fasta"
  log     = tmp_path / "log.tsv"
  fasta.write_text(fasta_text)
  tblout.write_text(tblout_text)
  result = subprocess.run(
    [sys.executable, str(SCRIPT),
     "--tblout",  str(tblout),
     "--fasta",   str(fasta),
     "--clean",   str(clean),
     "--flagged", str(flagged),
     "--log",     str(log)],
    capture_output=True, text=True,
  )
  assert result.returncode == 0, result.stderr
  return clean, flagged, log

def tsv_rows(path):
  """Return TSV as list of stripped row-lists, skipping the header."""
  lines = path.read_text().splitlines()
  return [line.rstrip('\t').split('\t') for line in lines[1:] if line.strip()]


# ── end-to-end integration tests ─────────────────────────────────────────────

class TestEndToEnd:
  """Run parse_cmsearch.py on the test tblout + FASTA subset and verify outputs."""

  def test_verified_fasta_and_log_match_expected(self, tmp_path):
    tblout_fp      = TEST_DIR / "cmsearch_log_ssu_bacteria_test.tblout"
    fasta_fp       = TEST_DIR / "silva_ssu_dom_bacteria_test.fasta"
    expected_tsv   = TEST_DIR / "cmsearch_log_ssu_bacteria_test.tsv"
    expected_fasta = TEST_DIR / "verified_ssu_bacteria_test.fasta"

    out_fasta   = tmp_path / "verified_ssu_bacteria_test_local.fasta"
    out_flagged = tmp_path / "flagged_local.fasta"
    out_log     = tmp_path / "log_local.tsv"

    result = subprocess.run(
      [
        sys.executable, str(SCRIPT),
        "--tblout",  str(tblout_fp),
        "--fasta",   str(fasta_fp),
        "--clean",   str(out_fasta),
        "--flagged", str(out_flagged),
        "--log",     str(out_log),
      ],
      capture_output=True,
      text=True,
    )
    assert result.returncode == 0, result.stderr

    # Log TSV: compare rows ignoring trailing whitespace on note column
    assert tsv_rows(out_log) == tsv_rows(expected_tsv)

    # Verified FASTA: content match (header + sequence), ignoring line-wrapping differences
    out_records = [(hdr, str(seq)) for hdr, seq in read_fasta(str(out_fasta))]
    exp_records = [(hdr, str(seq)) for hdr, seq in read_fasta(str(expected_fasta))]
    assert out_records == exp_records

  def test_multi_hit_high_coverage_trims_to_merged_span(self, tmp_path):
    # seq1: 1000 bp; hits at [1,600] (best, score 900) and [550,900]
    # merged span [1,900], coverage = 900/1000 = 0.9 >= COVERAGE_THRESHOLD
    seq = "A" * 1000
    fasta_text  = f">seq1\n{seq}\n"
    tblout_text = (
      TBLOUT_HEADER +
      make_tblout_line("seq1", 1,   600, 900.0, "0") +
      make_tblout_line("seq1", 550, 900, 800.0, "1e-200")
    )
    clean, _, log = run_script(tmp_path, fasta_text, tblout_text)

    rows = tsv_rows(log)
    assert len(rows) == 1
    _, _, _, sf, st, _, n_hits, note = rows[0]
    assert sf     == "1"
    assert st     == "900"
    assert n_hits == "2"
    assert note   == "merged_2_hits"

    out_records = [(hdr, str(seq)) for hdr, seq in read_fasta(str(clean))]
    assert out_records[0][1] == "A" * 900

  def test_multi_hit_low_coverage_falls_back_to_best_hit(self, tmp_path):
    # seq1: 1000 bp; best hit at [100,450] (score 900), second at [700,800]
    # merged span [100,800], coverage = 701/1000 = 0.701 < COVERAGE_THRESHOLD
    seq = "A" * 1000
    fasta_text  = f">seq1\n{seq}\n"
    tblout_text = (
      TBLOUT_HEADER +
      make_tblout_line("seq1", 100, 450, 900.0, "0") +
      make_tblout_line("seq1", 700, 800, 800.0, "1e-200")
    )
    clean, _, log = run_script(tmp_path, fasta_text, tblout_text)

    rows = tsv_rows(log)
    assert len(rows) == 1
    _, _, _, sf, st, _, n_hits, note = rows[0]
    assert sf     == "100"   # best hit sf, not merged_from
    assert st     == "450"   # best hit st, not merged_to
    assert n_hits == "2"
    assert "low_coverage" in note
    assert "REVIEW"        in note

    out_records = [(hdr, str(seq)) for hdr, seq in read_fasta(str(clean))]
    assert out_records[0][1] == "A" * 351  # positions 100-450 inclusive
