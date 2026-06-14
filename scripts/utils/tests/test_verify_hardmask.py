import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from verify_hardmask import check, iter_records


def write_fasta(path, records):
    """Write [(id, seq)] to a plain or .gz FASTA file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as f:
        for seq_id, seq in records:
            f.write(f">{seq_id}\n{seq}\n")


class TestIterRecords:
    def test_plain_fasta(self, tmp_path):
        p = tmp_path / "in.fasta"
        write_fasta(p, [("seq1", "ACGTacgt")])
        records = list(iter_records(str(p)))
        assert records == [("seq1", "ACGTacgt")]

    def test_gzipped_fasta(self, tmp_path):
        p = tmp_path / "in.fasta.gz"
        write_fasta(p, [("seq1", "ACGTacgt")])
        records = list(iter_records(str(p)))
        assert records == [("seq1", "ACGTacgt")]

    def test_preserves_lowercase(self, tmp_path):
        p = tmp_path / "in.fasta"
        write_fasta(p, [("seq1", "aaaaTTTT")])
        records = list(iter_records(str(p)))
        assert records[0][1] == "aaaaTTTT"

    def test_id_stops_at_first_space(self, tmp_path):
        p = tmp_path / "in.fasta"
        p.write_text(">seq1 some description\nACGT\n")
        records = list(iter_records(str(p)))
        assert records[0][0] == "seq1"

    def test_multiline_sequence_joined(self, tmp_path):
        p = tmp_path / "in.fasta"
        p.write_text(">seq1\nACGT\nacgt\n")
        records = list(iter_records(str(p)))
        assert records[0][1] == "ACGTacgt"

    def test_multiple_records(self, tmp_path):
        p = tmp_path / "in.fasta"
        write_fasta(p, [("s1", "ACGT"), ("s2", "acgt")])
        records = list(iter_records(str(p)))
        assert len(records) == 2
        assert records[1] == ("s2", "acgt")

    def test_empty_fasta(self, tmp_path):
        p = tmp_path / "in.fasta"
        p.write_text("")
        assert list(iter_records(str(p))) == []


class TestCheck:
    def test_clean_pass(self, tmp_path):
        soft = tmp_path / "soft.fasta"
        hard = tmp_path / "hard.fasta"
        write_fasta(soft, [("s1", "ACGTacgtACGT")])
        write_fasta(hard, [("s1", "ACGTNNNNACGT")])
        assert check(str(soft), str(hard)) == 0

    def test_no_masking_needed(self, tmp_path):
        soft = tmp_path / "soft.fasta"
        hard = tmp_path / "hard.fasta"
        write_fasta(soft, [("s1", "ACGTACGT")])
        write_fasta(hard, [("s1", "ACGTACGT")])
        assert check(str(soft), str(hard)) == 0

    def test_all_lowercase_masked(self, tmp_path):
        soft = tmp_path / "soft.fasta"
        hard = tmp_path / "hard.fasta"
        write_fasta(soft, [("s1", "acgtacgt")])
        write_fasta(hard, [("s1", "NNNNNNNN")])
        assert check(str(soft), str(hard)) == 0

    def test_lowercase_not_masked_in_hard(self, tmp_path):
        soft = tmp_path / "soft.fasta"
        hard = tmp_path / "hard.fasta"
        write_fasta(soft, [("s1", "ACGTacgt")])
        write_fasta(hard, [("s1", "ACGTacgt")])  # lowercase not converted
        assert check(str(soft), str(hard)) == 4

    def test_uppercase_wrongly_masked(self, tmp_path):
        # Mirrors the real bug: vsearch masking uppercase A -> N
        soft = tmp_path / "soft.fasta"
        hard = tmp_path / "hard.fasta"
        write_fasta(soft, [("s1", "AAAAaaaa")])
        write_fasta(hard, [("s1", "NNNNNNNN")])  # first 4 A's wrongly masked
        assert check(str(soft), str(hard)) == 4

    def test_id_mismatch(self, tmp_path):
        soft = tmp_path / "soft.fasta"
        hard = tmp_path / "hard.fasta"
        write_fasta(soft, [("s1", "ACGT")])
        write_fasta(hard, [("s2", "ACGT")])
        assert check(str(soft), str(hard)) == 1

    def test_length_mismatch(self, tmp_path):
        soft = tmp_path / "soft.fasta"
        hard = tmp_path / "hard.fasta"
        write_fasta(soft, [("s1", "ACGT")])
        write_fasta(hard, [("s1", "ACG")])
        assert check(str(soft), str(hard)) == 1

    def test_sequence_count_mismatch(self, tmp_path):
        soft = tmp_path / "soft.fasta"
        hard = tmp_path / "hard.fasta"
        write_fasta(soft, [("s1", "ACGT"), ("s2", "ACGT")])
        write_fasta(hard, [("s1", "ACGT")])
        assert check(str(soft), str(hard)) >= 1

    def test_gzipped_files(self, tmp_path):
        soft = tmp_path / "soft.fasta.gz"
        hard = tmp_path / "hard.fasta.gz"
        write_fasta(soft, [("s1", "ACGTacgtACGT")])
        write_fasta(hard, [("s1", "ACGTNNNNACGT")])
        assert check(str(soft), str(hard)) == 0

    def test_existing_N_in_soft_preserved(self, tmp_path):
        # N's already in the soft-masked file (from SILVA) must pass through unchanged
        soft = tmp_path / "soft.fasta"
        hard = tmp_path / "hard.fasta"
        write_fasta(soft, [("s1", "ACNTacNT")])
        write_fasta(hard, [("s1", "ACNTNNNT")])
        assert check(str(soft), str(hard)) == 0

    def test_multiple_sequences_all_pass(self, tmp_path):
        soft = tmp_path / "soft.fasta"
        hard = tmp_path / "hard.fasta"
        write_fasta(soft, [("s1", "ACGTacgt"), ("s2", "GGGGcccc")])
        write_fasta(hard, [("s1", "ACGTNNNN"), ("s2", "GGGGNNNN")])
        assert check(str(soft), str(hard)) == 0

    def test_multiple_sequences_one_fails(self, tmp_path):
        soft = tmp_path / "soft.fasta"
        hard = tmp_path / "hard.fasta"
        write_fasta(soft, [("s1", "ACGTacgt"), ("s2", "GGGG")])
        write_fasta(hard, [("s1", "ACGTNNNN"), ("s2", "NNNN")])  # s2 wrongly masked
        assert check(str(soft), str(hard)) == 4
