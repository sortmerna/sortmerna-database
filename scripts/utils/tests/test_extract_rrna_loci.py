import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from extract_rrna_loci import extract_rrna_loci

SCRIPT = Path(__file__).parent.parent / "extract_rrna_loci.py"


def write_gbff(path, content):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as f:
        f.write(content)


def read_bed(path):
    lines = Path(path).read_text().splitlines()
    return [tuple(l.split("\t")) for l in lines if l]


# ── extract_rrna_loci ─────────────────────────────────────────────────────────

class TestExtractRrnaLoci:
    def test_simple_forward(self, tmp_path):
        gbff = tmp_path / "test.gbff"
        bed  = tmp_path / "out.bed"
        write_gbff(gbff, (
            "LOCUS       NC_000001\n"
            "FEATURES             Location/Qualifiers\n"
            "     rRNA            100..500\n"
            "                     /product=\"18S ribosomal RNA\"\n"
            "ORIGIN\n"
            "//\n"
        ))
        count = extract_rrna_loci(gbff, bed)
        assert count == 1
        assert read_bed(bed) == [("NC_000001", "99", "500")]

    def test_complement_location(self, tmp_path):
        gbff = tmp_path / "test.gbff"
        bed  = tmp_path / "out.bed"
        write_gbff(gbff, (
            "LOCUS       NC_000002\n"
            "FEATURES             Location/Qualifiers\n"
            "     rRNA            complement(200..800)\n"
            "                     /product=\"23S ribosomal RNA\"\n"
            "ORIGIN\n"
            "//\n"
        ))
        count = extract_rrna_loci(gbff, bed)
        assert count == 1
        assert read_bed(bed) == [("NC_000002", "199", "800")]

    def test_multiline_join(self, tmp_path):
        gbff = tmp_path / "test.gbff"
        bed  = tmp_path / "out.bed"
        write_gbff(gbff, (
            "LOCUS       NC_000003\n"
            "FEATURES             Location/Qualifiers\n"
            "     rRNA            join(100..300,\n"
            "                     400..600)\n"
            "                     /product=\"rRNA\"\n"
            "ORIGIN\n"
            "//\n"
        ))
        count = extract_rrna_loci(gbff, bed)
        assert count == 2
        assert read_bed(bed) == [("NC_000003", "99", "300"), ("NC_000003", "399", "600")]

    def test_non_rrna_features_ignored(self, tmp_path):
        gbff = tmp_path / "test.gbff"
        bed  = tmp_path / "out.bed"
        write_gbff(gbff, (
            "LOCUS       NC_000004\n"
            "FEATURES             Location/Qualifiers\n"
            "     gene            1..100\n"
            "                     /gene=\"ACTB\"\n"
            "     CDS             1..100\n"
            "                     /product=\"actin\"\n"
            "ORIGIN\n"
            "//\n"
        ))
        count = extract_rrna_loci(gbff, bed)
        assert count == 0
        assert read_bed(bed) == []

    def test_multiple_chromosomes(self, tmp_path):
        gbff = tmp_path / "test.gbff"
        bed  = tmp_path / "out.bed"
        write_gbff(gbff, (
            "LOCUS       NC_000001\n"
            "FEATURES             Location/Qualifiers\n"
            "     rRNA            10..50\n"
            "ORIGIN\n"
            "//\n"
            "LOCUS       NC_000002\n"
            "FEATURES             Location/Qualifiers\n"
            "     rRNA            200..400\n"
            "ORIGIN\n"
            "//\n"
        ))
        count = extract_rrna_loci(gbff, bed)
        assert count == 2
        assert read_bed(bed) == [("NC_000001", "9", "50"), ("NC_000002", "199", "400")]

    def test_rrna_at_end_of_file(self, tmp_path):
        gbff = tmp_path / "test.gbff"
        bed  = tmp_path / "out.bed"
        write_gbff(gbff, (
            "LOCUS       NC_000001\n"
            "FEATURES             Location/Qualifiers\n"
            "     rRNA            1..100\n"
        ))
        count = extract_rrna_loci(gbff, bed)
        assert count == 1
        assert read_bed(bed) == [("NC_000001", "0", "100")]

    def test_gzipped_input(self, tmp_path):
        gbff = tmp_path / "test.gbff.gz"
        bed  = tmp_path / "out.bed"
        write_gbff(gbff, (
            "LOCUS       NC_000001\n"
            "FEATURES             Location/Qualifiers\n"
            "     rRNA            1..100\n"
            "ORIGIN\n"
            "//\n"
        ))
        count = extract_rrna_loci(gbff, bed)
        assert count == 1
        assert read_bed(bed) == [("NC_000001", "0", "100")]

    def test_margin_basic(self, tmp_path):
        gbff = tmp_path / "test.gbff"
        bed  = tmp_path / "out.bed"
        write_gbff(gbff, (
            "LOCUS       NC_000001              10000 bp\n"
            "FEATURES             Location/Qualifiers\n"
            "     rRNA            500..1000\n"
            "ORIGIN\n"
            "//\n"
        ))
        count = extract_rrna_loci(gbff, bed, margin=100)
        assert count == 1
        assert read_bed(bed) == [("NC_000001", "399", "1100")]

    def test_margin_clamps_start_at_zero(self, tmp_path):
        gbff = tmp_path / "test.gbff"
        bed  = tmp_path / "out.bed"
        write_gbff(gbff, (
            "LOCUS       NC_000001              10000 bp\n"
            "FEATURES             Location/Qualifiers\n"
            "     rRNA            50..200\n"
            "ORIGIN\n"
            "//\n"
        ))
        count = extract_rrna_loci(gbff, bed, margin=100)
        assert count == 1
        assert read_bed(bed) == [("NC_000001", "0", "300")]

    def test_margin_clamps_end_at_chrom_length(self, tmp_path):
        gbff = tmp_path / "test.gbff"
        bed  = tmp_path / "out.bed"
        write_gbff(gbff, (
            "LOCUS       NC_000001              1000 bp\n"
            "FEATURES             Location/Qualifiers\n"
            "     rRNA            800..950\n"
            "ORIGIN\n"
            "//\n"
        ))
        count = extract_rrna_loci(gbff, bed, margin=100)
        assert count == 1
        assert read_bed(bed) == [("NC_000001", "699", "1000")]

    def test_negative_margin_raises(self, tmp_path):
        gbff = tmp_path / "test.gbff"
        bed  = tmp_path / "out.bed"
        write_gbff(gbff, "LOCUS       NC_000001\n")
        with pytest.raises(ValueError, match="margin must be >= 0"):
            extract_rrna_loci(gbff, bed, margin=-1)

    def test_fuzzy_coordinates(self, tmp_path):
        gbff = tmp_path / "test.gbff"
        bed  = tmp_path / "out.bed"
        write_gbff(gbff, (
            "LOCUS       NC_000001\n"
            "FEATURES             Location/Qualifiers\n"
            "     rRNA            <1..500\n"
            "                     /product=\"18S\"\n"
            "     rRNA            100..>900\n"
            "                     /product=\"28S\"\n"
            "ORIGIN\n"
            "//\n"
        ))
        count = extract_rrna_loci(gbff, bed)
        assert count == 2
        assert read_bed(bed) == [("NC_000001", "0", "500"), ("NC_000001", "99", "900")]

    def test_multiple_rrna_features_same_chrom(self, tmp_path):
        gbff = tmp_path / "test.gbff"
        bed  = tmp_path / "out.bed"
        write_gbff(gbff, (
            "LOCUS       NC_000001\n"
            "FEATURES             Location/Qualifiers\n"
            "     rRNA            100..500\n"
            "                     /product=\"18S\"\n"
            "     rRNA            1000..3000\n"
            "                     /product=\"28S\"\n"
            "ORIGIN\n"
            "//\n"
        ))
        count = extract_rrna_loci(gbff, bed)
        assert count == 2
        assert read_bed(bed) == [("NC_000001", "99", "500"), ("NC_000001", "999", "3000")]
