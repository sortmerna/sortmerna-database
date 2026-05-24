import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from extract_rrna_loci import extract_rrna_loci, load_name_map


def write_gff(path, content):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as f:
        f.write(content)


def write_report(path, rows):
    """Write a minimal NCBI assembly report with given (refseq, genbank) pairs."""
    header = (
        "# Assembly name: test\n"
        "# comment\n"
        "# seq-name\trole\tmol\ttype\tgenbank\trel\trefseq\tunit\tlen\tucsc\n"
    )
    with open(path, "w") as f:
        f.write(header)
        for refseq, genbank in rows:
            f.write(f"x\tx\tx\tx\t{genbank}\tx\t{refseq}\tx\t0\tx\n")


def read_bed(path):
    lines = Path(path).read_text().splitlines()
    return [tuple(l.split("\t")) for l in lines if l]


def gff_row(seqname, start, end, feature="rRNA"):
    return f"{seqname}\tRefSeq\t{feature}\t{start}\t{end}\t.\t+\t.\tID=test\n"


# -- extract_rrna_loci --------------------------------------------------------

class TestExtractRrnaLoci:
    def test_simple(self, tmp_path):
        gff = tmp_path / "test.gff"
        bed = tmp_path / "out.bed"
        write_gff(gff, gff_row("chr1", 100, 500))
        count = extract_rrna_loci(gff, bed)
        assert count == 1
        assert read_bed(bed) == [("chr1", "99", "500")]

    def test_non_rrna_features_ignored(self, tmp_path):
        gff = tmp_path / "test.gff"
        bed = tmp_path / "out.bed"
        write_gff(gff,
            gff_row("chr1", 1, 100, feature="gene") +
            gff_row("chr1", 1, 100, feature="CDS") +
            gff_row("chr1", 1, 100, feature="tRNA")
        )
        count = extract_rrna_loci(gff, bed)
        assert count == 0
        assert read_bed(bed) == []

    def test_comment_lines_skipped(self, tmp_path):
        gff = tmp_path / "test.gff"
        bed = tmp_path / "out.bed"
        write_gff(gff,
            "##gff-version 3\n"
            "# this is a comment\n" +
            gff_row("chr1", 100, 500)
        )
        count = extract_rrna_loci(gff, bed)
        assert count == 1
        assert read_bed(bed) == [("chr1", "99", "500")]

    def test_multiple_chromosomes(self, tmp_path):
        gff = tmp_path / "test.gff"
        bed = tmp_path / "out.bed"
        write_gff(gff,
            gff_row("chr1", 10, 50) +
            gff_row("chr2", 200, 400)
        )
        count = extract_rrna_loci(gff, bed)
        assert count == 2
        assert read_bed(bed) == [("chr1", "9", "50"), ("chr2", "199", "400")]

    def test_multiple_rrna_same_chrom(self, tmp_path):
        gff = tmp_path / "test.gff"
        bed = tmp_path / "out.bed"
        write_gff(gff,
            gff_row("chr1", 100, 500) +
            gff_row("chr1", 1000, 3000)
        )
        count = extract_rrna_loci(gff, bed)
        assert count == 2
        assert read_bed(bed) == [("chr1", "99", "500"), ("chr1", "999", "3000")]

    def test_gzipped_input(self, tmp_path):
        gff = tmp_path / "test.gff.gz"
        bed = tmp_path / "out.bed"
        write_gff(gff, gff_row("chr1", 1, 100))
        count = extract_rrna_loci(gff, bed)
        assert count == 1
        assert read_bed(bed) == [("chr1", "0", "100")]

    def test_margin_basic(self, tmp_path):
        gff = tmp_path / "test.gff"
        bed = tmp_path / "out.bed"
        write_gff(gff, gff_row("chr1", 500, 1000))
        count = extract_rrna_loci(gff, bed, margin=100)
        assert count == 1
        assert read_bed(bed) == [("chr1", "399", "1100")]

    def test_margin_clamps_start_at_zero(self, tmp_path):
        gff = tmp_path / "test.gff"
        bed = tmp_path / "out.bed"
        write_gff(gff, gff_row("chr1", 50, 200))
        count = extract_rrna_loci(gff, bed, margin=100)
        assert count == 1
        assert read_bed(bed) == [("chr1", "0", "300")]

    def test_negative_margin_raises(self, tmp_path):
        gff = tmp_path / "test.gff"
        bed = tmp_path / "out.bed"
        write_gff(gff, "##gff-version 3\n")
        with pytest.raises(ValueError, match="margin must be >= 0"):
            extract_rrna_loci(gff, bed, margin=-1)

    def test_name_map_translates_refseq(self, tmp_path):
        gff    = tmp_path / "test.gff"
        bed    = tmp_path / "out.bed"
        report = tmp_path / "report.txt"
        write_gff(gff, gff_row("NC_060925.1", 100, 500))
        write_report(report, [("NC_060925.1", "chr1")])
        count = extract_rrna_loci(gff, bed, name_map_path=report)
        assert count == 1
        assert read_bed(bed) == [("chr1", "99", "500")]

    def test_name_map_unknown_seqname_kept(self, tmp_path):
        gff    = tmp_path / "test.gff"
        bed    = tmp_path / "out.bed"
        report = tmp_path / "report.txt"
        write_gff(gff, gff_row("NC_UNKNOWN.1", 100, 500))
        write_report(report, [("NC_060925.1", "chr1")])
        count = extract_rrna_loci(gff, bed, name_map_path=report)
        assert count == 1
        assert read_bed(bed) == [("NC_UNKNOWN.1", "99", "500")]

    def test_name_map_multiple_chroms(self, tmp_path):
        gff    = tmp_path / "test.gff"
        bed    = tmp_path / "out.bed"
        report = tmp_path / "report.txt"
        write_gff(gff,
            gff_row("NC_060925.1", 100, 500) +
            gff_row("NC_060926.1", 200, 800)
        )
        write_report(report, [("NC_060925.1", "chr1"), ("NC_060926.1", "chr2")])
        count = extract_rrna_loci(gff, bed, name_map_path=report)
        assert count == 2
        assert read_bed(bed) == [("chr1", "99", "500"), ("chr2", "199", "800")]


# -- load_name_map ------------------------------------------------------------

class TestLoadNameMap:
    def test_basic_mapping(self, tmp_path):
        report = tmp_path / "report.txt"
        write_report(report, [("NC_060925.1", "chr1"), ("NC_060926.1", "chr2")])
        m = load_name_map(report)
        assert m["NC_060925.1"] == "chr1"
        assert m["NC_060926.1"] == "chr2"

    def test_na_entries_excluded(self, tmp_path):
        report = tmp_path / "report.txt"
        with open(report, "w") as f:
            f.write("# comment\n")
            f.write("x\tx\tx\tx\tna\tx\tna\tx\t0\tx\n")
        m = load_name_map(report)
        assert len(m) == 0
