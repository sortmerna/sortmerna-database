import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from plot_scalability import parse_log, parse_blast, write_top_hits, plot_summary, plot_evalue_dist, plot_identity_dist


# Matches the real SortMeRNA aligned.log format from summary.cpp
LOG_TEMPLATE = """\
 Command:
    sortmerna --ref db.fasta --reads reads.fasta

 Process pid =

 Parameters summary:
    Reference file: db.fasta
        Seed length = 18
        Pass 1 = 18, Pass 2 = 9, Pass 3 = 3
        Gumbel lambda = {lambda_}
        Gumbel K = {gumbel_k}
        Minimal SW score based on E-value = {s_min}
    Number of seeds = 2
    Edges = 4
    SW match = 2
    SW mismatch = -3
    SW gap open penalty = 5
    SW gap extend penalty = 2
    SW ambiguous nucleotide = -3
    SQ tags are not output
    Number of alignment processing threads = 4
    Reads file: reads.fasta
    Total reads = {total_reads}

 Results:
    Total reads passing E-value threshold = {aligned} (1.50)
    Total reads failing E-value threshold = {failing} (98.50)
    Minimum read length = {min_len}
    Maximum read length = {max_len}
    Mean read length    = {mean_len}

 Coverage by database:
    db.fasta\t\t1.50

 Mon May 25 02:56:29 2026
"""


def make_log(tmp_path, total=1000000, aligned=14955, min_len=151, max_len=151,
             mean_len=151, s_min=59, lambda_=0.609747, gumbel_k=0.334118):
    p = tmp_path / 'aligned.log'
    p.write_text(LOG_TEMPLATE.format(
        total_reads=total, aligned=aligned, failing=total - aligned,
        min_len=min_len, max_len=max_len, mean_len=mean_len,
        s_min=s_min, lambda_=lambda_, gumbel_k=gumbel_k,
    ))
    return p


def make_blast(tmp_path, rows):
    """Write minimal BLAST tabular rows as (pident, evalue) pairs."""
    p = tmp_path / 'aligned.blast'
    with open(p, 'w') as f:
        for i, (pident, evalue) in enumerate(rows):
            f.write(f'read{i}\tref1\t{pident}\t151\t0\t0\t1\t151\t1\t151\t{evalue}\t100\n')
    return p


# -- parse_log ----------------------------------------------------------------

class TestParseLog:
    def test_basic_fields(self, tmp_path):
        p = make_log(tmp_path, total=1000000, aligned=14955)
        result = parse_log(p)
        assert result['total_reads'] == 1000000
        assert result['aligned']     == 14955

    def test_s_min_integer(self, tmp_path):
        p = make_log(tmp_path, s_min=59)
        assert parse_log(p)['s_min'] == pytest.approx(59.0)

    def test_gumbel_params(self, tmp_path):
        p = make_log(tmp_path, lambda_=0.609747, gumbel_k=0.334118)
        result = parse_log(p)
        assert result['lambda_']  == pytest.approx(0.609747)
        assert result['gumbel_k'] == pytest.approx(0.334118)

    def test_read_length_fields(self, tmp_path):
        p = make_log(tmp_path, min_len=151, max_len=151, mean_len=151)
        result = parse_log(p)
        assert result['min_len']  == 151
        assert result['max_len']  == 151
        assert result['mean_len'] == 151

    def test_missing_field_returns_none(self, tmp_path):
        p = tmp_path / 'aligned.log'
        p.write_text('    Total reads = 500\n')
        result = parse_log(p)
        assert result['total_reads'] == 500
        assert result['aligned']     is None
        assert result['s_min']       is None


# -- parse_blast --------------------------------------------------------------

class TestParseBlast:
    def test_basic_columns(self, tmp_path):
        p = make_blast(tmp_path, [(99.5, 1e-10), (85.0, 1e-5)])
        df = parse_blast(p)
        assert list(df.columns) == ['pident', 'evalue']
        assert len(df) == 2

    def test_pident_values(self, tmp_path):
        p = make_blast(tmp_path, [(100.0, 0.0), (75.3, 1e-3)])
        df = parse_blast(p)
        assert df['pident'].tolist() == pytest.approx([100.0, 75.3])

    def test_evalue_values(self, tmp_path):
        p = make_blast(tmp_path, [(99.0, 1e-20), (80.0, 5e-4)])
        df = parse_blast(p)
        assert df['evalue'].tolist() == pytest.approx([1e-20, 5e-4])

    def test_empty_file_returns_empty_dataframe(self, tmp_path):
        p = tmp_path / 'aligned.blast'
        p.write_text('')
        df = parse_blast(p)
        assert df.empty
        assert set(df.columns) == {'pident', 'evalue'}

    def test_missing_file_returns_empty_dataframe(self, tmp_path):
        df = parse_blast(tmp_path / 'nonexistent.blast')
        assert df.empty

    def test_chunked_concat(self, tmp_path, monkeypatch):
        """Rows spanning multiple chunks are all returned when chunksize is small."""
        import plot_scalability
        monkeypatch.setattr(plot_scalability, '_BLAST_CHUNKSIZE', 3)
        rows = [(90.0 + i % 10, 1e-10) for i in range(10)]
        p = make_blast(tmp_path, rows)
        df = parse_blast(p)
        assert len(df) == 10


# -- smoke tests for plotting functions --------------------------------------

def make_stats(tmp_path):
    """Build a minimal stats list matching the structure produced by main()."""
    stats = []
    for i, n in enumerate([10000, 100000]):
        d = tmp_path / f'scale_{n}'
        d.mkdir()
        blast_rows = [(95.0 + j % 5, 1e-10) for j in range(20)]
        stats.append({
            'n':       n,
            'runtime': 10 * (i + 1),
            'log': {
                'total_reads': n,
                'aligned':     n // 100,
                's_min':       40.0 + i * 10,
            },
            'blast': parse_blast(make_blast(d, blast_rows)),
        })
    return stats


# -- write_top_hits -----------------------------------------------------------

class TestWriteTopHits:
    def test_counts_and_order(self, tmp_path):
        p = tmp_path / 'aligned.blast'
        p.write_text(
            'r0\tsilva_A\t99.0\t151\t0\t0\t1\t151\t1\t151\t1e-10\t100\n'
            'r1\tsilva_B\t85.0\t151\t0\t0\t1\t151\t1\t151\t1e-5\t80\n'
            'r2\tsilva_A\t95.0\t151\t0\t0\t1\t151\t1\t151\t1e-8\t90\n'
        )
        out = tmp_path / 'top_hits.txt'
        write_top_hits(p, out)
        lines = out.read_text().splitlines()
        assert lines[0].strip().startswith('Count')
        data_lines = [l for l in lines[2:] if l.strip()]
        assert 'silva_A' in data_lines[0]   # most frequent first
        assert 'silva_B' in data_lines[1]

    def test_top_n_limit(self, tmp_path):
        rows = [(90.0, 1e-5)] * 30
        p = tmp_path / 'aligned.blast'
        with open(p, 'w') as f:
            for i in range(30):
                f.write(f'r{i}\tref_{i:02d}\t90.0\t151\t0\t0\t1\t151\t1\t151\t1e-5\t80\n')
        out = tmp_path / 'top_hits.txt'
        write_top_hits(p, out, top_n=5)
        data_lines = [l for l in out.read_text().splitlines()[2:] if l.strip()]
        assert len(data_lines) == 5

    def test_empty_file_writes_nothing(self, tmp_path):
        p = tmp_path / 'aligned.blast'
        p.write_text('')
        out = tmp_path / 'top_hits.txt'
        write_top_hits(p, out)
        assert not out.exists()

    def test_missing_file_writes_nothing(self, tmp_path):
        out = tmp_path / 'top_hits.txt'
        write_top_hits(tmp_path / 'nonexistent.blast', out)
        assert not out.exists()


try:
    import matplotlib  # noqa: F401
    _has_matplotlib = True
except ImportError:
    _has_matplotlib = False


@pytest.mark.skipif(not _has_matplotlib, reason='matplotlib not installed')
class TestPlotFunctions:
    def test_plot_summary_creates_file(self, tmp_path):
        plot_summary(make_stats(tmp_path), 'test', tmp_path)
        assert (tmp_path / 'test_summary.png').exists()

    def test_plot_evalue_dist_creates_file(self, tmp_path):
        plot_evalue_dist(make_stats(tmp_path), 'test', tmp_path)
        assert (tmp_path / 'test_evalue_dist.png').exists()

    def test_plot_identity_dist_creates_file(self, tmp_path):
        plot_identity_dist(make_stats(tmp_path), 'test', tmp_path)
        assert (tmp_path / 'test_identity_dist.png').exists()

    def test_plot_evalue_dist_empty_blast(self, tmp_path):
        stats = make_stats(tmp_path)
        for s in stats:
            s['blast'] = pd.DataFrame(columns=['pident', 'evalue'])
        plot_evalue_dist(stats, 'empty', tmp_path)
        assert (tmp_path / 'empty_evalue_dist.png').exists()
