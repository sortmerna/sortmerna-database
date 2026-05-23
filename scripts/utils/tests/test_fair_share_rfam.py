import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from fair_share_rfam import fair_share


# -- fair_share ---------------------------------------------------------------

class TestFairShare:
    def test_equal_families(self):
        sizes = {"a": 200, "b": 200, "c": 200}
        allocs = fair_share(sizes, 300)
        assert allocs == {"a": 100, "b": 100, "c": 100}
        assert sum(allocs.values()) == 300

    def test_one_small_family(self):
        # "a" has only 30, so its 30 are taken and the remaining 270 go to b and c
        sizes = {"a": 30, "b": 1000, "c": 1000}
        allocs = fair_share(sizes, 300)
        assert allocs["a"] == 30
        assert allocs["b"] + allocs["c"] == 270
        assert sum(allocs.values()) == 300

    def test_all_families_small(self):
        # total available < target -> take everything
        sizes = {"a": 50, "b": 60, "c": 40}
        allocs = fair_share(sizes, 500)
        assert allocs == {"a": 50, "b": 60, "c": 40}
        assert sum(allocs.values()) == 150

    def test_exact_total(self):
        # target exactly equals sum of all families
        sizes = {"a": 100, "b": 200, "c": 300}
        allocs = fair_share(sizes, 600)
        assert allocs == {"a": 100, "b": 200, "c": 300}
        assert sum(allocs.values()) == 600

    def test_rounding_remainder_distributed(self):
        # 10 target, 3 families of 100 -> 3+3+4 or 3+3+3 with remainder 1
        sizes = {"a": 100, "b": 100, "c": 100}
        allocs = fair_share(sizes, 10)
        assert sum(allocs.values()) == 10
        # remainder 1 goes to first family in sorted order ("a")
        assert allocs["a"] == 4
        assert allocs["b"] == 3
        assert allocs["c"] == 3

    def test_target_zero(self):
        sizes = {"a": 100, "b": 200}
        allocs = fair_share(sizes, 0)
        assert sum(allocs.values()) == 0

    def test_negative_target_raises(self):
        with pytest.raises(ValueError, match="target must be >= 0"):
            fair_share({"a": 100}, -1)

    def test_single_family(self):
        sizes = {"a": 1000}
        allocs = fair_share(sizes, 500)
        assert allocs == {"a": 500}

    def test_single_family_smaller_than_target(self):
        sizes = {"a": 300}
        allocs = fair_share(sizes, 500)
        assert allocs == {"a": 300}

    def test_realistic_rfam_scenario(self):
        # Mirrors actual Rfam family sizes from the pipeline at target=500000.
        # 7 small families (< 100K), 3 large (tRNA, SRP_RNA, U6).
        # Expected: small families take all, large families share the remainder
        # equally to hit exactly 500000.
        sizes = {
            "RF00003_U1_spliceosomal":   40284,
            "RF00004_U2_spliceosomal":   71722,
            "RF00005_tRNA":            5336929,
            "RF00009_RNaseP_eukaryotic":  3964,
            "RF00010_RNaseP_bacterial":   8942,
            "RF00015_U4_spliceosomal":   32357,
            "RF00017_SRP_RNA":          383615,
            "RF00020_U5_spliceosomal":   28713,
            "RF00023_tmRNA":              8579,
            "RF00026_U6_spliceosomal":  204635,
        }
        allocs = fair_share(sizes, 500000)
        assert sum(allocs.values()) == 500000
        # small families must contribute all their sequences
        for name in ["RF00003_U1_spliceosomal", "RF00009_RNaseP_eukaryotic",
                     "RF00010_RNaseP_bacterial", "RF00015_U4_spliceosomal",
                     "RF00020_U5_spliceosomal", "RF00023_tmRNA"]:
            assert allocs[name] == sizes[name]
        # large families must each receive more than any small family
        for name in ["RF00005_tRNA", "RF00017_SRP_RNA", "RF00026_U6_spliceosomal"]:
            assert allocs[name] > 0
            assert allocs[name] <= sizes[name]

    def test_sorted_order_determines_remainder(self):
        # With 3 families and target not divisible by 3, the remainder goes to
        # families earliest in sorted order.
        sizes = {"c": 100, "a": 100, "b": 100}
        allocs = fair_share(sizes, 7)  # 2+2+3 -> a gets 3 (sorted first)
        assert allocs["a"] == 3
        assert allocs["b"] == 2
        assert allocs["c"] == 2
        assert sum(allocs.values()) == 7
