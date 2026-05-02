"""Unit tests for cytoband loading and acrocentric helpers."""

from __future__ import annotations

import molamola as mm


def test_load_cytobands_returns_all_24_chromosomes(bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    assert set(cyto.keys()) == set(mm.CHROM_SET)


def test_load_cytobands_each_chr_has_bands(bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    for c, bands in cyto.items():
        assert len(bands) > 0, f"{c} has no bands"
        for (start, end, name, stain) in bands:
            assert isinstance(start, int) and isinstance(end, int)
            assert end > start
            assert isinstance(name, str)
            assert isinstance(stain, str)


def test_load_cytobands_has_acen_for_each_chrom(bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    for c, bands in cyto.items():
        stains = {s for (_, _, _, s) in bands}
        assert "acen" in stains, f"{c} missing acen band"


def test_acrocentric_p_arm_ends_returns_all_five(bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    ends = mm.acrocentric_p_arm_ends(cyto)
    assert set(ends.keys()) == set(mm.ACROCENTRIC_CHROMS)


def test_acrocentric_p_arm_ends_within_known_ranges(bundled_cytoband):
    """Sanity-check the p/q boundary against published hg38 centromere positions."""
    cyto = mm.load_cytobands(bundled_cytoband)
    ends = mm.acrocentric_p_arm_ends(cyto)
    # hg38 centromere starts (UCSC): chr13~16Mb, chr14~16Mb, chr15~17Mb,
    # chr21~12Mb, chr22~14Mb. Allow generous bounds.
    assert 10_000_000 < ends["chr13"] < 20_000_000
    assert 10_000_000 < ends["chr14"] < 20_000_000
    assert 14_000_000 < ends["chr15"] < 20_000_000
    assert  8_000_000 < ends["chr21"] < 16_000_000
    assert 10_000_000 < ends["chr22"] < 18_000_000
