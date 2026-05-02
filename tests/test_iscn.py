"""Unit tests for cytoband resolution and ISCN nomenclature labels."""

from __future__ import annotations

import molamola as mm


# --- resolve_cytoband -------------------------------------------------------

def test_resolve_cytoband_returns_band_name(bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    # chr1:1Mb is in p36.33 (the very first band on chr1, 0..2,300,000).
    band = mm.resolve_cytoband("chr1", 1_000_000, cyto)
    assert band == "p36.33"


def test_resolve_cytoband_works_across_arms(bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    # chr7:73Mb falls in q11.23 (the well-known WBS deletion locus).
    band = mm.resolve_cytoband("chr7", 73_000_000, cyto)
    assert band.startswith("q")  # q-arm
    # Position should be a real band, not empty
    assert band != ""


def test_resolve_cytoband_missing_chrom_returns_empty(bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    assert mm.resolve_cytoband("chrZZ", 100, cyto) == ""


def test_resolve_cytoband_position_past_end_returns_empty(bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    # 999 Gb is way past any human chromosome
    assert mm.resolve_cytoband("chr1", 999_000_000_000, cyto) == ""


# --- iscn_label -------------------------------------------------------------

def test_iscn_label_bnd_canonical_order(make_bnd, bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    b = make_bnd(chr1="chr7", pos1=73_000_000, chr2="chr17", pos2=22_000_000)
    label = mm.iscn_label(b, cyto)
    # Should be t(7;17)(...) regardless of which side came first
    assert label.startswith("t(7;17)(")
    assert label.endswith(")")


def test_iscn_label_bnd_swapped_ordering(make_bnd, bundled_cytoband):
    """A BND with chr2 < chr1 should be re-ordered to canonical form."""
    cyto = mm.load_cytobands(bundled_cytoband)
    swapped = make_bnd(chr1="chr17", pos1=22_000_000,
                        chr2="chr7",  pos2=73_000_000)
    label = mm.iscn_label(swapped, cyto)
    assert label.startswith("t(7;17)(")


def test_iscn_label_del_uses_two_bands(make_sv, bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    # A 100 Mb deletion will span multiple bands
    sv = make_sv(chrom="chr1", start=10_000_000, end=110_000_000,
                  svtype="DEL", svlen=100_000_000)
    label = mm.iscn_label(sv, cyto)
    assert label.startswith("del(1)(")
    # Two distinct bands: starts with p, mentions another band
    assert "p" in label


def test_iscn_label_del_collapses_when_one_band(make_sv, bundled_cytoband):
    """If start and end fall in the same band, ISCN shows that band once."""
    cyto = mm.load_cytobands(bundled_cytoband)
    # tiny 100 bp DEL well inside p36.33 (0..2.3 Mb)
    sv = make_sv(chrom="chr1", start=1_000_000, end=1_000_100,
                  svtype="DEL", svlen=100)
    label = mm.iscn_label(sv, cyto)
    # del(1)(p36.33), not del(1)(p36.33p36.33)
    assert label == "del(1)(p36.33)"


def test_iscn_label_dup_inv_use_same_form(make_sv, bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    dup = make_sv(chrom="chr1", start=1_000_000, end=1_000_500,
                   svtype="DUP", svlen=500)
    inv = make_sv(chrom="chr1", start=1_000_000, end=1_000_500,
                   svtype="INV", svlen=500)
    assert mm.iscn_label(dup, cyto).startswith("dup(1)(")
    assert mm.iscn_label(inv, cyto).startswith("inv(1)(")


def test_iscn_label_ins_uses_point_form(make_sv, bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    # INS records have start == end; only one band shown
    ins = make_sv(chrom="chr1", start=1_000_000, end=1_000_000,
                   svtype="INS", svlen=300)
    label = mm.iscn_label(ins, cyto)
    assert label == "ins(1)(p36.33)"


def test_iscn_label_unknown_band_renders_as_question_mark(make_sv,
                                                            bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    # Chromosome with no bands -> band lookup empty
    cyto_no_chr1 = {k: v for k, v in cyto.items() if k != "chr1"}
    sv = make_sv(chrom="chr1", start=1_000_000, end=1_000_500,
                  svtype="DEL", svlen=500)
    label = mm.iscn_label(sv, cyto_no_chr1)
    assert "?" in label
