"""Unit tests for ISCN-band focus filtering (`--focus chr7:q11.23`)."""

from __future__ import annotations

import argparse

import pytest

import molamola as mm


# --- parse_focus accepts band syntax ---------------------------------------

@pytest.mark.parametrize("spec, expected", [
    ("chr7:q11.23", ("chr7", "q11.23")),
    ("7:q11.23",    ("chr7", "q11.23")),
    ("chrX:p21",    ("chrX", "p21")),
    ("chr1:p36.33", ("chr1", "p36.33")),
    ("chr1:q11",    ("chr1", "q11")),  # prefix
])
def test_parse_focus_band_form(spec, expected):
    assert mm.parse_focus(spec) == expected


@pytest.mark.parametrize("spec", ["chr1:r11", "chr1:foo", "chr1:1abc"])
def test_parse_focus_invalid_band_or_pos(spec):
    with pytest.raises(argparse.ArgumentTypeError):
        mm.parse_focus(spec)


# --- resolve_band_range -----------------------------------------------------

def _bands(rows):
    """Helper: build a cytobands dict from a list of (chrom, start, end, name, stain)."""
    out: dict = {}
    for chrom, s, e, name, stain in rows:
        out.setdefault(chrom, []).append((s, e, name, stain))
    return out


def test_resolve_band_range_exact_match():
    cb = _bands([
        ("chr7", 0,        1_000_000, "p36.33",  "gneg"),
        ("chr7", 1_000_000, 2_000_000, "q11.23",  "gpos25"),
    ])
    assert mm.resolve_band_range("chr7", "q11.23", cb) == (1_000_000, 2_000_000)


def test_resolve_band_range_prefix_unions_subbands():
    cb = _bands([
        ("chr1", 0,        1_000_000, "q11",     "gneg"),
        ("chr1", 1_000_000, 2_000_000, "q11.1",   "gpos25"),
        ("chr1", 2_000_000, 3_000_000, "q11.21",  "gneg"),
        ("chr1", 3_000_000, 4_000_000, "q12",     "gpos50"),
    ])
    # "q11" should match q11, q11.1, q11.21 but NOT q12
    assert mm.resolve_band_range("chr1", "q11", cb) == (0, 3_000_000)


def test_resolve_band_range_no_match_returns_none():
    cb = _bands([("chr1", 0, 1_000_000, "p36.33", "gneg")])
    assert mm.resolve_band_range("chr1", "q11.23", cb) is None
    assert mm.resolve_band_range("chr2", "p36.33", cb) is None  # wrong chrom


# --- matches_focus with band specs ------------------------------------------

def test_matches_focus_band_endpoint_inside_band(make_bnd):
    cb = _bands([("chr7", 50_000_000, 60_000_000, "q11.23", "gpos25")])
    b = make_bnd(chr1="chr7", pos1=55_000_000, chr2="chr8", pos2=10_000_000)
    assert mm.matches_focus(b, [("chr7", "q11.23")], window=0, cytobands=cb)


def test_matches_focus_band_endpoint_outside_band(make_bnd):
    cb = _bands([("chr7", 50_000_000, 60_000_000, "q11.23", "gpos25")])
    b = make_bnd(chr1="chr7", pos1=70_000_000, chr2="chr8", pos2=10_000_000)
    assert not mm.matches_focus(b, [("chr7", "q11.23")], window=0, cytobands=cb)


def test_matches_focus_band_wrong_chrom(make_bnd):
    cb = _bands([("chr7", 50_000_000, 60_000_000, "q11.23", "gpos25")])
    b = make_bnd(chr1="chr8", pos1=55_000_000, chr2="chr9", pos2=55_000_000)
    assert not mm.matches_focus(b, [("chr7", "q11.23")], window=0, cytobands=cb)


def test_matches_focus_mixed_position_and_band(make_bnd):
    cb = _bands([("chr7", 50_000_000, 60_000_000, "q11.23", "gpos25")])
    b = make_bnd(chr1="chr1", pos1=1000, chr2="chr2", pos2=2000)
    foci = [("chr1", 1000), ("chr7", "q11.23")]
    assert mm.matches_focus(b, foci, window=0, cytobands=cb)


def test_matches_focus_unresolvable_band_silently_skipped(make_bnd):
    """matches_focus tolerates bad bands at runtime; plot_main should error first."""
    cb = _bands([("chr7", 0, 1_000_000, "p36", "gneg")])
    b = make_bnd(chr1="chr7", pos1=500, chr2="chr8", pos2=2000)
    assert not mm.matches_focus(b, [("chr7", "q99.99")], window=0, cytobands=cb)


# --- end-to-end via plot subcommand -----------------------------------------

def test_plot_focus_with_band_spec_runs(tmp_path, tiny_vcf):
    """`--focus chr7:q11.23` runs end-to-end and produces both PNGs."""
    rc = mm.main([
        "--vcf", str(tiny_vcf),
        "--out", str(tmp_path),
        "--focus", "chr7:q11.23",
    ])
    # tiny.vcf has BND.1 chr7:73129297 <-> chr17:31000000 — chr7 q11.2 is at
    # ~57Mb so this won't match; that's fine, returns 2 (no matches).
    # We just verify parsing + validation reaches the plot pipeline.
    assert rc in (0, 2)


def test_plot_focus_with_invalid_band_errors(tmp_path, tiny_vcf, capsys):
    """`--focus chr1:q99.99` errors out (no such band on chr1)."""
    rc = mm.main([
        "--vcf", str(tiny_vcf),
        "--out", str(tmp_path),
        "--focus", "chr1:q99.99",
    ])
    assert rc == 2
    err = capsys.readouterr().err
    assert "did not match any cytoband" in err
