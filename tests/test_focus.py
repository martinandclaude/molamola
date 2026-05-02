"""Unit tests for the --focus chr:pos parser and matcher."""

from __future__ import annotations

import argparse

import pytest

import molamola as mm


# --- parse_focus ------------------------------------------------------------

@pytest.mark.parametrize(
    "spec, expected",
    [
        ("chr1:73129297",     ("chr1", 73129297)),
        ("1:73129297",        ("chr1", 73129297)),
        ("1:73,129,297",      ("chr1", 73129297)),
        ("1:73_129_297",      ("chr1", 73129297)),
        ("chrX:1000000",      ("chrX", 1000000)),
        ("chrY:1",            ("chrY", 1)),
    ],
)
def test_parse_focus_valid(spec, expected):
    assert mm.parse_focus(spec) == expected


@pytest.mark.parametrize("spec", ["chr1", "chr1:abc", "no-colon", "chr1:"])
def test_parse_focus_invalid_raises(spec):
    with pytest.raises(argparse.ArgumentTypeError):
        mm.parse_focus(spec)


# --- matches_focus ----------------------------------------------------------

def test_matches_focus_matches_chr1_endpoint(make_bnd):
    b = make_bnd(chr1="chr1", pos1=1000, chr2="chr2", pos2=2000)
    assert mm.matches_focus(b, [("chr1", 1000)], window=0) is True


def test_matches_focus_matches_chr2_endpoint(make_bnd):
    b = make_bnd(chr1="chr1", pos1=1000, chr2="chr2", pos2=2000)
    assert mm.matches_focus(b, [("chr2", 2000)], window=0) is True


def test_matches_focus_respects_window(make_bnd):
    b = make_bnd(chr1="chr1", pos1=1000)
    assert mm.matches_focus(b, [("chr1", 1500)], window=1000) is True
    assert mm.matches_focus(b, [("chr1", 1500)], window=400) is False


def test_matches_focus_wrong_chrom_no_match(make_bnd):
    b = make_bnd(chr1="chr1", pos1=1000, chr2="chr2", pos2=2000)
    assert mm.matches_focus(b, [("chr3", 1000)], window=10000) is False


def test_matches_focus_multiple_specs_any_match(make_bnd):
    b = make_bnd(chr1="chr1", pos1=1000, chr2="chr2", pos2=2000)
    foci = [("chr5", 100), ("chr2", 2000)]
    assert mm.matches_focus(b, foci, window=0) is True


def test_matches_focus_no_specs_returns_false(make_bnd):
    b = make_bnd()
    assert mm.matches_focus(b, [], window=1000) is False
