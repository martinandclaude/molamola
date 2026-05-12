"""Tests for karyotype-mode BAF parsing (pure-Python, plain text)."""

from __future__ import annotations

import gzip

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest

import molamola as mm


def test_read_baf_keeps_het_pass_only(tiny_baf):
    """Hom, non-PASS, low-DP, multi-allelic, non-canonical chrom all dropped."""
    df = mm.read_baf_vcf(tiny_baf, min_dp=10)
    assert len(df) == 5
    assert list(df.columns) == ["chrom", "pos", "baf"]
    # Expected rows: chr1:1000, chr1:2000, chr2:100, chrX:200, chrY:300
    by_chrom = df.groupby("chrom", observed=True)["pos"].apply(list).to_dict()
    assert by_chrom["chr1"] == [1000, 2000]
    assert by_chrom["chr2"] == [100]
    assert by_chrom["chrX"] == [200]
    assert by_chrom["chrY"] == [300]


def test_read_baf_uses_format_af(tiny_baf):
    df = mm.read_baf_vcf(tiny_baf, min_dp=10)
    chr1 = df[df["chrom"] == "chr1"].reset_index(drop=True)
    np.testing.assert_allclose(chr1["baf"].to_numpy(), [0.50, 0.45], rtol=1e-6)


def test_read_baf_min_dp_filter(tiny_baf):
    """min_dp=30 drops the DP=25 records."""
    df = mm.read_baf_vcf(tiny_baf, min_dp=30)
    assert len(df) == 2  # chr1:2000 (DP=30) and chr2:100 (DP=50)


def test_read_baf_falls_back_to_ad(tiny_baf_no_af):
    df = mm.read_baf_vcf(tiny_baf_no_af, min_dp=10)
    # chr1:1000 -> 15/30=0.5; chr1:2000 -> 22/40=0.55; chr1:3000 -> 0/50=0.0
    # chr2:100 has malformed AD -> dropped (3 kept)
    assert len(df) == 3
    np.testing.assert_allclose(df["baf"].to_numpy(), [0.5, 0.55, 0.0], rtol=1e-6)


def test_read_baf_reads_gzipped(tmp_path, tiny_baf):
    gz = tmp_path / "tiny_baf.vcf.gz"
    with open(tiny_baf, "rb") as src, gzip.open(gz, "wb") as dst:
        dst.write(src.read())
    df = mm.read_baf_vcf(gz, min_dp=10)
    assert len(df) == 5


def test_read_baf_empty_returns_empty_df(tmp_path):
    vcf = tmp_path / "empty.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )
    df = mm.read_baf_vcf(vcf, min_dp=10)
    assert len(df) == 0
    assert list(df.columns) == ["chrom", "pos", "baf"]


def test_read_baf_chrom_is_canonical_categorical(tiny_baf):
    df = mm.read_baf_vcf(tiny_baf, min_dp=10)
    assert df["chrom"].cat.ordered
    # canonical chroms only -> chrM was filtered out
    assert "chrM" not in df["chrom"].astype(str).unique()


@pytest.mark.parametrize("af,expected", [
    ("0.45", 0.45),
    ("0.5,0.3", 0.5),       # multi-value AF: take first
    (".", None),
    ("", None),
    ("not_a_number", None),
])
def test_parse_baf_handles_af_variants(af, expected):
    out = mm._parse_baf({"AF": af})
    if expected is None:
        assert out is None
    else:
        assert out == pytest.approx(expected)


@pytest.mark.parametrize("ad,expected", [
    ("15,15", 0.5),
    ("18,22", 22 / 40),
    ("50,0", 0.0),
    ("0,0", None),            # divide-by-zero guard
    (".", None),
    ("not,numbers", None),
])
def test_parse_baf_handles_ad_variants(ad, expected):
    out = mm._parse_baf({"AD": ad})
    if expected is None:
        assert out is None
    else:
        assert out == pytest.approx(expected)


def test_parse_baf_prefers_af_over_ad():
    out = mm._parse_baf({"AF": "0.42", "AD": "10,90"})
    assert out == pytest.approx(0.42)
