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


def test_read_baf_refuses_sv_vcf(tmp_path):
    """Header has ##INFO=<ID=SVTYPE,...> -> refuse with a helpful error."""
    vcf = tmp_path / "sv.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">\n'
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t1000\t.\tA\t<DEL>\t30\tPASS\tSVTYPE=DEL\tGT\t0/1\n"
    )
    with pytest.raises(ValueError, match="SVTYPE"):
        mm.read_baf_vcf(vcf, min_dp=10)


def test_read_baf_refusal_message_points_at_small_variant_input(tmp_path):
    """The refusal text tells the user what kind of VCF the BAF needs."""
    vcf = tmp_path / "cnv.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of CNV">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )
    with pytest.raises(ValueError) as excinfo:
        mm.read_baf_vcf(vcf, min_dp=10)
    msg = str(excinfo.value)
    assert "small-variant" in msg
    assert "Clair3" in msg or "DeepVariant" in msg


def test_read_baf_skips_symbolic_alt(tmp_path):
    """ALT starting with '<' is silently skipped (no error, no row)."""
    vcf = tmp_path / "mixed.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
        "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t1000\t.\tA\tT\t30\tPASS\t.\tGT:DP:AF\t0/1:30:0.5\n"
        "chr1\t2000\t.\tA\t<INS>\t30\tPASS\t.\tGT:DP:AF\t0/1:30:0.5\n"
        "chr1\t3000\t.\tA\t<DEL>\t30\tPASS\t.\tGT:DP:AF\t0/1:30:0.5\n"
    )
    df = mm.read_baf_vcf(vcf, min_dp=10)
    assert len(df) == 1
    assert df["pos"].iloc[0] == 1000


def test_read_baf_skips_oversize_alt(tmp_path):
    """ALT longer than _KARY_BAF_MAX_ALT_LEN bp is silently skipped."""
    long_alt = "A" * (mm._KARY_BAF_MAX_ALT_LEN + 1)
    vcf = tmp_path / "with_oversize.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
        "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t1000\t.\tA\tT\t30\tPASS\t.\tGT:DP:AF\t0/1:30:0.5\n"
        f"chr1\t2000\t.\tA\t{long_alt}\t30\tPASS\t.\tGT:DP:AF\t0/1:30:0.5\n"
        # 50 bp is the cap; exactly 50 should be kept
        f"chr1\t3000\t.\tA\t{'C' * mm._KARY_BAF_MAX_ALT_LEN}\t30\tPASS\t.\tGT:DP:AF\t0/1:30:0.45\n"
    )
    df = mm.read_baf_vcf(vcf, min_dp=10)
    assert len(df) == 2
    assert list(df["pos"]) == [1000, 3000]
