"""Tests for the top-level CLI: --vcf is the only required arg, plot
mode is auto-detected from the VCF header.
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import pytest

import molamola as mm


def test_bare_vcf_dispatches_to_sv(tmp_path, tiny_vcf):
    """A VCF with ##INFO=<ID=SVTYPE> auto-routes to the SV report."""
    rc = mm.main(["--vcf", str(tiny_vcf), "--out", str(tmp_path)])
    assert rc == 0
    assert (tmp_path / "tiny.report.html").exists()


def test_bare_vcf_dispatches_to_compound_het(
    tmp_path, tiny_phased_vcf, tiny_canonical_exons, tiny_clinvar,
):
    """A VCF with CSQ + PS auto-routes to per-gene phased panels."""
    rc = mm.main([
        "--vcf", str(tiny_phased_vcf),
        "--canonical-exons", str(tiny_canonical_exons),
        "--clinvar", str(tiny_clinvar),
        "--gene", "GENE_A",
        "--out", str(tmp_path),
    ])
    assert rc == 0
    assert (tmp_path / "tiny_phased.compound_het.report.html").exists()


def test_unrecognised_vcf_shape_refuses(tmp_path):
    """A VCF with neither SVTYPE nor CSQ+PS is refused, not silently rendered."""
    vcf = tmp_path / "boring.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    rc = mm.main(["--vcf", str(vcf), "--out", str(tmp_path)])
    assert rc == 1


def test_help_shows_both_flag_groups(capsys):
    """--help surfaces SV-mode flags and compound-het flags side by side."""
    with pytest.raises(SystemExit):
        mm.main(["--help"])
    out = capsys.readouterr().out
    assert "SV-mode flags" in out
    assert "Compound-het mode flags" in out
    assert "--vcf" in out
    assert "--gene" in out
    assert "--filter" in out


def test_compound_het_refuses_t2t(tmp_path, tiny_phased_vcf):
    """Compound-het is hg38-only; --reference t2t is rejected."""
    rc = mm.main([
        "--vcf", str(tiny_phased_vcf),
        "--reference", "t2t",
        "--out", str(tmp_path),
    ])
    assert rc == 1


def test_detect_vcf_mode_returns_sv(tiny_vcf):
    assert mm.detect_vcf_mode(tiny_vcf) == "sv"


def test_detect_vcf_mode_returns_compound_het(tiny_phased_vcf):
    assert mm.detect_vcf_mode(tiny_phased_vcf) == "compound-het"


def test_detect_vcf_mode_refuses_unknown(tmp_path):
    vcf = tmp_path / "no_signal.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    with pytest.raises(ValueError, match="doesn't match either"):
        mm.detect_vcf_mode(vcf)


def test_detect_vcf_mode_csq_alone_is_not_compound_het(tmp_path):
    """CSQ without PS isn't enough — compound-het needs both."""
    vcf = tmp_path / "csq_only.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: A">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    with pytest.raises(ValueError, match="doesn't match either"):
        mm.detect_vcf_mode(vcf)
