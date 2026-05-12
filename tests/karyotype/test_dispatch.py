"""Top-level CLI dispatch tests for karyotype coverage mode.

Karyotype mode is triggered by ``--mosdepth PATH``. ``--vcf`` becomes
optional; either ``--vcf`` or ``--mosdepth`` is required, and when
both are supplied karyotype mode wins (the VCF is consumed as the
BAF source, not run through SV / compound-het header dispatch).
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import pytest

import molamola as mm


def test_mosdepth_alone_dispatches_to_karyotype(tmp_path, capsys):
    """``--mosdepth foo.bed.gz`` routes to karyotype mode (stub returns 0)."""
    bed = tmp_path / "fake.regions.bed.gz"
    rc = mm.main(["--mosdepth", str(bed), "--out", str(tmp_path)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "karyotype coverage mode" in out


def test_mosdepth_plus_vcf_dispatches_to_karyotype(tmp_path, capsys):
    """Both flags: karyotype mode wins; VCF is treated as the BAF source.

    The VCF passed here has no SV/CSQ/PS headers; if dispatch tried to
    run :func:`detect_vcf_mode` on it the CLI would exit 1. A clean
    exit-0 proves karyotype dispatch intercepted before VCF header
    inspection.
    """
    bed = tmp_path / "fake.regions.bed.gz"
    vcf = tmp_path / "no_signal.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    rc = mm.main([
        "--mosdepth", str(bed),
        "--vcf", str(vcf),
        "--out", str(tmp_path),
    ])
    assert rc == 0
    out = capsys.readouterr().out
    assert "karyotype coverage mode" in out
    assert f"BAF source: {vcf}" in out


def test_neither_input_flag_returns_2(capsys):
    """Neither ``--vcf`` nor ``--mosdepth`` is given: hard refusal."""
    rc = mm.main([])
    assert rc == 2
    err = capsys.readouterr().err
    assert "either --vcf or --mosdepth is required" in err


def test_vcf_alone_still_dispatches_to_sv(tmp_path, tiny_vcf):
    """Existing SV-mode dispatch is unchanged when ``--mosdepth`` is absent."""
    rc = mm.main(["--vcf", str(tiny_vcf), "--out", str(tmp_path)])
    assert rc == 0
    assert (tmp_path / "tiny.report.html").exists()


def test_help_surfaces_mosdepth_flag(capsys):
    """``--help`` lists the new ``--mosdepth`` flag alongside ``--vcf``."""
    with pytest.raises(SystemExit):
        mm.main(["--help"])
    out = capsys.readouterr().out
    assert "--vcf" in out
    assert "--mosdepth" in out
