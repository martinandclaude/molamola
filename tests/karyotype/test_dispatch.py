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


def test_mosdepth_alone_dispatches_to_karyotype(tmp_path, tiny_regions):
    """``--mosdepth foo.bed.gz`` runs karyotype mode end-to-end.

    Uses ``--no-mask --no-gc`` because the bundled exclusion mask and
    GC table are calibrated for ~500 bp mosdepth runs; the synthetic
    fixture is at 1 Mb where the mask would over-apply. Mask /
    GC-correction logic is exercised separately in test_normalise.py.
    """
    rc = mm.main([
        "--mosdepth", str(tiny_regions),
        "--reference", "hg38",
        "--no-mask", "--no-gc",
        "--out", str(tmp_path),
    ])
    assert rc == 0
    out_html = tmp_path / "tiny_regions.karyotype.report.html"
    assert out_html.exists()
    body = out_html.read_text()
    # Two embedded karyotype figures + one header-fish PNG = 3
    assert body.count("data:image/png;base64,") == 3
    assert 'id="fig-karyotype-genome"' in body
    assert 'id="fig-karyotype-per-chrom"' in body


def test_mosdepth_plus_vcf_dispatches_to_karyotype(tmp_path, tiny_regions):
    """Both flags: karyotype mode wins; VCF goes into the BAF panel.

    Pass a VCF with no SV/CSQ/PS headers. If the dispatch tried to
    run :func:`detect_vcf_mode` on it the CLI would exit 1; a clean
    exit-0 proves karyotype dispatch intercepted before VCF header
    inspection.
    """
    vcf = tmp_path / "no_signal.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )
    rc = mm.main([
        "--mosdepth", str(tiny_regions),
        "--vcf", str(vcf),
        "--reference", "hg38",
        "--no-mask", "--no-gc",
        "--out", str(tmp_path),
    ])
    assert rc == 0
    body = (tmp_path / "tiny_regions.karyotype.report.html").read_text()
    assert vcf.name in body


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
