"""Unit tests for T2T-CHM13v2.0 reference support (`--reference t2t`)."""

from __future__ import annotations

import pytest

import molamola as mm


# --- find_cytoband_file: bundled cytobands ---------------------------------

def test_find_cytoband_file_t2t_uses_bundled(package_data_dir):
    """T2T default returns the bundled cytoBand.t2t.txt.gz."""
    found = mm.find_cytoband_file(reference="t2t")
    assert found == package_data_dir / "cytoBand.t2t.txt.gz"


def test_find_cytoband_file_hg38_default_unchanged(package_data_dir):
    """hg38 default returns the bundled cytoBand.txt.gz."""
    found = mm.find_cytoband_file()
    assert found == package_data_dir / "cytoBand.txt.gz"


# --- end-to-end plot --reference t2t ----------------------------------------

def test_plot_reference_t2t_uses_bundled_cytoband(repo_root, tmp_path):
    """`--reference t2t` runs end-to-end with the bundled cytoband."""
    rc = mm.main([
        "--vcf", str(repo_root / "tests" / "data" / "tiny_t2t.vcf"),
        "--out", str(tmp_path),
        "--reference", "t2t",
    ])
    assert rc == 0
    assert (tmp_path / "tiny_t2t.report.html").exists()


def test_acrocentric_default_off_for_t2t(repo_root, tmp_path, capsys):
    """Without --mark-acrocentric, T2T should NOT flag the chr21<->chr22 BND."""
    rc = mm.main([
        "--vcf", str(repo_root / "tests" / "data" / "tiny_t2t.vcf"),
        "--out", str(tmp_path),
        "--reference", "t2t",
    ])
    assert rc == 0
    out = capsys.readouterr().out
    assert "acrocentric=0" in out


def test_acrocentric_default_on_for_hg38(tiny_vcf, tmp_path, capsys):
    """Without --mark-acrocentric, hg38 still flags acrocentric BNDs."""
    rc = mm.main([
        "--vcf", str(tiny_vcf),
        "--out", str(tmp_path),
    ])
    assert rc == 0
    out = capsys.readouterr().out
    assert "acrocentric=1" in out


def test_acrocentric_explicit_override_wins_for_t2t(repo_root, tmp_path, capsys):
    """Explicit --mark-acrocentric on T2T should still flag."""
    rc = mm.main([
        "--vcf", str(repo_root / "tests" / "data" / "tiny_t2t.vcf"),
        "--out", str(tmp_path),
        "--reference", "t2t",
        "--mark-acrocentric",
    ])
    assert rc == 0
    out = capsys.readouterr().out
    assert "acrocentric=1" in out


# --- detect_reference_hint -------------------------------------------------

@pytest.mark.parametrize("name, expected", [
    ("sample.hg38.sniffles.vcf",   "hg38"),
    ("sample.GRCh38.sniffles.vcf", "hg38"),
    ("sample.t2t.sniffles.vcf",    "t2t"),
    ("sample.T2T.sniffles.vcf",    "t2t"),
    ("sample.chm13.sniffles.vcf",  "t2t"),
    ("sample.CHM13v2.sniffles.vcf", "t2t"),
    ("sample.sniffles.vcf",        None),
    ("PATIENT123.sniffles.vcf",    None),
    # Both hints present (e.g. liftover): no clear winner
    ("sample.hg38_to_t2t.vcf",     None),
])
def test_detect_reference_hint(name, expected):
    assert mm.detect_reference_hint(name) == expected


# --- filename / --reference mismatch sanity check --------------------------

def test_reference_mismatch_t2t_filename_with_hg38_reference_errors(
    repo_root, tmp_path, capsys
):
    """A '*.t2t.*' filename with --reference hg38 should error out."""
    vcf_src = repo_root / "tests" / "data" / "tiny_t2t.vcf"
    vcf_renamed = tmp_path / "sample.t2t.sniffles.vcf"
    vcf_renamed.write_bytes(vcf_src.read_bytes())
    rc = mm.main([
        "--vcf", str(vcf_renamed),
        "--out", str(tmp_path / "out"),
    ])
    assert rc == 2
    err = capsys.readouterr().err
    assert "suggests --reference t2t" in err
    assert "Use --force to override" in err


def test_reference_mismatch_hg38_filename_with_t2t_reference_errors(
    tiny_vcf, tmp_path, capsys
):
    """A '*.hg38.*' filename with --reference t2t should error out."""
    vcf_renamed = tmp_path / "sample.hg38.sniffles.vcf"
    vcf_renamed.write_bytes(tiny_vcf.read_bytes())
    rc = mm.main([
        "--vcf", str(vcf_renamed),
        "--out", str(tmp_path / "out"),
        "--reference", "t2t",
    ])
    assert rc == 2
    err = capsys.readouterr().err
    assert "suggests --reference hg38" in err


def test_reference_mismatch_force_continues_with_warning(
    tiny_vcf, tmp_path, capsys
):
    """--force lets a mismatch through (warning to stderr only)."""
    # hg38-content VCF renamed to look like t2t. With --reference hg38 the
    # downstream pipeline runs cleanly (data + cytoband match); the filename
    # hint of t2t triggers the gate, which --force bypasses with a warning.
    vcf_renamed = tmp_path / "sample.t2t.sniffles.vcf"
    vcf_renamed.write_bytes(tiny_vcf.read_bytes())
    rc = mm.main([
        "--vcf", str(vcf_renamed),
        "--out", str(tmp_path / "out"),
        "--reference", "hg38",
        "--force",
    ])
    assert rc == 0
    err = capsys.readouterr().err
    assert "Continuing because of --force" in err
    assert "suggests --reference t2t" in err


def test_reference_match_no_warning(repo_root, tmp_path, capsys):
    """Filename hint matching --reference: no warning, no error."""
    vcf_src = repo_root / "tests" / "data" / "tiny_t2t.vcf"
    vcf_renamed = tmp_path / "sample.t2t.sniffles.vcf"
    vcf_renamed.write_bytes(vcf_src.read_bytes())
    rc = mm.main([
        "--vcf", str(vcf_renamed),
        "--out", str(tmp_path / "out"),
        "--reference", "t2t",
    ])
    assert rc == 0
    err = capsys.readouterr().err
    assert "suggests --reference" not in err
    assert "Continuing because of --force" not in err


def test_no_reference_hint_no_check(tiny_vcf, tmp_path, capsys):
    """Filename without t2t/hg38/chm13/grch38: check is silent."""
    rc = mm.main([
        "--vcf", str(tiny_vcf),
        "--out", str(tmp_path),
    ])
    assert rc == 0
    err = capsys.readouterr().err
    assert "suggests --reference" not in err
