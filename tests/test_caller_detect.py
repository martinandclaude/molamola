"""Unit tests for detect_caller (Phase A: warn-only fingerprinting)."""

from __future__ import annotations

import molamola as mm


# --- by ##source= ----------------------------------------------------------

def test_detect_caller_recognises_sniffles2_via_source(tiny_vcf):
    caller, basis = mm.detect_caller(tiny_vcf)
    assert caller == "sniffles2"
    assert "Sniffles2" in basis or "STDEV_POS" in basis


def test_detect_caller_recognises_sniffles1_via_source(repo_root):
    caller, basis = mm.detect_caller(repo_root / "tests" / "data" / "tiny_sniffles1.vcf")
    # Source is "Sniffles" without 2; we route through INFO-field fallback
    # which sees SUPPORT but no COVERAGE -> sniffles1.
    assert caller == "sniffles1"
    assert "SUPPORT" in basis or "Sniffles" in basis


def test_detect_caller_recognises_cutesv_via_source(repo_root):
    caller, _ = mm.detect_caller(repo_root / "tests" / "data" / "tiny_cutesv.vcf")
    assert caller == "cutesv"


def test_detect_caller_recognises_svim_via_source(repo_root):
    caller, _ = mm.detect_caller(repo_root / "tests" / "data" / "tiny_svim.vcf")
    assert caller == "svim"


# --- INFO-field fingerprinting (no ##source available) ---------------------

def _write_vcf(path, source_line: str, info_examples: list[str]):
    """Write a minimal VCF with one INS record per INFO example."""
    header = [
        "##fileformat=VCFv4.2",
        source_line,
        "##contig=<ID=chr1,length=248956422>",
        "##FILTER=<ID=PASS,Description=\"\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    rows = [
        f"chr1\t{1_000_000 + i}\trec.{i}\tN\t<DEL>\t60\tPASS\t{info}\tGT\t0/1"
        for i, info in enumerate(info_examples)
    ]
    path.write_text("\n".join(header + rows) + "\n")


def test_detect_caller_uses_info_when_source_blank(tmp_path):
    """A bcftools-merged VCF often loses ##source -- fall back to INFO."""
    vcf = tmp_path / "merged.vcf"
    _write_vcf(
        vcf,
        source_line="##source=",  # blank
        info_examples=[
            "SVTYPE=DEL;SVLEN=-100;END=1000100;SUPPORT=10;COVERAGE=30,30,15,30,30;STDEV_POS=0.5",
        ],
    )
    caller, basis = mm.detect_caller(vcf)
    assert caller == "sniffles2"
    assert "SUPPORT" in basis


def test_detect_caller_unknown_when_no_signature(tmp_path):
    vcf = tmp_path / "weird.vcf"
    _write_vcf(
        vcf,
        source_line="##source=",
        info_examples=["SVTYPE=DEL;SVLEN=-100;END=1000100"],
    )
    caller, _ = mm.detect_caller(vcf)
    assert caller == "unknown"


def test_detect_caller_handles_empty_data_section(tmp_path):
    """A header-only VCF with no records still resolves cleanly."""
    vcf = tmp_path / "empty.vcf"
    _write_vcf(vcf, source_line="##source=Sniffles2_2.7.5", info_examples=[])
    caller, _ = mm.detect_caller(vcf)
    assert caller == "sniffles2"


# --- plot subcommand wiring -------------------------------------------------

def test_plot_auto_dispatches_cutesv_parser(repo_root, tmp_path, capsys):
    """Running plot on a cuteSV VCF auto-routes to the cuteSV parser."""
    rc = mm.main([
        "--vcf", str(repo_root / "tests" / "data" / "tiny_cutesv.vcf"),
        "--out", str(tmp_path),
    ])
    captured = capsys.readouterr()
    assert "Detected caller: cutesv" in captured.out
    # Phase B replaced the warn-only behaviour with real dispatch, so
    # there is no longer a stderr warning for non-Sniffles2 input.
    assert "WARNING" not in captured.err
    assert rc == 0


def test_plot_explicit_caller_mismatch_warns(repo_root, tmp_path, capsys):
    """--caller sniffles2 on a cuteSV VCF warns but proceeds."""
    rc = mm.main([
        "--vcf", str(repo_root / "tests" / "data" / "tiny_cutesv.vcf"),
        "--out", str(tmp_path),
        "--caller", "sniffles2",
    ])
    captured = capsys.readouterr()
    # The fingerprint disagrees with the explicit override; we warn.
    assert "WARNING" in captured.err
    assert "fingerprint suggests cutesv" in captured.err
    assert rc == 0


def test_plot_explicit_caller_matches_silently(repo_root, tmp_path, capsys):
    """--caller cutesv on a cuteSV VCF passes through with no warning."""
    rc = mm.main([
        "--vcf", str(repo_root / "tests" / "data" / "tiny_cutesv.vcf"),
        "--out", str(tmp_path),
        "--caller", "cutesv",
    ])
    captured = capsys.readouterr()
    assert "WARNING" not in captured.err
    assert "Caller: cutesv" in captured.out
    assert rc == 0
