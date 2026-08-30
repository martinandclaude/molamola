"""Header-driven mode dispatch, and the compound-het removal refusal.

This coverage used to live in ``tests/compound_het/test_cli_dispatch.py``
and came out with that directory when compound-het mode was removed.
:func:`detect_vcf_mode` still gates every VCF molamola accepts, so it is
restored here rather than lost.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import molamola as mm


def _write(tmp_path: Path, name: str, header: str) -> Path:
    vcf = tmp_path / name
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        f"{header}"
        "##contig=<ID=chr1,length=248956422>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )
    return vcf


# --- detect_vcf_mode -------------------------------------------------------

def test_sv_vcf_detected(tiny_vcf):
    assert mm.detect_vcf_mode(tiny_vcf) == "sv"


def test_phased_vep_vcf_reports_the_removal(tiny_phased_vcf):
    """The CSQ + PS shape is still recognised, so a returning user can be
    told the mode was removed rather than that their file is malformed."""
    assert mm.detect_vcf_mode(tiny_phased_vcf) == "compound-het-removed"


def test_svtype_wins_over_csq_and_ps(tmp_path):
    """A VEP-annotated phased SV VCF is rare but plottable, and plotting
    beats refusing. This precedence flipped when compound-het was removed."""
    vcf = _write(
        tmp_path, "both.vcf",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"t\">\n"
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Format: Allele\">\n"
        "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"p\">\n",
    )
    assert mm.detect_vcf_mode(vcf) == "sv"


def test_unrecognised_shape_refuses(tmp_path):
    vcf = _write(tmp_path, "plain.vcf",
                 "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n")
    with pytest.raises(ValueError):
        mm.detect_vcf_mode(vcf)


def test_csq_without_ps_is_not_the_removed_mode(tmp_path):
    vcf = _write(tmp_path, "csq.vcf",
                 "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Format: A\">\n")
    with pytest.raises(ValueError):
        mm.detect_vcf_mode(vcf)


def test_ps_without_csq_is_not_the_removed_mode(tmp_path):
    vcf = _write(tmp_path, "ps.vcf",
                 "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"p\">\n")
    with pytest.raises(ValueError):
        mm.detect_vcf_mode(vcf)


def test_refusal_message_does_not_mention_compound_het_as_available(tmp_path):
    vcf = _write(tmp_path, "plain.vcf",
                 "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n")
    with pytest.raises(ValueError) as exc:
        mm.detect_vcf_mode(vcf)
    assert "CSQ" not in str(exc.value)


# --- CLI -------------------------------------------------------------------

def test_sv_vcf_dispatches_and_writes_a_report(tmp_path, tiny_vcf):
    out = tmp_path / "out"
    assert mm.main(["--vcf", str(tiny_vcf), "--out", str(out)]) == 0
    assert list(out.glob("*.report.html"))


def test_phased_vep_vcf_exits_one_and_names_the_removal(
    tmp_path, tiny_phased_vcf, capsys,
):
    out = tmp_path / "out"
    rc = mm.main(["--vcf", str(tiny_phased_vcf), "--out", str(out)])
    err = capsys.readouterr().err
    assert rc == 1
    assert "0.5.1" in err
    assert not list(out.glob("*.html"))


def test_compound_het_flags_are_gone():
    parser = mm.build_argparser()
    for flag in ("--gene", "--clinvar", "--min-pair-count", "--max-genes"):
        with pytest.raises(SystemExit):
            parser.parse_args(["--vcf", "x.vcf", "--out", "o", flag, "V"])


def test_help_offers_only_sv_and_karyotype_groups(capsys):
    parser = mm.build_argparser()
    parser.print_help()
    text = capsys.readouterr().out.lower()
    assert "sv-mode flags" in text
    assert "karyotype-mode flags" in text
    assert "compound-het mode flags" not in text
