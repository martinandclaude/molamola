"""Per-caller VCF parsing round-trip tests (Phase B)."""

from __future__ import annotations

import pytest

import molamola as mm


# --- Sniffles2 (the canonical case) ----------------------------------------

def test_caller_sniffles2_round_trip(tiny_vcf):
    contigs, bnds, svs, _ = mm.read_vcf(tiny_vcf, caller="sniffles2")
    assert contigs["chr1"] == 248956422
    assert len(bnds) == 6
    # First INS in tiny.vcf has VAF 0.500 in INFO
    ins_records = [s for s in svs if s.svtype == "INS"]
    assert ins_records[0].vaf == 0.500


# --- Sniffles1 -------------------------------------------------------------

def test_caller_sniffles1_round_trip(repo_root):
    vcf = repo_root / "tests" / "data" / "tiny_sniffles1.vcf"
    contigs, _bnds, svs, _ = mm.read_vcf(vcf, caller="sniffles1")
    # Two records in the fixture (DEL + INS); SUPPORT comes from INFO
    assert {s.svtype for s in svs} == {"DEL", "INS"}
    # VAF defaults to 0 when no DR/DV available (the fixture has no FORMAT data)
    assert all(s.vaf == 0.0 for s in svs)


def test_caller_sniffles1_derives_vaf_from_dr_dv(tmp_path):
    """Build a Sniffles1-style VCF with DR/DV in FORMAT and check VAF."""
    vcf = tmp_path / "s1.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##source=Sniffles\n"
        "##contig=<ID=chr1,length=248956422>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">\n"
        "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"\">\n"
        "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"\">\n"
        "##FILTER=<ID=PASS,Description=\"\">\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">\n"
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"\">\n"
        "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t1000000\tS1.1\tN\t<DEL>\t60\tPASS\t"
        "SVTYPE=DEL;SVLEN=-100;END=1000100;SUPPORT=20\tGT:DR:DV\t0/1:10:30\n"
    )
    _, _, svs, _ = mm.read_vcf(vcf, caller="sniffles1")
    assert len(svs) == 1
    # DR=10, DV=30 -> VAF = 30/40 = 0.75
    assert abs(svs[0].vaf - 0.75) < 1e-6
    assert svs[0].support == 20  # from INFO/SUPPORT


# --- cuteSV ----------------------------------------------------------------

def test_caller_cutesv_uses_re_and_af(repo_root):
    vcf = repo_root / "tests" / "data" / "tiny_cutesv.vcf"
    _, _, svs, _ = mm.read_vcf(vcf, caller="cutesv")
    assert len(svs) == 2
    by_type = {s.svtype: s for s in svs}
    # cuteSV.DEL.1: RE=10, AF=0.5
    assert by_type["DEL"].support == 10
    assert by_type["DEL"].vaf == 0.5
    # cuteSV.INS.1: RE=12, AF=0.45
    assert by_type["INS"].support == 12
    assert by_type["INS"].vaf == 0.45


# --- SVIM ------------------------------------------------------------------

def test_caller_svim_uses_support_and_af(repo_root):
    vcf = repo_root / "tests" / "data" / "tiny_svim.vcf"
    _, _, svs, _ = mm.read_vcf(vcf, caller="svim")
    assert len(svs) == 1
    s = svs[0]
    assert s.svtype == "DEL"
    assert s.support == 10
    # SVIM doesn't have AF in our fixture -> defaults to 0
    assert s.vaf == 0.0


# --- pbsv ------------------------------------------------------------------

def test_caller_pbsv_derives_from_format_ad(repo_root):
    vcf = repo_root / "tests" / "data" / "tiny_pbsv.vcf"
    _, _, svs, _ = mm.read_vcf(vcf, caller="pbsv")
    assert len(svs) == 2
    by_type = {s.svtype: s for s in svs}
    # DEL: AD=10,8 -> support=8, VAF=8/18 ~ 0.444
    assert by_type["DEL"].support == 8
    assert abs(by_type["DEL"].vaf - 8 / 18) < 1e-6
    # INS: AD=5,15 -> support=15, VAF=15/20=0.75
    assert by_type["INS"].support == 15
    assert abs(by_type["INS"].vaf - 0.75) < 1e-6


# --- NanoVar ---------------------------------------------------------------

def test_caller_nanovar_uses_sr_and_af(repo_root):
    vcf = repo_root / "tests" / "data" / "tiny_nanovar.vcf"
    _, _, svs, _ = mm.read_vcf(vcf, caller="nanovar")
    assert len(svs) == 2
    by_type = {s.svtype: s for s in svs}
    assert by_type["DEL"].support == 12
    assert by_type["DEL"].vaf == 0.4
    assert by_type["INS"].support == 15
    assert by_type["INS"].vaf == 0.5


# --- auto-dispatch agreement -----------------------------------------------

@pytest.mark.parametrize("fixture, expected_caller", [
    ("tiny.vcf",          "sniffles2"),
    ("tiny_sniffles1.vcf", "sniffles1"),
    ("tiny_cutesv.vcf",    "cutesv"),
    ("tiny_svim.vcf",      "svim"),
    ("tiny_pbsv.vcf",      "pbsv"),
    ("tiny_nanovar.vcf",   "nanovar"),
])
def test_caller_auto_matches_explicit(repo_root, fixture, expected_caller):
    """`caller='auto'` and `caller=<expected>` produce identical objects."""
    path = repo_root / "tests" / "data" / fixture
    auto = mm.read_vcf(path, caller="auto")
    explicit = mm.read_vcf(path, caller=expected_caller)
    # Same lengths, same first-record fields
    assert len(auto[1]) == len(explicit[1])
    assert len(auto[2]) == len(explicit[2])


# --- error paths -----------------------------------------------------------

def test_caller_unknown_raises(tiny_vcf):
    with pytest.raises(ValueError, match="unknown caller"):
        mm.read_vcf(tiny_vcf, caller="not-a-caller")
