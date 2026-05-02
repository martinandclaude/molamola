"""Unit tests for VCF parsing helpers."""

from __future__ import annotations

import pytest

import molamola as mm


# --- parse_alt_for_mate -----------------------------------------------------

@pytest.mark.parametrize(
    "alt, expected",
    [
        ("N[chr2:50000[",   ("chr2", 50000, "++")),
        ("N]chr2:50000]",   ("chr2", 50000, "+-")),
        ("[chr2:50000[N",   ("chr2", 50000, "-+")),
        ("]chr2:50000]N",   ("chr2", 50000, "--")),
        ("N[chrX:1[",       ("chrX", 1, "++")),
    ],
)
def test_parse_alt_for_mate_valid(alt, expected):
    assert mm.parse_alt_for_mate(alt) == expected


@pytest.mark.parametrize("alt", ["A", "GATC", "N", "<DEL>", "N<chr2:50000>"])
def test_parse_alt_for_mate_invalid_raises(alt):
    with pytest.raises(ValueError):
        mm.parse_alt_for_mate(alt)


# --- parse_info -------------------------------------------------------------

def test_parse_info_basic():
    assert mm.parse_info("SVTYPE=BND;SUPPORT=10;VAF=0.5") == {
        "SVTYPE": "BND", "SUPPORT": "10", "VAF": "0.5",
    }


def test_parse_info_flag_only():
    out = mm.parse_info("PRECISE;SVTYPE=INS;SUPPORT=20")
    assert out["PRECISE"] is True
    assert out["SVTYPE"] == "INS"
    assert out["SUPPORT"] == "20"


def test_parse_info_value_with_equals():
    # values may legitimately contain '=' (e.g. embedded URLs) - we split
    # only on the first '='
    out = mm.parse_info("KEY=a=b=c")
    assert out == {"KEY": "a=b=c"}


# --- parse_coverage ---------------------------------------------------------

def test_parse_coverage_five_values():
    assert mm.parse_coverage("30,30,30,30,30") == [30.0]*5


def test_parse_coverage_with_nulls():
    out = mm.parse_coverage("30,null,30,.,30")
    assert out == [30.0, None, 30.0, None, 30.0]


def test_parse_coverage_empty():
    assert mm.parse_coverage("") == []
    assert mm.parse_coverage(None) == []


def test_parse_coverage_invalid_token_becomes_none():
    out = mm.parse_coverage("30,abc,30")
    assert out == [30.0, None, 30.0]


# --- BND canonical key + dedup ---------------------------------------------

def test_bnd_canonical_key_swaps_to_same_value(make_bnd):
    a = make_bnd(chr1="chr1", pos1=100, chr2="chr2", pos2=200)
    b = make_bnd(chr1="chr2", pos1=200, chr2="chr1", pos2=100)
    assert a.canonical_key == b.canonical_key


def test_deduplicate_reciprocal_collapses_pair(make_bnd):
    a = make_bnd(chr1="chr1", pos1=100, chr2="chr2", pos2=200, sv_id="A")
    b = make_bnd(chr1="chr2", pos1=200, chr2="chr1", pos2=100, sv_id="B")
    c = make_bnd(chr1="chr3", pos1=300, chr2="chr5", pos2=400, sv_id="C")
    out = mm.deduplicate_reciprocal([a, b, c])
    assert len(out) == 2
    # First record per canonical key wins
    assert {x.sv_id for x in out} == {"A", "C"}


# --- read_vcf end-to-end on tiny.vcf ----------------------------------------

def test_read_vcf_counts_match_tiny_fixture(tiny_vcf):
    contigs, bnds, svs, median_cov = mm.read_vcf(tiny_vcf)

    # Contigs from header
    assert contigs["chr1"] == 248956422
    assert contigs["chrY"] == 57227415
    assert len(contigs) == 9

    # BND count: 6 records in tiny.vcf
    assert len(bnds) == 6

    # Non-BND SVs: 5 INS + 4 DEL + 2 DUP + 1 INV = 12
    by_type = {}
    for s in svs:
        by_type[s.svtype] = by_type.get(s.svtype, 0) + 1
    assert by_type == {"INS": 5, "DEL": 4, "DUP": 2, "INV": 1}

    # Baseline median cov: nonbnd_center_cov collected from INS+DEL center
    # values. INS centers all 30; DEL centers are 15, 12, 15, 15.
    # The list is [30,30,30,30,30, 15,12,200,15] -> 9 values, median = 30.
    # (DEL.3_NOISE has center=200; its center contributes too because
    # we don't pre-filter on noise when building the baseline.)
    assert 25.0 <= median_cov <= 35.0


def test_read_vcf_bnd_alt_orientation_round_trip(tiny_vcf):
    _, bnds, _, _ = mm.read_vcf(tiny_vcf)
    by_id = {b.sv_id: b for b in bnds}
    # Picked deliberately to cover all 4 bracket forms
    assert by_id["BND.1"].orientation == "++"   # N[chr2:50000000[
    assert by_id["BND.2"].orientation == "--"   # ]chr1:73000000]N
    assert by_id["BND.3_NOISE"].orientation == "+-"  # N]chr14:30000000]
    assert by_id["BND.5_GT"].orientation == "-+"  # [chrX:50000000[N
    assert by_id["BND.5_GT"].filter_ == "GT"


def test_read_vcf_sv_lengths_are_absolute(tiny_vcf):
    _, _, svs, _ = mm.read_vcf(tiny_vcf)
    for s in svs:
        assert s.svlen >= 0, f"{s.sv_id}: SVLEN must be stored as |value|"


# --- _normalize_chrom ------------------------------------------------------

@pytest.mark.parametrize("name, expected", [
    # Canonical bare → prepend chr
    ("1", "chr1"), ("2", "chr2"), ("22", "chr22"), ("X", "chrX"), ("Y", "chrY"),
    # Already prefixed → leave alone
    ("chr1", "chr1"), ("chr22", "chr22"), ("chrX", "chrX"),
    # Non-canonical → leave alone
    ("chrM", "chrM"), ("MT", "MT"), ("GL000220.1", "GL000220.1"),
    ("chrEBV", "chrEBV"), ("23", "23"),
])
def test_normalize_chrom(name, expected):
    assert mm._normalize_chrom(name) == expected


def test_read_vcf_normalises_bare_contig_names(tmp_path):
    """A VCF that uses '1' / '2' / 'X' instead of 'chr1' / ... still parses."""
    vcf = tmp_path / "ensembl_style.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##source=Sniffles2_2.7.5\n"
        "##contig=<ID=1,length=248956422>\n"
        "##contig=<ID=X,length=156040895>\n"
        "##FILTER=<ID=PASS,Description=\"\">\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">\n"
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">\n"
        "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"\">\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"\">\n"
        "##INFO=<ID=VAF,Number=1,Type=Float,Description=\"\">\n"
        "##INFO=<ID=COVERAGE,Number=.,Type=Float,Description=\"\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t1000\tINS.1\tN\tGATC\t60\tPASS\t"
        "SVTYPE=INS;SVLEN=4;END=1000;SUPPORT=10;VAF=0.5;COVERAGE=30,30,30,30,30\n"
    )
    contigs, _bnds, svs, _ = mm.read_vcf(vcf, caller="sniffles2")
    assert "chr1" in contigs and "chrX" in contigs
    assert contigs["chr1"] == 248956422
    assert svs[0].chrom == "chr1"
