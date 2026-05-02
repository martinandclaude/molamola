"""Tests for the VEP CSQ format parser and entry picker."""

from __future__ import annotations

import molamola as mm


def test_parse_csq_format_extracts_field_list():
    line = (
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence '
        'annotations from Ensembl VEP. Format: '
        'Allele|Consequence|SYMBOL|Feature|HGVSp|CANONICAL">'
    )
    fields = mm.parse_csq_format(line)
    assert fields == [
        "Allele", "Consequence", "SYMBOL", "Feature", "HGVSp", "CANONICAL",
    ]


def test_parse_csq_format_returns_none_for_non_csq_line():
    assert mm.parse_csq_format(
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">'
    ) is None
    assert mm.parse_csq_format("#CHROM\tPOS\tID\tREF\tALT") is None


def test_parse_csq_format_resilient_to_field_reordering():
    """Plugin order varies between VEP runs; parser stays index-free."""
    line = (
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: '
        'CANONICAL|HGVSp|Feature|SYMBOL|Consequence|Allele">'
    )
    fields = mm.parse_csq_format(line)
    assert fields == [
        "CANONICAL", "HGVSp", "Feature", "SYMBOL", "Consequence", "Allele",
    ]


def test_pick_canonical_csq_selects_canonical_yes():
    fields = ["Consequence", "SYMBOL", "Feature", "CANONICAL"]
    csq = (
        "missense_variant|GENE_A|NM_001.1|YES,"
        "missense_variant|GENE_A|NM_999.1|"
    )
    entry, is_canonical = mm._pick_canonical_csq(csq, fields)
    assert is_canonical is True
    assert entry["Feature"] == "NM_001.1"


def test_pick_canonical_csq_falls_back_to_first_when_none_canonical():
    fields = ["Consequence", "SYMBOL", "Feature", "CANONICAL"]
    csq = (
        "missense_variant|GENE_A|NM_111.1|,"
        "missense_variant|GENE_A|NM_222.1|"
    )
    entry, is_canonical = mm._pick_canonical_csq(csq, fields)
    assert is_canonical is False
    assert entry["Feature"] == "NM_111.1"


def test_pick_canonical_csq_handles_missing_trailing_fields():
    """VEP sometimes truncates trailing empty fields; parser tolerates it."""
    fields = ["Consequence", "SYMBOL", "Feature", "HGVSp", "CANONICAL"]
    csq = "missense_variant|GENE_A|NM_001.1"  # truncated, no CANONICAL
    entry, is_canonical = mm._pick_canonical_csq(csq, fields)
    assert is_canonical is False
    assert entry["HGVSp"] == ""
    assert entry["CANONICAL"] == ""


def test_is_missense_consequence_handles_compound_terms():
    assert mm.is_missense_consequence("missense_variant") is True
    assert mm.is_missense_consequence(
        "missense_variant&splice_region_variant"
    ) is True
    assert mm.is_missense_consequence(
        "splice_region_variant&missense_variant"
    ) is True
    assert mm.is_missense_consequence("synonymous_variant") is False
    assert mm.is_missense_consequence("") is False
