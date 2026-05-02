"""Tests for the phased VCF reader."""

from __future__ import annotations

import pytest

import molamola as mm


def test_parses_kept_records(tiny_phased_vcf):
    variants, summary = mm.read_phased_vcf(tiny_phased_vcf)
    # tiny fixture: 8 records, 5 should pass per-record filters.
    assert summary["kept"] == 5
    assert len(variants) == 5
    assert summary["unphased"] == 1
    assert summary["non_het"] == 1
    assert summary["multi_allelic"] == 1
    assert summary["no_canonical_csq"] == 0
    assert summary["has_ps_in_header"] is True


def test_csq_fields_round_trip(tiny_phased_vcf):
    _, summary = mm.read_phased_vcf(tiny_phased_vcf)
    assert summary["csq_fields"] == [
        "Allele", "Consequence", "SYMBOL", "Feature", "HGVSc", "HGVSp",
        "CANONICAL",
    ]


def test_variant_hap_assignment(tiny_phased_vcf):
    variants, _ = mm.read_phased_vcf(tiny_phased_vcf)
    by_pos = {v.pos: v for v in variants}
    # gt=0|1 -> hap 2; gt=1|0 -> hap 1.
    assert by_pos[1100].variant_hap == 2  # 0|1
    assert by_pos[1500].variant_hap == 1  # 1|0
    assert by_pos[2200].variant_hap == 1  # 1|0
    assert by_pos[5000].variant_hap == 2  # 0|1


def test_canonical_transcript_picks_first_yes(tiny_phased_vcf):
    """v1 has two CSQ entries; the YES one (NM_001.1) wins."""
    variants, _ = mm.read_phased_vcf(tiny_phased_vcf)
    v = next(v for v in variants if v.pos == 1100)
    assert v.feature == "NM_001.1"
    assert v.is_canonical_transcript is True


def test_chrom_normalised_to_ucsc(tiny_phased_vcf):
    variants, _ = mm.read_phased_vcf(tiny_phased_vcf)
    assert all(v.chrom.startswith("chr") for v in variants)


def test_refuses_unphased_vcf(tiny_phased_unphased_vcf):
    with pytest.raises(ValueError, match="unphased"):
        mm.read_phased_vcf(tiny_phased_unphased_vcf)


def test_refuses_no_csq_vcf(tiny_phased_no_csq_vcf):
    with pytest.raises(ValueError, match="CSQ"):
        mm.read_phased_vcf(tiny_phased_no_csq_vcf)


def test_clnsig_starts_unset(tiny_phased_vcf):
    """Parser leaves clnsig=None; ClinVar annotation runs as a 2nd pass."""
    variants, _ = mm.read_phased_vcf(tiny_phased_vcf)
    assert all(v.clnsig is None for v in variants)
