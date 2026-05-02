"""Tests for find_compound_het_candidates and classify_pairs."""

from __future__ import annotations

import molamola as mm


def test_classify_pairs_counts_trans_and_cis(make_phased_variant):
    a = make_phased_variant(pos=1100, gt="0|1", ps=12345)  # hap 2
    b = make_phased_variant(pos=1500, gt="1|0", ps=12345)  # hap 1
    c = make_phased_variant(pos=1700, gt="1|0", ps=12345)  # hap 1
    n_trans, n_cis = mm.classify_pairs([a, b, c])
    # a-b trans, a-c trans, b-c cis
    assert n_trans == 2
    assert n_cis == 1


def test_classify_pairs_ignores_different_ps(make_phased_variant):
    a = make_phased_variant(pos=1100, gt="0|1", ps=12345)
    b = make_phased_variant(pos=2200, gt="1|0", ps=99999)
    n_trans, n_cis = mm.classify_pairs([a, b])
    assert n_trans == 0
    assert n_cis == 0


def test_classify_pairs_ignores_non_canonical(make_phased_variant):
    """Non-canonical or non-missense don't contribute to title counts."""
    a = make_phased_variant(pos=1100, gt="0|1", ps=12345,
                             is_canonical_transcript=False)
    b = make_phased_variant(pos=1500, gt="1|0", ps=12345)
    c = make_phased_variant(pos=1700, gt="0|1", ps=12345,
                             consequence="synonymous_variant")
    n_trans, n_cis = mm.classify_pairs([a, b, c])
    assert n_trans == 0
    assert n_cis == 0


def test_find_candidates_picks_strict_pair(make_phased_variant, make_gene):
    """Both variants in {p_or_lp, vus} -> strict bucket."""
    genes = {"GENE_A": make_gene()}
    a = make_phased_variant(pos=1100, gt="0|1", ps=12345, clnsig="vus")
    b = make_phased_variant(pos=1500, gt="1|0", ps=12345, clnsig="p_or_lp")
    strict, extended_only = mm.find_compound_het_candidates(
        [a, b], genes, min_pair_count=1,
    )
    assert strict == ["GENE_A"]
    assert extended_only == []


def test_find_candidates_extended_only(make_phased_variant, make_gene):
    """P/LP + conflicting -> extended (anchor + non-benign partner) only."""
    genes = {"GENE_A": make_gene()}
    a = make_phased_variant(pos=1100, gt="0|1", ps=12345, clnsig="p_or_lp")
    b = make_phased_variant(pos=1500, gt="1|0", ps=12345, clnsig="conflicting")
    strict, extended_only = mm.find_compound_het_candidates(
        [a, b], genes, min_pair_count=1,
    )
    assert strict == []
    assert extended_only == ["GENE_A"]


def test_find_candidates_no_clinvar_partner_is_extended(
    make_phased_variant, make_gene,
):
    """P/LP + no-ClinVar -> extended (partner is not 'benign')."""
    genes = {"GENE_A": make_gene()}
    a = make_phased_variant(pos=1100, gt="0|1", ps=12345, clnsig="p_or_lp")
    b = make_phased_variant(pos=1500, gt="1|0", ps=12345, clnsig=None)
    strict, extended_only = mm.find_compound_het_candidates(
        [a, b], genes, min_pair_count=1,
    )
    assert strict == []
    assert extended_only == ["GENE_A"]


def test_find_candidates_excludes_both_benign(
    make_phased_variant, make_gene,
):
    """Benign+benign trans pair: not auto-selectable at all."""
    genes = {"GENE_A": make_gene()}
    a = make_phased_variant(pos=1100, gt="0|1", ps=12345, clnsig="benign")
    b = make_phased_variant(pos=1500, gt="1|0", ps=12345, clnsig="benign")
    strict, extended_only = mm.find_compound_het_candidates(
        [a, b], genes, min_pair_count=1,
    )
    assert strict == []
    assert extended_only == []


def test_find_candidates_excludes_benign_partner(
    make_phased_variant, make_gene,
):
    """P/LP + benign: excluded entirely (partner must be non-benign)."""
    genes = {"GENE_A": make_gene()}
    a = make_phased_variant(pos=1100, gt="0|1", ps=12345, clnsig="p_or_lp")
    b = make_phased_variant(pos=1500, gt="1|0", ps=12345, clnsig="benign")
    strict, extended_only = mm.find_compound_het_candidates(
        [a, b], genes, min_pair_count=1,
    )
    assert strict == []
    assert extended_only == []


def test_find_candidates_excludes_pure_no_clinvar(
    make_phased_variant, make_gene,
):
    """Two no-ClinVar variants in trans: no anchor, excluded."""
    genes = {"GENE_A": make_gene()}
    a = make_phased_variant(pos=1100, gt="0|1", ps=12345, clnsig=None)
    b = make_phased_variant(pos=1500, gt="1|0", ps=12345, clnsig=None)
    strict, extended_only = mm.find_compound_het_candidates(
        [a, b], genes, min_pair_count=1,
    )
    assert strict == []
    assert extended_only == []


def test_find_candidates_excludes_pure_conflicting(
    make_phased_variant, make_gene,
):
    """Two conflicting variants in trans: no anchor, excluded."""
    genes = {"GENE_A": make_gene()}
    a = make_phased_variant(pos=1100, gt="0|1", ps=12345, clnsig="conflicting")
    b = make_phased_variant(pos=1500, gt="1|0", ps=12345, clnsig="conflicting")
    strict, extended_only = mm.find_compound_het_candidates(
        [a, b], genes, min_pair_count=1,
    )
    assert strict == []
    assert extended_only == []


def test_find_candidates_strict_takes_priority(
    make_phased_variant, make_gene,
):
    """Gene with mixed pairs lands in 'strict' if any strict pair exists."""
    genes = {"GENE_A": make_gene()}
    # Strict pair in PS=11
    a = make_phased_variant(pos=1100, gt="0|1", ps=11, clnsig="p_or_lp")
    b = make_phased_variant(pos=1500, gt="1|0", ps=11, clnsig="vus")
    # Extended-only pair in PS=22
    c = make_phased_variant(pos=1700, gt="0|1", ps=22, clnsig="p_or_lp")
    d = make_phased_variant(pos=1900, gt="1|0", ps=22, clnsig=None)
    strict, extended_only = mm.find_compound_het_candidates(
        [a, b, c, d], genes, min_pair_count=1,
    )
    assert strict == ["GENE_A"]
    assert extended_only == []


def test_find_candidates_min_pair_count_threshold(
    make_phased_variant, make_gene,
):
    """min_pair_count=2 rejects a single-qualifying-pair gene."""
    genes = {"GENE_A": make_gene()}
    a = make_phased_variant(pos=1100, gt="0|1", ps=12345, clnsig="p_or_lp")
    b = make_phased_variant(pos=1500, gt="1|0", ps=12345, clnsig="vus")
    strict, extended_only = mm.find_compound_het_candidates(
        [a, b], genes, min_pair_count=2,
    )
    assert strict == []
    assert extended_only == []


def test_find_candidates_excludes_variants_outside_window(
    make_phased_variant, make_gene,
):
    """Variants past gene.end aren't paired even if same PS / opposite hap."""
    genes = {"GENE_A": make_gene(start=1000, end=2000)}
    a = make_phased_variant(pos=1100, gt="0|1", ps=12345, clnsig="p_or_lp")
    b = make_phased_variant(pos=2500, gt="1|0", ps=12345, clnsig="vus")
    strict, extended_only = mm.find_compound_het_candidates(
        [a, b], genes, min_pair_count=1,
    )
    assert strict == []
    assert extended_only == []


def test_find_candidates_sorted_alphabetically(
    make_phased_variant, make_gene,
):
    """Both lists are alphabetically sorted; their union spans both."""
    genes = {
        "ZZZ_GENE": make_gene(symbol="ZZZ_GENE"),
        "AAA_GENE": make_gene(symbol="AAA_GENE"),
        "MMM_GENE": make_gene(symbol="MMM_GENE"),
    }
    # AAA: strict pair
    a1 = make_phased_variant(pos=1100, gt="0|1", ps=11,
                              gene_symbol="AAA_GENE", clnsig="p_or_lp")
    a2 = make_phased_variant(pos=1500, gt="1|0", ps=11,
                              gene_symbol="AAA_GENE", clnsig="vus")
    # ZZZ: strict pair
    z1 = make_phased_variant(pos=1100, gt="0|1", ps=22,
                              gene_symbol="ZZZ_GENE", clnsig="vus")
    z2 = make_phased_variant(pos=1500, gt="1|0", ps=22,
                              gene_symbol="ZZZ_GENE", clnsig="p_or_lp")
    # MMM: extended-only pair (anchor + no-ClinVar)
    m1 = make_phased_variant(pos=1100, gt="0|1", ps=33,
                              gene_symbol="MMM_GENE", clnsig="p_or_lp")
    m2 = make_phased_variant(pos=1500, gt="1|0", ps=33,
                              gene_symbol="MMM_GENE", clnsig=None)
    strict, extended_only = mm.find_compound_het_candidates(
        [z1, z2, m1, m2, a1, a2], genes, min_pair_count=1,
    )
    assert strict == ["AAA_GENE", "ZZZ_GENE"]
    assert extended_only == ["MMM_GENE"]
