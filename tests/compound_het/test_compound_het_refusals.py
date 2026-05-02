"""Refusal-path tests for compound-het mode."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import molamola as mm


def _ch_args(vcf, exons, clinvar, out, *extra):
    return [
        
        "--vcf", str(vcf),
        "--reference", "hg38",
        "--canonical-exons", str(exons),
        "--clinvar", str(clinvar),
        "--out", str(out),
        *extra,
    ]


def test_refuses_unphased_vcf(
    capsys, tmp_path, tiny_phased_unphased_vcf,
    tiny_canonical_exons, tiny_clinvar,
):
    rc = mm.main(_ch_args(
        tiny_phased_unphased_vcf, tiny_canonical_exons, tiny_clinvar, tmp_path,
    ))
    assert rc == 1
    err = capsys.readouterr().err
    assert "unphased" in err


def test_refuses_no_csq_vcf(
    capsys, tmp_path, tiny_phased_no_csq_vcf,
    tiny_canonical_exons, tiny_clinvar,
):
    rc = mm.main(_ch_args(
        tiny_phased_no_csq_vcf, tiny_canonical_exons, tiny_clinvar, tmp_path,
    ))
    assert rc == 1
    err = capsys.readouterr().err
    assert "CSQ" in err


def test_refuses_unknown_gene(
    capsys, tmp_path, tiny_phased_vcf,
    tiny_canonical_exons, tiny_clinvar,
):
    rc = mm.main(_ch_args(
        tiny_phased_vcf, tiny_canonical_exons, tiny_clinvar, tmp_path,
        "--gene", "NOPE",
    ))
    assert rc == 1
    err = capsys.readouterr().err
    assert "NOPE" in err
    assert "not found" in err


def test_empty_gene_exits_zero_with_message(
    capsys, tmp_path, tiny_phased_vcf,
    tiny_canonical_exons, tiny_clinvar,
):
    """--gene GENE_C has no variants in window: exit 0 + placeholder + msg."""
    rc = mm.main(_ch_args(
        tiny_phased_vcf, tiny_canonical_exons, tiny_clinvar, tmp_path,
        "--gene", "GENE_C",
    ))
    assert rc == 0
    err = capsys.readouterr().err
    assert "GENE_C: 0 phased hets" in err


def test_auto_select_zero_candidates_exits_zero(
    capsys, tmp_path, tiny_phased_vcf,
    tiny_canonical_exons, tiny_clinvar,
):
    """When the threshold isn't met, auto-select returns 0 genes (exit 0).

    Synthetic GENE_A has exactly one qualifying trans pair (P/LP +
    VUS in PS=12345); pushing min_pair_count to 2 falsifies the
    rule and we expect a clean empty result.
    """
    rc = mm.main(_ch_args(
        tiny_phased_vcf, tiny_canonical_exons, tiny_clinvar, tmp_path,
        "--min-pair-count", "2",
    ))
    assert rc == 0
    err = capsys.readouterr().err
    assert "0 candidate genes" in err


def test_max_genes_caps_auto_selection(
    capsys, tmp_path, tiny_phased_vcf,
    tiny_canonical_exons, tiny_clinvar,
):
    """With --min-pair-count=1, GENE_A passes; --max-genes 0 caps everything."""
    rc = mm.main(_ch_args(
        tiny_phased_vcf, tiny_canonical_exons, tiny_clinvar, tmp_path,
        "--min-pair-count", "1",
        "--max-genes", "0",
    ))
    assert rc == 0
    err = capsys.readouterr().err
    assert "capping at 0" in err
