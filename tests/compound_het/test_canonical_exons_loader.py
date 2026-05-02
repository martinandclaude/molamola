"""Tests for the canonical-exon TSV loader."""

from __future__ import annotations

import gzip

import pytest

import molamola as mm


def test_load_canonical_exons_round_trip(tiny_canonical_exons):
    genes = mm.load_canonical_exons(tiny_canonical_exons)
    assert set(genes.keys()) == {"GENE_A", "GENE_B", "GENE_C"}
    g = genes["GENE_A"]
    assert g.chrom == "chr1"
    assert g.start == 1000
    assert g.end == 3000
    assert g.strand == "+"
    assert g.transcript_id == "NM_001.1"
    assert g.canonical_exons == (
        (1100, 1200), (1500, 1600), (1800, 1900), (2200, 2300),
    )


def test_chrom_normalised(tmp_path):
    """Bare contig names (Ensembl style) get the chr prefix on load."""
    p = tmp_path / "ensembl_style.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(
            "gene_symbol\tchrom\tstart\tend\tstrand\ttranscript_id\t"
            "exon_starts\texon_ends\n"
        )
        fh.write("GENE_X\t1\t100\t200\t+\tNM_X.1\t110\t190\n")
    genes = mm.load_canonical_exons(p)
    assert genes["GENE_X"].chrom == "chr1"


def test_missing_required_column_raises(tmp_path):
    p = tmp_path / "broken.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("gene_symbol\tchrom\tstart\tend\n")  # missing strand+transcript+exons
        fh.write("GENE_X\tchr1\t100\t200\n")
    with pytest.raises(ValueError, match="missing columns"):
        mm.load_canonical_exons(p)


def test_empty_exon_lists_allowed(tmp_path):
    """Genes with no canonical exons get () — panel still renders."""
    p = tmp_path / "no_exons.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(
            "gene_symbol\tchrom\tstart\tend\tstrand\ttranscript_id\t"
            "exon_starts\texon_ends\n"
        )
        fh.write("GENE_NOEXONS\tchr1\t100\t200\t+\tNM_NX.1\t\t\n")
    genes = mm.load_canonical_exons(p)
    assert genes["GENE_NOEXONS"].canonical_exons == ()


def test_duplicate_symbols_keep_first(tmp_path):
    p = tmp_path / "dup.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(
            "gene_symbol\tchrom\tstart\tend\tstrand\ttranscript_id\t"
            "exon_starts\texon_ends\n"
        )
        fh.write("GENE_X\tchr1\t100\t200\t+\tNM_FIRST.1\t110\t190\n")
        fh.write("GENE_X\tchr1\t300\t400\t+\tNM_SECOND.1\t310\t390\n")
    genes = mm.load_canonical_exons(p)
    assert genes["GENE_X"].transcript_id == "NM_FIRST.1"


def test_find_canonical_exon_file_missing_bundled(monkeypatch, tmp_path):
    """Without a bundled file, the loader raises a clear error."""
    monkeypatch.setattr(mm, "__file__", str(tmp_path / "molamola.py"))
    with pytest.raises(FileNotFoundError, match="bundled"):
        mm.find_canonical_exon_file()
