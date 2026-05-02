"""Shared fixtures for compound-het tests.

Mirrors the factory pattern in ``tests/conftest.py`` (``make_bnd`` /
``make_sv``) for the compound-het dataclasses, plus session-scoped
fixtures pointing at the synthetic VCF / canonical-exon / ClinVar
files in ``tests/compound_het/data/``.
"""

from __future__ import annotations

from pathlib import Path

import pytest


CH_DATA = Path(__file__).resolve().parent / "data"


@pytest.fixture(scope="session")
def ch_data_dir() -> Path:
    """Absolute path to ``tests/compound_het/data/``."""
    return CH_DATA


@pytest.fixture(scope="session")
def tiny_phased_vcf() -> Path:
    """Synthetic 8-record phased VCF.

    Contents (one PATIENT123 sample column, all variants on chr1):
    - v1, v2: missense canonical in PS=12345, opposite haps (trans pair)
    - v3: synonymous canonical in PS=12345 (non-missense tick coverage)
    - v4: missense canonical in PS=99999 (singleton, no pair)
    - v5: missense canonical in GENE_B PS=55555 (single hit)
    - v6: unphased GT (refused per-record)
    - v7: multi-allelic ALT (refused per-record)
    - v8: homozygous 1|1 (refused per-record)
    """
    return CH_DATA / "tiny_phased.vcf"


@pytest.fixture(scope="session")
def tiny_phased_unphased_vcf() -> Path:
    """Phased-VCF stand-in with no PS in header (whole-file refusal)."""
    return CH_DATA / "tiny_phased_unphased.vcf"


@pytest.fixture(scope="session")
def tiny_phased_no_csq_vcf() -> Path:
    """Phased-VCF stand-in with no CSQ in header (whole-file refusal)."""
    return CH_DATA / "tiny_phased_no_csq.vcf"


@pytest.fixture(scope="session")
def tiny_canonical_exons() -> Path:
    """Synthetic 3-gene canonical-exon TSV (gzipped).

    GENE_A chr1:1000-3000 (+) with 4 exons
    GENE_B chr1:5000-7000 (-) with 2 exons
    GENE_C chr1:10000-12000 (+) with 3 exons
    """
    return CH_DATA / "tiny_canonical_exons.tsv.gz"


@pytest.fixture(scope="session")
def tiny_clinvar() -> Path:
    """Synthetic ClinVar in molamola's reduced TSV format (5 cols).

    Mirrors the bundled-data schema (chrom, pos, ref, alt, bucket).
    Pre-computed buckets:

    chr1:1100 A>T -> p_or_lp (matches v1 in tiny_phased.vcf)
    chr1:1500 G>C -> vus     (matches v2; together with v1 = strict
                              trans pair under the auto-select rule)
    chr1:2200 T>G -> benign  (matches v4)
    chr1:5000 A>G -> vus     (matches v5)

    Other phased hets in the test VCF stay unmatched (no-ClinVar grey).
    """
    return CH_DATA / "tiny_clinvar.tsv.gz"


@pytest.fixture(scope="session")
def tiny_clinvar_vcf() -> Path:
    """Synthetic ClinVar in NCBI raw VCF format.

    Same data as ``tiny_clinvar`` but in upstream VCF shape. Used by
    the dispatch test that exercises the VCF code path of the
    auto-dispatching loader.
    """
    return CH_DATA / "tiny_clinvar.vcf"


@pytest.fixture
def make_phased_variant():
    """Factory for :class:`PhasedVariant` instances with sensible defaults.

    Mirrors ``make_bnd`` / ``make_sv`` in ``tests/conftest.py``.
    """
    import molamola as mm

    def _make(
        chrom: str = "chr1",
        pos: int = 1000,
        ref: str = "A",
        alt: str = "T",
        gt: str = "0|1",
        ps: int = 12345,
        consequence: str = "missense_variant",
        is_canonical_transcript: bool = True,
        hgvs_p: str | None = "p.Asn4Ile",
        hgvs_c: str | None = "NM_001.1:c.10A>T",
        gene_symbol: str | None = "GENE_A",
        feature: str | None = "NM_001.1",
        clnsig: str | None = None,
    ) -> "mm.PhasedVariant":
        variant_hap = 1 if gt == "1|0" else 2
        return mm.PhasedVariant(
            chrom=chrom, pos=pos, ref=ref, alt=alt, gt=gt, ps=ps,
            variant_hap=variant_hap, consequence=consequence,
            is_canonical_transcript=is_canonical_transcript,
            hgvs_p=hgvs_p, hgvs_c=hgvs_c,
            gene_symbol=gene_symbol, feature=feature, clnsig=clnsig,
        )

    return _make


@pytest.fixture
def make_gene():
    """Factory for :class:`Gene` instances with sensible defaults."""
    import molamola as mm

    def _make(
        symbol: str = "GENE_A",
        chrom: str = "chr1",
        start: int = 1000,
        end: int = 3000,
        strand: str = "+",
        transcript_id: str = "NM_001.1",
        canonical_exons: tuple = ((1100, 1200), (1500, 1600)),
    ) -> "mm.Gene":
        return mm.Gene(
            symbol=symbol, chrom=chrom, start=start, end=end,
            strand=strand, transcript_id=transcript_id,
            canonical_exons=canonical_exons,
        )

    return _make
