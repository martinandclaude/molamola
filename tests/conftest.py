"""Shared pytest fixtures.

The repo root is on ``sys.path`` (configured in ``pytest.ini``), so
tests can ``import molamola as mm`` directly.
"""

from __future__ import annotations

from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent
TEST_DATA = Path(__file__).resolve().parent / "data"
PACKAGE_DATA = REPO_ROOT / "molamola" / "data"


@pytest.fixture(scope="session")
def repo_root() -> Path:
    """Absolute path to the repository root (one level up from tests/)."""
    return REPO_ROOT


@pytest.fixture(scope="session")
def package_data_dir() -> Path:
    """Absolute path to the bundled-refs directory inside the package."""
    return PACKAGE_DATA


@pytest.fixture(scope="session")
def tiny_vcf() -> Path:
    """Synthetic 18-record Sniffles2 VCF used by parser and regression tests."""
    return TEST_DATA / "tiny.vcf"


@pytest.fixture(scope="session")
def tiny_phased_vcf() -> Path:
    """Minimal phased + VEP-annotated VCF.

    Kept after compound-het mode was removed so the dispatch tests can
    still assert that this header shape is recognised and refused with
    a message naming the removal.
    """
    return TEST_DATA / "tiny_phased.vcf"


@pytest.fixture(scope="session")
def bundled_cytoband() -> Path:
    """The hg38 cytoBand.txt.gz shipped with the package."""
    return PACKAGE_DATA / "cytoBand.txt.gz"


@pytest.fixture
def make_bnd():
    """Factory for constructing minimal :class:`BND` instances in unit tests.

    Returns a function that builds a BND with sensible defaults; the
    caller can override any attribute by keyword argument.
    """
    import molamola as mm

    def _make(
        chr1: str = "chr1",
        pos1: int = 1000,
        chr2: str = "chr2",
        pos2: int = 2000,
        orientation: str = "++",
        support: int = 10,
        vaf: float = 0.5,
        filter_: str = "PASS",
        sv_id: str = "BND.test",
        coverage: list | None = None,
    ) -> "mm.BND":
        return mm.BND(
            chr1=chr1, pos1=pos1, chr2=chr2, pos2=pos2,
            orientation=orientation, support=support, vaf=vaf,
            filter_=filter_, sv_id=sv_id,
            coverage=coverage if coverage is not None else [30.0]*5,
        )

    return _make


@pytest.fixture
def make_sv():
    """Factory for constructing minimal :class:`SV` instances in unit tests."""
    import molamola as mm

    def _make(
        chrom: str = "chr1",
        start: int = 1000,
        end: int | None = None,
        svtype: str = "DEL",
        svlen: int = 500,
        filter_: str = "PASS",
        support: int = 10,
        vaf: float = 0.5,
        sv_id: str = "SV.test",
        coverage: list | None = None,
    ) -> "mm.SV":
        if end is None:
            end = start + svlen
        return mm.SV(
            chrom=chrom, start=start, end=end,
            svtype=svtype, svlen=svlen,
            filter_=filter_, support=support, vaf=vaf, sv_id=sv_id,
            coverage=coverage if coverage is not None else [30.0]*5,
        )

    return _make
