"""Shared fixtures for karyotype-mode tests.

The karyotype-mode test fixtures are all synthetic — depths come from
``numpy.random.default_rng(seed=...)``, cytobands and masks are
hand-written for the four test chromosomes. No real-sample-derived
data ever lands here.
"""

from __future__ import annotations

from pathlib import Path

import pytest


TEST_DATA = Path(__file__).resolve().parent / "data"


@pytest.fixture(scope="session")
def karyotype_data_dir() -> Path:
    """Absolute path to ``tests/karyotype/data/``."""
    return TEST_DATA


@pytest.fixture(scope="session")
def tiny_regions() -> Path:
    """Synthetic mosdepth regions.bed.gz (4 chroms x 50 1 Mb bins)."""
    return TEST_DATA / "tiny_regions.bed.gz"


@pytest.fixture(scope="session")
def tiny_cytoband() -> Path:
    """Synthetic UCSC cytoband (4 chroms x 3 bands: p / acen / q)."""
    return TEST_DATA / "tiny_cytoband.txt.gz"


@pytest.fixture(scope="session")
def tiny_mask() -> Path:
    """Synthetic 3-col mask BED with a handful of intervals."""
    return TEST_DATA / "tiny_mask.bed.gz"


@pytest.fixture(scope="session")
def tiny_gc() -> Path:
    """Synthetic 10 kb GC table at mosdepth-bin midpoints only."""
    return TEST_DATA / "tiny_gc.bed.gz"
