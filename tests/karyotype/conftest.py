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
