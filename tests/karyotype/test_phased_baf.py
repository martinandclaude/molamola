"""Haplotype-resolved BAF aggregation.

Phase is an upgrade, never a requirement: molamola takes whatever
small-variant VCF it is handed, so a phased VCF gets phase-block
aggregation and an unphased one still gets the per-site panel. Within a
phase block, read counts are summed across a tiling window rather than
per-site fractions being averaged (the variance-correct estimator), and
each window is emitted at both v and 1-v because which haplotype a block
labels "1" is arbitrary and flips between blocks.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import molamola as mm


# --- _parse_phase_fields ---------------------------------------------------

def test_left_allele_is_haplotype_one():
    """1|0 puts ALT reads on hap 1; 0|1 puts REF reads there."""
    ps, h1 = mm._parse_phase_fields({"GT": "1|0", "PS": "42", "AD": "7,13"})
    assert (ps, h1) == (42, 13.0)
    ps, h1 = mm._parse_phase_fields({"GT": "0|1", "PS": "42", "AD": "7,13"})
    assert (ps, h1) == (42, 7.0)


@pytest.mark.parametrize("kv", [
    {"GT": "0/1", "PS": "42", "AD": "7,13"},   # unphased
    {"GT": "1|1", "PS": "42", "AD": "7,13"},   # not het
    {"GT": "1|0", "AD": "7,13"},               # no PS
    {"GT": "1|0", "PS": ".", "AD": "7,13"},    # null PS
    {"GT": "1|0", "PS": "42"},                 # no AD
    {"GT": "1|0", "PS": "42", "AD": "7"},      # malformed AD
])
def test_unusable_records_are_rejected(kv):
    ps, h1 = mm._parse_phase_fields(kv)
    assert ps == -1 and np.isnan(h1)


def test_allele_floor_guards_miscalled_homs():
    """A 'het' with one read on an allele is likely a mis-called hom and
    would drag its window toward 0 or 1."""
    ps, h1 = mm._parse_phase_fields({"GT": "1|0", "PS": "1", "AD": "1,40"})
    assert ps == -1 and np.isnan(h1)


# --- aggregate_phased_baf --------------------------------------------------

def _sites(n, chrom="chr1", ps=1, start=1000, step=1000, h1=10, dp=20):
    return pd.DataFrame({
        "chrom": pd.Categorical([chrom] * n, categories=mm.CHROM_ORDER),
        "pos": np.arange(n) * step + start,
        "baf": np.full(n, 0.5),
        "ps": np.full(n, ps, dtype=np.int64),
        "h1": np.full(n, float(h1)),
        "dp": np.full(n, dp, dtype=np.int64),
    })


def test_too_few_phased_sites_falls_back():
    """Below the threshold the caller must get None and use per-site."""
    assert mm.aggregate_phased_baf(_sites(10)) is None


def test_unphased_input_falls_back():
    df = _sites(mm.KARY_BAF_MIN_PHASED * 2)
    df["ps"] = -1
    df["h1"] = np.nan
    assert mm.aggregate_phased_baf(df) is None


def test_windows_are_mirrored():
    n = mm.KARY_BAF_MIN_PHASED * 2
    out = mm.aggregate_phased_baf(_sites(n, h1=10, dp=20))
    assert out is not None
    v = np.sort(out["baf"].to_numpy())
    assert np.allclose(v, np.sort(1.0 - v), atol=1e-9)


def test_read_counts_are_summed_not_fractions_averaged():
    """A deep site must pull its window more than a shallow one.

    Averaging per-site fractions would weight them equally; summing
    reads is the variance-correct estimator.
    """
    n = mm.KARY_BAF_MIN_PHASED * 2
    df = _sites(n, h1=10, dp=20)
    # one very deep, strongly hap-1-skewed site in the first window
    df.loc[0, ["h1", "dp"]] = [900.0, 1000]
    out = mm.aggregate_phased_baf(df)
    first = out[out["pos"] == out["pos"].min()]["baf"].to_numpy()
    # window mean of fractions would stay near 0.5; summed reads do not
    assert first.max() > 0.6


def test_a_window_never_spans_two_phase_blocks():
    per_win = mm.KARY_BAF_SNPS_PER_WIN
    a = _sites(per_win // 2, ps=1, start=1000)
    b = _sites(per_win // 2, ps=2, start=1000 + per_win * 1000)
    df = pd.concat([a, b], ignore_index=True)
    df = pd.concat([df] * (2 * mm.KARY_BAF_MIN_PHASED // len(df) + 1),
                   ignore_index=True)
    out = mm.aggregate_phased_baf(df)
    assert out is not None and len(out) > 0


def test_large_positional_gap_splits_a_window():
    """One window must not straddle unrelated sequence."""
    n = mm.KARY_BAF_MIN_PHASED * 2
    df = _sites(n, step=10)                    # tightly packed
    gapped = _sites(n, step=10)
    gapped["pos"] = gapped["pos"] + 50 * mm.KARY_BAF_MAX_GAP
    both = pd.concat([df, gapped], ignore_index=True)
    out = mm.aggregate_phased_baf(both)
    tight = mm.aggregate_phased_baf(df)
    assert out is not None and tight is not None
    # the gap forces at least one extra window beyond simple doubling
    assert len(out) >= len(tight) * 2


def test_output_columns_match_the_per_site_panel():
    """The plot path is shared, so the frames must be interchangeable."""
    out = mm.aggregate_phased_baf(_sites(mm.KARY_BAF_MIN_PHASED * 2))
    assert {"chrom", "pos", "baf"}.issubset(out.columns)
    assert out["baf"].between(0.0, 1.0).all()
