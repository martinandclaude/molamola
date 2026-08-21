"""Exclusion-mask overlap semantics.

The bundled masks are built from 500 bp runs. An any-overlap test makes
the excluded fraction depend on the caller's mosdepth bin size rather
than on the mask: on a 1 kb run it excluded 57.8 % of bins for a mask
covering 39.8 % of hg38, and those bins are dropped from the CN
normalisation anchor as well as from the plot. Bins that are exactly
half masked measure like clean sequence (SD 0.230 log2 vs 0.225 clean,
against 1.072 for fully-masked bins), so "majority masked" is both
bin-size-robust and the empirically correct cut.
"""

from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
import pytest

import molamola as mm


def _cov(chrom: str, n_bins: int, bin_size: int) -> pd.DataFrame:
    starts = np.arange(n_bins, dtype=np.int64) * bin_size
    return pd.DataFrame({
        "chrom": pd.Categorical([chrom] * n_bins, categories=mm.CHROM_ORDER),
        "start": starts,
        "end": starts + bin_size,
        "depth": np.ones(n_bins, dtype=float),
    })


def _mask(tmp_path, intervals, name="m.bed.gz"):
    path = tmp_path / name
    with gzip.open(path, "wt") as fh:
        for s, e in intervals:
            fh.write(f"chr1\t{s}\t{e}\n")
    return path


def test_half_masked_bin_is_kept(tmp_path):
    """A 500 bp mask hit on a 1 kb bin must not drop the whole bin."""
    cov = _cov("chr1", 4, 1000)
    mask = _mask(tmp_path, [(0, 500)])
    out = mm.annotate_mask(cov, mask, 1000, {"chr1": 4000})
    assert out.tolist() == [False, False, False, False]


def test_fully_masked_bin_is_dropped(tmp_path):
    cov = _cov("chr1", 4, 1000)
    mask = _mask(tmp_path, [(0, 1000)])
    out = mm.annotate_mask(cov, mask, 1000, {"chr1": 4000})
    assert out.tolist() == [True, False, False, False]


def test_just_over_half_is_dropped(tmp_path):
    cov = _cov("chr1", 2, 1000)
    mask = _mask(tmp_path, [(0, 501)])
    out = mm.annotate_mask(cov, mask, 1000, {"chr1": 2000})
    assert bool(out[0]) is True


def test_overlap_is_summed_across_intervals(tmp_path):
    """Two separate 300 bp hits on one bin total 600 bp -> majority."""
    cov = _cov("chr1", 1, 1000)
    mask = _mask(tmp_path, [(0, 300), (400, 700)])
    out = mm.annotate_mask(cov, mask, 1000, {"chr1": 1000})
    assert bool(out[0]) is True


def test_long_interval_masks_every_bin_it_spans(tmp_path):
    """The difference-array path for interior bins must not lose bins."""
    cov = _cov("chr1", 10, 1000)
    mask = _mask(tmp_path, [(0, 10_000)])
    out = mm.annotate_mask(cov, mask, 1000, {"chr1": 10_000})
    assert out.all()


def test_at_mask_resolution_behaviour_is_unchanged(tmp_path):
    """At 500 bp bins every hit is total, so the old and new rules agree."""
    cov = _cov("chr1", 4, 500)
    mask = _mask(tmp_path, [(500, 1000)])
    out = mm.annotate_mask(cov, mask, 500, {"chr1": 2000})
    assert out.tolist() == [False, True, False, False]


@pytest.mark.parametrize("bin_size", [500, 1000, 5000])
def test_excluded_fraction_tracks_the_mask_not_the_bin_size(bin_size, tmp_path):
    """The whole point: mask half the contig, get ~half excluded, at any
    bin size. The old any-overlap rule inflated this as bins got coarser."""
    span = 200_000
    cov = _cov("chr1", span // bin_size, bin_size)
    mask = _mask(tmp_path, [(0, span // 2)])
    out = mm.annotate_mask(cov, mask, bin_size, {"chr1": span})
    assert out.mean() == pytest.approx(0.5, abs=0.02)


def test_zero_threshold_restores_any_overlap(tmp_path):
    """--mask-overlap 0 is documented as the pre-v0.5.0 escape hatch."""
    cov = _cov("chr1", 4, 1000)
    mask = _mask(tmp_path, [(0, 100)])
    out = mm.annotate_mask(cov, mask, 1000, {"chr1": 4000},
                           min_overlap_frac=0.0)
    assert bool(out[0]) is True
