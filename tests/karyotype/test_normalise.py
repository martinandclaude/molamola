"""Tests for karyotype-mode normalisation + analytics functions."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

import molamola as mm


@pytest.fixture(scope="module")
def cov_and_bin(tiny_regions):
    return mm.read_mosdepth(tiny_regions)


@pytest.fixture(scope="module")
def cb(tiny_cytoband):
    return mm.read_cytoband_df(tiny_cytoband)


@pytest.fixture(scope="module")
def lengths(cb):
    return mm.chrom_lengths_from_cb(cb)


def test_annotate_mask_marks_listed_intervals(cov_and_bin, tiny_mask, lengths):
    cov, bin_size = cov_and_bin
    masked = mm.annotate_mask(cov, tiny_mask, bin_size, lengths)
    assert masked.dtype == bool
    assert len(masked) == len(cov)
    # chr1: mask is 5-7 Mb -> bins 5,6 of chr1 are True
    chr1 = (cov["chrom"] == "chr1").to_numpy()
    chr1_idx = np.where(chr1)[0]
    chr1_masked_indices = np.where(masked & chr1)[0] - chr1_idx[0]
    assert set(chr1_masked_indices.tolist()) == {5, 6}
    # chr2: mask 10-12 Mb -> bins 10,11 of chr2 are True
    chr2 = (cov["chrom"] == "chr2").to_numpy()
    chr2_idx = np.where(chr2)[0]
    chr2_masked_indices = np.where(masked & chr2)[0] - chr2_idx[0]
    assert set(chr2_masked_indices.tolist()) == {10, 11}


def test_annotate_centromere_pad_within_pad(cov_and_bin, cb, lengths):
    cov, _ = cov_and_bin
    # cytoband fixture has acen at 20-22 Mb on every test chrom; pad 1 Mb
    # should cover bins 19..22 (the acen overlap + 1 Mb pad on each side)
    padded = mm.annotate_centromere_pad(cov, cb, lengths, pad_kb=1000.0)
    chr1 = (cov["chrom"] == "chr1").to_numpy()
    chr1_idx = np.where(chr1)[0]
    chr1_padded = np.where(padded & chr1)[0] - chr1_idx[0]
    assert set(chr1_padded.tolist()) == {19, 20, 21, 22}


def test_annotate_centromere_pad_zero_returns_all_false(
    cov_and_bin, cb, lengths,
):
    cov, _ = cov_and_bin
    padded = mm.annotate_centromere_pad(cov, cb, lengths, pad_kb=0.0)
    assert not padded.any()


def test_annotate_gc_returns_int16_with_few_missing(cov_and_bin, tiny_gc):
    cov, _ = cov_and_bin
    gc = mm.annotate_gc(cov, tiny_gc, gc_bin_size=10_000)
    assert gc.dtype == np.int16
    assert len(gc) == len(cov)
    # The fixture emits one GC bin per mosdepth bin midpoint -> every
    # row should resolve to a real value (not GC_MISSING).
    assert (gc != mm.GC_MISSING).all()
    assert gc.min() >= 0
    assert gc.max() <= 100


def test_fit_gc_correction_returns_factors(cov_and_bin, tiny_gc):
    cov, _ = cov_and_bin
    gc = mm.annotate_gc(cov, tiny_gc, gc_bin_size=10_000)
    mask_pass = np.ones(len(cov), dtype=bool)
    factors = mm.fit_gc_correction(
        cov["depth"].to_numpy(), gc, mask_pass, cov["chrom"],
        min_bin_count=1,
    )
    assert isinstance(factors, dict)
    assert len(factors) > 0
    assert all(f > 0 for f in factors.values())


def test_apply_gc_correction_preserves_shape(cov_and_bin, tiny_gc):
    cov, _ = cov_and_bin
    gc = mm.annotate_gc(cov, tiny_gc, gc_bin_size=10_000)
    factors = {int(v): 2.0 for v in np.unique(gc) if v != mm.GC_MISSING}
    corrected = mm.apply_gc_correction(cov["depth"].to_numpy(), gc, factors)
    assert corrected.shape == cov["depth"].shape
    np.testing.assert_allclose(
        corrected, cov["depth"].to_numpy().astype(np.float64) * 2.0,
        rtol=1e-6,
    )


def test_apply_gc_correction_empty_factors_is_identity(cov_and_bin, tiny_gc):
    cov, _ = cov_and_bin
    gc = mm.annotate_gc(cov, tiny_gc, gc_bin_size=10_000)
    out = mm.apply_gc_correction(cov["depth"].to_numpy(), gc, factors={})
    np.testing.assert_allclose(out, cov["depth"].to_numpy().astype(np.float64))


def test_detect_sex_from_cn_male_when_chry_high():
    chroms = pd.Series(["chr1"] * 50 + ["chrY"] * 50)
    cn = np.concatenate([np.full(50, 2.0), np.full(50, 1.0)])
    mask = np.ones(100, dtype=bool)
    assert mm.detect_sex_from_cn(cn, chroms, mask) == "male"


def test_detect_sex_from_cn_female_when_chry_low():
    chroms = pd.Series(["chr1"] * 50 + ["chrY"] * 50)
    cn = np.concatenate([np.full(50, 2.0), np.full(50, 0.05)])
    mask = np.ones(100, dtype=bool)
    assert mm.detect_sex_from_cn(cn, chroms, mask) == "female"


def test_detect_sex_from_cn_no_y_defaults_female():
    chroms = pd.Series(["chr1"] * 50)
    cn = np.full(50, 2.0)
    mask = np.ones(50, dtype=bool)
    assert mm.detect_sex_from_cn(cn, chroms, mask) == "female"


@pytest.mark.parametrize("chrom,sex,expected", [
    ("chr1", "male", 2.0),
    ("chr15", "female", 2.0),
    ("chrX", "male", 1.0),
    ("chrX", "female", 2.0),
    ("chrY", "male", 1.0),
    ("chrY", "female", 0.0),
])
def test_expected_copy_number(chrom, sex, expected):
    assert mm.expected_copy_number(chrom, sex) == expected


def test_rolling_median_per_chrom_smooths_deletion(cov_and_bin):
    cov, bin_size = cov_and_bin
    # Convert depth -> CN using a global autosomal median anchor of 30
    cn = cov["depth"].to_numpy().astype(np.float64) * 2.0 / 30.0
    cov2 = cov.copy()
    cov2["cn"] = cn
    mask_pass = np.ones(len(cov2), dtype=bool)
    out = mm.rolling_median_per_chrom(
        cov2, bin_size=bin_size, window_mb=3.0, mask_pass=mask_pass,
    )
    assert "smooth" in out.columns
    assert len(out) == len(cov2)
    # chr1 has a 5 Mb het deletion centred ~22 Mb (bins 20-24)
    chr1_mask = (out["chrom"] == "chr1") & out["start"].between(
        20_000_000, 24_000_000,
    )
    in_del = out.loc[chr1_mask, "smooth"].dropna()
    assert (in_del < 1.6).any()


def test_aggregate_for_scatter_reduces_rows(cov_and_bin):
    cov, bin_size = cov_and_bin
    cov2 = cov.copy()
    cov2["cn"] = cov2["depth"].to_numpy() * 2.0 / 30.0
    cov2["xpos"] = cov2["start"].to_numpy().astype(np.float64)
    cov2["mask_pass"] = True
    agg = mm.aggregate_for_scatter(cov2, factor=5)
    # 4 chroms x (50 // 5 == 10) groups = ~40 rows
    assert len(agg) == 40
    assert {"chrom", "xpos", "cn", "mask_pass"} <= set(agg.columns)
    assert agg["mask_pass"].all()


def test_aggregate_for_scatter_factor_one_returns_unmasked(cov_and_bin):
    cov, bin_size = cov_and_bin
    cov2 = cov.copy()
    cov2["cn"] = cov2["depth"].to_numpy() * 2.0 / 30.0
    cov2["xpos"] = cov2["start"].to_numpy().astype(np.float64)
    cov2["mask_pass"] = True
    out = mm.aggregate_for_scatter(cov2, factor=1)
    assert len(out) == len(cov2)


def test_downsample_systematic_caps_rows():
    df = pd.DataFrame({"x": np.arange(1000), "y": np.arange(1000)})
    out = mm.downsample_systematic(df, max_points=100)
    assert len(out) <= 100
    assert out["x"].iloc[0] == 0


def test_downsample_systematic_passthrough_below_cap():
    df = pd.DataFrame({"x": np.arange(50)})
    out = mm.downsample_systematic(df, max_points=100)
    assert len(out) == 50
