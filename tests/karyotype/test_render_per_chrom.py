"""Tests for the per-chromosome 3 x 8 A4-portrait karyotype renderer."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

import molamola as mm


def _make_args(tiny_regions: Path, **overrides) -> argparse.Namespace:
    base = dict(
        mosdepth=tiny_regions,
        reference="hg38",
        scatter_bin_kb=50.0,
        max_points=200_000,
        max_baf_points=80_000,
        smooth_window_mb=0.5,
        ymax=5.0,
        no_mask=False,
        no_gc=False,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def _prepare_cov(cov, bin_size, window_mb=0.5):
    cov2 = cov.copy()
    autosomes = cov2["chrom"].astype(str).isin(mm.KARY_AUTOSOMES)
    median = float(cov2.loc[autosomes, "depth"].median())
    cov2["cn"] = cov2["depth"].to_numpy() * 2.0 / median
    cov2["mask_pass"] = True
    cov2 = mm.rolling_median_per_chrom(
        cov2, bin_size=bin_size, window_mb=window_mb,
        mask_pass=cov2["mask_pass"].to_numpy(),
    )
    return cov2


def test_render_per_chrom_returns_png_bytes(tiny_regions, tiny_cytoband):
    cov, bin_size = mm.read_mosdepth(tiny_regions)
    cov = _prepare_cov(cov, bin_size)
    cb = mm.read_cytoband_df(tiny_cytoband)
    lengths = mm.chrom_lengths_from_cb(cb)
    args = _make_args(tiny_regions)
    png, label = mm.render_karyotype_per_chrom_png(
        cov, cb, lengths, sex="male", bin_size=bin_size, args=args,
    )
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    assert label.endswith(("kb", "Mb", "bp"))


def test_render_per_chrom_closes_figure(tiny_regions, tiny_cytoband):
    plt.close("all")
    cov, bin_size = mm.read_mosdepth(tiny_regions)
    cov = _prepare_cov(cov, bin_size)
    cb = mm.read_cytoband_df(tiny_cytoband)
    lengths = mm.chrom_lengths_from_cb(cb)
    args = _make_args(tiny_regions)
    mm.render_karyotype_per_chrom_png(
        cov, cb, lengths, sex="male", bin_size=bin_size, args=args,
    )
    assert len(plt.get_fignums()) == 0


def test_render_per_chrom_does_not_leak_rcparams(tiny_regions, tiny_cytoband):
    before = (
        plt.rcParams["font.family"],
        plt.rcParams["axes.edgecolor"],
        plt.rcParams["grid.color"],
    )
    cov, bin_size = mm.read_mosdepth(tiny_regions)
    cov = _prepare_cov(cov, bin_size)
    cb = mm.read_cytoband_df(tiny_cytoband)
    lengths = mm.chrom_lengths_from_cb(cb)
    args = _make_args(tiny_regions)
    mm.render_karyotype_per_chrom_png(
        cov, cb, lengths, sex="female", bin_size=bin_size, args=args,
    )
    after = (
        plt.rcParams["font.family"],
        plt.rcParams["axes.edgecolor"],
        plt.rcParams["grid.color"],
    )
    assert before == after


def test_render_per_chrom_refuses_empty_lengths(tiny_regions, tiny_cytoband):
    cov, bin_size = mm.read_mosdepth(tiny_regions)
    cov = _prepare_cov(cov, bin_size)
    cb = mm.read_cytoband_df(tiny_cytoband)
    args = _make_args(tiny_regions)
    with pytest.raises(ValueError, match="no chromosomes"):
        mm.render_karyotype_per_chrom_png(
            cov, cb, lengths={}, sex="male", bin_size=bin_size, args=args,
        )


def test_render_per_chrom_layout_constants():
    assert mm._KARY_PER_CHROM_COLS == 3
    assert mm._KARY_PER_CHROM_ROWS == 8


def test_render_region_panels_returns_bin_label(tiny_regions, tiny_cytoband):
    cov, bin_size = mm.read_mosdepth(tiny_regions)
    cov = _prepare_cov(cov, bin_size)
    cb = mm.read_cytoband_df(tiny_cytoband)
    args = _make_args(tiny_regions)
    fig, axes = plt.subplots(2, 1)
    label = mm._render_kary_region_panels(
        cov, cb, sex="male", bin_size=bin_size,
        chrom="chr1", start=0, end=50_000_000,
        ax_band=axes[0], ax_cov=axes[1], args=args,
    )
    assert label.endswith(("kb", "Mb", "bp"))
    plt.close(fig)
