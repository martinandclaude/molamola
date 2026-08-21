"""Tests for the genome-wide karyotype renderer."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd

import molamola as mm


def _make_args(tiny_regions: Path, **overrides) -> argparse.Namespace:
    """Stub Namespace with the karyotype-mode attrs the renderer reads."""
    base = dict(
        mosdepth=tiny_regions,
        reference="hg38",
        scatter_bin_kb=50.0,
        max_points=200_000,
        max_baf_points=80_000,
        smooth_window_mb=0.5,
        ymax=1.5,
        ymin=-2.0,
        no_mask=False,
        no_gc=False,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def _prepare_cov(cov, bin_size, window_mb=0.5):
    """Add `cn`, `mask_pass`, and `smooth` columns to a tiny-fixture cov frame.

    Uses the autosome-only median as the CN-anchor and assumes nothing
    is masked. Sufficient for renderer smoke tests; the real
    karyotype_main wires this through annotate_mask + GC correction.
    """
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


def test_render_genome_returns_png_bytes(tiny_regions, tiny_cytoband):
    cov, bin_size = mm.read_mosdepth(tiny_regions)
    cov = _prepare_cov(cov, bin_size)
    cb = mm.read_cytoband_df(tiny_cytoband)
    lengths = mm.chrom_lengths_from_cb(cb)
    args = _make_args(tiny_regions)
    png, label = mm.render_karyotype_genome_png(
        cov, cb, lengths, baf_df=None, sex="male",
        bin_size=bin_size, args=args,
    )
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    assert label.endswith(("kb", "Mb", "bp"))


def test_render_genome_with_baf(tiny_regions, tiny_cytoband, tiny_baf):
    cov, bin_size = mm.read_mosdepth(tiny_regions)
    cov = _prepare_cov(cov, bin_size)
    cb = mm.read_cytoband_df(tiny_cytoband)
    lengths = mm.chrom_lengths_from_cb(cb)
    baf = mm.read_baf_vcf(tiny_baf, min_dp=10)
    args = _make_args(tiny_regions)
    png, _ = mm.render_karyotype_genome_png(
        cov, cb, lengths, baf_df=baf, sex="male",
        bin_size=bin_size, args=args,
    )
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    # BAF-present render is taller (FIG_H_GENOME_BAF > FIG_H_GENOME_ONLY)
    # so the PNG should be larger than the BAF-less version
    png_no_baf, _ = mm.render_karyotype_genome_png(
        cov, cb, lengths, baf_df=None, sex="male",
        bin_size=bin_size, args=args,
    )
    assert len(png) > len(png_no_baf)


def test_render_genome_closes_figure(tiny_regions, tiny_cytoband):
    """No leaked matplotlib figures after rendering."""
    plt.close("all")
    cov, bin_size = mm.read_mosdepth(tiny_regions)
    cov = _prepare_cov(cov, bin_size)
    cb = mm.read_cytoband_df(tiny_cytoband)
    lengths = mm.chrom_lengths_from_cb(cb)
    args = _make_args(tiny_regions)
    mm.render_karyotype_genome_png(
        cov, cb, lengths, baf_df=None, sex="female",
        bin_size=bin_size, args=args,
    )
    assert len(plt.get_fignums()) == 0


def test_render_genome_does_not_leak_rcparams(tiny_regions, tiny_cytoband):
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
    mm.render_karyotype_genome_png(
        cov, cb, lengths, baf_df=None, sex="male",
        bin_size=bin_size, args=args,
    )
    after = (
        plt.rcParams["font.family"],
        plt.rcParams["axes.edgecolor"],
        plt.rcParams["grid.color"],
    )
    assert before == after




def test_kary_attach_xpos_aligns_with_offsets():
    df = pd.DataFrame({"chrom": ["chr1", "chr2"], "start": [100, 200]})
    offsets = {"chr1": 0, "chr2": 1_000_000}
    out = mm._kary_attach_xpos(df, "start", offsets)
    assert list(out["xpos"]) == [100, 1_000_200]
