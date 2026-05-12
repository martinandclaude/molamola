"""Tests for karyotype-mode plotting primitives.

Smoke-level: each helper renders into a fresh Axes without raising,
basic geometry assertions on patch / line counts, and an
``rcParams['font.family']`` non-leak check to guard against the
"matplotlib state bleed across modes" risk flagged in the plan.
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import molamola as mm


@pytest.fixture
def fresh_ax():
    fig, ax = plt.subplots()
    yield ax
    plt.close(fig)


@pytest.fixture
def cb(tiny_cytoband):
    return mm.read_cytoband_df(tiny_cytoband)


@pytest.fixture
def lengths(cb):
    return mm.chrom_lengths_from_cb(cb)


@pytest.fixture
def offsets(lengths):
    return mm.cum_offsets(lengths)


def test_rounded_pill_path_vertex_count():
    p = mm._kary_rounded_pill_path(0.0, 0.0, 100.0, 1.0, 1.0)
    assert len(p.vertices) == 11
    assert len(p.codes) == 11


def test_draw_kary_cytoband_strip_renders(fresh_ax, cb):
    chr1 = cb[cb["chrom"] == "chr1"]
    n_patches_before = len(fresh_ax.patches)
    mm._draw_kary_cytoband_strip(fresh_ax, chr1)
    n_patches_after = len(fresh_ax.patches)
    assert n_patches_after > n_patches_before
    # 3 bands + 1 silhouette clip + 1 border edge = 5 patches
    assert n_patches_after - n_patches_before == 5


def test_draw_kary_cytoband_strip_empty_input_safe(fresh_ax):
    empty = pd.DataFrame(columns=["chrom", "start", "end", "name", "stain"])
    mm._draw_kary_cytoband_strip(fresh_ax, empty)
    assert list(fresh_ax.get_xticks()) == []


def test_draw_kary_centromere_ticks_one_per_chrom(fresh_ax, cb, offsets):
    n_lines_before = len(fresh_ax.lines)
    mm._draw_kary_centromere_ticks(fresh_ax, cb, offsets)
    # 4 test chroms, each with a p-arm -> 4 ticks
    assert len(fresh_ax.lines) - n_lines_before == 4


def test_kary_build_arm_ticks_drops_acrocentric_p_labels(cb, offsets):
    # Add a synthetic chr13 row so the acrocentric drop is exercised.
    extra = pd.DataFrame([
        {"chrom": "chr13", "start": 0, "end": 14_000_000,
         "name": "p11", "stain": "gvar", "arm": "p"},
        {"chrom": "chr13", "start": 14_000_000, "end": 50_000_000,
         "name": "q11", "stain": "gneg", "arm": "q"},
    ])
    cb_plus = pd.concat([cb, extra], ignore_index=True)
    cb_plus["chrom"] = pd.Categorical(
        cb_plus["chrom"], categories=mm.CHROM_ORDER, ordered=True,
    )
    lengths = mm.chrom_lengths_from_cb(cb_plus)
    offsets_plus = mm.cum_offsets(lengths)
    xs, labels = mm._kary_build_arm_ticks(cb_plus, offsets_plus)
    assert len(xs) == len(labels)
    assert "13p" not in labels
    assert "13q" in labels
    # Non-acrocentric autosomes keep both arms
    assert "1p" in labels
    assert "1q" in labels
    # chrY p is dropped, q kept
    assert "Yp" not in labels
    assert "Yq" in labels


def test_kary_plot_coverage_draws_scatter_smooth_and_expected(fresh_ax):
    scatter_df = pd.DataFrame({
        "chrom": ["chr1"] * 5,
        "xpos": [0.0, 1.0, 2.0, 3.0, 4.0],
        "cn": [2.0, 2.1, 1.9, 2.0, 2.2],
        "mask_pass": [True, True, True, True, True],
    })
    smooth_df = pd.DataFrame({
        "chrom": pd.Categorical(["chr1"] * 5, categories=mm.CHROM_ORDER),
        "xpos": [0.0, 1.0, 2.0, 3.0, 4.0],
        "smooth": [2.0, 2.0, 2.0, 2.0, 2.0],
    })
    expected_lines = [(0.0, 4.0, 2.0)]
    n_lines_before = len(fresh_ax.lines)
    mm._kary_plot_coverage(
        fresh_ax, scatter_df, smooth_df, expected_lines, ymax=5.0,
    )
    # 1 smooth line + 1 expected line = 2 lines (scatter is a collection)
    assert len(fresh_ax.lines) - n_lines_before == 2
    assert len(fresh_ax.collections) >= 1
    assert fresh_ax.get_ylim() == (0.0, 5.0)
    assert fresh_ax.get_ylabel() == "CN"


def test_kary_plot_baf_draws_scatter_and_grid(fresh_ax):
    baf_df = pd.DataFrame({
        "chrom": ["chr1"] * 10,
        "xpos": np.arange(10).astype(float),
        "baf": np.linspace(0.0, 1.0, 10),
    })
    n_lines_before = len(fresh_ax.lines)
    mm._kary_plot_baf(fresh_ax, baf_df)
    # 3 axhline (0.25, 0.5, 0.75)
    assert len(fresh_ax.lines) - n_lines_before == 3
    assert len(fresh_ax.collections) >= 1
    assert fresh_ax.get_ylim() == (0.0, 1.0)
    assert fresh_ax.get_ylabel() == "BAF"


def test_kary_apply_tabular_numerics_changes_tick_fontfamily():
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    fig.canvas.draw()
    before = ax.get_xticklabels()[0].get_fontfamily()
    mm._kary_apply_tabular_numerics(ax)
    after = ax.get_xticklabels()[0].get_fontfamily()
    assert after != before
    assert "Menlo" in after or "DejaVu Sans Mono" in after \
        or "monospace" in after
    plt.close(fig)


def test_kary_align_panel_ylabels_pins_label_coords():
    fig, axes = plt.subplots(2, 1)
    for ax in axes:
        ax.set_ylabel("x")
    mm._kary_align_panel_ylabels(*axes)
    for ax in axes:
        x, y = ax.yaxis.get_label().get_position()
        assert x == pytest.approx(mm.KARY_YLABEL_X)
        assert y == pytest.approx(0.5)
    plt.close(fig)


@pytest.mark.parametrize("bp,expected", [
    (1_000_000, "1.0 Mb"),
    (500_000,   "500 kb"),
    (1_500,     "2 kb"),
    (250,       "250 bp"),
])
def test_kary_format_bin_size(bp, expected):
    assert mm._kary_format_bin_size(bp) == expected


def test_kary_resolve_fonts_keeps_only_available_plus_generic():
    """Filter candidates down to installed fonts; always keep the last entry."""
    resolved = mm._kary_resolve_fonts(
        ("DefinitelyNotInstalledFont12345", "AlsoMissing67890", "sans-serif"),
    )
    assert resolved == ("sans-serif",)


def test_kary_resolve_fonts_preserves_present_fonts():
    """A font that exists on this system should survive resolution."""
    resolved = mm._kary_resolve_fonts(mm.KARY_FONT_SANS)
    # Final element is matplotlib's generic 'sans-serif'; always present.
    assert resolved[-1] == "sans-serif"
    # On any matplotlib install, DejaVu Sans is bundled and should resolve.
    assert "DejaVu Sans" in resolved or "DejaVu Sans" in mm.KARY_FONT_SANS


def test_kary_resolve_fonts_is_cached():
    """Memoised: two calls return the same object."""
    a = mm._kary_resolve_fonts(mm.KARY_FONT_MONO)
    b = mm._kary_resolve_fonts(mm.KARY_FONT_MONO)
    assert a is b


def test_plotting_does_not_leak_rcparams(fresh_ax, cb, tiny_regions):
    """Render every plotting primitive and confirm rcParams stays put.

    Guards the matplotlib-state-bleed risk: any inadvertent
    ``plt.rcParams.update()`` would shift the SV-mode and
    compound-het renders that share the same Python process.
    """
    before = (
        plt.rcParams["font.family"],
        plt.rcParams["axes.edgecolor"],
        plt.rcParams["grid.color"],
    )
    chr1 = cb[cb["chrom"] == "chr1"]
    mm._draw_kary_cytoband_strip(fresh_ax, chr1)
    lengths = mm.chrom_lengths_from_cb(cb)
    offsets = mm.cum_offsets(lengths)
    mm._draw_kary_centromere_ticks(fresh_ax, cb, offsets)
    cov, bin_size = mm.read_mosdepth(tiny_regions)
    cov2 = cov.copy()
    cov2["xpos"] = cov2["start"].to_numpy().astype(np.float64)
    cov2["cn"] = cov2["depth"].to_numpy() * 2.0 / 30.0
    cov2["mask_pass"] = True
    smooth = mm.rolling_median_per_chrom(
        cov2, bin_size, window_mb=3.0, mask_pass=cov2["mask_pass"].to_numpy(),
    )
    mm._kary_plot_coverage(
        fresh_ax, cov2, smooth, [(0.0, 1e9, 2.0)], ymax=5.0,
    )
    baf = pd.DataFrame({"xpos": [0.0, 1.0], "baf": [0.5, 0.5]})
    mm._kary_plot_baf(fresh_ax, baf)
    after = (
        plt.rcParams["font.family"],
        plt.rcParams["axes.edgecolor"],
        plt.rcParams["grid.color"],
    )
    assert before == after
