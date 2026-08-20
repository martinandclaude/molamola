"""Tests for the discrete VAF colour classes.

VAF is binned rather than ramped: in a typical ONT call set ~90 % of
BND arcs sit below VAF 0.6, so a linear [0, 1] ramp crowds nearly all
of them into one narrow stretch of the colormap. These tests pin the
class boundaries, the class colours, and the invariant that the arcs
and the colorbar they are read against share one norm.
"""

from __future__ import annotations

import pytest
from matplotlib.colors import to_rgb

import molamola as mm


def _relative_luminance(rgb) -> float:
    def chan(c):
        return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (chan(c) for c in rgb)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _contrast_on_paper(color) -> float:
    """WCAG contrast ratio of a colour against the paper background."""
    lo, hi = sorted((_relative_luminance(to_rgb(color)),
                     _relative_luminance(to_rgb(mm.PAPER_BG))))
    return (hi + 0.05) / (lo + 0.05)


def test_edges_and_colors_are_consistent():
    assert len(mm.VAF_CLASS_EDGES) == len(mm.VAF_CLASS_COLORS) + 1
    assert len(mm.VAF_CLASS_LABELS) == len(mm.VAF_CLASS_COLORS)
    assert mm.VAF_CLASS_EDGES[0] == 0.0
    assert mm.VAF_CLASS_EDGES[-1] == 1.0
    assert list(mm.VAF_CLASS_EDGES) == sorted(mm.VAF_CLASS_EDGES)


@pytest.mark.parametrize("color", mm.VAF_CLASS_COLORS)
def test_every_class_is_legible_on_the_paper_background(color):
    """Thin arcs need real contrast; plasma's yellow end gave ~1.6:1,
    which is what made the highest-VAF arcs hardest to see."""
    assert _contrast_on_paper(color) >= 3.0


@pytest.mark.parametrize(
    "vaf, expected_index",
    [
        (0.00, 0), (0.10, 0), (0.32, 0), (0.329, 0),
        (0.33, 1), (0.50, 1), (0.659, 1),
        (0.66, 2), (0.80, 2), (1.00, 2),
    ],
)
def test_vaf_maps_to_expected_class(vaf, expected_index):
    expected = to_rgb(mm.VAF_CLASS_COLORS[expected_index])
    assert to_rgb(mm.vaf_to_color(vaf)[:3]) == pytest.approx(expected, abs=0.01)


def test_boundaries_are_left_inclusive():
    """A value exactly on a boundary belongs to the class above it."""
    for edge in mm.VAF_CLASS_EDGES[1:-1]:
        below = mm.vaf_to_color(edge - 1e-6)
        at = mm.vaf_to_color(edge)
        assert below != at


def test_out_of_range_vaf_is_clamped_not_crashed():
    assert mm.vaf_to_color(-0.5) == mm.vaf_to_color(0.0)
    assert mm.vaf_to_color(1.5) == mm.vaf_to_color(1.0)


def test_arcs_and_colorbar_share_one_norm():
    """The colorbar is built from VAF_NORM; if vaf_to_color used a
    different mapping the scale would silently mislabel every arc."""
    for vaf in (0.05, 0.2, 0.35, 0.5, 0.8):
        assert mm.VAF_CMAP(mm.VAF_NORM(vaf)) == mm.vaf_to_color(vaf)


def test_noise_color_is_distinct_from_every_class():
    """Noise-flagged BNDs must never be mistaken for a VAF class."""
    assert to_rgb(mm.NOISE_COLOR) not in [to_rgb(c) for c in mm.VAF_CLASS_COLORS]


# --- VAF classes are thirds ------------------------------------------------

def test_classes_are_even_thirds():
    """Three classes at 0-33 / 33-66 / 66-100 %. Uneven boundaries were
    harder to hold in the head than the ramp they replaced."""
    assert mm.VAF_CLASS_EDGES == (0.0, 0.33, 0.66, 1.0)
    assert len(mm.VAF_CLASS_COLORS) == 3


def test_colorbar_ticks_are_labelled_as_percentages():
    """House style is percentages, not fractions, in user-facing output."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    sm = plt.cm.ScalarMappable(cmap=mm.VAF_CMAP, norm=mm.VAF_NORM)
    cb = fig.colorbar(sm, ax=ax)
    mm._style_vaf_colorbar(cb)
    labels = [t.get_text() for t in cb.ax.get_yticklabels()]
    plt.close(fig)
    assert labels == ["0 %", "33 %", "66 %", "100 %"]


# --- --plotvaf -------------------------------------------------------------

def test_plotvaf_flag_defaults_off():
    """WGS runs must not start emitting hundreds of labels by default."""
    parser = mm.build_argparser()
    args = parser.parse_args(["--vcf", "x.vcf", "--out", "o"])
    assert args.plotvaf is False


def test_plotvaf_flag_can_be_enabled():
    parser = mm.build_argparser()
    args = parser.parse_args(["--vcf", "x.vcf", "--out", "o", "--plotvaf"])
    assert args.plotvaf is True


def test_plot_circos_accepts_plot_vaf_keyword():
    """Guards the plot_main -> plot_circos wiring for the flag."""
    import inspect
    sig = inspect.signature(mm.plot_circos)
    assert "plot_vaf" in sig.parameters
    assert sig.parameters["plot_vaf"].default is False


# --- --plotvaf label staggering -------------------------------------------

class _FakeSector:
    def __init__(self, sink):
        self.sink = sink

    def text(self, text, x=None, r=None, **kwargs):
        self.sink.append({"text": text, "x": x, "r": r})


class _FakeCircos:
    """Minimal stand-in: _draw_vaf_labels only needs get_sector().text()."""

    def __init__(self):
        self.drawn: list[dict] = []

    def get_sector(self, name):
        return _FakeSector(self.drawn)


def test_vaf_labels_are_percentages():
    circos = _FakeCircos()
    n_drawn, n_crowded = mm._draw_vaf_labels(
        circos, [("chr1", 1_000_000, 0.35, "#000", "bnd_a")],
        {"chr1": 248_956_422},
    )
    assert (n_drawn, n_crowded) == (1, 0)
    assert circos.drawn[0]["text"] == "35 %"


def test_well_separated_labels_share_the_innermost_ring():
    """Nothing should be pushed outward when there is room."""
    contigs = {"chr1": 248_956_422}
    span = contigs["chr1"]
    labels = [("chr1", int(span * f), 0.5, "#000", f"b{i}")
              for i, f in enumerate((0.0, 0.45, 0.9))]
    circos = _FakeCircos()
    n_drawn, n_crowded = mm._draw_vaf_labels(circos, labels, contigs)
    assert n_drawn == 3
    assert n_crowded == 0
    assert {d["r"] for d in circos.drawn} == {mm.VAF_LABEL_RADII[0]}


def test_clustered_labels_are_staggered_across_rings():
    """Co-located breakpoints must not stack on one ring."""
    contigs = {"chr1": 248_956_422}
    labels = [("chr1", 1_000_000 + i * 1000, 0.5, "#000", f"b{i}")
              for i in range(len(mm.VAF_LABEL_RADII))]
    circos = _FakeCircos()
    n_drawn, n_crowded = mm._draw_vaf_labels(circos, labels, contigs)
    assert n_drawn == len(mm.VAF_LABEL_RADII)
    assert n_crowded == 0
    assert len({d["r"] for d in circos.drawn}) == len(mm.VAF_LABEL_RADII)


def test_overcrowded_labels_are_drawn_and_counted_not_dropped():
    """More clustered breakpoints than rings: still drawn, but reported."""
    contigs = {"chr1": 248_956_422}
    n = len(mm.VAF_LABEL_RADII) + 3
    labels = [("chr1", 1_000_000 + i * 1000, 0.5, "#000", f"b{i}")
              for i in range(n)]
    circos = _FakeCircos()
    n_drawn, n_crowded = mm._draw_vaf_labels(circos, labels, contigs)
    assert n_drawn == n
    assert n_crowded == 3


def test_no_labels_is_not_an_error():
    assert mm._draw_vaf_labels(_FakeCircos(), [], {"chr1": 1000}) == (0, 0)
