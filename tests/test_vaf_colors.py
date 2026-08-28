"""Tests for the discrete VAF colour classes.

VAF is binned rather than ramped: in a typical ONT call set ~90 % of
BND arcs sit below VAF 0.6, so a linear [0, 1] ramp crowds nearly all
of them into one narrow stretch of the colormap. These tests pin the
class boundaries, the class colours, and the invariant that the arcs
and the colorbar they are read against share one norm.
"""

from __future__ import annotations

import numpy as np
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


# --- colour-vision deficiency ---------------------------------------------
#
# The VAF classes were picked under simulated CVD, not by eye. Two things
# had gone wrong with an eye-picked palette and neither was visible to
# normal vision: `het` and `hom` sat at the same lightness so they merged
# under tritanopia, and `het` desaturated to near-grey under deuteranopia
# -- the same grey molamola uses for noise-flagged BNDs, so for ~6 % of
# men real het breakends read as artefacts. These tests pin both.

# Machado et al. (2009) severity-1.0 matrices, linear sRGB.
_CVD_MATRICES = {
    "protan": ((0.152286, 1.052583, -0.204868),
               (0.114503, 0.786281, 0.099216),
               (-0.003882, -0.048116, 1.051998)),
    "deutan": ((0.367322, 0.860646, -0.227968),
               (0.280085, 0.672501, 0.047413),
               (-0.011820, 0.042940, 0.968881)),
    "tritan": ((1.255528, -0.076749, -0.178779),
               (-0.078411, 0.930809, 0.147602),
               (0.004733, 0.691367, 0.303900)),
}


def _to_linear(c):
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _simulate(color, kind):
    """Return `color` as seen with the given dichromacy."""
    rgb = _to_linear(np.array(to_rgb(color)))
    if kind != "normal":
        rgb = np.array(_CVD_MATRICES[kind]) @ rgb
    return np.clip(rgb, 0.0, 1.0)


def _lab(linear_rgb):
    m = np.array([[0.4124, 0.3576, 0.1805],
                  [0.2126, 0.7152, 0.0722],
                  [0.0193, 0.1192, 0.9505]])
    xyz = m @ linear_rgb
    wp = np.array([0.95047, 1.0, 1.08883])

    def f(t):
        return t ** (1 / 3) if t > 0.008856 else 7.787 * t + 16 / 116

    fx, fy, fz = (f(v) for v in xyz / wp)
    return np.array([116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz)])


def _blend(color, alpha):
    """Colour as actually drawn: alpha-composited over the page.

    Arcs render at alpha 0.70 and noise-flagged arcs at 0.18, which
    lifts both toward the background. Judging the raw hex values
    overstates how distinct they are on the plot.
    """
    fg, bg = np.array(to_rgb(color)), np.array(to_rgb(mm.PAPER_BG))
    return tuple(alpha * fg + (1 - alpha) * bg)


ARC_ALPHA = 0.70
NOISE_ALPHA = 0.18


def _delta_e(a, b, kind):
    return float(np.linalg.norm(_lab(_simulate(a, kind)) - _lab(_simulate(b, kind))))


CVD_KINDS = ["normal", "protan", "deutan", "tritan"]


@pytest.mark.parametrize("kind", CVD_KINDS)
def test_classes_stay_distinct_under_colour_blindness(kind):
    drawn = [_blend(c, ARC_ALPHA) for c in mm.VAF_CLASS_COLORS]
    worst = min(_delta_e(drawn[i], drawn[j], kind)
                for i in range(len(drawn)) for j in range(i + 1, len(drawn)))
    assert worst >= 30, f"{kind}: closest VAF classes are only dE {worst:.1f} apart"


@pytest.mark.parametrize("kind", CVD_KINDS)
def test_no_class_looks_like_a_noise_flagged_bnd(kind):
    """A VAF class that lands on the noise grey makes real events read as
    artefacts -- worse than merely being hard to tell apart."""
    noise = _blend(mm.NOISE_COLOR, NOISE_ALPHA)
    worst = min(_delta_e(_blend(c, ARC_ALPHA), noise, kind)
                for c in mm.VAF_CLASS_COLORS)
    assert worst >= 25, f"{kind}: a VAF class is only dE {worst:.1f} from noise"


def test_every_class_is_legible_as_drawn():
    """The 3:1 floor has to hold on the composited colour, not the raw
    hex -- alpha 0.70 costs about a third of the nominal contrast."""
    for color in mm.VAF_CLASS_COLORS:
        assert _contrast_on_paper(_blend(color, ARC_ALPHA)) >= 3.0, color


def test_lightness_is_monotonic_with_vaf():
    """Ordering must survive even total loss of hue discrimination, so
    lightness carries it as well as hue."""
    lightness = [_lab(_simulate(c, "normal"))[0] for c in mm.VAF_CLASS_COLORS]
    assert lightness == sorted(lightness, reverse=True), lightness
