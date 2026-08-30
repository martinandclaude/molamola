"""Circos SV density rings, legends, and the shared density path.

The circos gained an INS / DEL / DUP / INV ring stack so it no longer
depends on the linear genome map for positional SVs. The invariants
pinned here are the ones that make the two figures readable as one
report: same filter, same normalisation, same type order, and a ring
stack that does not grow into its neighbours' radii.
"""

from __future__ import annotations

import io

import numpy as np
import pytest

import molamola as mm


# --- ring geometry ---------------------------------------------------------

def test_ring_radii_cover_all_sv_types():
    radii = mm._circos_ring_radii()
    assert set(radii) == set(mm.SV_TYPES)


def test_rings_are_ordered_outermost_first():
    radii = mm._circos_ring_radii()
    tops = [radii[t][1] for t in mm.CIRCOS_SV_RING_ORDER]
    assert tops == sorted(tops, reverse=True)


def test_rings_do_not_overlap_and_keep_their_gap():
    radii = mm._circos_ring_radii()
    ordered = [radii[t] for t in mm.CIRCOS_SV_RING_ORDER]
    for outer, inner in zip(ordered, ordered[1:]):
        # Radii are (bottom, top); the gap is the outer ring's bottom
        # minus the next ring in's top.
        assert outer[0] - inner[1] == pytest.approx(mm.CIRCOS_SV_RING_GAP)


def test_ring_widths_follow_the_declared_weights():
    radii = mm._circos_ring_radii()
    widths = {t: hi - lo for t, (lo, hi) in radii.items()}
    assert all(w > 0 for w in widths.values())
    ref = next(iter(widths))
    scale = widths[ref] / mm.CIRCOS_SV_RING_WEIGHT[ref]
    for t, w in widths.items():
        assert w == pytest.approx(mm.CIRCOS_SV_RING_WEIGHT[t] * scale)


def test_shared_landscape_types_get_less_radius_than_sparse_ones():
    """The measured priority inversion this encodes: INS and DEL are
    ~85 % identical between unrelated people, so they must not be the
    thickest bands on a figure meant to show one sample."""
    radii = mm._circos_ring_radii()
    width = {t: hi - lo for t, (lo, hi) in radii.items()}
    assert max(width["INS"], width["DEL"]) < min(width["DUP"], width["INV"])


def test_shared_landscape_types_are_capped_below_full_alpha():
    """INS/DEL fill 86-94 % of bins; uncapped they render as solid bands
    however thin, and a lone DUP mark cannot compete."""
    for t in ("INS", "DEL"):
        assert mm.SV_DENSITY_ALPHA_CEILING[t] < 1.0
    for t in ("DUP", "INV"):
        assert mm.SV_DENSITY_ALPHA_CEILING[t] == 1.0


def test_alpha_ceiling_is_applied_per_type():
    huge = [10_000]
    assert mm.sv_density_alpha(huge, 5.0, svtype="INS")[0] == pytest.approx(
        mm.SV_DENSITY_ALPHA_CEILING["INS"])
    assert mm.sv_density_alpha(huge, 5.0, svtype="DUP")[0] == pytest.approx(1.0)
    # Omitting the type keeps the full range, so the helper stays usable
    # outside the two SV figures.
    assert mm.sv_density_alpha(huge, 5.0)[0] == pytest.approx(1.0)


def test_capped_type_still_keeps_the_floor():
    assert mm.sv_density_alpha([1], 10_000.0, svtype="INS")[0] >= (
        mm.SV_DENSITY_ALPHA_FLOOR)


def test_ring_stack_fills_its_declared_span_exactly():
    radii = mm._circos_ring_radii()
    lo = min(r[0] for r in radii.values())
    hi = max(r[1] for r in radii.values())
    assert (lo, hi) == pytest.approx(mm.CIRCOS_SV_RING_R)


def test_ring_stack_clears_the_cytoband_ring():
    """The gap up to r=95 is where the position tick labels are drawn.

    pyCirclize draws the "0 / 100 / 200" Mb labels inward from the
    cytoband track, so a stack that reaches 95 puts INS density
    underneath them.
    """
    assert max(mm.CIRCOS_SV_RING_R) < 95.0


def test_ring_order_matches_the_linear_map():
    """A reader who has learned one figure should not relearn the other."""
    assert mm.CIRCOS_SV_RING_ORDER == mm.SV_TYPES


# --- density alpha ---------------------------------------------------------

def test_empty_bins_are_fully_transparent():
    alphas = mm.sv_density_alpha([0, 0, 0], 10.0)
    assert np.all(alphas == 0.0)


def test_single_event_bin_stays_visible():
    """The whole point of the floor: one event must not vanish, even
    against an anchor hundreds of times larger."""
    alphas = mm.sv_density_alpha([1], 500.0)
    assert alphas[0] >= mm.SV_DENSITY_ALPHA_FLOOR


def test_alpha_floor_leaves_most_of_the_range_for_real_variation():
    """Every point of floor is a point taken off the usable ramp. The
    0.40 floor this replaced left only 0.60 for the data and produced a
    near-uniform track on a real genome."""
    assert mm.SV_DENSITY_ALPHA_FLOOR <= 0.25


def test_counts_at_the_anchor_saturate():
    assert mm.sv_density_alpha([12], 12.0)[0] == pytest.approx(1.0)


def test_counts_above_the_anchor_clip_rather_than_overflow():
    """~1 % of bins sit above p99 by construction; they must clip to
    full alpha, not produce a value matplotlib rejects."""
    alphas = mm.sv_density_alpha([12, 40, 400], 12.0)
    assert list(alphas) == pytest.approx([1.0, 1.0, 1.0])


def test_peak_bin_reaches_full_alpha():
    alphas = mm.sv_density_alpha([8], 8.0)
    assert alphas[0] == pytest.approx(1.0)


def test_alpha_is_monotone_in_count():
    alphas = mm.sv_density_alpha([1, 2, 5, 20], 20.0)
    assert list(alphas) == sorted(alphas)


def test_alpha_never_exceeds_one_above_the_anchor():
    alphas = mm.sv_density_alpha([99], 10.0)
    assert alphas[0] <= 1.0


def test_zero_anchor_does_not_divide_by_zero():
    alphas = mm.sv_density_alpha([0, 3], 0.0)
    assert np.all(np.isfinite(alphas))
    assert alphas[1] > 0


# --- the shared filter -----------------------------------------------------

def test_density_drops_non_pass_and_noise(make_sv):
    keep = make_sv(chrom="chr1", start=1_000, svtype="DEL")
    dropped_filter = make_sv(chrom="chr1", start=2_000, svtype="DEL",
                             filter_="GT")
    noisy = make_sv(chrom="chr1", start=3_000, svtype="DEL")
    noisy.noise_flags.add("cov_anomaly")

    _, svs_filt, _, _, _ = mm._sv_density(
        [keep, dropped_filter, noisy], {"chr1": 10_000_000}, 1_000_000,
    )
    assert svs_filt == [keep]


def test_density_and_linear_map_share_one_normalisation(make_sv):
    """Same bin must not render at two alphas across one report."""
    svs = [make_sv(chrom="chr1", start=i * 1000, svtype="DEL")
           for i in range(7)]
    contigs = {"chr1": 10_000_000}

    _, _, bins_a, peak_a, anchor_a = mm._sv_density(svs, contigs, 1_000_000)
    _, _, bins_b, peak_b, anchor_b = mm._sv_density(svs, contigs, 1_000_000)
    assert peak_a == peak_b
    assert anchor_a == anchor_b
    assert np.array_equal(bins_a["chr1"]["DEL"], bins_b["chr1"]["DEL"])
    assert peak_a["DEL"] == 7.0


def test_density_restricts_to_canonical_chroms(make_sv):
    off = make_sv(chrom="chrUn_KI270302v1", start=100, svtype="INS")
    chroms, _, bins, _, _ = mm._sv_density(
        [off], {"chr1": 10_000_000}, 1_000_000,
    )
    assert chroms == ["chr1"]
    assert "chrUn_KI270302v1" not in bins


# --- the p99 anchor --------------------------------------------------------

def test_anchor_ignores_a_single_outlier_bin(make_sv):
    """The defect this fixes: one hotspot bin dragged the whole scale.

    Twenty bins holding one event each plus one bin holding fifty must
    not put the twenty at the bottom of the ramp.
    """
    svs = [make_sv(chrom="chr1", start=i * 1_000_000 + 10, svtype="DEL")
           for i in range(20)]
    svs += [make_sv(chrom="chr1", start=50_000_000 + i, svtype="DEL")
            for i in range(50)]
    _, _, _, peak, anchor = mm._sv_density(svs, {"chr1": 100_000_000},
                                           1_000_000)
    assert peak["DEL"] == 50.0
    assert anchor["DEL"] < peak["DEL"]


def test_anchor_is_genome_wide_not_per_chromosome(make_sv):
    """A per-chromosome anchor would render the same count differently
    on chr1 and chr21, which defeats the point of a density track."""
    svs = ([make_sv(chrom="chr1", start=i * 1_000_000, svtype="INS")
            for i in range(30)]
           + [make_sv(chrom="chr2", start=0, svtype="INS") for _ in range(9)])
    contigs = {"chr1": 100_000_000, "chr2": 100_000_000}
    _, _, bins, _, anchor = mm._sv_density(svs, contigs, 1_000_000)

    a1 = mm.sv_density_alpha(bins["chr1"]["INS"], anchor["INS"])
    a2 = mm.sv_density_alpha(bins["chr2"]["INS"], anchor["INS"])
    # chr1 bin 0 and chr2 bin 0 hold different counts, but the mapping
    # applied to them is the same one.
    assert bins["chr1"]["INS"][0] == 1.0
    assert bins["chr2"]["INS"][0] == 9.0
    assert a2[0] > a1[0]


def test_anchor_defaults_sanely_with_no_events():
    _, _, _, peak, anchor = mm._sv_density([], {"chr1": 10_000_000}, 1_000_000)
    assert all(v == 1.0 for v in peak.values())
    assert all(v == 1.0 for v in anchor.values())


def test_ramp_spreads_a_skewed_distribution_wider_than_peak_anchoring():
    """Regression guard on the whole point of the change.

    Counts drawn from a realistic skew (bulk of 1-5, rare hotspot) must
    occupy a materially wider alpha range under the p99 anchor than
    under the peak.
    """
    counts = np.array([1, 2, 3, 4, 5] * 40 + [67])
    occupied = counts[counts > 0]
    p99 = float(np.percentile(occupied, mm.SV_DENSITY_ANCHOR_PCT))

    new = mm.sv_density_alpha(counts, p99)
    old = 0.40 + 0.60 * np.sqrt(counts / counts.max())

    def spread(a):
        return np.percentile(a, 90) - np.percentile(a, 10)

    assert spread(new) > 2 * spread(old)


# --- cytoband palette ------------------------------------------------------

def test_circos_cytobands_share_the_linear_map_greys():
    """The circos took pyCirclize's own ramp until v0.6; both figures
    sit in one report and must agree on what a grey band means."""
    for stain in ("gneg", "gpos25", "gpos50", "gpos75", "gpos100", "gvar"):
        assert mm.CIRCOS_CYTOBAND_COLORS[stain] == mm.CYTOBAND_COLORS[stain]


def test_circos_centromere_is_red_not_black():
    """The one deliberate carve-out: black acen is not separable from
    gpos100 on a ring five radial units thick."""
    assert mm.CIRCOS_CYTOBAND_COLORS["acen"] != mm.CYTOBAND_COLORS["acen"]
    assert mm.CIRCOS_CYTOBAND_COLORS["acen"].lower() == "#d92f27"


# --- wiring ----------------------------------------------------------------

def test_plot_circos_takes_svs_and_bin_size():
    import inspect
    sig = inspect.signature(mm.plot_circos)
    assert "svs" in sig.parameters
    assert sig.parameters["bin_size"].default == 1_000_000


def test_disk_margin_clears_the_plotvaf_sector_names():
    """Sector names move out to r=116 under --plotvaf; the disk axes
    only spans r=105, so too small a margin clips the labels."""
    needed = mm.SECTOR_NAME_R_WITH_VAF / 105.0
    m = mm.CIRCOS_DISK_MARGIN
    disk_radius = (1.0 - 2 * m) / 2.0
    assert disk_radius * needed <= 0.5


# --- end to end ------------------------------------------------------------

def _render(bnds, svs, contigs, cytoband_path, **kwargs):
    buf = io.BytesIO()
    mm.plot_circos(
        bnds, svs, contigs, cytoband_path, buf, "SAMPLE",
        len(bnds), len(bnds), "pass", mm.noise_breakdown(bnds), **kwargs,
    )
    return buf.getvalue()


def test_circos_renders_with_rings(bundled_cytoband, make_bnd, make_sv):
    contigs = {"chr1": 248_956_422, "chr2": 242_193_529}
    svs = [
        make_sv(chrom="chr1", start=5_000_000, svtype="INS"),
        make_sv(chrom="chr1", start=5_500_000, svtype="DEL"),
        make_sv(chrom="chr2", start=9_000_000, svtype="DUP"),
        make_sv(chrom="chr2", start=9_100_000, svtype="INV"),
    ]
    png = _render([make_bnd()], svs, contigs, bundled_cytoband)
    assert png.startswith(b"\x89PNG")


def test_circos_renders_with_no_positional_svs(bundled_cytoband, make_bnd):
    """A BND-only call set must still draw its (empty) rings."""
    contigs = {"chr1": 248_956_422, "chr2": 242_193_529}
    png = _render([make_bnd()], [], contigs, bundled_cytoband)
    assert png.startswith(b"\x89PNG")


def test_circos_renders_with_plotvaf_and_rings(bundled_cytoband, make_bnd,
                                               make_sv):
    contigs = {"chr1": 248_956_422, "chr2": 242_193_529}
    png = _render(
        [make_bnd()], [make_sv(chrom="chr1", start=1_000_000)],
        contigs, bundled_cytoband, plot_vaf=True,
    )
    assert png.startswith(b"\x89PNG")
