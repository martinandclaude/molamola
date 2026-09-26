"""Circos drawing driven by rearrangement tiers.

Each classified event is one arc styled by its tier; BNDs outside every
event are single junctions; candidates are numbered with badges that do
not overprint. The colour guarantees for the tiers live in
``test_vaf_colors.py``.
"""

from __future__ import annotations

import io

import pytest

import molamola as mm


P8, P21 = 92_050_000, 34_900_000
CONTIGS = {"chr8": 145_138_636, "chr16": 90_338_345, "chr21": 46_709_983}


@pytest.fixture(scope="module")
def cyto(bundled_cytoband):
    return mm.load_cytobands(bundled_cytoband)


def _t821(make_bnd, **kw):
    return [
        make_bnd(chr1="chr8", pos1=P8, chr2="chr21", pos2=P21,
                 orientation="++", sv_id="der8", **kw),
        make_bnd(chr1="chr8", pos1=P8 + 10, chr2="chr21", pos2=P21 - 10,
                 orientation="--", sv_id="der21", **kw),
    ]


def _classify(bnds, svs, cyto):
    return mm.classify_rearrangements(mm.collect_junctions(bnds, svs),
                                      cyto, mask=None)


# --- arcs ------------------------------------------------------------------

def test_a_paired_event_is_one_arc(make_bnd, cyto):
    bnds = _t821(make_bnd)
    arcs = mm._circos_arcs(_classify(bnds, [], cyto), [], CONTIGS, cyto)
    assert len(arcs) == 1
    a = arcs[0]
    assert a.tier == "candidate"
    assert a.name == "t(8;21)(q21.3;q22.12)"
    assert a.region1[0] == "chr8" and a.region1[1] < P8 < a.region1[2]
    assert set(a.ident.split("+")) == {"der8", "der21"}


def test_an_insertion_ribbon_spans_the_donor_segment(make_bnd, cyto):
    bnds = [
        make_bnd(chr1="chr8", pos1=P8, chr2="chr21", pos2=P21,
                 orientation="++"),
        make_bnd(chr1="chr8", pos1=P8 + 20, chr2="chr21",
                 pos2=P21 + 3_000_000, orientation="--"),
    ]
    arcs = mm._circos_arcs(_classify(bnds, [], cyto), [], CONTIGS, cyto)
    donor = arcs[0].region2
    assert (donor[1], donor[2]) == (P21, P21 + 3_000_000)


def test_bnds_outside_every_event_are_single_arcs(make_bnd, cyto):
    loose = make_bnd(chr1="chr16", pos1=10_000_000, chr2="chr16",
                     pos2=10_200_000, orientation="+-", sv_id="short")
    arcs = mm._circos_arcs([], [loose], CONTIGS, cyto)
    assert [(a.tier, a.ident) for a in arcs] == [("single", "short")]


def test_arcs_off_the_plotted_chromosomes_are_skipped(make_bnd, cyto):
    arcs = mm._circos_arcs(_classify(_t821(make_bnd), [], cyto), [],
                           {"chr8": CONTIGS["chr8"]}, cyto)
    assert arcs == []


# --- styling ---------------------------------------------------------------

def _arc(**kw):
    base = dict(region1=("chr8", 0, 1), region2=("chr21", 0, 1), vaf=0.5,
                support=10, tier="single", noisy=False, non_pass=False,
                ident="x")
    base.update(kw)
    return mm._Arc(**base)


def test_tiers_rank_by_emphasis():
    alpha = {t: mm._arc_style(_arc(tier=t))[1] for t in mm.ARC_TIER_STYLE}
    width = {t: mm._arc_style(_arc(tier=t))[2] for t in mm.ARC_TIER_STYLE}
    assert alpha["candidate"] >= alpha["repeat"] > alpha["single"]
    assert width["candidate"] > width["repeat"] > width["single"]


def test_repeat_tier_is_as_readable_as_a_candidate():
    """A real fusion demoted by the mask rule is drawn in this tier, and
    the VAF palette's 3:1 contrast holds only at the candidates' alpha -
    so the two share it and differ by width."""
    assert mm.ARC_TIER_STYLE["repeat"][0] == mm.ARC_TIER_STYLE["candidate"][0]


def test_noise_flag_greys_an_arc_whatever_its_tier():
    for tier in mm.ARC_TIER_STYLE:
        color, alpha, _ = mm._arc_style(_arc(tier=tier, noisy=True))
        assert (color, alpha) == (mm.NOISE_COLOR, mm.NOISE_ARC_ALPHA)


def test_non_pass_is_fainter_than_pass():
    for tier in mm.ARC_TIER_STYLE:
        assert (mm._arc_style(_arc(tier=tier, non_pass=True))[1]
                < mm._arc_style(_arc(tier=tier))[1])


def test_candidates_draw_last():
    arcs = [_arc(tier="candidate"), _arc(tier="single"),
            _arc(tier="repeat"), _arc(tier="candidate", noisy=True)]
    order = [(a.tier, a.noisy) for a in sorted(arcs, key=mm._arc_draw_key)]
    assert order == [("candidate", True), ("single", False),
                     ("repeat", False), ("candidate", False)]


# --- badges ----------------------------------------------------------------

def test_close_badges_step_inward():
    placed = mm._place_badges(
        [("chr21", P21, 1), ("chr21", P21 + 100_000, 2)], CONTIGS)
    radii = {n: r for _, _, r, n in placed}
    assert radii[1] != radii[2]


def test_both_ends_of_a_short_event_share_one_badge():
    placed = mm._place_badges(
        [("chr16", 20_000_000, 1), ("chr16", 20_500_000, 1)], CONTIGS)
    assert len(placed) == 1


def test_far_apart_badges_stay_on_the_outer_radius():
    placed = mm._place_badges(
        [("chr8", P8, 1), ("chr21", P21, 1)], CONTIGS)
    assert {r for _, _, r, _ in placed} == {mm.CIRCOS_BADGE_RADII[0]}


# --- end to end ------------------------------------------------------------

def test_plot_circos_draws_events(make_bnd, make_sv, bundled_cytoband, cyto):
    bnds = _t821(make_bnd) + [make_bnd(chr1="chr8", pos1=10_000_000,
                                       chr2="chr16", pos2=5_000_000,
                                       sv_id="lone")]
    inv = make_sv(chrom="chr16", start=15_800_000, end=67_070_000,
                  svtype="INV", svlen=51_270_000)
    buf = io.BytesIO()
    mm.plot_circos(bnds, [inv], CONTIGS, bundled_cytoband, buf, "S",
                   events=_classify(bnds, [inv], cyto))
    assert buf.getvalue().startswith(b"\x89PNG")


def _vcf(tmp_path):
    header = (
        "##fileformat=VCFv4.2\n##source=Sniffles2_2.8.0\n"
        + "".join(f"##contig=<ID={c},length={L}>\n"
                  for c, L in {"chr1": 248_956_422, **CONTIGS}.items())
        + "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"t\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )
    info = "SUPPORT=20;COVERAGE=30,30,30,30,30;VAF=0.400"
    recs = [
        f"chr8\t{P8}\tB1\tG\tG[chr21:{P21}[\t60\tPASS\tSVTYPE=BND;{info}\tGT\t0/1",
        f"chr8\t{P8 + 10}\tB2\tG\t]chr21:{P21 - 10}]G\t60\tPASS\tSVTYPE=BND;{info}\tGT\t0/1",
        f"chr1\t5000000\tD1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=-500;END=5000500;{info}\tGT\t0/1",
    ]
    p = tmp_path / "t821.vcf"
    p.write_text(header + "\n".join(recs) + "\n")
    return p


def test_only_sv_chroms_draws_just_the_rearranged_chromosomes(tmp_path,
                                                              capsys):
    rc = mm.main(["--vcf", str(_vcf(tmp_path)), "--out", str(tmp_path / "o"),
                  "--only-sv-chroms"])
    assert rc == 0
    assert "--only-sv-chroms: drawing 2 of 4 chromosomes" in capsys.readouterr().out


def test_report_states_the_rearrangement_summary(tmp_path):
    mm.main(["--vcf", str(_vcf(tmp_path)), "--out", str(tmp_path / "o")])
    html = (tmp_path / "o" / "t821.report.html").read_text()
    assert ("<strong>Rearrangements:</strong> 1 candidate (1 translocation)"
            in html)
