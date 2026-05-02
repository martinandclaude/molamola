"""Unit tests for the coverage and acrocentric noise filters."""

from __future__ import annotations

import molamola as mm


# --- annotate_noise (BNDs) --------------------------------------------------

def test_cov_anomaly_flagged_when_high_cov_and_low_vaf(make_bnd):
    b = make_bnd(coverage=[200.0]*5, vaf=0.20)
    mm.annotate_noise(
        [b], cytobands={}, median_cov=30.0,
        mark_acrocentric=True, cov_ratio_thr=2.5, cov_vaf_max=0.35,
    )
    assert "cov_anomaly" in b.noise_flags


def test_cov_anomaly_not_flagged_when_vaf_high(make_bnd):
    # High coverage but VAF above threshold -> NOT noise (could be a real
    # tandem-dup-style event with elevated coverage)
    b = make_bnd(coverage=[200.0]*5, vaf=0.55)
    mm.annotate_noise(
        [b], cytobands={}, median_cov=30.0,
        mark_acrocentric=True, cov_ratio_thr=2.5, cov_vaf_max=0.35,
    )
    assert "cov_anomaly" not in b.noise_flags


def test_cov_anomaly_not_flagged_when_cov_normal(make_bnd):
    b = make_bnd(coverage=[30.0]*5, vaf=0.20)
    mm.annotate_noise(
        [b], cytobands={}, median_cov=30.0,
        mark_acrocentric=True, cov_ratio_thr=2.5, cov_vaf_max=0.35,
    )
    assert "cov_anomaly" not in b.noise_flags


def test_acrocentric_flag_requires_both_ends_in_p_arm(make_bnd, bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    # chr21 acen ~ 12 Mb, chr22 acen ~ 14.4 Mb in hg38
    both_in = make_bnd(chr1="chr21", pos1=5_000_000, chr2="chr22", pos2=10_000_000)
    one_in  = make_bnd(chr1="chr21", pos1=5_000_000, chr2="chr22", pos2=30_000_000)
    none_in = make_bnd(chr1="chr1",  pos1=5_000_000, chr2="chr2",  pos2=10_000_000)
    mm.annotate_noise(
        [both_in, one_in, none_in], cytobands=cyto, median_cov=30.0,
        mark_acrocentric=True, cov_ratio_thr=2.5, cov_vaf_max=0.35,
    )
    assert "acrocentric" in both_in.noise_flags
    assert "acrocentric" not in one_in.noise_flags
    assert "acrocentric" not in none_in.noise_flags


def test_acrocentric_flag_disabled_when_mark_acrocentric_false(make_bnd, bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    b = make_bnd(chr1="chr21", pos1=5_000_000, chr2="chr22", pos2=10_000_000)
    mm.annotate_noise(
        [b], cytobands=cyto, median_cov=30.0,
        mark_acrocentric=False, cov_ratio_thr=2.5, cov_vaf_max=0.35,
    )
    assert "acrocentric" not in b.noise_flags


def test_in_acrocentric_p_returns_false_for_non_acrocentric_chrom():
    # chr1 is not acrocentric, so any pos returns False even with
    # an artificial p_arm_ends mapping
    assert mm.in_acrocentric_p("chr1", 1_000_000, {"chr21": 12_000_000}) is False


def test_in_acrocentric_p_returns_true_inside_pa(bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    p_arm_ends = mm.acrocentric_p_arm_ends(cyto)
    assert mm.in_acrocentric_p("chr21", 5_000_000, p_arm_ends) is True


def test_in_acrocentric_p_returns_false_past_acen(bundled_cytoband):
    cyto = mm.load_cytobands(bundled_cytoband)
    p_arm_ends = mm.acrocentric_p_arm_ends(cyto)
    # 30 Mb on chr21 is well past the centromere
    assert mm.in_acrocentric_p("chr21", 30_000_000, p_arm_ends) is False


# --- annotate_sv_noise (DEL / DUP) ------------------------------------------

def test_sv_noise_flags_del_with_high_cov_low_vaf(make_sv):
    s = make_sv(svtype="DEL", coverage=[200.0]*5, vaf=0.20)
    mm.annotate_sv_noise([s], median_cov=30.0,
                         cov_ratio_thr=2.5, cov_vaf_max=0.35)
    assert s.is_noise is True


def test_sv_noise_flags_dup_with_high_cov_low_vaf(make_sv):
    s = make_sv(svtype="DUP", coverage=[200.0]*5, vaf=0.20)
    mm.annotate_sv_noise([s], median_cov=30.0,
                         cov_ratio_thr=2.5, cov_vaf_max=0.35)
    assert s.is_noise is True


def test_sv_noise_does_not_touch_ins(make_sv):
    s = make_sv(svtype="INS", coverage=[200.0]*5, vaf=0.20)
    mm.annotate_sv_noise([s], median_cov=30.0,
                         cov_ratio_thr=2.5, cov_vaf_max=0.35)
    assert s.is_noise is False


def test_sv_noise_does_not_touch_inv(make_sv):
    s = make_sv(svtype="INV", coverage=[200.0]*5, vaf=0.20)
    mm.annotate_sv_noise([s], median_cov=30.0,
                         cov_ratio_thr=2.5, cov_vaf_max=0.35)
    assert s.is_noise is False


def test_sv_noise_real_dup_at_2x_coverage_not_flagged(make_sv):
    # Real homozygous DUP coverage is ~2x baseline - threshold is 2.5x,
    # so this should NOT be flagged
    s = make_sv(svtype="DUP", coverage=[60.0]*5, vaf=0.95)
    mm.annotate_sv_noise([s], median_cov=30.0,
                         cov_ratio_thr=2.5, cov_vaf_max=0.35)
    assert s.is_noise is False


def test_sv_noise_no_op_when_median_zero(make_sv):
    s = make_sv(svtype="DEL", coverage=[200.0]*5, vaf=0.20)
    mm.annotate_sv_noise([s], median_cov=0.0,
                         cov_ratio_thr=2.5, cov_vaf_max=0.35)
    assert s.is_noise is False


# --- breakdowns -------------------------------------------------------------

def test_noise_breakdown_counts(make_bnd):
    a = make_bnd(sv_id="clean")
    b_acen = make_bnd(sv_id="acen")
    b_acen.noise_flags = {"acrocentric"}
    b_cov = make_bnd(sv_id="cov")
    b_cov.noise_flags = {"cov_anomaly"}
    b_both = make_bnd(sv_id="both")
    b_both.noise_flags = {"acrocentric", "cov_anomaly"}

    bd = mm.noise_breakdown([a, b_acen, b_cov, b_both])
    assert bd["n"] == 4
    assert bd["clean"] == 1
    assert bd["acrocentric"] == 2
    assert bd["cov_anomaly"] == 2
    assert bd["any_noise"] == 3


def test_sv_noise_breakdown_per_type(make_sv):
    svs = [
        make_sv(svtype="DEL"),
        make_sv(svtype="DEL"),
        make_sv(svtype="INS"),
    ]
    svs[0].noise_flags.add("cov_anomaly")  # one noisy DEL

    bd = mm.sv_noise_breakdown(svs)
    assert bd["DEL"]["pass"] == 2
    assert bd["DEL"]["clean"] == 1
    assert bd["DEL"]["noise"] == 1
    assert bd["DEL"]["cov"] == 1
    assert bd["INS"]["pass"] == 1
    assert bd["INS"]["clean"] == 1
    assert bd["INS"]["noise"] == 0
    assert bd["DUP"]["pass"] == 0
    assert bd["DUP"]["clean"] == 0
    assert bd["DUP"]["noise"] == 0


# --- _parse_cov_ratio -------------------------------------------------------

import argparse  # noqa: E402

import pytest  # noqa: E402


@pytest.mark.parametrize("s, expected", [
    ("2.5",  2.5),
    ("3",    3.0),
    ("auto", "auto"),
    ("AUTO", "auto"),
])
def test_parse_cov_ratio_valid(s, expected):
    assert mm._parse_cov_ratio(s) == expected


@pytest.mark.parametrize("s", ["abc", "-1", "0", ""])
def test_parse_cov_ratio_invalid(s):
    with pytest.raises(argparse.ArgumentTypeError):
        mm._parse_cov_ratio(s)


# --- auto_cov_ratio_threshold -----------------------------------------------

def test_auto_cov_ratio_floor_when_no_outliers(make_bnd):
    """All events at baseline coverage -> threshold pinned to floor (2.0)."""
    bnds = [make_bnd(coverage=[30.0]*5) for _ in range(20)]
    thr, n = mm.auto_cov_ratio_threshold(bnds, [], median_cov=30.0)
    assert thr == 2.0
    assert n == 20


def test_auto_cov_ratio_picks_p99_when_outliers_present(make_bnd):
    """High-coverage outliers should push p99 above the floor."""
    # 90 events at 30x (ratio 1.0), 10 events at 300x (ratio 10.0).
    # p99 of this 100-element distribution lands inside the outlier tail.
    baseline = [make_bnd(coverage=[30.0]*5, sv_id=f"BND.{i}") for i in range(90)]
    outliers = [make_bnd(coverage=[300.0]*5, sv_id=f"BND.out{i}") for i in range(10)]
    thr, n = mm.auto_cov_ratio_threshold(baseline + outliers, [], median_cov=30.0)
    assert thr > 5.0
    assert n == 100


def test_auto_cov_ratio_zero_median_returns_floor(make_bnd):
    bnds = [make_bnd(coverage=[30.0]*5)]
    thr, n = mm.auto_cov_ratio_threshold(bnds, [], median_cov=0.0)
    assert thr == 2.0
    assert n == 0


def test_auto_cov_ratio_skips_non_pass_events(make_bnd):
    """Filtered-out events must not contribute to the empirical distribution."""
    flagged = make_bnd(coverage=[300.0]*5, filter_="GT", sv_id="BND.gt")
    pass_b = make_bnd(coverage=[30.0]*5, filter_="PASS", sv_id="BND.pass")
    thr, n = mm.auto_cov_ratio_threshold([flagged, pass_b], [], median_cov=30.0)
    assert n == 1  # only the PASS event was considered


def test_auto_cov_ratio_end_to_end(tmp_path, tiny_vcf, capsys):
    """`--cov-ratio auto` runs end-to-end and reports the chosen threshold."""
    rc = mm.main([
        "--vcf", str(tiny_vcf),
        "--out", str(tmp_path),
        "--cov-ratio", "auto",
    ])
    assert rc == 0
    out = capsys.readouterr().out
    assert "--cov-ratio auto:" in out
    assert "PASS events" in out
