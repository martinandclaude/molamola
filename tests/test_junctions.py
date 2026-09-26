"""Rearrangement junction classifier: sides, pairing, tiers, labels.

The classifier groups BND junctions (and large INV records) into
reciprocal translocations, inversions, insertions and single junctions,
then assigns display tiers. The invariants pinned here are the ones the
tiers rest on: the side each orientation keeps, what counts as the same
junction versus its reciprocal partner, and that nothing is labelled
``t(A;B)`` unless both junctions were seen.
"""

from __future__ import annotations

import gzip

import numpy as np
import pytest

import molamola as mm


# Breakpoints inside RUNX1T1 (8q21.3) and RUNX1 (21q22.12), hg38 - far
# from any centromere, so the pericentromeric rule stays out of the way.
P8, P21 = 92_050_000, 34_900_000


@pytest.fixture(scope="module")
def cyto(bundled_cytoband):
    return mm.load_cytobands(bundled_cytoband)


def _t821(make_bnd, offset8=10, offset21=10, **kw):
    """The two junctions of a balanced t(8;21), anchored on chr8.

    der(8) keeps 8pter->P8 and 21:P21->qter: near side L, mate side R.
    der(21) keeps 21pter->P21 and 8:P8->qter: near side R, mate side L.
    """
    return [
        make_bnd(chr1="chr8", pos1=P8, chr2="chr21", pos2=P21,
                 orientation="++", sv_id="der8", **kw),
        make_bnd(chr1="chr8", pos1=P8 + offset8, chr2="chr21",
                 pos2=P21 - offset21, orientation="--", sv_id="der21", **kw),
    ]


# --- sides -----------------------------------------------------------------

@pytest.mark.parametrize("orientation, near, mate", [
    ("++", "L", "R"),   # t[p[
    ("+-", "L", "L"),   # t]p]
    ("-+", "R", "R"),   # [p[t
    ("--", "R", "L"),   # ]p]t
])
def test_bnd_orientation_maps_to_kept_sides(make_bnd, orientation, near, mate):
    j = mm.junction_from_bnd(make_bnd(chr1="chr1", pos1=100, chr2="chr2",
                                      pos2=200, orientation=orientation))
    assert (j.chr_a, j.side_a, j.chr_b, j.side_b) == ("chr1", near, "chr2", mate)


def test_junction_ends_are_ordered_whichever_record_anchors_it(make_bnd):
    a = mm.junction_from_bnd(make_bnd(chr1="chr1", pos1=100, chr2="chr2",
                                      pos2=200, orientation="+-"))
    b = mm.junction_from_bnd(make_bnd(chr1="chr2", pos1=200, chr2="chr1",
                                      pos2=100, orientation="+-"))
    assert (a.chr_a, a.pos_a, a.side_a, a.chr_b, a.pos_b, a.side_b) == \
           (b.chr_a, b.pos_a, b.side_a, b.chr_b, b.pos_b, b.side_b)


# --- collecting ------------------------------------------------------------

def test_same_junction_reported_twice_merges(make_bnd):
    """cuteSV writes each mate as its own record, a few bp apart."""
    r1 = make_bnd(chr1="chr8", pos1=P8, chr2="chr21", pos2=P21,
                  orientation="++", sv_id="m1", support=12)
    r2 = make_bnd(chr1="chr21", pos1=P21 + 3, chr2="chr8", pos2=P8 - 2,
                  orientation="--", sv_id="m2", support=9)
    js = mm.collect_junctions([r1, r2], [])
    assert len(js) == 1
    assert set(js[0].sv_ids) == {"m1", "m2"}
    assert js[0].support == 12


def test_short_intrachromosomal_events_are_left_out(make_bnd, make_sv):
    small_bnd = make_bnd(chr1="chr1", pos1=50_000_000, chr2="chr1",
                         pos2=50_400_000, orientation="+-")
    small_inv = make_sv(chrom="chr1", start=60_000_000, end=60_500_000,
                        svtype="INV", svlen=500_000)
    big_inv = make_sv(chrom="chr1", start=60_000_000, end=75_000_000,
                      svtype="INV", svlen=15_000_000)
    js = mm.collect_junctions([small_bnd], [small_inv, big_inv])
    assert len(js) == 2
    assert all(j.source == "INV" for j in js)


def test_non_canonical_contigs_are_ignored(make_bnd):
    b = make_bnd(chr1="chrUn_KI270302v1", pos1=100, chr2="chr1", pos2=5_000)
    assert mm.collect_junctions([b], []) == []


# --- pairing ---------------------------------------------------------------

def test_balanced_translocation_is_one_candidate(make_bnd, cyto):
    ev = mm.classify_rearrangements(
        mm.collect_junctions(_t821(make_bnd), []), cyto, mask=None)
    assert [(e.kind, e.tier) for e in ev] == [("translocation", "candidate")]


def test_same_sides_do_not_pair(make_bnd, cyto):
    """Two junctions keeping the same sides 5 kb apart are two artefacts
    (or one junction called twice loosely), not a reciprocal pair."""
    bnds = [
        make_bnd(chr1="chr8", pos1=P8, chr2="chr21", pos2=P21,
                 orientation="++", sv_id="a"),
        make_bnd(chr1="chr8", pos1=P8 + 5_000, chr2="chr21", pos2=P21 + 5_000,
                 orientation="++", sv_id="b"),
    ]
    ev = mm.classify_rearrangements(mm.collect_junctions(bnds, []), cyto, None)
    assert sorted(e.kind for e in ev) == ["single", "single"]


def test_junctions_further_apart_than_the_window_do_not_pair(make_bnd, cyto):
    far = mm.JUNCTION_PAIR_WINDOW + 1
    ev = mm.classify_rearrangements(
        mm.collect_junctions(_t821(make_bnd, offset8=far), []), cyto, None)
    assert all(e.kind == "single" for e in ev)


def test_inv_record_is_an_inversion(make_sv, cyto):
    inv = make_sv(chrom="chr16", start=15_800_000, end=67_070_000,
                  svtype="INV", svlen=51_270_000)
    ev = mm.classify_rearrangements(mm.collect_junctions([], [inv]), cyto, None)
    assert [(e.kind, e.tier) for e in ev] == [("inversion", "candidate")]


def test_intrachromosomal_bnd_pair_with_inversion_sides_is_an_inversion(
        make_bnd, cyto):
    a, b = 15_800_000, 67_070_000
    bnds = [
        make_bnd(chr1="chr16", pos1=a, chr2="chr16", pos2=b, orientation="+-"),
        make_bnd(chr1="chr16", pos1=a + 5, chr2="chr16", pos2=b + 5,
                 orientation="-+"),
    ]
    ev = mm.classify_rearrangements(mm.collect_junctions(bnds, []), cyto, None)
    assert [e.kind for e in ev] == ["inversion"]


def test_megabase_insertion_pairs_at_one_end(make_bnd, cyto):
    """ins(8;21)-like: two junctions share the chr8 site; on chr21 they
    bound a 3 Mb donor segment that both keep."""
    bnds = [
        make_bnd(chr1="chr8", pos1=P8, chr2="chr21", pos2=P21,
                 orientation="++", sv_id="left"),       # 8 L :: 21 R
        make_bnd(chr1="chr8", pos1=P8 + 20, chr2="chr21", pos2=P21 + 3_000_000,
                 orientation="--", sv_id="right"),      # 8 R :: 21 L
    ]
    ev = mm.classify_rearrangements(mm.collect_junctions(bnds, []), cyto, None)
    assert [(e.kind, e.tier) for e in ev] == [("insertion", "candidate")]
    assert mm.rearrangement_label(ev[0], cyto).startswith("ins(8;21)(q21.3;")


# --- tiers -----------------------------------------------------------------

def test_mobile_element_shape_is_demoted(make_bnd, cyto):
    """A reference L1 copy on chr21 and a 15 bp target-site duplication
    at the chr8 insertion site: looks reciprocal, is an insertion."""
    bnds = [
        make_bnd(chr1="chr8", pos1=P8 + 15, chr2="chr21", pos2=P21,
                 orientation="++", sv_id="l1_5p"),      # 8 L :: 21 R
        make_bnd(chr1="chr8", pos1=P8, chr2="chr21", pos2=P21 + 6_000,
                 orientation="--", sv_id="l1_3p"),      # 8 R :: 21 L
    ]
    ev = mm.classify_rearrangements(mm.collect_junctions(bnds, []), cyto, None)
    assert [(e.kind, e.tier) for e in ev] == [("insertion", "repeat")]
    assert "mobile-element-sized insertion" in ev[0].reasons


def test_kilobase_offsets_at_both_ends_stay_a_translocation(make_bnd, cyto):
    ev = mm.classify_rearrangements(
        mm.collect_junctions(_t821(make_bnd, offset8=3_000, offset21=2_000),
                             []), cyto, None)
    assert [(e.kind, e.tier) for e in ev] == [("translocation", "candidate")]


def test_pericentromeric_pair_is_demoted(make_bnd, cyto):
    acen_start = min(s for s, e, n, st in cyto["chr8"] if st == "acen")
    bnds = [
        make_bnd(chr1="chr8", pos1=acen_start - 500_000, chr2="chr21",
                 pos2=P21, orientation="++"),
        make_bnd(chr1="chr8", pos1=acen_start - 499_990, chr2="chr21",
                 pos2=P21 - 10, orientation="--"),
    ]
    ev = mm.classify_rearrangements(mm.collect_junctions(bnds, []), cyto, None)
    assert ev[0].tier == "repeat"
    assert ev[0].reasons == ("pericentromeric",)


def _mask(intervals):
    return {c: (np.array([s for s, e in iv]), np.array([e for s, e in iv]))
            for c, iv in intervals.items()}


def test_pair_is_demoted_only_when_every_breakpoint_is_masked(make_bnd, cyto):
    js = mm.collect_junctions(_t821(make_bnd), [])
    both = _mask({"chr8": [(P8 - 100, P8 + 100)],
                  "chr21": [(P21 - 100, P21 + 100)]})
    one = _mask({"chr8": [(P8 - 100, P8 + 100)]})
    assert mm.classify_rearrangements(js, cyto, both)[0].reasons == \
        ("all breakpoints masked",)
    assert mm.classify_rearrangements(js, cyto, one)[0].tier == "candidate"


def test_noise_flag_demotes(make_bnd, cyto):
    bnds = _t821(make_bnd)
    bnds[0].noise_flags.add("cov_anomaly")
    ev = mm.classify_rearrangements(mm.collect_junctions(bnds, []), cyto, None)
    assert ev[0].reasons == ("noise-flagged",)


def test_mask_pad_stays_zero():
    """The mask is fragmented 500 bp runs: a 1 kb pad takes it from 40 %
    of hg38 to 84 % and would demote most real fusions."""
    assert mm.JUNCTION_MASK_PAD == 0


# --- mask intervals --------------------------------------------------------

def test_mask_loader_merges_touching_runs(tmp_path):
    p = tmp_path / "m.bed.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("chr1\t100\t200\nchr1\t200\t300\nchr1\t500\t600\n1\t700\t800\n")
    m = mm.load_mask_intervals(p)
    starts, ends = m["chr1"]
    assert list(starts) == [100, 500, 700]
    assert list(ends) == [300, 600, 800]
    assert mm.in_mask(m, "chr1", 250)
    assert not mm.in_mask(m, "chr1", 400)
    assert mm.in_mask(m, "chr1", 450, pad=50)
    assert not mm.in_mask(m, "chr2", 150)


# --- labels ----------------------------------------------------------------

def test_labels(make_bnd, make_sv, cyto):
    t = mm.classify_rearrangements(
        mm.collect_junctions(_t821(make_bnd), []), cyto, None)[0]
    assert mm.rearrangement_label(t, cyto) == "t(8;21)(q21.3;q22.12)"

    inv = make_sv(chrom="chr16", start=15_800_000, end=67_070_000,
                  svtype="INV", svlen=51_270_000)
    i = mm.classify_rearrangements(mm.collect_junctions([], [inv]), cyto, None)[0]
    assert mm.rearrangement_label(i, cyto) == "inv(16)(p13.11q22.1)"


def test_single_junction_is_never_called_a_translocation(make_bnd, cyto):
    one = _t821(make_bnd)[:1]
    ev = mm.classify_rearrangements(mm.collect_junctions(one, []), cyto, None)
    label = mm.rearrangement_label(ev[0], cyto)
    assert not label.startswith("t(")
    assert label == "8q21.3::21q22.12 (single junction)"


def test_summary_counts_tiers(make_bnd, cyto):
    bnds = _t821(make_bnd) + [make_bnd(chr1="chr1", pos1=50_000_000,
                                       chr2="chr2", pos2=50_000_000)]
    ev = mm.classify_rearrangements(mm.collect_junctions(bnds, []), cyto, None)
    assert mm.rearrangement_summary(ev) == (
        "1 candidate (1 translocation), 0 paired in repeats, "
        "1 single junctions")


# --- end to end through the parser ------------------------------------------

def test_real_format_bnd_records_classify_as_a_translocation(tmp_path, cyto):
    """Sniffles2 >= 2.8 writes the real reference base in the ALT
    (``G[chr21:...[``, ``]chr21:...]G``). Through the parser, the two
    junctions of a t(8;21) and a 51 Mb INV record must come out named
    correctly - this is the path the side mapping has to survive."""
    header = (
        "##fileformat=VCFv4.2\n"
        "##source=Sniffles2_2.8.0\n"
        "##contig=<ID=chr8,length=145138636>\n"
        "##contig=<ID=chr16,length=90338345>\n"
        "##contig=<ID=chr21,length=46709983>\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"t\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )
    info = "SUPPORT=20;COVERAGE=30,30,30,30,30;VAF=0.400"
    recs = [
        f"chr8\t{P8}\tBND.1\tG\tG[chr21:{P21}[\t60\tPASS\tSVTYPE=BND;{info}\tGT\t0/1",
        f"chr8\t{P8 + 10}\tBND.2\tG\t]chr21:{P21 - 10}]G\t60\tPASS\tSVTYPE=BND;{info}\tGT\t0/1",
        f"chr16\t20000000\tINV.1\tN\t<INV>\t60\tPASS\t"
        f"SVTYPE=INV;SVLEN=51000000;END=71000000;{info}\tGT\t0/1",
    ]
    p = tmp_path / "t821.vcf"
    p.write_text(header + "\n".join(recs) + "\n")
    _contigs, bnds, svs, _cov = mm.read_vcf(p, caller="sniffles2")
    ev = mm.classify_rearrangements(mm.collect_junctions(bnds, svs), cyto, None)
    labels = sorted(mm.rearrangement_label(e, cyto) for e in ev)
    assert [e.tier for e in ev] == ["candidate", "candidate"]
    assert labels[1] == "t(8;21)(q21.3;q22.12)"
    assert labels[0].startswith("inv(16)(p")
