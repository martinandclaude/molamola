"""End-to-end regression tests: run main() on tiny.vcf, validate outputs."""

from __future__ import annotations

import matplotlib

# Force a headless backend so matplotlib doesn't try to open windows in CI
matplotlib.use("Agg")

import molamola as mm


# --- defaults: PASS-only run ------------------------------------------------

def test_default_run_writes_only_html(tiny_vcf, tmp_path):
    """A default run produces exactly one artefact: the self-contained HTML."""
    rc = mm.main(["--vcf", str(tiny_vcf), "--out", str(tmp_path)])
    assert rc == 0

    html = tmp_path / "tiny.report.html"
    assert html.exists() and html.stat().st_size > 50_000

    # No PNGs written to disk; figures are embedded in the HTML.
    assert list(tmp_path.glob("*.png")) == []


# --- noise filter behaviour -------------------------------------------------

def test_default_run_flags_expected_noise(tiny_vcf, tmp_path, capsys):
    """The synthetic VCF has 1 acrocentric BND and 1 cov-anomaly BND.

    The cov-anomaly assertion uses an explicit `--cov-ratio 2.5` so the
    test asserts the flagging logic rather than the auto-default
    (which on a 17-event fixture lands well above 2.5x).
    """
    mm.main(["--vcf", str(tiny_vcf), "--out", str(tmp_path),
             "--cov-ratio", "2.5"])
    out = capsys.readouterr().out
    # BND.4_ACEN -> acrocentric flag; BND.3_NOISE -> cov_anomaly
    assert "acrocentric=1" in out
    assert "cov_anomaly=1" in out
    # DEL.3_NOISE and DUP.2_NOISE flagged
    assert "cov-noise=1" in out  # DEL row
    # DUP row also has 1 noisy
    assert "DUP=2 (clean=1, cov-noise=1)" in out


def test_cov_filter_none_keeps_everything(tiny_vcf, tmp_path, capsys):
    mm.main(["--vcf", str(tiny_vcf), "--out", str(tmp_path), "--cov-filter", "none"])
    out = capsys.readouterr().out
    assert "cov_anomaly=0" in out


def test_cov_filter_drop_removes_noise_bnds(tiny_vcf, tmp_path, capsys):
    mm.main(["--vcf", str(tiny_vcf), "--out", str(tmp_path), "--cov-filter", "drop"])
    out = capsys.readouterr().out
    # The cov-anomaly BND was removed entirely, so cov_anomaly=0 in the
    # final breakdown
    assert "cov_anomaly=0" in out


def test_no_mark_acrocentric_disables_flag(tiny_vcf, tmp_path, capsys):
    mm.main(["--vcf", str(tiny_vcf), "--out", str(tmp_path), "--no-mark-acrocentric"])
    out = capsys.readouterr().out
    assert "acrocentric=0" in out


# --- focus mode -------------------------------------------------------------

def test_focus_filters_to_matching_bnds(tiny_vcf, tmp_path, capsys):
    rc = mm.main([
        "--vcf", str(tiny_vcf), "--out", str(tmp_path),
        "--focus", "chr1:73000000",
    ])
    assert rc == 0
    out = capsys.readouterr().out
    # The reciprocal pair (BND.1, BND.2) collapses to one event before
    # focus filtering, so we should see "1 -> 1" or similar
    assert "--focus filtered" in out
    # Output filename carries the focus tag
    assert (tmp_path / "tiny.focus_chr1_73000000.report.html").exists()


def test_focus_no_match_exits_nonzero(tiny_vcf, tmp_path):
    rc = mm.main([
        "--vcf", str(tiny_vcf), "--out", str(tmp_path),
        "--focus", "chr1:999999999",
    ])
    assert rc == 2


# --- --filter all -----------------------------------------------------------

def test_filter_all_includes_gt_filtered_bnd(tiny_vcf, tmp_path, capsys):
    mm.main(["--vcf", str(tiny_vcf), "--out", str(tmp_path), "--filter", "all"])
    out = capsys.readouterr().out
    # tiny.vcf has 6 BNDs total (5 PASS + 1 GT). After dedup of the
    # reciprocal pair (BND.1/BND.2), 5 unique events remain in --filter all.
    assert "After --filter=all: 6 kept" in out
    assert "Unique BND events after reciprocal dedupe: 5" in out


# --- --min-svlen hard filter ------------------------------------------------

def test_min_svlen_default_drops_subthreshold_events(tiny_vcf, tmp_path, capsys):
    """tiny.vcf contains a 30 bp INS (INS.4); default --min-svlen 50 drops it."""
    rc = mm.main(["--vcf", str(tiny_vcf), "--out", str(tmp_path)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "--min-svlen 50: dropped 1 non-BND SVs shorter than 50 bp" in out


def test_min_svlen_zero_keeps_subthreshold_events(tiny_vcf, tmp_path, capsys):
    """`--min-svlen 0` keeps every SV regardless of size."""
    rc = mm.main([
        "--vcf", str(tiny_vcf), "--out", str(tmp_path),
        "--min-svlen", "0",
    ])
    assert rc == 0
    out = capsys.readouterr().out
    assert "dropped" not in out


def test_min_svlen_zero_does_not_drop_subthreshold(tiny_vcf, tmp_path, capsys):
    """`--min-svlen 0` runs cleanly without printing any 'dropped' line."""
    rc = mm.main([
        "--vcf", str(tiny_vcf), "--out", str(tmp_path),
        "--min-svlen", "0",
    ])
    assert rc == 0
    out = capsys.readouterr().out
    assert "dropped" not in out


# --- assembly-mismatch refuse-to-plot ---------------------------------------

def test_wrong_assembly_errors_with_helpful_message(tmp_path, capsys):
    """A VCF with hg19 contig lengths must error out, not render a misleading plot."""
    vcf = tmp_path / "hg19_sample.vcf"
    # Use bare contig names ('1', '2') so name-normalization succeeds — then
    # the length check is what should catch the assembly mismatch.
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##source=Sniffles2_2.7.5\n"
        # hg19 chr1 + chr2 lengths (vs hg38: 248956422 / 242193529)
        "##contig=<ID=1,length=249250621>\n"
        "##contig=<ID=2,length=243199373>\n"
        "##contig=<ID=3,length=198022430>\n"
        "##contig=<ID=4,length=191154276>\n"
        "##FILTER=<ID=PASS,Description=\"\">\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">\n"
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">\n"
        "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"\">\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"\">\n"
        "##INFO=<ID=VAF,Number=1,Type=Float,Description=\"\">\n"
        "##INFO=<ID=COVERAGE,Number=.,Type=Float,Description=\"\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t1000000\tINS.1\tN\tGATC\t60\tPASS\t"
        "SVTYPE=INS;SVLEN=100;END=1000000;SUPPORT=10;VAF=0.5;COVERAGE=30,30,30,30,30\n"
    )
    rc = mm.main(["--vcf", str(vcf), "--out", str(tmp_path)])
    assert rc == 1
    err = capsys.readouterr().err
    assert "ERROR: contig lengths" in err
    assert "different assembly" in err
    assert "hg38" in err
