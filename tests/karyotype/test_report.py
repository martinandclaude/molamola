"""End-to-end tests for the karyotype HTML report."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import molamola as mm


def test_karyotype_report_full_run_renders_html(tmp_path, tiny_regions):
    rc = mm.main([
        "--mosdepth", str(tiny_regions),
        "--reference", "hg38",
        "--no-mask", "--no-gc",
        "--out", str(tmp_path),
    ])
    assert rc == 0
    out_html = tmp_path / "tiny_regions.karyotype.report.html"
    assert out_html.exists()
    body = out_html.read_text()
    assert "<title>molamola karyotype" in body
    assert "inferred genomic sex" in body
    assert "scatter bin" in body
    assert "smooth: 0.5 Mb" in body
    assert "mask: off" in body
    assert "GC: off" in body


def test_karyotype_report_with_baf_includes_baf_chip(
    tmp_path, tiny_regions, tiny_baf,
):
    rc = mm.main([
        "--mosdepth", str(tiny_regions),
        "--vcf", str(tiny_baf),
        "--reference", "hg38",
        "--no-mask", "--no-gc",
        "--out", str(tmp_path),
    ])
    assert rc == 0
    body = (tmp_path / "tiny_regions.karyotype.report.html").read_text()
    assert tiny_baf.name in body
    assert "het sites" in body


def test_karyotype_report_png_flag_writes_standalone_png(
    tmp_path, tiny_regions,
):
    rc = mm.main([
        "--mosdepth", str(tiny_regions),
        "--reference", "hg38",
        "--no-mask", "--no-gc",
        "--png",
        "--out", str(tmp_path),
    ])
    assert rc == 0
    assert (tmp_path / "tiny_regions.karyotype.report.html").exists()
    genome = tmp_path / "tiny_regions.karyotype.genome.png"
    assert genome.exists()
    assert genome.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    # The previous v0.3.0 also produced a *.karyotype.per_chrom.png;
    # confirm we no longer emit it.
    assert not (tmp_path / "tiny_regions.karyotype.per_chrom.png").exists()


def test_karyotype_report_sample_override(tmp_path, tiny_regions):
    rc = mm.main([
        "--mosdepth", str(tiny_regions),
        "--reference", "hg38",
        "--no-mask", "--no-gc",
        "--sample", "MYSAMPLE",
        "--out", str(tmp_path),
    ])
    assert rc == 0
    out_html = tmp_path / "MYSAMPLE.karyotype.report.html"
    assert out_html.exists()


def test_karyotype_report_sex_override(tmp_path, tiny_regions):
    rc = mm.main([
        "--mosdepth", str(tiny_regions),
        "--reference", "hg38",
        "--no-mask", "--no-gc",
        "--sex", "female",
        "--out", str(tmp_path),
    ])
    assert rc == 0
    body = (tmp_path / "tiny_regions.karyotype.report.html").read_text()
    assert "inferred genomic sex: female" in body


def test_karyotype_report_reference_mismatch_refused(tmp_path, tmp_path_factory):
    """Filename hints t2t but --reference hg38 -> exit 2 without --force."""
    import gzip
    bed = tmp_path / "sample.t2t.regions.bed.gz"
    with gzip.open(bed, "wt") as fh:
        fh.write("chr1\t0\t1000000\t30\n")
    rc = mm.main([
        "--mosdepth", str(bed),
        "--reference", "hg38",
        "--no-mask", "--no-gc",
        "--out", str(tmp_path),
    ])
    assert rc == 2


def test_karyotype_report_html_no_section_headers(tmp_path, tiny_regions):
    """No 'Genome-wide coverage' / 'Per-chromosome coverage' headers anywhere."""
    rc = mm.main([
        "--mosdepth", str(tiny_regions),
        "--reference", "hg38",
        "--no-mask", "--no-gc",
        "--out", str(tmp_path),
    ])
    assert rc == 0
    body = (tmp_path / "tiny_regions.karyotype.report.html").read_text()
    assert "Genome-wide coverage" not in body
    assert "Per-chromosome coverage" not in body


def test_karyotype_refuses_sv_vcf_as_baf_source(tmp_path, tiny_regions, capsys):
    """``--vcf`` pointing at an SV VCF is rejected with a clear error."""
    sv_vcf = tmp_path / "sample.sv.vcf"
    sv_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t1000\t.\tA\t<DEL>\t30\tPASS\tSVTYPE=DEL\tGT\t0/1\n"
    )
    rc = mm.main([
        "--mosdepth", str(tiny_regions),
        "--vcf", str(sv_vcf),
        "--reference", "hg38",
        "--no-mask", "--no-gc",
        "--out", str(tmp_path),
    ])
    assert rc == 1
    err = capsys.readouterr().err
    assert "SVTYPE" in err
    assert "small-variant" in err
