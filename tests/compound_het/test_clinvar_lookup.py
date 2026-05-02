"""Tests for the ClinVar lookup (TSV + raw-VCF dispatch, hg38)."""

from __future__ import annotations

import gzip
import lzma

import pytest

import molamola as mm


def test_canon_clnsig_buckets():
    assert mm.canon_clnsig("Pathogenic") == "p_or_lp"
    assert mm.canon_clnsig("Likely_pathogenic") == "p_or_lp"
    assert mm.canon_clnsig("Pathogenic/Likely_pathogenic") == "p_or_lp"
    assert mm.canon_clnsig("Uncertain_significance") == "vus"
    assert (
        mm.canon_clnsig("Conflicting_classifications_of_pathogenicity")
        == "conflicting"
    )
    assert mm.canon_clnsig("Benign") == "benign"
    assert mm.canon_clnsig("Likely_benign") == "benign"
    assert mm.canon_clnsig("Benign/Likely_benign") == "benign"
    assert mm.canon_clnsig(None) is None
    assert mm.canon_clnsig("") is None
    assert mm.canon_clnsig("not_provided") == "other"


def test_p_lp_share_color():
    """Locked-spec invariant: P and LP must collapse to one bucket."""
    assert (
        mm.canon_clnsig("Pathogenic")
        == mm.canon_clnsig("Likely_pathogenic")
    )
    assert mm.CLNSIG_COLOR[mm.canon_clnsig("Pathogenic")] == "#c0143c"


def test_load_clinvar_lookup_tsv_returns_buckets(tiny_clinvar):
    """Bundled-format TSV: lookup yields pre-computed buckets."""
    keys = {
        ("chr1", 1100, "A", "T"),
        ("chr1", 2200, "T", "G"),
        ("chr1", 5000, "A", "G"),
    }
    lookup = mm.load_clinvar_lookup(tiny_clinvar, keys_of_interest=keys)
    assert lookup[("chr1", 1100, "A", "T")] == "p_or_lp"
    assert lookup[("chr1", 2200, "T", "G")] == "benign"
    assert lookup[("chr1", 5000, "A", "G")] == "vus"


def test_load_clinvar_lookup_vcf_yields_same_buckets(tiny_clinvar_vcf):
    """Override path: NCBI raw VCF should produce identical buckets."""
    keys = {
        ("chr1", 1100, "A", "T"),
        ("chr1", 2200, "T", "G"),
        ("chr1", 5000, "A", "G"),
    }
    lookup = mm.load_clinvar_lookup(tiny_clinvar_vcf, keys_of_interest=keys)
    assert lookup[("chr1", 1100, "A", "T")] == "p_or_lp"
    assert lookup[("chr1", 2200, "T", "G")] == "benign"
    assert lookup[("chr1", 5000, "A", "G")] == "vus"


def test_dispatch_by_extension(tiny_clinvar, tiny_clinvar_vcf):
    """Same data, two formats, identical return values."""
    keys = {("chr1", 1100, "A", "T")}
    tsv_lookup = mm.load_clinvar_lookup(tiny_clinvar, keys_of_interest=keys)
    vcf_lookup = mm.load_clinvar_lookup(tiny_clinvar_vcf, keys_of_interest=keys)
    assert tsv_lookup == vcf_lookup


def test_load_clinvar_lookup_normalises_chrom(tiny_clinvar_vcf):
    """ClinVar VCFs use bare contig names; loader normalises to UCSC."""
    keys = {("chr1", 1100, "A", "T")}
    lookup = mm.load_clinvar_lookup(tiny_clinvar_vcf, keys_of_interest=keys)
    assert ("chr1", 1100, "A", "T") in lookup
    assert ("1", 1100, "A", "T") not in lookup


def test_load_clinvar_lookup_filters_by_keys(tiny_clinvar):
    """Only requested keys land in the dict."""
    keys = {("chr1", 1100, "A", "T")}
    lookup = mm.load_clinvar_lookup(tiny_clinvar, keys_of_interest=keys)
    assert len(lookup) == 1


def test_xz_tsv_load(tmp_path):
    """The bundled file is xz-compressed; loader handles .tsv.xz."""
    p = tmp_path / "clinvar.tsv.xz"
    with lzma.open(p, "wt", preset=6) as fh:
        fh.write("chrom\tpos\tref\talt\tbucket\n")
        fh.write("chr1\t1100\tA\tT\tp_or_lp\n")
        fh.write("chr1\t1500\tG\tC\tvus\n")
    lookup = mm.load_clinvar_lookup(p)
    assert lookup[("chr1", 1100, "A", "T")] == "p_or_lp"
    assert lookup[("chr1", 1500, "G", "C")] == "vus"


def test_tsv_loader_drops_other_bucket(tmp_path):
    """Records with bucket='other' are dropped (not colour-mappable)."""
    p = tmp_path / "with_other.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("chrom\tpos\tref\talt\tbucket\n")
        fh.write("chr1\t100\tA\tT\tother\n")
        fh.write("chr1\t200\tA\tT\tp_or_lp\n")
    lookup = mm.load_clinvar_lookup(p)
    assert ("chr1", 100, "A", "T") not in lookup
    assert ("chr1", 200, "A", "T") in lookup


def test_tsv_loader_validates_required_columns(tmp_path):
    """A truncated TSV header errors out with a clear message."""
    p = tmp_path / "broken.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("chrom\tpos\tref\talt\n")  # missing bucket
        fh.write("chr1\t100\tA\tT\n")
    with pytest.raises(ValueError, match="missing"):
        mm.load_clinvar_lookup(p)
