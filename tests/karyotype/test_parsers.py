"""Tests for karyotype-mode I/O parsers."""

from __future__ import annotations

import gzip

import matplotlib
matplotlib.use("Agg")

import pytest

import molamola as mm


def test_read_mosdepth_returns_df_and_bin_size(tiny_regions):
    cov, bin_size = mm.read_mosdepth(tiny_regions)
    assert bin_size == 1_000_000
    assert len(cov) == 200
    assert list(cov.columns) == ["chrom", "start", "end", "depth"]


def test_read_mosdepth_canonical_chrom_categorical(tiny_regions):
    cov, _ = mm.read_mosdepth(tiny_regions)
    assert cov["chrom"].cat.ordered
    assert list(cov["chrom"].cat.categories)[:3] == ["chr1", "chr2", "chr3"]


def test_read_mosdepth_refuses_zero_canonical_rows(tmp_path):
    bed = tmp_path / "nope.bed.gz"
    with gzip.open(bed, "wt") as fh:
        fh.write("1\t0\t1000000\t30\n")
        fh.write("2\t0\t1000000\t30\n")
    with pytest.raises(ValueError, match="zero rows on canonical"):
        mm.read_mosdepth(bed)


def test_read_mosdepth_refuses_non_uniform_bins(tmp_path):
    bed = tmp_path / "uneven.bed.gz"
    rows = []
    for i in range(10):
        rows.append(f"chr1\t{i * 100}\t{(i + 1) * 100}\t30\n")
    for i in range(10):
        rows.append(f"chr2\t{i * 5_000_000}\t{(i + 1) * 5_000_000}\t30\n")
    with gzip.open(bed, "wt") as fh:
        fh.writelines(rows)
    with pytest.raises(ValueError, match="non-uniform bins"):
        mm.read_mosdepth(bed)


def test_read_cytoband_df_basic(tiny_cytoband):
    cb = mm.read_cytoband_df(tiny_cytoband)
    assert len(cb) == 12
    assert set(cb.columns) >= {"chrom", "start", "end", "name", "stain", "arm"}
    assert (cb["arm"].isin({"p", "q"})).all()


def test_chrom_lengths_from_cb(tiny_cytoband):
    cb = mm.read_cytoband_df(tiny_cytoband)
    lengths = mm.chrom_lengths_from_cb(cb)
    assert lengths == {
        "chr1": 50_000_000,
        "chr2": 50_000_000,
        "chrX": 50_000_000,
        "chrY": 50_000_000,
    }


def test_cum_offsets_orders_per_chrom_order(tiny_cytoband):
    cb = mm.read_cytoband_df(tiny_cytoband)
    lengths = mm.chrom_lengths_from_cb(cb)
    offsets = mm.cum_offsets(lengths)
    # In CHROM_ORDER chr1 < chr2 < ... < chrX < chrY. With only those
    # four chroms present (chr3..chr22 absent -> length 0), offsets
    # accumulate as 0, 50M, 50M (chr3..22 contribute 0), 100M, 150M.
    assert offsets["chr1"] == 0
    assert offsets["chr2"] == 50_000_000
    assert offsets["chrX"] == 100_000_000
    assert offsets["chrY"] == 150_000_000
