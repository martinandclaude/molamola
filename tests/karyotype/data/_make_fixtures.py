"""Regenerate the karyotype-mode test fixtures.

One-shot — run this script once when the test suite first lands or
when you want to refresh the synthetic inputs:

::

    python tests/karyotype/data/_make_fixtures.py

All fixtures here are entirely synthetic. Depths and BAF allele
fractions come from ``numpy.random.default_rng(seed=42)`` so the
output is bit-stable. The four test "chromosomes" (chr1, chr2, chrX,
chrY) are each treated as 50 Mb long with 1 Mb bins, regardless of
real reference lengths — this keeps the fixtures tiny and exercises
the autosome / chrX / chrY code paths without dragging in a real
reference.

Embedded simulated signals:

- chr1, bins 20-24 (Mb): het deletion (depth x 0.5)
- chr2, bins 35-39 (Mb): het duplication (depth x 1.5)
- chrY: male-sample depth (half autosomal median)
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent

CHROMS = ("chr1", "chr2", "chrX", "chrY")
BINS_PER_CHROM = 50
BIN_SIZE_BP = 1_000_000
GC_BIN_SIZE_BP = 10_000


def write_tiny_regions(path: Path, rng: np.random.Generator) -> None:
    """Synthetic mosdepth ``regions.bed.gz``: 4 chroms x 50 1 Mb bins."""
    rows = []
    for chrom in CHROMS:
        for i in range(BINS_PER_CHROM):
            base = rng.normal(loc=15.0, scale=1.5) if chrom in ("chrX", "chrY") \
                else rng.normal(loc=30.0, scale=2.0)
            if chrom == "chr1" and 20 <= i < 25:
                base *= 0.5
            if chrom == "chr2" and 35 <= i < 40:
                base *= 1.5
            base = max(0.0, base)
            start = i * BIN_SIZE_BP
            end = start + BIN_SIZE_BP
            rows.append(f"{chrom}\t{start}\t{end}\t{base:.2f}\n")
    with gzip.open(path, "wt") as fh:
        fh.writelines(rows)


def write_tiny_cytoband(path: Path) -> None:
    """Synthetic UCSC-style cytoband: 3 bands per chrom (p, acen, q)."""
    rows = []
    for chrom in CHROMS:
        # p-arm: 0-20 Mb (gpos50)
        rows.append(f"{chrom}\t0\t20000000\tp15\tgpos50\n")
        # centromere: 20-22 Mb (acen)
        rows.append(f"{chrom}\t20000000\t22000000\tp11.1\tacen\n")
        # q-arm: 22-50 Mb (gneg)
        rows.append(f"{chrom}\t22000000\t50000000\tq22\tgneg\n")
    with gzip.open(path, "wt") as fh:
        fh.writelines(rows)


def write_tiny_mask(path: Path) -> None:
    """Synthetic 3-col BED: a couple of mask intervals per autosome."""
    rows = [
        "chr1\t5000000\t7000000\n",
        "chr2\t10000000\t12000000\n",
        "chrX\t30000000\t31000000\n",
    ]
    with gzip.open(path, "wt") as fh:
        fh.writelines(rows)


def write_tiny_gc(path: Path) -> None:
    """Synthetic GC table at 10 kb bins covering the mosdepth midpoints.

    For every mosdepth bin's midpoint (500 kb, 1.5 Mb, ...), emit one
    10 kb GC row at the floored midpoint. GC% is a smooth function of
    bin index so the test can verify both the GC fit and the apply.
    """
    rows = []
    for chrom in CHROMS:
        for i in range(BINS_PER_CHROM):
            midpoint = i * BIN_SIZE_BP + BIN_SIZE_BP // 2
            start = (midpoint // GC_BIN_SIZE_BP) * GC_BIN_SIZE_BP
            end = start + GC_BIN_SIZE_BP
            gc = int(38 + (i % 25))  # 38..62, repeats every 25 bins
            rows.append(f"{chrom}\t{start}\t{end}\t{gc}\n")
    with gzip.open(path, "wt") as fh:
        fh.writelines(rows)


def main() -> int:
    rng = np.random.default_rng(seed=42)
    write_tiny_regions(HERE / "tiny_regions.bed.gz", rng)
    write_tiny_cytoband(HERE / "tiny_cytoband.txt.gz")
    write_tiny_mask(HERE / "tiny_mask.bed.gz")
    write_tiny_gc(HERE / "tiny_gc.bed.gz")
    for name in ("tiny_regions.bed.gz", "tiny_cytoband.txt.gz",
                 "tiny_mask.bed.gz", "tiny_gc.bed.gz"):
        size = (HERE / name).stat().st_size
        print(f"  wrote {name}  ({size:,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
