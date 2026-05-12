#!/usr/bin/env python3
"""Derive molamola's karyotype-mode bundled refs from karyotype-cov-plot.

karyotype-mode needs two new bundled references per build:

- ``exclusion.<build>.bed.gz`` — the per-build mappability + polyTR
  exclusion mask (no transformation; we just rename and copy from the
  upstream tool).
- ``gc_10kb.<build>.bed.gz`` — per-bin GC%, downsampled from the
  upstream 500 bp table to 10 kb bins. The coarser bin saves ~30 MB
  per build while staying well within the resolution needed for a
  karyotype-overview GC correction. All-N bins are tracked through
  aggregation as ``255``.

Upstream:
    https://github.com/martinandclaude/karyotype-cov-plot
    (refs/hg38/{exclusion.bed.gz,gc_500bp.bed.gz}
     refs/t2t/{exclusion.bed.gz,gc_500bp.bed.gz})

Output schemas
--------------

``exclusion.<build>.bed.gz``: 3-column BED ``chrom\\tstart\\tend``
(merged excluded runs; gzipped).

``gc_10kb.<build>.bed.gz``: 4-column BED ``chrom\\tstart\\tend\\tgc_pct``
where ``gc_pct`` is an int 0-100, or ``255`` if every 500 bp source bin
in that 10 kb window was all-N. Bins are uniform 10 kb (final bin per
chromosome may be shorter).

Usage
-----

::

    python scripts/derive_karyotype_refs.py \\
        --upstream-refs /path/to/karyotype-cov-plot/refs \\
        --out-dir molamola/data \\
        --build hg38

    python scripts/derive_karyotype_refs.py \\
        --upstream-refs /path/to/karyotype-cov-plot/refs \\
        --out-dir molamola/data \\
        --build t2t

Re-run when the upstream tool publishes refreshed accessibility
masks, polyTR catalogs, or new reference FASTAs.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
from pathlib import Path


GC_SOURCE_BIN_BP = 500
GC_TARGET_BIN_BP = 10_000
GC_AGG_FACTOR = GC_TARGET_BIN_BP // GC_SOURCE_BIN_BP  # 20

SENTINEL_ALL_N = 255


def copy_exclusion(src: Path, dst: Path) -> None:
    """Copy ``exclusion.bed.gz`` verbatim under the molamola-style name."""
    if not src.exists():
        raise FileNotFoundError(f"missing upstream exclusion: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)
    print(f"  wrote {dst}  ({dst.stat().st_size:,} bytes)")


def aggregate_gc(src: Path, dst: Path) -> None:
    """Aggregate a 500 bp GC table down to 10 kb bins.

    Within each 10 kb window, the output GC% is the rounded mean of
    the source bins' GC% values, excluding any ``SENTINEL_ALL_N``
    source bins from the mean. If every source bin in the window is
    sentinel, the output is also sentinel.
    """
    if not src.exists():
        raise FileNotFoundError(f"missing upstream GC table: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)

    rows_in = 0
    rows_out = 0
    cur_chrom: str | None = None
    buf: list[tuple[int, int, int]] = []  # (start, end, gc)

    def flush(out_fh, chrom: str, buf: list[tuple[int, int, int]]) -> int:
        """Aggregate ``buf`` into 10 kb output rows. Returns rows emitted."""
        emitted = 0
        i = 0
        n = len(buf)
        while i < n:
            window = buf[i:i + GC_AGG_FACTOR]
            start = window[0][0]
            end = window[-1][1]
            vals = [g for (_s, _e, g) in window if g != SENTINEL_ALL_N]
            if vals:
                gc = int(round(sum(vals) / len(vals)))
                gc = max(0, min(100, gc))
            else:
                gc = SENTINEL_ALL_N
            out_fh.write(f"{chrom}\t{start}\t{end}\t{gc}\n")
            emitted += 1
            i += GC_AGG_FACTOR
        return emitted

    with gzip.open(src, "rt") as fh, gzip.open(dst, "wt") as out_fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            chrom, start, end, gc = f[0], int(f[1]), int(f[2]), int(f[3])
            rows_in += 1
            if cur_chrom is None:
                cur_chrom = chrom
            if chrom != cur_chrom:
                rows_out += flush(out_fh, cur_chrom, buf)
                cur_chrom = chrom
                buf = []
            buf.append((start, end, gc))
        if cur_chrom is not None and buf:
            rows_out += flush(out_fh, cur_chrom, buf)

    print(
        f"  wrote {dst}  ({dst.stat().st_size:,} bytes; "
        f"{rows_in:,} -> {rows_out:,} rows)",
    )


def derive_one_build(
    build: str,
    upstream_refs: Path,
    out_dir: Path,
) -> None:
    """Copy + aggregate refs for one build (hg38 or t2t)."""
    subdir = upstream_refs / build
    print(f"[{build}] reading from {subdir}")

    copy_exclusion(
        subdir / "exclusion.bed.gz",
        out_dir / f"exclusion.{build}.bed.gz",
    )
    aggregate_gc(
        subdir / "gc_500bp.bed.gz",
        out_dir / f"gc_10kb.{build}.bed.gz",
    )


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="derive_karyotype_refs",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--upstream-refs",
        type=Path,
        required=True,
        help="path to karyotype-cov-plot's refs/ directory (contains "
             "hg38/ and t2t/ subdirs)",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="output directory; usually molamola/data/",
    )
    p.add_argument(
        "--build",
        choices=["hg38", "t2t", "both"],
        default="both",
        help="reference build to derive (default: both)",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_argparser().parse_args(argv)
    builds = ["hg38", "t2t"] if args.build == "both" else [args.build]
    for build in builds:
        derive_one_build(build, args.upstream_refs, args.out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
