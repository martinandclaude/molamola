#!/usr/bin/env python3
"""Derive molamola's bundled gene tables from NCBI RefSeq feature tables.

SV mode labels rearrangement breakpoints with the gene they fall in. The
tables come from NCBI's RefSeq annotation, which annotates GRCh38 and
T2T-CHM13v2.0 natively in the same release - so a gene has the same
symbol on both builds, and the T2T coordinates are the annotation's own,
not a liftover.

Upstream (annotation release RS_2025_08, 2025-08-06)::

    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/
        GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_feature_table.txt.gz
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/
        GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_feature_table.txt.gz

What is kept
------------

- Every ``gene`` feature of class ``protein_coding`` on the primary
  assembly's chromosomes 1-22, X and Y, one row per annotated copy.
- The immunoglobulin and T-cell receptor loci as one row each - IGH,
  IGK, IGL, TRA/TRD, TRB, TRG - spanning their V / D / J / C segment
  genes on the locus's home chromosome. Rearrangements into these loci
  are named after the locus, not after whichever V segment a breakpoint
  happens to land in. TRD sits inside TRA, so the two are one locus.
  Orphon segments (copies far from the locus) are left out: a segment
  more than ``ORPHON_GAP`` from the locus median does not widen it.

Output
------

``genes.<build>.bed.gz``: BED6, ``chrom start end name 0 strand``, 0-based
half-open, sorted. IG/TR loci carry strand ``.``.

Usage
-----

::

    python scripts/derive_gene_tables.py \\
        --feature-table ~/Documents/refs/ncbi_refseq/GCF_000001405.40_GRCh38.p14_feature_table.txt.gz \\
        --out molamola/data/genes.hg38.bed.gz

    python scripts/derive_gene_tables.py \\
        --feature-table ~/Documents/refs/ncbi_refseq/GCF_009914755.1_T2T-CHM13v2.0_feature_table.txt.gz \\
        --out molamola/data/genes.t2t.bed.gz
"""

from __future__ import annotations

import argparse
import gzip
import re
import statistics
import sys
from pathlib import Path

CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]

#: Segment-gene classes that make up the IG / TR loci.
SEGMENT_CLASSES = {
    "V_segment", "D_segment", "J_segment", "C_region",
    "V_segment_pseudogene", "D_segment_pseudogene",
    "J_segment_pseudogene", "C_region_pseudogene",
}

#: Symbol prefix -> locus name. TRD genes sit inside TRA.
LOCI = {
    "IGH": "IGH", "IGK": "IGK", "IGL": "IGL",
    "TRA": "TRA/TRD", "TRD": "TRA/TRD", "TRB": "TRB", "TRG": "TRG",
}

#: A segment gene further than this from its locus median is an orphon.
ORPHON_GAP = 5_000_000


#: Segment and constant-region gene symbols: V / D / J segments
#: (IGHV3-21, TRBD1, IGHVII-1-1), constant regions with or without a
#: number (IGKC, TRBC1, IGLC7), and the IGH isotypes (IGHM, IGHD, IGHG1,
#: IGHA2, IGHE, IGHGP). Anchored at both ends so ordinary genes sharing
#: the prefix - IGHMBP2, TRAF1, TRDN, IGLL1 - never match.
SEGMENT_SYMBOL = re.compile(
    r"^(IG[HKL]|TR[ABDG])([VDJ][0-9IVX].*|C\d*|[MDGAE]\d?|GP)$")


def locus_of(symbol: str) -> str | None:
    """IG/TR locus a segment or constant-region gene belongs to, if any:
    ``IGHV3-21`` -> IGH, ``TRBC1`` -> TRB, ``IGHMBP2`` -> None."""
    m = SEGMENT_SYMBOL.match(symbol)
    return LOCI[m.group(1)] if m else None


def derive(feature_table: Path) -> list[tuple]:
    genes: list[tuple] = []
    segments: dict[str, list[tuple]] = {}
    with gzip.open(feature_table, "rt") as fh:
        header = fh.readline().lstrip("# ").rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[col["feature"]] != "gene":
                continue
            # GRCh38's alt scaffolds and patches have their own seq_type,
            # so this also keeps the primary assembly only.
            if f[col["seq_type"]] != "chromosome":
                continue
            chrom = f[col["chromosome"]]
            if chrom not in CHROMS:
                continue
            cls, sym = f[col["class"]], f[col["symbol"]]
            start, end = int(f[col["start"]]) - 1, int(f[col["end"]])
            if cls == "protein_coding":
                genes.append((f"chr{chrom}", start, end, sym, 0,
                              f[col["strand"]]))
            elif cls in SEGMENT_CLASSES and locus_of(sym):
                segments.setdefault(locus_of(sym), []).append(
                    (f"chr{chrom}", start, end))

    for locus, segs in sorted(segments.items()):
        home = statistics.mode(c for c, _, _ in segs)
        on_home = [(s, e) for c, s, e in segs if c == home]
        mid = statistics.median((s + e) / 2 for s, e in on_home)
        kept = [(s, e) for s, e in on_home if abs((s + e) / 2 - mid) <= ORPHON_GAP]
        dropped = len(segs) - len(kept)
        lo, hi = min(s for s, _ in kept), max(e for _, e in kept)
        print(f"  {locus:8s} {home}:{lo:,}-{hi:,} ({(hi - lo) / 1e6:.2f} Mb) "
              f"from {len(kept)} segment genes"
              + (f", {dropped} orphons / off-chromosome left out" if dropped else ""))
        genes.append((home, lo, hi, locus, 0, "."))

    order = {f"chr{c}": i for i, c in enumerate(CHROMS)}
    genes.sort(key=lambda g: (order[g[0]], g[1], g[2], g[3]))
    return genes


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--feature-table", type=Path, required=True,
                   help="NCBI RefSeq *_feature_table.txt.gz for one assembly")
    p.add_argument("--out", type=Path, required=True,
                   help="output genes.<build>.bed.gz")
    args = p.parse_args(argv)

    print(f"reading {args.feature_table.name}")
    rows = derive(args.feature_table)
    n_pc = sum(1 for r in rows if r[5] != ".")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    # mtime=0 keeps the gzip byte-identical across re-runs.
    with open(args.out, "wb") as raw, \
            gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
        for r in rows:
            gz.write(("\t".join(map(str, r)) + "\n").encode())
    print(f"wrote {args.out}: {n_pc:,} protein-coding genes + "
          f"{len(rows) - n_pc} IG/TR loci, {args.out.stat().st_size / 1e3:.0f} kB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
