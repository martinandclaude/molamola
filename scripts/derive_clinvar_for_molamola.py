#!/usr/bin/env python3
"""Derive molamola's reduced ClinVar TSV from NCBI's raw ClinVar VCF.

molamola's compound-het mode only needs five things from each ClinVar
record: chrom, pos, ref, alt, and the canonicalised significance
bucket (``p_or_lp`` / ``vus`` / ``conflicting`` / ``benign``).
Everything else NCBI ships in the VCF INFO field — disease names,
submitter IDs, allele frequencies, SO terms, HGVS expressions, etc.
— is unused.

Stripping to those five columns and bucketing CLNSIG up front cuts a
typical weekly ClinVar release from ~190 MB compressed (full VCF) to
~13 MB (xz-compressed TSV). The bundled file at
``molamola/data/clinvar.hg38.tsv.xz`` is produced by this script;
re-run when you want to refresh against a newer NCBI release.

Source:
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    (or the dated archive at ``vcf_GRCh38/archive/`` for a pinned
    snapshot)

Usage:
    python scripts/derive_clinvar_for_molamola.py \\
        --in /path/to/clinvar.vcf.gz \\
        --out molamola/data/clinvar.hg38.tsv.xz

Output format: 5-column TSV ``chrom\\tpos\\tref\\talt\\tbucket``,
xz-compressed via Python stdlib ``lzma`` (preset 6 — good ratio,
fast enough for re-derivation, no extra dependencies). The file is
consumed by ``mm.load_clinvar_lookup`` directly.
"""

from __future__ import annotations

import argparse
import gzip
import lzma
import re
import sys
import time
from pathlib import Path


CLNSIG_RE = re.compile(rb"(?:^|;)CLNSIG=([^;]+)")
CANONICAL_CONTIGS: set[bytes] = (
    {str(i).encode() for i in range(1, 23)} | {b"X", b"Y"}
)


def bucket_of(raw: bytes) -> str:
    """Same logic as ``molamola.canon_clnsig`` but on raw bytes.

    Kept inline so this script does not import molamola — it must
    work even from a fresh clone before molamola's deps are
    installed.
    """
    s = raw.decode("ascii", errors="replace")
    if "Conflicting" in s:
        return "conflicting"
    if "Pathogenic" in s or "Likely_pathogenic" in s:
        return "p_or_lp"
    if "Uncertain" in s:
        return "vus"
    if "Benign" in s or "Likely_benign" in s:
        return "benign"
    return "other"


COLOUR_MAPPABLE: frozenset[str] = frozenset(
    {"p_or_lp", "vus", "conflicting", "benign"}
)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--in", dest="in_path", required=True, type=Path,
                   help="NCBI ClinVar VCF (vcf_GRCh38/clinvar.vcf.gz)")
    p.add_argument("--out", dest="out_path", required=True, type=Path,
                   help="output reduced TSV path "
                        "(molamola/data/clinvar.hg38.tsv.xz)")
    args = p.parse_args(argv)

    if not args.in_path.exists():
        print(f"ERROR: input file not found: {args.in_path}",
              file=sys.stderr)
        return 1
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    n_in = n_out = 0
    bucket_counts: dict[str, int] = {}
    file_date: str | None = None
    t0 = time.time()
    with gzip.open(args.in_path, "rb") as fhi, \
         lzma.open(args.out_path, "wt", preset=6,
                    encoding="ascii") as fho:
        fho.write("chrom\tpos\tref\talt\tbucket\n")
        for raw in fhi:
            if raw.startswith(b"#"):
                if raw.startswith(b"##fileDate="):
                    file_date = raw[len(b"##fileDate="):].rstrip().decode()
                continue
            n_in += 1
            cols = raw.split(b"\t", 8)
            if len(cols) < 8:
                continue
            chrom = cols[0]
            if chrom not in CANONICAL_CONTIGS:
                continue
            m_sig = CLNSIG_RE.search(cols[7])
            if m_sig is None:
                continue
            bucket = bucket_of(m_sig.group(1))
            if bucket not in COLOUR_MAPPABLE:
                continue
            chrom_norm = "chr" + chrom.decode()
            fho.write(
                f"{chrom_norm}\t{cols[1].decode()}\t"
                f"{cols[3].decode()}\t{cols[4].decode()}\t"
                f"{bucket}\n"
            )
            bucket_counts[bucket] = bucket_counts.get(bucket, 0) + 1
            n_out += 1

    in_size = args.in_path.stat().st_size
    out_size = args.out_path.stat().st_size
    print(f"ClinVar release date: {file_date or '(no ##fileDate header found)'}")
    print(f"Input records:        {n_in:>12,}  "
          f"({in_size / 1e6:6.1f} MB compressed)")
    print(f"Output records:       {n_out:>12,}  "
          f"({out_size / 1e6:6.1f} MB compressed) "
          f"= {100 * out_size / in_size:.1f}% of input")
    print()
    print("Bucket distribution:")
    for b in ("p_or_lp", "vus", "conflicting", "benign"):
        c = bucket_counts.get(b, 0)
        pct = 100 * c / n_out if n_out else 0
        print(f"  {c:>10,} ({pct:5.1f}%)  {b}")
    print()
    print(f"Wrote {args.out_path}  ({time.time() - t0:.0f}s)")
    print()
    print("Next: confirm the file_date is what you expect, then "
          "commit the new TSV.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
