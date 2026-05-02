#!/usr/bin/env python3
"""Derive molamola's canonical-transcript exon TSV from a GFF3.

molamola's compound-het mode needs one canonical transcript per
protein-coding gene, with chrom + gene span + per-exon coordinates.
This script extracts that from a MANE Select GFF (hg38) or an
Ensembl GFF (T2T-CHM13v2.0).

Tag preference (auto-detected by scanning the file's first ~5k
records):

- ``tag=MANE_Select``: NCBI's gold-standard canonical-transcript
  picker. Used for hg38 — fetched from
  https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/
- ``tag=Ensembl_canonical``: Ensembl's per-build canonical pick.
  Used for T2T-CHM13v2.0 because MANE doesn't directly publish T2T
  coordinates. For protein-coding genes Ensembl_canonical mostly
  tracks MANE Select on hg38, but is build-specific and may pick a
  different transcript in T2T-resolved regions.

Output schema (gzipped TSV with header)::

    gene_symbol  chrom  start  end  strand  transcript_id
                                     exon_starts  exon_ends

Coordinates are 0-based half-open, matching the BED convention used
by ``mm.load_canonical_exons``. Multi-comma-separated ``exon_starts``
/ ``exon_ends`` are sorted ascending (regardless of strand) so the
renderer can iterate them in genomic order.

Usage::

    # hg38 (MANE Select)
    python scripts/derive_canonical_exons.py \\
        --gff /path/to/MANE.GRCh38.v1.4.ensembl_genomic.gff.gz \\
        --out molamola/data/canonical_exons.hg38.tsv.gz

    # T2T (Ensembl canonical)
    python scripts/derive_canonical_exons.py \\
        --gff /path/to/t2t.ensembl.sorted.gff3.gz \\
        --out molamola/data/canonical_exons.t2t.tsv.gz
"""

from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path


# Bare and chr-prefixed forms of canonical contigs. MANE GFF emits
# UCSC-style ``chr1`` / ``chrX``; Ensembl emits bare ``1`` / ``X``.
_CANONICAL_BARE: set[str] = {str(i) for i in range(1, 23)} | {"X", "Y"}
_CANONICAL_UCSC: set[str] = {f"chr{c}" for c in _CANONICAL_BARE}
CANONICAL_CHROMS: set[str] = _CANONICAL_BARE | _CANONICAL_UCSC


def _normalise_chrom(c: str) -> str:
    """Return the UCSC-style ``chr*`` form regardless of input convention."""
    return c if c.startswith("chr") else "chr" + c


def _get_biotype(attrs: dict[str, str]) -> str | None:
    """Extract gene/transcript biotype handling both conventions."""
    return attrs.get("biotype") or attrs.get("gene_type") or attrs.get("transcript_type")


def _get_refseq_id(attrs: dict[str, str]) -> str | None:
    """Pull the RefSeq transcript ID from Ensembl-style ``refseq_id``
    or MANE-style ``Dbxref=...,RefSeq:NM_*,...``.
    """
    rid = attrs.get("refseq_id")
    if rid:
        return rid
    dbxref = attrs.get("Dbxref", "")
    for ref in dbxref.split(","):
        if ref.startswith("RefSeq:"):
            return ref[len("RefSeq:"):]
    return None


def _parse_attrs(s: str) -> dict[str, str]:
    """GFF3 attribute column -> dict. Values are kept as raw strings."""
    out: dict[str, str] = {}
    for kv in s.rstrip(";").split(";"):
        kv = kv.strip()
        if not kv or "=" not in kv:
            continue
        k, v = kv.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def _detect_tag(gff_path: Path, sample: int = 5000) -> str:
    """Pick the canonical-transcript tag by scanning the first records.

    Returns ``"MANE_Select"`` if present (any non-zero count), else
    ``"Ensembl_canonical"``.  Raises if neither is present in the
    sample window — molamola will not silently fall back to
    "any transcript" because the canonical-pick affects HGVS lookup
    against ClinVar.
    """
    have_mane = have_ec = False
    n = 0
    with gzip.open(gff_path, "rt") if str(gff_path).endswith(".gz") \
            else open(gff_path) as fh:
        for raw in fh:
            if raw.startswith("#"):
                continue
            if "tag=MANE_Select" in raw:
                have_mane = True
            if "tag=Ensembl_canonical" in raw:
                have_ec = True
            n += 1
            if n >= sample and (have_mane or have_ec):
                break
    if have_mane:
        return "MANE_Select"
    if have_ec:
        return "Ensembl_canonical"
    raise ValueError(
        f"GFF at {gff_path} has neither tag=MANE_Select nor "
        f"tag=Ensembl_canonical in the first {n} records — refusing "
        f"to guess which transcripts are canonical"
    )


def _liftover_genes(
    genes: dict[str, dict],
    chain_path: Path,
) -> tuple[dict[str, dict], dict[str, int]]:
    """Project per-gene exon coords across builds using a UCSC chain.

    Used to take MANE Select transcripts derived from the hg38 GFF
    and produce a T2T-coordinate TSV with the *same canonical
    transcripts* (rather than Ensembl's per-build canonical pick,
    which can diverge from MANE Select in T2T-resolved regions).

    Returns ``(lifted_genes, drop_counters)``. A gene is dropped
    when at least one of its exon endpoints fails to lift, when
    exons land on multiple T2T contigs, or when the resulting span
    is degenerate. The caller uses ``drop_counters`` for the
    run summary.
    """
    from pyliftover import LiftOver  # script-time dep; molamola
                                       # itself never imports this
    lo = LiftOver(str(chain_path))
    lifted: dict[str, dict] = {}
    drops = {"no_lift": 0, "split_contig": 0, "degenerate": 0}
    for gid, g in genes.items():
        new_exons: list[tuple[int, int]] = []
        new_chrom: str | None = None
        ok = True
        for s, e in g["exons"]:
            s_lift = lo.convert_coordinate(g["chrom"], s)
            e_lift = lo.convert_coordinate(g["chrom"], e)
            if not s_lift or not e_lift:
                ok = False
                drops["no_lift"] += 1
                break
            s_chrom, s_new, *_ = s_lift[0]
            e_chrom, e_new, *_ = e_lift[0]
            if s_chrom != e_chrom:
                ok = False
                drops["split_contig"] += 1
                break
            if new_chrom is None:
                new_chrom = s_chrom
            elif new_chrom != s_chrom:
                ok = False
                drops["split_contig"] += 1
                break
            new_exons.append((min(s_new, e_new), max(s_new, e_new)))
        if not ok:
            continue
        if not new_exons or new_chrom is None:
            drops["degenerate"] += 1
            continue
        lifted[gid] = {
            "symbol": g["symbol"],
            "chrom": new_chrom,
            "start": min(s for s, _ in new_exons),
            "end":   max(e for _, e in new_exons),
            "strand": g["strand"],
            "transcript_id": g["transcript_id"],
            "exons": sorted(new_exons),
        }
    return lifted, drops


def derive(
    gff_path: Path,
    out_path: Path,
    *,
    biotype_filter: str | None = "protein_coding",
    canonical_tag: str | None = None,
    lift_chain: Path | None = None,
) -> dict[str, int]:
    """Stream the GFF, write the canonical-exon TSV, return counters."""
    if canonical_tag is None:
        canonical_tag = _detect_tag(gff_path)

    # First pass: collect gene records (so we know symbol/strand/span
    # before we encounter the transcript and its child exons).
    genes: dict[str, dict] = {}
    n_genes = 0
    n_tx = 0
    n_exons_total = 0

    with gzip.open(gff_path, "rt") if str(gff_path).endswith(".gz") \
            else open(gff_path) as fh:
        for raw in fh:
            if raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").rstrip("\r").split("\t")
            if len(cols) < 9:
                continue
            seqid, _src, ftype, start, end, _score, strand, _phase, attrs = cols[:9]
            if ftype != "gene":
                continue
            if seqid not in CANONICAL_CHROMS:
                continue
            attr_map = _parse_attrs(attrs)
            biotype = _get_biotype(attr_map)
            if biotype_filter is not None and biotype != biotype_filter:
                continue
            # Ensembl GFF prefixes IDs with ``gene:`` / ``transcript:``;
            # MANE GFF uses bare IDs (with version suffix). Normalise.
            gid = attr_map.get("ID", "")
            if gid.startswith("gene:"):
                gid = gid[len("gene:"):]
            symbol = (attr_map.get("Name")
                      or attr_map.get("gene_name") or gid)
            chrom = _normalise_chrom(seqid)
            try:
                gstart = int(start) - 1   # GFF is 1-based inclusive
                gend = int(end)
            except ValueError:
                continue
            genes[gid] = {
                "symbol": symbol, "chrom": chrom, "start": gstart,
                "end": gend, "strand": strand,
                "transcript_id": None, "exons": [],
            }
            n_genes += 1

    # Second pass: pick the one canonical transcript per gene that
    # carries the chosen tag. MANE GFF is 1 mRNA per gene already;
    # Ensembl GFF has many transcripts per gene but only one with
    # tag=Ensembl_canonical.
    tx_to_gene: dict[str, str] = {}
    canonical_tag_marker = f"tag={canonical_tag}"
    with gzip.open(gff_path, "rt") if str(gff_path).endswith(".gz") \
            else open(gff_path) as fh:
        for raw in fh:
            if raw.startswith("#"):
                continue
            if canonical_tag_marker not in raw:
                continue
            cols = raw.rstrip("\n").rstrip("\r").split("\t")
            if len(cols) < 9:
                continue
            ftype = cols[2]
            if ftype not in ("mRNA", "transcript"):
                continue
            attrs = _parse_attrs(cols[8])
            tx_id = attrs.get("ID", "")
            if tx_id.startswith("transcript:"):
                tx_id = tx_id[len("transcript:"):]
            parent = attrs.get("Parent", "")
            if parent.startswith("gene:"):
                parent = parent[len("gene:"):]
            biotype = _get_biotype(attrs)
            if biotype_filter is not None and biotype != biotype_filter:
                continue
            if parent in genes and genes[parent]["transcript_id"] is None:
                # Prefer RefSeq ID when present (matches ClinVar's
                # CLNHGVS convention); else Ensembl ID.
                refseq_id = _get_refseq_id(attrs)
                genes[parent]["transcript_id"] = refseq_id or tx_id
                tx_to_gene[tx_id] = parent
                n_tx += 1

    # Third pass: gather exons of the picked transcripts.
    with gzip.open(gff_path, "rt") if str(gff_path).endswith(".gz") \
            else open(gff_path) as fh:
        for raw in fh:
            if raw.startswith("#"):
                continue
            if "exon" not in raw:
                continue
            cols = raw.rstrip("\n").rstrip("\r").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            attrs = _parse_attrs(cols[8])
            parent = attrs.get("Parent", "")
            if parent.startswith("transcript:"):
                parent = parent[len("transcript:"):]
            if parent not in tx_to_gene:
                continue
            try:
                estart = int(cols[3]) - 1
                eend = int(cols[4])
            except ValueError:
                continue
            genes[tx_to_gene[parent]]["exons"].append((estart, eend))
            n_exons_total += 1

    # Filter to genes with both a transcript and at least one exon.
    derived = {
        gid: g for gid, g in genes.items()
        if g["transcript_id"] is not None and g["exons"]
    }

    drops: dict[str, int] = {}
    if lift_chain is not None:
        derived, drops = _liftover_genes(derived, lift_chain)

    # Write TSV.
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    with gzip.open(out_path, "wt", compresslevel=9) as fh:
        fh.write(
            "gene_symbol\tchrom\tstart\tend\tstrand\ttranscript_id\t"
            "exon_starts\texon_ends\n"
        )
        for gid, g in sorted(derived.items(), key=lambda kv: kv[1]["symbol"]):
            exons = sorted(g["exons"])
            starts = ",".join(str(s) for s, _ in exons)
            ends   = ",".join(str(e) for _, e in exons)
            fh.write(
                f"{g['symbol']}\t{g['chrom']}\t{g['start']}\t{g['end']}\t"
                f"{g['strand']}\t{g['transcript_id']}\t{starts}\t{ends}\n"
            )
            n_written += 1

    out: dict[str, int] = {
        "canonical_tag": canonical_tag,
        "genes_seen": n_genes,
        "transcripts_picked": n_tx,
        "exons_collected": n_exons_total,
        "rows_written": n_written,
    }
    out.update({f"lift_drop_{k}": v for k, v in drops.items()})
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--gff", required=True, type=Path,
                   help="input GFF3 (gzipped OK). MANE Select for hg38, "
                        "Ensembl GFF for T2T.")
    p.add_argument("--out", required=True, type=Path,
                   help="output TSV path (gzipped, e.g. "
                        "molamola/data/canonical_exons.hg38.tsv.gz)")
    p.add_argument("--biotype", default="protein_coding",
                   help="biotype filter on gene records (default "
                        "'protein_coding'); pass empty string to disable.")
    p.add_argument("--canonical-tag", default=None,
                   help="canonical-transcript tag to filter on. Default "
                        "auto-detect (MANE_Select if present, else "
                        "Ensembl_canonical).")
    p.add_argument("--lift-to", type=Path, default=None,
                   help="Optional UCSC chain file. When given, the "
                        "derived gene+exon coordinates are projected "
                        "across builds via pyliftover before writing. "
                        "Used to produce a T2T canonical-exons TSV by "
                        "lifting MANE Select hg38 coordinates onto "
                        "CHM13v2.0 (chain: hg38ToHs1.over.chain.gz).")
    args = p.parse_args(argv)

    if not args.gff.exists():
        print(f"ERROR: input GFF not found: {args.gff}", file=sys.stderr)
        return 1
    if args.lift_to is not None and not args.lift_to.exists():
        print(f"ERROR: chain file not found: {args.lift_to}", file=sys.stderr)
        return 1
    biotype = args.biotype if args.biotype else None
    t0 = time.time()
    stats = derive(
        args.gff, args.out,
        biotype_filter=biotype,
        canonical_tag=args.canonical_tag,
        lift_chain=args.lift_to,
    )
    print(f"Source GFF:         {args.gff}")
    if args.lift_to is not None:
        print(f"Lift chain:         {args.lift_to}")
    print(f"Output TSV:         {args.out}  "
          f"({args.out.stat().st_size / 1e6:.2f} MB)")
    print(f"Canonical tag used: {stats['canonical_tag']}")
    print(f"Genes seen:         {stats['genes_seen']:,}")
    print(f"Transcripts picked: {stats['transcripts_picked']:,}")
    print(f"Exons collected:    {stats['exons_collected']:,}")
    print(f"Rows written:       {stats['rows_written']:,}")
    if args.lift_to is not None:
        for k in ("no_lift", "split_contig", "degenerate"):
            v = stats.get(f"lift_drop_{k}", 0)
            print(f"Lift drop ({k:<13}): {v:,}")
    print(f"Time:               {time.time() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
