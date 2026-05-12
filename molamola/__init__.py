#!/usr/bin/env python3
"""molamola — plot Oxford Nanopore variation as self-contained HTML reports.

A single VCF in, one self-contained HTML out. molamola inspects the
VCF header and picks one of two plot types:

- **SV / cytogenetics report** — a circos plot plus a linear cytoband
  ideogram with per-type density tracks (INS / DEL / DUP / INV) and
  BND arcs. Selected when the VCF carries ``##INFO=<ID=SVTYPE,...>``
  (Sniffles2 / cuteSV / SVIM / pbsv / NanoVar). Supports hg38 and
  T2T-CHM13v2.0 via bundled cytobands.

- **Compound-het panels** — one per-gene phased-haplotype panel per
  candidate gene: canonical-transcript exon track, H1 / H2 hap
  lines, mint phase blocks, ClinVar-coloured missense lollipops and
  synonymous-variant ticks. Selected when the VCF carries
  ``##INFO=<ID=CSQ,...>`` AND ``##FORMAT=<ID=PS,...>`` (a phased
  small-variant VCF with VEP CSQ annotation). hg38-only.

Both modes embed figures as base64 PNG data URIs; no separate image
files are written. The report lands next to the input VCF by default
(or in ``--out DIR`` if specified).

Bundled references ship in ``molamola/data/``: cytobands for both
SV-mode references, MANE Select v1.x canonical-exon coordinates for
compound-het, and a reduced ClinVar TSV (chrom + pos + ref + alt +
significance bucket; xz-compressed). No auto-download, no online
lookups — molamola is self-contained, offline-friendly.

VCFs that don't match either shape are refused with a clear error
rather than silently producing a default plot.
"""

from __future__ import annotations

__version__ = "0.2.0"

import argparse
import base64
import gzip
import html as html_mod
import io
import lzma
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.path import Path as MplPath


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CHROM_ORDER: list[str] = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CHROM_SET: set[str] = set(CHROM_ORDER)

#: Bare canonical names without the ``chr`` prefix (Ensembl / GRCh* style).
#: Used to detect VCFs that ship ``1``, ``2``, ..., ``X``, ``Y`` instead of
#: the UCSC ``chr1``/``chr2``/... convention. We auto-prepend ``chr`` so
#: those VCFs match the bundled cytoband.
_BARE_CANONICAL: set[str] = {str(i) for i in range(1, 23)} | {"X", "Y"}


def _normalize_chrom(name: str) -> str:
    """Prepend ``chr`` to canonical chromosome names that lack it.

    Leaves anything else untouched (so non-canonical contigs like
    decoys, alt scaffolds, or ``chrM`` stay as-is and get filtered
    out downstream by :data:`CHROM_SET` membership).
    """
    if name in _BARE_CANONICAL:
        return "chr" + name
    return name

#: Cytoband palette - greyscale ISCN G-banding style. Centromere
#: (acen) is solid black so it still reads as the most condensed
#: region; gpos100 is very-dark-grey so it stays distinguishable.
CYTOBAND_COLORS: dict[str, str] = {
    "gneg":    "#FFFFFF",
    "gpos25":  "#CFCFCF",
    "gpos50":  "#9C9C9C",
    "gpos75":  "#666666",
    "gpos100": "#2A2A2A",
    "acen":    "#000000",
    "gvar":    "#B5B5B5",
    "stalk":   "#7A7A7A",
}

SV_TYPES: tuple[str, ...] = ("INS", "DEL", "DUP", "INV")

#: Density-track colours. Saturated, perceptually distinct, and not on
#: the greyscale cytoband palette. Density per bin is encoded as alpha
#: (transparent = empty bin, solid = hotspot) so single-event bins are
#: still visible against the white background.
SV_TYPE_COLOR: dict[str, str] = {
    "INS": "#1f77b4",  # blue
    "DEL": "#d62728",  # red
    "DUP": "#2ca02c",  # green
    "INV": "#9467bd",  # purple
}

#: SV types to which the coverage-based noise filter is applied. A real
#: het DEL halves coverage; a real DUP sits at 1.5-2x. Anything where
#: max coverage is far above baseline AND VAF is low looks like
#: repeat-collapse rather than a true event.
SV_COV_FILTER_TYPES: tuple[str, ...] = ("DEL", "DUP")

#: VAF colorbar limits. Maps the full [0, 1] VAF range across the
#: plasma colormap. Combined with PASS filtering, acrocentric flagging,
#: and cov-anomaly detection, the surviving BNDs are mostly real
#: events whose true VAF (subgermline -> mosaic -> het -> hom) is
#: information worth surfacing directly.
VAF_VMIN: float = 0.0
VAF_VMAX: float = 1.0
VAF_CMAP = plt.get_cmap("plasma")

NOISE_COLOR: str = "#888888"

#: Acrocentric chromosomes. Their short arms (chr13/14/15/21/22 p) are
#: highly repetitive in hg38 and a common source of long-read
#: mismapping that produces phantom BND calls.
ACROCENTRIC_CHROMS: tuple[str, ...] = ("chr13", "chr14", "chr15", "chr21", "chr22")


# ---------------------------------------------------------------------------
# Karyotype-mode constants
# ---------------------------------------------------------------------------

#: Canonical autosomes used by karyotype-mode normalisation (the
#: autosomal median anchors CN = 2.0).
KARY_AUTOSOMES: tuple[str, ...] = tuple(f"chr{i}" for i in range(1, 23))

#: Chromosomes whose p-arm label is dropped when rendering arm-tick
#: labels in karyotype mode. Acrocentric short arms (chr13/14/15/21/22)
#: are too short to label without colliding with the q-label;
#: chrY's p-arm is similarly tiny.
KARY_ACROCENTRIC_FOR_LABELS: tuple[str, ...] = (
    "chr13", "chr14", "chr15", "chr21", "chr22", "chrY",
)

#: Sentinel value used by the bundled 10 kb GC tables to mark bins
#: where every source base was N (no GC% defined). Joined-on-bin
#: lookups for cov rows that fall outside every GC bin also use this
#: sentinel.
GC_MISSING: int = 255

#: Karyotype-mode palette (Direction A "Paper"). Distinct from the
#: SV-mode ISCN greyscale; applied per-call so matplotlib rcParams
#: stay unchanged between modes.
KARY_PAPER: str   = "#FAF8F4"
KARY_INK: str     = "#1F2024"
KARY_INK_2: str   = "#5C5D63"
KARY_INK_3: str   = "#9A9892"
KARY_RULE: str    = "#E4E0D8"
KARY_SCATTER: str = "#3D3D45"
KARY_ROSE: str    = "#DC5A99"
KARY_OXFORD: str  = "#1F4F7A"
KARY_MIST: str    = "#DDE6D8"

#: Warm-palette cytoband stain map used by the karyotype cytoband
#: silhouette. Distinct from the SV-mode greyscale ``CYTOBAND_COLORS``
#: because the karyotype strip is the chromosome's primary visual
#: identity (rather than a backdrop under density bars).
KARY_CYTO_STAIN_COLOR: dict[str, str] = {
    "gneg":    KARY_PAPER,
    "gpos25":  "#D9D5CC",
    "gpos50":  KARY_INK_3,
    "gpos75":  KARY_INK_2,
    "gpos100": KARY_INK,
    "acen":    "#A33A3A",
    "gvar":    "#B89A6F",
    "stalk":   KARY_ROSE,
}

#: Font stacks for karyotype-mode tick / label text. Applied per-call;
#: matplotlib walks the list and picks the first installed font.
KARY_FONT_SANS: tuple[str, ...] = (
    "IBM Plex Sans", "Helvetica", "Arial", "DejaVu Sans", "sans-serif",
)
KARY_FONT_MONO: tuple[str, ...] = (
    "IBM Plex Mono", "Menlo", "Consolas", "DejaVu Sans Mono", "monospace",
)

#: Karyotype-mode figure / layout constants. Locked at port time.
KARY_FIG_W: float                          = 18.0
KARY_FIG_H_GENOME_BAF: float               = 6.4
KARY_FIG_H_GENOME_ONLY: float              = 4.8
KARY_FIG_H_REGION_BAF: float               = 5.2
KARY_FIG_H_REGION_ONLY: float              = 3.8
KARY_YLABEL_X: float                       = -0.030
KARY_HEIGHT_RATIOS_CN_BAF: tuple[float, float]        = (2.4, 1.0)
KARY_HEIGHT_RATIOS_BAND_CN_BAF: tuple[float, float, float] = (0.35, 2.4, 1.0)
KARY_HEIGHT_RATIOS_BAND_CN: tuple[float, float]       = (0.35, 2.4)


# ---------------------------------------------------------------------------
# Compound-het constants (locked 2026-05-02; do not relitigate)
# ---------------------------------------------------------------------------
# Visual spec ported verbatim from the laptop-only prototype at
# /Users/martin/Documents/compound-het-plotter/plot_gene.py (draw_panel,
# lines 266-415). Three styling variants were tried and these were
# locked. P/LP are merged to a single colour per the v0.3 handoff.

#: Lollipop fill colours by canonicalised ClinVar significance. Every
#: missense canonical-transcript variant gets one of these. Pathogenic
#: and Likely_pathogenic share ``"p_or_lp"`` so the dot-pair-in-trans
#: read of "two red dots on opposite haps" works whether ClinVar calls
#: them P or LP. Variants with no ClinVar match get the default grey.
CLNSIG_COLOR: dict[str, str] = {
    "p_or_lp":     "#c0143c",
    "vus":         "#f4a013",
    "conflicting": "#d9c200",
    "benign":      "#5fa860",
}
DEFAULT_CLNSIG_COLOR: str = "#888888"

#: Phase-block rectangle aesthetic. Mint fill + dark-green border,
#: linewidth 1.2. Each rectangle spans one WhatsHap PS group across
#: both haplotypes, with off-edge arrows when the block stretches past
#: the gene window.
BLOCK_FILL: str = "#e2f0e3"
BLOCK_EDGE: str = "#3f6e44"
BLOCK_EDGE_LW: float = 1.2

#: Non-missense markers: ``x`` on the hap line, dark grey, edge
#: width 1.4. Used for synonymous / intronic / 5'/3' UTR variants
#: that landed in the same phase block as a missense lollipop. The
#: ``x`` glyph (vs the original tick ``|``) keeps these visually
#: distinct from the dark-grey rectangles of the exon track at the
#: top of the panel.
NON_MISSENSE_COLOR: str = "#222222"
NON_MISSENSE_MARKER: str = "x"
NON_MISSENSE_SIZE: float = 7.5
NON_MISSENSE_EDGE_WIDTH: float = 1.4

#: Missense lollipop dimensions. Stem extends from the hap line at
#: ``HAP_Y[hap-1]`` *down* to ``LOLLIPOP_TIP_Y[hap-1]``; the filled
#: circle sits at the tip with a thin black edge. Both H1 and H2
#: lollipops point downward so the visual reading is consistent
#: ("dots below their hap"), with H2 nudged up so the four bands
#: (exon track, H1, H1 dots, H2, H2 dots) are evenly spaced.
LOLLIPOP_MS: float = 8.0
LOLLIPOP_EDGE_LW: float = 0.4
HAP_Y: tuple[float, float] = (0.85, -0.4)            # H1, H2
LOLLIPOP_TIP_Y: tuple[float, float] = (0.25, -1.0)   # H1, H2

#: Exon-track aesthetic. Thin grey connector + IGV-style blue solid
#: rectangles per exon, drawn at the top of each gene panel. The
#: blue (vs the previous dark grey) reads as "gene track" in the
#: same visual register cytogeneticists already see in IGV.
EXON_Y: float = 1.75
EXON_H: float = 0.14
EXON_CONNECTOR_COLOR: str = "#888888"
EXON_FILL: str = "#1E5BA8"

#: Output PNG width cap. Anthropic's many-image upload tops out around
#: 2000 px; cap at 1800 px to leave a margin. ``target_dpi`` is
#: computed as ``min(150, MAX_PX_WIDTH / fig.get_figwidth())`` at
#: render time so wider figures step down dpi rather than producing
#: an oversize PNG.
MAX_PX_WIDTH: int = 1800

#: Consequence terms treated as missense. The renderer only draws
#: lollipops for variants whose VEP Consequence is in this set AND
#: whose CSQ entry is the canonical transcript. Everything else
#: renders as a non-missense tick.
MISSENSE_CONSEQUENCES: frozenset[str] = frozenset({"missense_variant"})


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class BND:
    """A single Sniffles BND record (one side of a translocation).

    Sniffles emits two BND records per inter-chromosomal rearrangement
    (one anchored on each chromosome). Reciprocal pairs are typically
    collapsed via :func:`deduplicate_reciprocal` before plotting.

    Attributes
    ----------
    chr1, pos1 : str, int
        The "near" breakend, taken from CHROM/POS of the VCF row.
    chr2, pos2 : str, int
        The "mate" breakend, parsed from the ALT bracket notation.
    orientation : str
        One of ``'++'``, ``'+-'``, ``'-+'``, ``'--'`` encoding the
        join orientation (see :func:`parse_alt_for_mate`).
    support : int
        Number of supporting reads (``INFO/SUPPORT``).
    vaf : float
        Variant allele fraction (``INFO/VAF``).
    filter_ : str
        ``FILTER`` column value (``"PASS"`` or a sniffles filter tag).
    sv_id : str
        VCF ID column (e.g. ``Sniffles2.BND.17A9S6``).
    coverage : list[float | None]
        5-vector ``[upstream, start, center, end, downstream]`` from
        ``INFO/COVERAGE``. Missing entries are ``None``.
    noise_flags : set[str]
        Populated by :func:`annotate_noise`. Possible members:
        ``"acrocentric"``, ``"cov_anomaly"``.
    """

    chr1: str
    pos1: int
    chr2: str
    pos2: int
    orientation: str
    support: int
    vaf: float
    filter_: str
    sv_id: str
    coverage: list = field(default_factory=list)
    noise_flags: set = field(default_factory=set)

    @property
    def is_pass(self) -> bool:
        """True iff this BND has FILTER == 'PASS'."""
        return self.filter_ == "PASS"

    @property
    def is_noise(self) -> bool:
        """True iff any noise flag has been set."""
        return bool(self.noise_flags)

    @property
    def max_coverage(self) -> float:
        """Maximum of the 5 COVERAGE values, ignoring missing entries."""
        vals = [c for c in self.coverage if c is not None]
        return max(vals) if vals else 0.0

    @property
    def canonical_key(self) -> tuple:
        """A pair-key invariant under the (chr1, chr2) swap.

        Used by :func:`deduplicate_reciprocal` to identify the two
        records describing the same translocation.
        """
        i1 = CHROM_ORDER.index(self.chr1) if self.chr1 in CHROM_SET else 99
        i2 = CHROM_ORDER.index(self.chr2) if self.chr2 in CHROM_SET else 99
        a, b = (i1, self.pos1), (i2, self.pos2)
        return tuple(sorted([a, b]))


@dataclass
class SV:
    """A non-BND structural variant (INS, DEL, DUP, or INV).

    Attributes
    ----------
    chrom, start, end : str, int, int
        Genomic coordinates of the event.
    svtype : str
        One of ``"INS"``, ``"DEL"``, ``"DUP"``, ``"INV"``.
    svlen : int
        Absolute event length in bp (always non-negative).
    filter_ : str
        ``FILTER`` column value.
    support : int, vaf : float, sv_id : str
        Same as :class:`BND`.
    coverage : list[float | None]
        5-vector from ``INFO/COVERAGE`` (see :class:`BND`).
    noise_flags : set[str]
        Populated by :func:`annotate_sv_noise` (``"cov_anomaly"``).
        :attr:`is_noise` is True iff the set is non-empty.
    """

    chrom: str
    start: int
    end: int
    svtype: str
    svlen: int
    filter_: str
    support: int
    vaf: float
    sv_id: str
    coverage: list = field(default_factory=list)
    noise_flags: set = field(default_factory=set)

    @property
    def is_pass(self) -> bool:
        """True iff this SV has FILTER == 'PASS'."""
        return self.filter_ == "PASS"

    @property
    def is_noise(self) -> bool:
        """True iff any noise flag has been set."""
        return bool(self.noise_flags)

    @property
    def max_coverage(self) -> float:
        """Maximum of the 5 COVERAGE values, ignoring missing entries."""
        vals = [c for c in self.coverage if c is not None]
        return max(vals) if vals else 0.0


@dataclass(frozen=True)
class Gene:
    """A canonical-transcript gene record for the compound-het mode.

    Loaded from the bundled per-build TSV (``data/canonical_exons.<ref>.tsv.gz``).
    One row per gene symbol; the canonical transcript is whichever the
    bundled table picks (MANE Select for hg38, NCBI "best refseq" for
    T2T). Coordinates are 0-based half-open BED-style.

    Attributes
    ----------
    symbol : str
        Gene symbol (e.g. ``"NEB"``, ``"APOB"``). Used as the primary
        key throughout compound-het mode.
    chrom : str
        UCSC-style chromosome name (``"chr1"`` ... ``"chr22"``,
        ``"chrX"``, ``"chrY"``). Already normalised on load.
    start, end : int
        Gene-body span (0-based half-open). Used as the panel x-limits.
    strand : str
        ``"+"`` or ``"-"``. Displayed in the panel title.
    transcript_id : str
        Canonical-transcript ID (e.g. ``"NM_004543.4"`` or
        ``"ENST00000397345.8"``); shown in run-metadata.
    canonical_exons : tuple[tuple[int, int], ...]
        Tuple of ``(start, end)`` exon spans, sorted by start. Empty
        tuple is allowed (gene panels still render with hap lines and
        any phased hets in range, just without an exon track).
    """
    symbol: str
    chrom: str
    start: int
    end: int
    strand: str
    transcript_id: str
    canonical_exons: tuple[tuple[int, int], ...] = ()


@dataclass(frozen=True)
class PhasedVariant:
    """A single phased het small variant for the compound-het mode.

    Produced by :func:`read_phased_vcf` (one record per kept VCF row).
    ``clnsig`` is filled in a second pass after ClinVar lookup.

    The hap convention matches the prototype's SQL: ``"1|0"`` puts the
    ALT on hap 1 (``variant_hap=1``); ``"0|1"`` puts the ALT on hap 2
    (``variant_hap=2``). Two variants in the same PS with the same
    ``variant_hap`` are in cis; opposite ``variant_hap`` are in trans.

    Attributes
    ----------
    chrom, pos, ref, alt : str, int, str, str
        Variant coordinates and alleles. ``chrom`` is normalised to
        UCSC style on load.
    gt : str
        Genotype string from FORMAT/GT, exactly as read. Always
        ``"0|1"`` or ``"1|0"`` (other patterns are filtered upstream).
    ps : int
        Phase set ID from FORMAT/PS. Variants with the same ``ps``
        and the same ``chrom`` form one phase block.
    variant_hap : int
        ``1`` if ``gt == "1|0"``, ``2`` if ``gt == "0|1"``.
    consequence : str
        VEP Consequence (canonical-transcript CSQ entry, e.g.
        ``"missense_variant"``, ``"synonymous_variant"``).
    is_canonical_transcript : bool
        True if a CSQ entry with ``CANONICAL == "YES"`` was found
        for this record. False means we fell back to the first CSQ
        entry; the renderer treats those as non-missense ticks
        regardless of ``consequence``.
    hgvs_p : str | None
        VEP HGVSp string for the canonical transcript, or ``None``.
    hgvs_c : str | None
        VEP HGVSc string for the canonical transcript, or ``None``.
        Surfaced in the run-metadata of the rendered report; not
        currently used as a lookup key.
    gene_symbol : str | None
        VEP SYMBOL for the canonical transcript, or ``None``.
    feature : str | None
        VEP Feature (transcript ID) for the canonical transcript, or
        ``None``.
    clnsig : str | None
        Canonicalised ClinVar significance: one of ``"p_or_lp"``,
        ``"vus"``, ``"conflicting"``, ``"benign"``, ``"other"``, or
        ``None`` (no ClinVar match). Filled by ClinVar lookup in
        :func:`compound_het_main`.
    """
    chrom: str
    pos: int
    ref: str
    alt: str
    gt: str
    ps: int
    variant_hap: int
    consequence: str
    is_canonical_transcript: bool
    hgvs_p: str | None = None
    hgvs_c: str | None = None
    gene_symbol: str | None = None
    feature: str | None = None
    clnsig: str | None = None

    def with_clnsig(self, clnsig: str | None) -> "PhasedVariant":
        """Return a copy of this variant with ``clnsig`` replaced.

        ``PhasedVariant`` is frozen for hashability (used as dict
        keys in pair-classification); annotation runs as a second
        pass that rebuilds the list rather than mutating in place.
        """
        return PhasedVariant(
            chrom=self.chrom, pos=self.pos, ref=self.ref, alt=self.alt,
            gt=self.gt, ps=self.ps, variant_hap=self.variant_hap,
            consequence=self.consequence,
            is_canonical_transcript=self.is_canonical_transcript,
            hgvs_p=self.hgvs_p, hgvs_c=self.hgvs_c,
            gene_symbol=self.gene_symbol, feature=self.feature,
            clnsig=clnsig,
        )


@dataclass(frozen=True)
class PhaseBlock:
    """One WhatsHap phase set (PS) within one gene window.

    Computed from the variants assigned to a gene: group by ``ps``,
    take ``min(pos)`` and ``max(pos)`` as the rendered extent. The
    renderer optionally extends ``start``/``end`` to the WhatsHap-
    reported full block extent (which may stretch past the gene
    window), with off-edge arrows when so.

    Attributes
    ----------
    ps : int
        Phase set ID.
    start, end : int
        Rendered extent (0-based half-open). Currently the min/max
        position of variants in this block within the gene window.
    n_phased : int
        Number of phased het variants in the block (across both haps).
    """
    ps: int
    start: int
    end: int
    n_phased: int


# ---------------------------------------------------------------------------
# VCF parsing
# ---------------------------------------------------------------------------

def parse_alt_for_mate(alt: str) -> tuple[str, int, str]:
    """Parse a VCF BND ALT bracket notation into mate position + orientation.

    The four possible forms encode different join orientations:

    ===================  ==================
    ALT                  orientation
    ===================  ==================
    ``N[chr:pos[``       ``'++'``
    ``N]chr:pos]``       ``'+-'``
    ``[chr:pos[N``       ``'-+'``
    ``]chr:pos]N``       ``'--'``
    ===================  ==================

    Parameters
    ----------
    alt : str
        The ALT field from a VCF BND record.

    Returns
    -------
    tuple of (str, int, str)
        ``(mate_chr, mate_pos, orientation)``.

    Raises
    ------
    ValueError
        If ``alt`` does not match any of the four expected forms.
    """
    s = alt.strip()
    if s.startswith("N["):
        chrm, pos = s[2:-1].split(":")
        return chrm, int(pos), "++"
    if s.startswith("N]"):
        chrm, pos = s[2:-1].split(":")
        return chrm, int(pos), "+-"
    if s.startswith("["):
        chrm, pos = s[1:-2].split(":")
        return chrm, int(pos), "-+"
    if s.startswith("]"):
        chrm, pos = s[1:-2].split(":")
        return chrm, int(pos), "--"
    raise ValueError(f"unrecognized BND ALT: {alt!r}")


def parse_info(info: str) -> dict:
    """Parse a VCF INFO field into a ``{key: value}`` dict.

    Flag-only entries (no ``=``) are stored as ``True``. Values are
    kept as raw strings - the caller is responsible for casting.
    """
    out: dict = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
        else:
            out[kv] = True
    return out


def parse_coverage(coverage_str: str | None) -> list:
    """Parse a comma-separated COVERAGE string to a list of floats.

    Sniffles emits coverage as five comma-separated values
    (``upstream,start,center,end,downstream``). Empty entries or the
    literal strings ``"null"`` / ``"."`` are converted to ``None``.

    Returns an empty list if ``coverage_str`` is falsy.
    """
    if not coverage_str:
        return []
    out = []
    for p in coverage_str.split(","):
        p = p.strip()
        if p in ("", "null", "."):
            out.append(None)
            continue
        try:
            out.append(float(p))
        except ValueError:
            out.append(None)
    return out


def open_text(path: Path):
    """Open a (possibly gzipped or xz-compressed) text file for reading."""
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(path, "rt")
    if p.endswith(".xz"):
        return lzma.open(path, "rt")
    return open(path, "rt")


def _format_field(format_col: str, sample_col: str, key: str) -> str | None:
    """Fetch one value from a VCF FORMAT/SAMPLE column pair by key name."""
    keys = format_col.split(":")
    vals = sample_col.split(":")
    try:
        i = keys.index(key)
    except ValueError:
        return None
    if i >= len(vals):
        return None
    return vals[i]


def _build_event(
    fields: list[str],
    info: dict,
    *,
    support: int,
    vaf: float,
    coverage: list,
) -> BND | SV | None:
    """Common BND / SV record construction once a parser has resolved
    the caller-specific support / VAF / coverage fields.
    """
    chrom, pos, vid, _ref, alt, _qual, flt, _info_str = fields[:8]
    chrom = _normalize_chrom(chrom)
    svt = info.get("SVTYPE")
    if not svt:
        return None
    try:
        pos_i = int(pos)
    except ValueError:
        return None
    if svt == "BND":
        try:
            mate_chr, mate_pos, ori = parse_alt_for_mate(alt)
        except ValueError:
            return None
        mate_chr = _normalize_chrom(mate_chr)
        return BND(
            chr1=chrom, pos1=pos_i,
            chr2=mate_chr, pos2=mate_pos,
            orientation=ori, support=support, vaf=vaf,
            filter_=flt, sv_id=vid, coverage=coverage,
        )
    if svt in SV_TYPES:
        try:
            end_i = int(info.get("END", pos_i))
        except (ValueError, TypeError):
            end_i = pos_i
        try:
            svlen = abs(int(info.get("SVLEN", 0)))
        except (ValueError, TypeError):
            svlen = 0
        return SV(
            chrom=chrom, start=pos_i, end=end_i,
            svtype=svt, svlen=svlen,
            filter_=flt, support=support, vaf=vaf, sv_id=vid,
            coverage=coverage,
        )
    return None


def _parse_record_sniffles2(fields: list[str], info: dict) -> BND | SV | None:
    """Sniffles2: INFO/SUPPORT, INFO/VAF, INFO/COVERAGE."""
    try:
        support = int(info.get("SUPPORT", 0))
    except (ValueError, TypeError):
        support = 0
    try:
        vaf = float(info.get("VAF", 0))
    except (ValueError, TypeError):
        vaf = 0.0
    return _build_event(fields, info, support=support, vaf=vaf,
                         coverage=parse_coverage(info.get("COVERAGE")))


def _parse_record_sniffles1(fields: list[str], info: dict) -> BND | SV | None:
    """Sniffles1: INFO/SUPPORT, derive VAF from FORMAT/DR + DV."""
    try:
        support = int(info.get("SUPPORT", 0))
    except (ValueError, TypeError):
        support = 0
    vaf = 0.0
    if len(fields) >= 10:
        dr_str = _format_field(fields[8], fields[9], "DR")
        dv_str = _format_field(fields[8], fields[9], "DV")
        try:
            dr = int(dr_str) if dr_str else 0
            dv = int(dv_str) if dv_str else 0
            if dr + dv > 0:
                vaf = dv / (dr + dv)
        except (ValueError, TypeError):
            pass
    return _build_event(fields, info, support=support, vaf=vaf, coverage=[])


def _parse_record_cutesv(fields: list[str], info: dict) -> BND | SV | None:
    """cuteSV: INFO/RE -> support, INFO/AF -> VAF, no COVERAGE."""
    try:
        support = int(info.get("RE", 0))
    except (ValueError, TypeError):
        support = 0
    try:
        vaf = float(info.get("AF", 0))
    except (ValueError, TypeError):
        vaf = 0.0
    return _build_event(fields, info, support=support, vaf=vaf, coverage=[])


def _parse_record_svim(fields: list[str], info: dict) -> BND | SV | None:
    """SVIM: INFO/SUPPORT, INFO/AF, no COVERAGE."""
    try:
        support = int(info.get("SUPPORT", 0))
    except (ValueError, TypeError):
        support = 0
    try:
        vaf = float(info.get("AF", 0))
    except (ValueError, TypeError):
        vaf = 0.0
    return _build_event(fields, info, support=support, vaf=vaf, coverage=[])


def _parse_record_pbsv(fields: list[str], info: dict) -> BND | SV | None:
    """pbsv: support and VAF derived from FORMAT/AD = (ref, alt)."""
    support = 0
    vaf = 0.0
    if len(fields) >= 10:
        ad_str = _format_field(fields[8], fields[9], "AD")
        if ad_str:
            try:
                parts = [int(p) for p in ad_str.split(",")]
                if len(parts) >= 2:
                    ref, alt_count = parts[0], parts[1]
                    support = alt_count
                    if ref + alt_count > 0:
                        vaf = alt_count / (ref + alt_count)
            except ValueError:
                pass
    return _build_event(fields, info, support=support, vaf=vaf, coverage=[])


def _parse_record_nanovar(fields: list[str], info: dict) -> BND | SV | None:
    """NanoVar: INFO/SR -> support, INFO/AF -> VAF."""
    try:
        support = int(info.get("SR", 0))
    except (ValueError, TypeError):
        support = 0
    try:
        vaf = float(info.get("AF", 0))
    except (ValueError, TypeError):
        vaf = 0.0
    return _build_event(fields, info, support=support, vaf=vaf, coverage=[])


CALLER_PARSERS: dict = {
    "sniffles2": _parse_record_sniffles2,
    "sniffles1": _parse_record_sniffles1,
    "cutesv":    _parse_record_cutesv,
    "svim":      _parse_record_svim,
    "pbsv":      _parse_record_pbsv,
    "nanovar":   _parse_record_nanovar,
}


def read_vcf(
    vcf_path: Path, caller: str = "auto",
) -> tuple[dict, list[BND], list[SV], float]:
    """Parse a long-read SV VCF using the specified caller's INFO map.

    Parameters
    ----------
    vcf_path : Path
        Path to the VCF (plain text or .vcf.gz).
    caller : str, optional
        One of ``"auto"`` (default), ``"sniffles2"``, ``"sniffles1"``,
        ``"cutesv"``, ``"svim"``, ``"pbsv"``, or ``"nanovar"``. ``auto``
        runs :func:`detect_caller` and falls back to ``"sniffles2"``
        if no signature matches.

    Returns
    -------
    contigs : dict[str, int]
        ``{chrom: length}`` from ``##contig=`` header lines.
    bnds : list[BND]
        All BND records.
    svs : list[SV]
        All non-BND SV records (INS / DEL / DUP / INV).
    median_coverage : float
        Median of ``INFO/COVERAGE[center]`` from INS+DEL records (only
        meaningful for Sniffles2; defaults to ``1.0`` for callers
        without COVERAGE).
    """
    if caller == "auto":
        detected, _ = detect_caller(vcf_path)
        caller = "sniffles2" if detected == "unknown" else detected
    if caller not in CALLER_PARSERS:
        raise ValueError(
            f"unknown caller {caller!r}; "
            f"choose from {sorted(CALLER_PARSERS)} or 'auto'",
        )
    parse_record = CALLER_PARSERS[caller]

    contigs: dict[str, int] = {}
    bnds: list[BND] = []
    svs: list[SV] = []
    nonbnd_center_cov: list[float] = []
    contig_re = re.compile(r"ID=([^,>]+).*length=(\d+)")
    with open_text(vcf_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##contig="):
                m = contig_re.search(line)
                if m:
                    contigs[_normalize_chrom(m.group(1))] = int(m.group(2))
                continue
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 8:
                continue
            inf = parse_info(f[7])
            record = parse_record(f, inf)
            if record is None:
                continue
            if isinstance(record, BND):
                bnds.append(record)
            else:
                svs.append(record)
                if (record.svtype in ("INS", "DEL")
                        and len(record.coverage) >= 3
                        and record.coverage[2] is not None
                        and record.coverage[2] > 0):
                    nonbnd_center_cov.append(record.coverage[2])
    median_cov = float(np.median(nonbnd_center_cov)) if nonbnd_center_cov else 1.0
    return contigs, bnds, svs, median_cov


_KARY_HET_GTS: frozenset[str] = frozenset({"0/1", "1/0", "0|1", "1|0"})


def _parse_baf(format_kv: dict[str, str]) -> float | None:
    """Extract a B-allele fraction from a parsed FORMAT/SAMPLE mapping.

    Prefers ``AF``; falls back to ``AD`` parsed as
    ``ref_count,alt_count`` and computed as ``alt / (ref + alt)``.
    Returns ``None`` if neither field is present or parseable.
    """
    af = format_kv.get("AF")
    if af and af not in (".", ""):
        try:
            return float(af.split(",")[0])
        except ValueError:
            pass
    ad = format_kv.get("AD")
    if ad and ad not in (".", ""):
        try:
            counts = [int(x) for x in ad.split(",")]
        except ValueError:
            return None
        if len(counts) >= 2 and counts[0] + counts[1] > 0:
            return counts[1] / (counts[0] + counts[1])
    return None


def read_baf_vcf(path: Path, min_dp: int) -> pd.DataFrame:
    """Stream het allele fractions from a small-variant VCF as plain text.

    No bcftools dependency; reads gzipped or plain VCF via the
    :func:`open_text` helper. Returns a DataFrame with columns
    ``chrom, pos, baf`` for the karyotype-mode BAF panel.

    Filters applied per record:

    - ``FILTER == "PASS"``
    - heterozygous biallelic GT (``0/1``, ``1/0``, ``0|1``, ``1|0``)
    - ``FORMAT/DP >= min_dp``
    - chrom in :data:`CHROM_ORDER` (canonical chr1-22, chrX, chrY)

    The BAF value is taken from ``FORMAT/AF`` if present, otherwise
    computed from ``FORMAT/AD`` as ``alt / (ref + alt)``. Records
    with neither field parseable are skipped.
    """
    rows: list[tuple[str, int, float]] = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom = f[0]
            if chrom not in CHROM_SET:
                continue
            if f[6] != "PASS":
                continue
            format_keys = f[8].split(":")
            sample_vals = f[9].split(":")
            kv = dict(zip(format_keys, sample_vals))
            if kv.get("GT", "") not in _KARY_HET_GTS:
                continue
            try:
                dp = int(kv.get("DP", "0"))
            except ValueError:
                continue
            if dp < min_dp:
                continue
            baf = _parse_baf(kv)
            if baf is None:
                continue
            try:
                pos = int(f[1])
            except ValueError:
                continue
            rows.append((chrom, pos, float(baf)))
    if not rows:
        return pd.DataFrame({
            "chrom": pd.Categorical(
                [], categories=CHROM_ORDER, ordered=True,
            ),
            "pos": np.array([], dtype=np.int64),
            "baf": np.array([], dtype=np.float64),
        })
    df = pd.DataFrame(rows, columns=["chrom", "pos", "baf"])
    df["chrom"] = pd.Categorical(
        df["chrom"], categories=CHROM_ORDER, ordered=True,
    )
    return df.sort_values(["chrom", "pos"]).reset_index(drop=True)


def detect_vcf_mode(vcf_path: Path) -> str:
    """Inspect the VCF header and pick the right plotting mode.

    Discriminators:

    - ``##INFO=<ID=CSQ,...>`` AND ``##FORMAT=<ID=PS,...>`` →
      ``"compound-het"`` (phased small-variant VCF with VEP CSQ
      annotation; the only shape that compound-het mode can plot).
    - ``##INFO=<ID=SVTYPE,...>`` (and not the CSQ+PS pair above) →
      ``"sv"`` (long-read SV VCF; Sniffles / cuteSV / SVIM / pbsv /
      NanoVar all emit SVTYPE).

    A VCF with both shapes (rare — would need to be VEP-annotated
    SVs that are also phased) prefers ``"compound-het"`` because the
    CSQ + PS pair is more specific.

    Raises ``ValueError`` if neither shape is present so molamola
    refuses rather than render a misleading default plot.
    """
    has_svtype = False
    has_csq = False
    has_ps = False
    with open_text(vcf_path) as fh:
        for line in fh:
            if not line.startswith("##"):
                break
            if line.startswith("##INFO=<ID=SVTYPE,"):
                has_svtype = True
            elif line.startswith("##INFO=<ID=CSQ,"):
                has_csq = True
            elif line.startswith("##FORMAT=<ID=PS,"):
                has_ps = True
    if has_csq and has_ps:
        return "compound-het"
    if has_svtype:
        return "sv"
    raise ValueError(
        f"VCF at {vcf_path} doesn't match either supported shape: "
        f"SV mode needs ##INFO=<ID=SVTYPE,...>; compound-het mode "
        f"needs ##INFO=<ID=CSQ,...> AND ##FORMAT=<ID=PS,...>."
    )


def detect_caller(vcf_path: Path, n_records: int = 10) -> tuple[str, str]:
    """Detect the SV caller that produced ``vcf_path``.

    Two-pass detection: first parse ``##source=`` for an explicit hint,
    then fingerprint the INFO fields seen in the first ``n_records``
    data rows. The fallback handles bcftools-merged or re-headered VCFs
    that lose ``##source`` entirely.

    Returns
    -------
    tuple[str, str]
        ``(caller, basis)`` where ``caller`` is one of
        ``"sniffles2"``, ``"sniffles1"``, ``"cutesv"``, ``"svim"``,
        ``"pbsv"``, ``"nanovar"``, or ``"unknown"``. ``basis`` is a
        short human-readable explanation of the decision (intended for
        logging at startup).

    Notes
    -----
    Phase A is warn-only - the parser still treats every caller as
    Sniffles2 in practice. Phase B will add per-caller dispatch.
    """
    source_line = ""
    info_keys: set[str] = set()
    n_seen = 0
    with open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("##source="):
                source_line = line.rstrip("\n")
                continue
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 8:
                continue
            inf = parse_info(f[7])
            info_keys.update(inf.keys())
            n_seen += 1
            if n_seen >= n_records:
                break

    src = source_line.lower()

    # Explicit ##source mentions
    if "sniffles2" in src or "sniffles_2" in src or "sniffles 2" in src:
        return ("sniffles2", "##source mentions Sniffles2")
    if "cutesv" in src:
        return ("cutesv", "##source mentions cuteSV")
    if "svim" in src:
        return ("svim", "##source mentions SVIM")
    if "pbsv" in src:
        return ("pbsv", "##source mentions pbsv")
    if "nanovar" in src:
        return ("nanovar", "##source mentions NanoVar")

    # Fingerprint INFO field signatures
    if {"SUPPORT", "COVERAGE", "STDEV_POS"} <= info_keys:
        return ("sniffles2", "INFO has SUPPORT + COVERAGE + STDEV_POS")
    if {"RE", "STRAND"} <= info_keys:
        return ("cutesv", "INFO has RE + STRAND")
    if "BREAKEND_TYPE" in info_keys:
        return ("nanovar", "INFO has BREAKEND_TYPE")
    if "STD_SPAN" in info_keys or "STD_POS" in info_keys:
        return ("svim", "INFO has STD_SPAN / STD_POS")
    if "SUPPORT" in info_keys and "COVERAGE" not in info_keys:
        return ("sniffles1", "INFO has SUPPORT without COVERAGE")
    if "sniffles" in src:
        return ("sniffles1", "##source mentions Sniffles, no v2 fields seen")

    return ("unknown", "no caller-specific signature matched")


def deduplicate_reciprocal(bnds: list[BND]) -> list[BND]:
    """Collapse the two records describing each translocation to one.

    Each inter-chromosomal rearrangement produces two BND records in
    the VCF (one anchored on each side). They share a
    :attr:`BND.canonical_key`, so we keep the first occurrence of
    each unique key.
    """
    seen: dict = {}
    for b in bnds:
        seen.setdefault(b.canonical_key, b)
    return list(seen.values())


# ---------------------------------------------------------------------------
# Compound-het loaders & selection
# ---------------------------------------------------------------------------

# Match the ``Format: A|B|C|...`` substring inside a VEP CSQ INFO line.
# VEP emits the format list inside the Description="..." string, so the
# field list ends at the closing double-quote.
_CSQ_FORMAT_RE = re.compile(r"Format:\s*([^\"]+?)\s*\"?\s*>")


def parse_csq_format(header_line: str) -> list[str] | None:
    """Extract the VEP CSQ field order from a VCF header line.

    Returns the field list (e.g. ``["Allele", "Consequence", "SYMBOL",
    ...]``) when ``header_line`` is the ``##INFO=<ID=CSQ,...>`` line,
    or ``None`` when the line is not a CSQ definition.

    Resilient to VEP plugin reordering: callers should look up fields
    by name, not by index, since plugin order varies between runs.
    """
    if not header_line.startswith("##INFO=<ID=CSQ,"):
        return None
    m = _CSQ_FORMAT_RE.search(header_line)
    if m is None:
        return None
    return [f.strip() for f in m.group(1).split("|") if f.strip()]


def _parse_csq_entry(entry: str, csq_fields: list[str]) -> dict[str, str]:
    """Split one ``|``-delimited CSQ entry into a ``{field: value}`` dict.

    Trailing fields beyond ``len(csq_fields)`` are dropped; missing
    trailing fields stay as empty strings. Field values are kept as
    raw strings (no casting).
    """
    parts = entry.split("|")
    out: dict[str, str] = {}
    for i, name in enumerate(csq_fields):
        out[name] = parts[i] if i < len(parts) else ""
    return out


def _pick_canonical_csq(
    csq_value: str,
    csq_fields: list[str],
) -> tuple[dict[str, str], bool]:
    """Pick the canonical-transcript CSQ entry from a multi-entry CSQ.

    VEP packs all transcript annotations for one variant into a single
    ``CSQ=A|...,B|...,C|...`` INFO value. We walk the entries and pick
    the first one with ``CANONICAL == "YES"``. If none has it, we fall
    back to the first entry and flag ``is_canonical=False`` so the
    renderer can treat it as a non-canonical record.

    Returns ``(entry_dict, is_canonical)``.
    """
    entries = csq_value.split(",")
    parsed = [_parse_csq_entry(e, csq_fields) for e in entries]
    for p in parsed:
        if p.get("CANONICAL", "") == "YES":
            return p, True
    if parsed:
        return parsed[0], False
    return {}, False


def is_missense_consequence(consequence: str) -> bool:
    """Return True iff ``consequence`` includes ``missense_variant``.

    VEP joins multi-consequence annotations with ``&`` (e.g.
    ``"missense_variant&splice_region_variant"``); we treat any string
    that contains ``missense_variant`` as a token as missense.
    """
    return "missense_variant" in consequence.split("&")


def read_phased_vcf(
    path: Path,
) -> tuple[list[PhasedVariant], dict[str, object]]:
    """Read a phased small-variant VCF with VEP CSQ annotation.

    Behaviour:

    - Hard refusal (raises ``ValueError``) when the header lacks
      ``##INFO=<ID=CSQ,...>`` (no annotation we can colour-code) or
      ``##FORMAT=<ID=PS,...>`` (no phase info we can plot).
    - Per-record skips (counted, never raised): unphased GT
      (no ``|``), homozygous (``0|0`` / ``1|1``), multi-allelic
      ALT (``,`` in ALT column).
    - Canonical transcript pick: prefer the CSQ entry with
      ``CANONICAL == "YES"``; fall back to the first entry with
      ``is_canonical_transcript=False`` so the renderer can treat
      it as a non-canonical record.

    Parameters
    ----------
    path : Path
        Phased VCF (``.vcf`` or ``.vcf.gz``).

    Returns
    -------
    variants : list[PhasedVariant]
        Phased het records, autosome filter NOT applied (caller
        decides; explicit ``--gene FOO`` on chrX must still render).
    summary : dict
        Run-metadata for the HTML report. Keys: ``csq_fields``,
        ``has_ps_in_header``, ``unphased``, ``non_het``,
        ``multi_allelic``, ``no_canonical_csq``, ``kept``.

    Raises
    ------
    ValueError
        On file-level CSQ or PS absence.
    """
    csq_fields: list[str] | None = None
    has_ps_in_header = False
    refusal: dict[str, int] = {
        "unphased": 0,
        "non_het": 0,
        "multi_allelic": 0,
        "no_canonical_csq": 0,
        "kept": 0,
    }
    variants: list[PhasedVariant] = []

    with open_text(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith("##"):
                if csq_fields is None:
                    parsed = parse_csq_format(line)
                    if parsed is not None:
                        csq_fields = parsed
                if line.startswith("##FORMAT=<ID=PS,"):
                    has_ps_in_header = True
                continue
            if line.startswith("#CHROM"):
                if csq_fields is None:
                    raise ValueError(
                        "VCF has no CSQ INFO field; molamola compound-het "
                        "requires VEP-annotated input"
                    )
                if not has_ps_in_header:
                    raise ValueError(
                        "VCF is unphased (no ##FORMAT=<ID=PS,...> in "
                        "header); supply a WhatsHap/HiPhase-phased VCF"
                    )
                continue
            if line.startswith("#"):
                continue

            # Data row.
            fields = line.split("\t")
            if len(fields) < 10:
                continue
            chrom_raw, pos, _vid, ref, alt, _qual, _flt, info_str, fmt, sample = fields[:10]

            if "," in alt:
                refusal["multi_allelic"] += 1
                continue

            gt = _format_field(fmt, sample, "GT")
            if gt is None or "|" not in gt:
                refusal["unphased"] += 1
                continue
            if gt not in ("0|1", "1|0"):
                refusal["non_het"] += 1
                continue
            ps_str = _format_field(fmt, sample, "PS")
            if ps_str is None or ps_str in (".", ""):
                refusal["unphased"] += 1
                continue
            try:
                ps = int(ps_str)
            except ValueError:
                refusal["unphased"] += 1
                continue
            try:
                pos_i = int(pos)
            except ValueError:
                continue

            info = parse_info(info_str)
            csq_raw = info.get("CSQ")
            if not isinstance(csq_raw, str):
                refusal["no_canonical_csq"] += 1
                continue
            csq_entry, is_canonical = _pick_canonical_csq(csq_raw, csq_fields)
            if not csq_entry:
                refusal["no_canonical_csq"] += 1
                continue

            variant_hap = 1 if gt == "1|0" else 2
            v = PhasedVariant(
                chrom=_normalize_chrom(chrom_raw),
                pos=pos_i,
                ref=ref,
                alt=alt,
                gt=gt,
                ps=ps,
                variant_hap=variant_hap,
                consequence=csq_entry.get("Consequence", ""),
                is_canonical_transcript=is_canonical,
                hgvs_p=csq_entry.get("HGVSp") or None,
                hgvs_c=csq_entry.get("HGVSc") or None,
                gene_symbol=csq_entry.get("SYMBOL") or None,
                feature=csq_entry.get("Feature") or None,
                clnsig=None,
            )
            variants.append(v)
            refusal["kept"] += 1

    if csq_fields is None:
        raise ValueError(
            "VCF has no CSQ INFO field; molamola compound-het requires "
            "VEP-annotated input"
        )
    if not has_ps_in_header:
        raise ValueError(
            "VCF is unphased (no ##FORMAT=<ID=PS,...> in header); supply a "
            "WhatsHap/HiPhase-phased VCF"
        )

    summary: dict[str, object] = dict(refusal)
    summary["csq_fields"] = csq_fields
    summary["has_ps_in_header"] = has_ps_in_header
    return variants, summary


def find_canonical_exon_file() -> Path:
    """Return the bundled hg38 canonical-exon TSV path.

    Raises ``FileNotFoundError`` when the bundled file is missing —
    molamola does not auto-download.
    """
    here = Path(__file__).resolve().parent
    p = here / "data" / "canonical_exons.hg38.tsv.gz"
    if not p.exists():
        raise FileNotFoundError(
            f"bundled canonical-exon table not found at {p}; molamola "
            "does not auto-download. Provide one via --canonical-exons."
        )
    return p


def load_canonical_exons(path: Path) -> dict[str, Gene]:
    """Load a canonical-exon TSV into a ``{symbol: Gene}`` dict.

    Schema (gzipped TSV with header)::

        gene_symbol  chrom  start  end  strand  transcript_id
                                         exon_starts  exon_ends

    ``exon_starts`` and ``exon_ends`` are comma-separated 0-based
    half-open coordinates. Empty strings produce an empty exon tuple.
    Duplicate symbols keep the first occurrence; the loader does not
    silently merge — bundled tables are produced one row per symbol.
    """
    out: dict[str, Gene] = {}
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").rstrip("\r").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        required = ("gene_symbol", "chrom", "start", "end", "strand",
                    "transcript_id", "exon_starts", "exon_ends")
        missing = [c for c in required if c not in idx]
        if missing:
            raise ValueError(
                f"canonical-exon TSV at {path} missing columns: {missing}"
            )
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < len(header):
                continue
            symbol = cols[idx["gene_symbol"]]
            if symbol in out:
                continue
            chrom = _normalize_chrom(cols[idx["chrom"]])
            try:
                start = int(cols[idx["start"]])
                end = int(cols[idx["end"]])
            except ValueError:
                continue
            strand = cols[idx["strand"]]
            transcript_id = cols[idx["transcript_id"]]
            starts_raw = cols[idx["exon_starts"]]
            ends_raw = cols[idx["exon_ends"]]
            exon_starts = [int(x) for x in starts_raw.split(",") if x]
            exon_ends = [int(x) for x in ends_raw.split(",") if x]
            exons = tuple(zip(exon_starts, exon_ends))
            out[symbol] = Gene(
                symbol=symbol, chrom=chrom, start=start, end=end,
                strand=strand, transcript_id=transcript_id,
                canonical_exons=exons,
            )
    return out


def canon_clnsig(raw: str | None) -> str | None:
    """Canonicalise a ClinVar CLNSIG string to a colour-key bucket.

    Returns one of ``"p_or_lp"``, ``"vus"``, ``"conflicting"``,
    ``"benign"``, ``"other"``, or ``None`` when ``raw`` is empty.

    P/LP merge into a single bucket so the dot-pair-in-trans read
    of "two red dots on opposite haps" works whether ClinVar calls
    a variant Pathogenic or Likely_pathogenic.
    """
    if not raw:
        return None
    if "Conflicting" in raw:
        return "conflicting"
    if "Pathogenic" in raw or "Likely_pathogenic" in raw:
        return "p_or_lp"
    if "Uncertain" in raw:
        return "vus"
    if "Benign" in raw or "Likely_benign" in raw:
        return "benign"
    return "other"


def find_clinvar_file() -> Path:
    """Return the bundled ClinVar reduced-TSV path.

    The bundled file is a derived 5-column TSV (not the raw NCBI
    VCF) — see :func:`load_clinvar_lookup` for the schema. xz-
    compressed reduced TSV vs raw VCF is roughly 13 MB vs 191 MB.

    ClinVar is hg38-coordinate; the compound-het mode supports hg38
    only as of v0.3.1 (T2T deferred until a non-coordinate-based
    matching path lands).

    Raises ``FileNotFoundError`` when missing — bundled-only refs.
    """
    here = Path(__file__).resolve().parent
    p = here / "data" / "clinvar.hg38.tsv.xz"
    if not p.exists():
        raise FileNotFoundError(
            f"bundled ClinVar TSV not found at {p}; molamola does not "
            "auto-download. Provide one via --clinvar (TSV or VCF)."
        )
    return p


def _load_clinvar_lookup_tsv(
    path: Path,
    keys_of_interest: set[tuple[str, int, str, str]] | None,
) -> dict[tuple[str, int, str, str], str]:
    """Read molamola's reduced ClinVar TSV (already bucketed)."""
    out: dict[tuple[str, int, str, str], str] = {}
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").rstrip("\r").split("\t")
        idx = {n: i for i, n in enumerate(header)}
        required = ("chrom", "pos", "ref", "alt", "bucket")
        missing = [c for c in required if c not in idx]
        if missing:
            raise ValueError(
                f"reduced ClinVar TSV at {path} missing columns: {missing}"
            )
        for raw in fh:
            line = raw.rstrip("\n").rstrip("\r")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) <= idx["bucket"]:
                continue
            chrom = _normalize_chrom(cols[idx["chrom"]])
            try:
                pos = int(cols[idx["pos"]])
            except ValueError:
                continue
            ref, alt = cols[idx["ref"]], cols[idx["alt"]]
            key = (chrom, pos, ref, alt)
            if keys_of_interest is not None and key not in keys_of_interest:
                continue
            bucket = cols[idx["bucket"]]
            if bucket in CLNSIG_COLOR:
                out[key] = bucket
    return out


def _load_clinvar_lookup_vcf(
    path: Path,
    keys_of_interest: set[tuple[str, int, str, str]] | None,
) -> dict[tuple[str, int, str, str], str]:
    """Read NCBI's raw ClinVar VCF, computing buckets via canon_clnsig."""
    out: dict[tuple[str, int, str, str], str] = {}
    with open_text(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n").rstrip("\r")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            chrom = _normalize_chrom(fields[0])
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            ref, alt = fields[3], fields[4]
            if "," in alt:
                continue
            key = (chrom, pos, ref, alt)
            if keys_of_interest is not None and key not in keys_of_interest:
                continue
            info = parse_info(fields[7])
            raw_clnsig = info.get("CLNSIG")
            if not isinstance(raw_clnsig, str) or not raw_clnsig:
                continue
            bucket = canon_clnsig(raw_clnsig)
            if bucket in CLNSIG_COLOR:
                out[key] = bucket
    return out


def load_clinvar_lookup(
    path: Path,
    keys_of_interest: set[tuple[str, int, str, str]] | None = None,
) -> dict[tuple[str, int, str, str], str]:
    """Stream a ClinVar source, return ``{(chrom, pos, ref, alt): bucket}``.

    Bucket is one of ``"p_or_lp"``, ``"vus"``, ``"conflicting"``, or
    ``"benign"`` — the colour-mappable subset (keys of
    :data:`CLNSIG_COLOR`). Records that would canonicalise to
    ``"other"`` or have no CLNSIG are dropped — they would render
    in the no-ClinVar grey anyway.

    Format is auto-detected from the file extension:

    - ``.tsv`` / ``.tsv.gz`` / ``.tsv.xz``: molamola's reduced TSV
      with header ``chrom\\tpos\\tref\\talt\\tbucket``. Bucket is
      read verbatim — no re-canonicalisation. This is what
      :func:`find_clinvar_file` returns. Produced by
      ``scripts/derive_clinvar_for_molamola.py`` from the raw NCBI
      VCF.
    - ``.vcf`` / ``.vcf.gz``: NCBI's raw ClinVar VCF release. CLNSIG
      is bucketed at load time via :func:`canon_clnsig`. Useful for
      one-off overrides via ``--clinvar`` before the bundled TSV
      has been re-derived against a fresh release.

    When ``keys_of_interest`` is given, only matching rows are kept —
    this keeps the dict small even on the full ClinVar release
    (~4M sites). Without it, every mappable record is kept.

    Both source formats use bare contig names internally; this
    loader normalises to UCSC ``chrN`` so dict keys match the
    convention used elsewhere in molamola.
    """
    p = str(path)
    if (p.endswith(".tsv") or p.endswith(".tsv.gz")
            or p.endswith(".tsv.xz")):
        return _load_clinvar_lookup_tsv(path, keys_of_interest)
    return _load_clinvar_lookup_vcf(path, keys_of_interest)


def classify_pairs(
    variants: list[PhasedVariant],
) -> tuple[int, int]:
    """Count trans and cis pairs of canonical-transcript missense variants.

    Two variants are paired iff they are in the same phase set
    (``ps``) and both are missense canonical-transcript variants.
    Same ``variant_hap`` ⇒ cis; opposite ⇒ trans. Pairs are counted
    in the unordered ``C(n, 2)`` sense (each pair contributes once).

    Used to populate the panel title's ``T trans, C cis`` suffix.
    Pure counting helper; the renderer does NOT draw arcs (the
    locked spec intentionally omits them — opposite-hap dots in the
    same block already read as trans).
    """
    pool = [
        v for v in variants
        if v.is_canonical_transcript and is_missense_consequence(v.consequence)
    ]
    n_trans = 0
    n_cis = 0
    for i, a in enumerate(pool):
        for b in pool[i + 1:]:
            if a.ps != b.ps:
                continue
            if a.variant_hap == b.variant_hap:
                n_cis += 1
            else:
                n_trans += 1
    return n_trans, n_cis


#: ClinVar bucket set treated as "informative for compound-het":
#: pathogenic / likely-pathogenic and uncertain-significance carry
#: enough signal that an anchor in this set, paired with a non-benign
#: partner, is worth surfacing in the auto-select sweep. The
#: stricter "true compound-het" rule additionally requires the
#: partner to also be in this set (strict-section in the report).
_AUTO_SELECT_BUCKETS: frozenset[str] = frozenset({"p_or_lp", "vus"})


def _is_strict_pair(a: "PhasedVariant", b: "PhasedVariant") -> bool:
    """Both variants are P/LP or VUS — the "true compound-het" filter."""
    return (a.clnsig in _AUTO_SELECT_BUCKETS
            and b.clnsig in _AUTO_SELECT_BUCKETS)


def _is_extended_pair(a: "PhasedVariant", b: "PhasedVariant") -> bool:
    """At least one anchor is P/LP-or-VUS and partner is not benign.

    Extended pairs are a superset of strict pairs.
    """
    return (
        (a.clnsig in _AUTO_SELECT_BUCKETS and b.clnsig != "benign")
        or (b.clnsig in _AUTO_SELECT_BUCKETS and a.clnsig != "benign")
    )


def find_compound_het_candidates(
    variants: list[PhasedVariant],
    genes: dict[str, Gene],
    *,
    min_pair_count: int = 1,
) -> tuple[list[str], list[str]]:
    """Auto-select candidate genes, split by strictness.

    Returns ``(strict_genes, extended_only_genes)``: two disjoint,
    alphabetically-sorted lists. Their union is the full
    auto-select set; each gene appears in exactly one list.

    A *trans pair* is two phased het missense canonical-transcript
    variants in the same phase set on opposite haps. Two pair
    classifications matter:

    - **strict** — both variants have ClinVar CLNSIG in
      :data:`_AUTO_SELECT_BUCKETS` (``p_or_lp`` or ``vus``). True
      compound-het filter; rare in real data outside actual
      recessive-disease cases.
    - **extended** (a strict superset) — at least one variant is in
      ``_AUTO_SELECT_BUCKETS`` AND the partner is not benign.

    Selection:

    - Gene → ``strict_genes`` iff ``>= min_pair_count`` strict pairs
      in a single PS.
    - Gene → ``extended_only_genes`` iff ``>= min_pair_count``
      extended pairs in a single PS, AND the gene does not already
      qualify as strict.

    Pairs where both variants are conflicting / no-ClinVar / both
    benign are excluded from auto-select entirely. Use ``--gene``
    to plot those explicitly.

    "In trans" pair semantics: the two variants share PS but have
    opposite ``variant_hap``. :func:`classify_pairs` uses the same
    definition for the panel's `T trans, C cis` summary so the
    auto-select rule and the rendered panel agree on what "trans"
    means.
    """
    by_gene: dict[str, list[PhasedVariant]] = {}
    for v in variants:
        if v.gene_symbol is None or v.gene_symbol not in genes:
            continue
        if not v.is_canonical_transcript:
            continue
        if not is_missense_consequence(v.consequence):
            continue
        gene = genes[v.gene_symbol]
        if v.chrom != gene.chrom or not (gene.start <= v.pos <= gene.end):
            continue
        by_gene.setdefault(v.gene_symbol, []).append(v)

    strict: list[str] = []
    extended_only: list[str] = []
    for symbol, vs in by_gene.items():
        by_ps: dict[int, list[PhasedVariant]] = {}
        for v in vs:
            by_ps.setdefault(v.ps, []).append(v)
        n_strict = 0
        n_extended = 0
        for ps_vs in by_ps.values():
            for i, a in enumerate(ps_vs):
                for b in ps_vs[i + 1:]:
                    if a.variant_hap == b.variant_hap:
                        continue
                    if _is_strict_pair(a, b):
                        n_strict += 1
                        n_extended += 1
                    elif _is_extended_pair(a, b):
                        n_extended += 1
        if n_strict >= min_pair_count:
            strict.append(symbol)
        elif n_extended >= min_pair_count:
            extended_only.append(symbol)
    strict.sort()
    extended_only.sort()
    return strict, extended_only


# ---------------------------------------------------------------------------
# Focus / CLI parsing helpers
# ---------------------------------------------------------------------------

def parse_focus(s: str) -> tuple[str, int]:
    """Parse a ``CHR:POS`` or ``CHR:BAND`` string for ``--focus``.

    The second component can be either:
    - A genomic position (digits, commas/underscores allowed) → matched
      with `--focus-window` bp tolerance.
    - An ISCN cytoband name like ``q11.23`` or a prefix like ``q11``
      that matches every band on the chromosome whose name starts with
      that prefix → matched as a range with no window.

    Returns ``(chrom, key)`` where ``key`` is either an int (position)
    or a str (band name / prefix, leading ``p``/``q`` preserved).

    Examples
    --------
    >>> parse_focus("chr1:73129297")
    ('chr1', 73129297)
    >>> parse_focus("1:73,129,297")
    ('chr1', 73129297)
    >>> parse_focus("chr7:q11.23")
    ('chr7', 'q11.23')
    >>> parse_focus("chrX:p21")
    ('chrX', 'p21')

    Raises
    ------
    argparse.ArgumentTypeError
        If the string does not parse as ``CHR:POS`` or ``CHR:BAND``.
    """
    if ":" not in s:
        raise argparse.ArgumentTypeError(
            f"--focus must be CHR:POS or CHR:BAND, got {s!r}",
        )
    chrom, rest = s.split(":", 1)
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    rest = rest.strip()
    if not rest:
        raise argparse.ArgumentTypeError(f"empty position/band in {s!r}")
    # Cytoband names always start with p or q (followed by digits/dots).
    if rest[0] in ("p", "q"):
        return chrom, rest
    # Otherwise treat as a position.
    try:
        pos_i = int(rest.replace(",", "").replace("_", ""))
    except ValueError as e:
        raise argparse.ArgumentTypeError(
            f"invalid position/band in {s!r}: {e}",
        )
    return chrom, pos_i


def resolve_band_range(
    chrom: str,
    band_prefix: str,
    cytobands: dict,
) -> tuple[int, int] | None:
    """Resolve a cytoband name or prefix to a ``(start, end)`` range.

    Scans ``cytobands[chrom]`` for bands whose name starts with
    ``band_prefix`` and returns the union range ``(min start, max end)``.
    Returns ``None`` if no band matches (e.g. typo or wrong reference).
    """
    bands = cytobands.get(chrom, [])
    matches = [(s, e) for (s, e, name, _stain) in bands
               if name.startswith(band_prefix)]
    if not matches:
        return None
    return (min(s for s, _ in matches), max(e for _, e in matches))


def matches_focus(
    b: BND,
    focus: list,
    window: int,
    cytobands: dict | None = None,
) -> bool:
    """Return True iff either endpoint of *b* matches any focus entry.

    Position foci match if the endpoint is within ``window`` bp;
    band foci match if the endpoint falls inside the band's range.
    """
    for fc, key in focus:
        if isinstance(key, int):
            if b.chr1 == fc and abs(b.pos1 - key) <= window:
                return True
            if b.chr2 == fc and abs(b.pos2 - key) <= window:
                return True
        else:
            if cytobands is None:
                continue
            rng = resolve_band_range(fc, key, cytobands)
            if rng is None:
                continue
            start, end = rng
            if b.chr1 == fc and start <= b.pos1 < end:
                return True
            if b.chr2 == fc and start <= b.pos2 < end:
                return True
    return False


# ---------------------------------------------------------------------------
# Cytobands & noise
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Reference-data lookup (cytoband)
# ---------------------------------------------------------------------------

#: Supported reference assemblies. hg38 is the long-standing default;
#: t2t targets T2T-CHM13v2.0 (the with-chrY release).
SUPPORTED_REFERENCES: tuple[str, ...] = ("hg38", "t2t")


def detect_reference_hint(filename: str) -> str | None:
    """Guess the reference from a VCF filename.

    Returns ``"hg38"``, ``"t2t"``, or ``None`` if the filename has no
    recognisable hint. Case-insensitive substring match on tokens
    {hg38, grch38, t2t, chm13}.
    """
    lower = filename.lower()
    has_t2t = "t2t" in lower or "chm13" in lower
    has_hg38 = "hg38" in lower or "grch38" in lower
    if has_t2t and not has_hg38:
        return "t2t"
    if has_hg38 and not has_t2t:
        return "hg38"
    return None


def find_cytoband_file(reference: str = "hg38") -> Path:
    """Return the bundled cytoband file path for the given reference.

    Both supported references ship in the repo:
    ``data/cytoBand.txt.gz`` (hg38) and ``data/cytoBand.t2t.txt.gz``
    (T2T-CHM13v2.0). No external lookup, no auto-download.
    """
    name = "cytoBand.t2t.txt.gz" if reference == "t2t" else "cytoBand.txt.gz"
    bundled = Path(__file__).resolve().parent / "data" / name
    if not bundled.exists():
        raise FileNotFoundError(
            f"bundled cytoband missing: {bundled} (expected to ship "
            f"with molamola)",
        )
    return bundled


def find_mask_file(reference: str = "hg38") -> Path:
    """Return the bundled karyotype-mode exclusion mask path.

    Both supported references ship in the repo:
    ``data/exclusion.hg38.bed.gz`` and ``data/exclusion.t2t.bed.gz``.
    Used by karyotype mode to drop low-mappability bins from the CN
    scatter, smooth, and normalisation.
    """
    build = "t2t" if reference == "t2t" else "hg38"
    bundled = Path(__file__).resolve().parent / "data" / f"exclusion.{build}.bed.gz"
    if not bundled.exists():
        raise FileNotFoundError(
            f"bundled exclusion mask missing: {bundled} (expected to "
            f"ship with molamola)",
        )
    return bundled


def find_gc_file(reference: str = "hg38") -> Path:
    """Return the bundled karyotype-mode GC table path (10 kb bins).

    Both supported references ship in the repo:
    ``data/gc_10kb.hg38.bed.gz`` and ``data/gc_10kb.t2t.bed.gz``.
    Used by karyotype mode for per-1 % GC bucket median ratio
    correction.
    """
    build = "t2t" if reference == "t2t" else "hg38"
    bundled = Path(__file__).resolve().parent / "data" / f"gc_10kb.{build}.bed.gz"
    if not bundled.exists():
        raise FileNotFoundError(
            f"bundled GC table missing: {bundled} (expected to ship "
            f"with molamola)",
        )
    return bundled


# ---------------------------------------------------------------------------
# Karyotype-mode data pipeline (mosdepth -> CN -> scatter / smooth)
# ---------------------------------------------------------------------------

#: Modal-bin-frequency threshold below which read_mosdepth refuses
#: rather than aggregating a non-uniform-bins file into meaningless
#: output. mosdepth ``--by <int>`` is the supported shape.
_KARY_MIN_MODAL_BIN_FRACTION: float = 0.95


def read_mosdepth(path: Path) -> tuple[pd.DataFrame, int]:
    """Load a mosdepth ``regions.bed.gz`` and return ``(cov, bin_size)``.

    The input must come from a uniform-bin mosdepth run
    (``mosdepth --by <int>``). Karyotype mode refuses non-uniform
    inputs because its scatter aggregation and rolling-median smooth
    both assume a constant bin size; aggregating variable bins would
    yield meaningless plots.

    ``cov`` has columns ``chrom, start, end, depth`` with ``chrom``
    ordered as a pandas ``Categorical`` over :data:`CHROM_ORDER`.
    ``bin_size`` is the modal ``end - start`` across all
    canonical-chrom rows.

    Raises ``ValueError`` if fewer than 95 % of rows share the modal
    bin size, or zero rows survive after filtering to canonical
    chromosomes (chr1-22, chrX, chrY).
    """
    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["chrom", "start", "end", "depth"],
        dtype={"chrom": str, "start": np.int64, "end": np.int64,
               "depth": np.float32},
    )
    df = df[df["chrom"].isin(CHROM_SET)].copy()
    if len(df) == 0:
        raise ValueError(
            f"mosdepth file {path} has zero rows on canonical "
            f"chromosomes (chr1-22, chrX, chrY) -- was it called "
            f"against a build using '1'/'2'/'X' chrom naming?",
        )
    spans = (df["end"] - df["start"]).to_numpy()
    modal_bin = int(np.bincount(spans).argmax())
    modal_share = float((spans == modal_bin).mean())
    if modal_share < _KARY_MIN_MODAL_BIN_FRACTION:
        raise ValueError(
            f"mosdepth file {path} has non-uniform bins "
            f"(modal {modal_bin} bp covers only {modal_share:.1%} of "
            f"rows). Re-run mosdepth with --by <int> for karyotype mode.",
        )
    df["chrom"] = pd.Categorical(
        df["chrom"], categories=CHROM_ORDER, ordered=True,
    )
    df = df.sort_values(["chrom", "start"]).reset_index(drop=True)
    return df, modal_bin


def read_cytoband_df(path: Path) -> pd.DataFrame:
    """Load a UCSC cytoband file as a pandas ``DataFrame``.

    Distinct from :func:`load_cytobands` (which returns a per-chrom
    dict of tuples and is read by SV mode): karyotype mode prefers
    the flat tabular form for groupby-based queries. ``chrom`` is
    ordered as a ``Categorical`` over :data:`CHROM_ORDER`; an
    ``arm`` column is derived from the band name's first character
    (``p`` or ``q``).
    """
    cb = pd.read_csv(
        path, sep="\t", header=None,
        names=["chrom", "start", "end", "name", "stain"],
        dtype={"chrom": str, "start": np.int64, "end": np.int64,
               "name": str, "stain": str},
    )
    cb = cb[cb["chrom"].isin(CHROM_SET)].copy()
    cb["chrom"] = pd.Categorical(
        cb["chrom"], categories=CHROM_ORDER, ordered=True,
    )
    cb["arm"] = cb["name"].str[0]
    return cb


def chrom_lengths_from_cb(cb: pd.DataFrame) -> dict[str, int]:
    """Return ``{chrom: max(end)}`` from a cytoband DataFrame."""
    return (
        cb.groupby("chrom", observed=True)["end"]
        .max()
        .astype(int)
        .to_dict()
    )


def cum_offsets(lengths: dict[str, int]) -> dict[str, int]:
    """Per-chrom cumulative bp offset for genome-wide x-axis coordinates."""
    out: dict[str, int] = {}
    running = 0
    for c in CHROM_ORDER:
        out[c] = running
        running += lengths.get(c, 0)
    return out


def annotate_mask(cov: pd.DataFrame, mask_bed: Path, bin_size: int,
                  lengths: dict[str, int]) -> np.ndarray:
    """Per-cov-row bool: ``True`` if the bin overlaps the exclusion mask.

    Aligned 1:1 with ``cov`` rows. Downstream code typically uses
    ``mask_pass = ~result`` as the keep-mask.
    """
    masked = {
        c: np.zeros((lengths[c] + bin_size - 1) // bin_size, dtype=bool)
        for c in CHROM_ORDER if c in lengths
    }
    with open_text(mask_bed) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            if chrom not in masked:
                continue
            s, e = int(parts[1]), int(parts[2])
            b0 = s // bin_size
            b1 = (e + bin_size - 1) // bin_size
            masked[chrom][b0:b1] = True
    out = np.zeros(len(cov), dtype=bool)
    for chrom, idx in cov.groupby("chrom", observed=True).indices.items():
        if chrom not in masked:
            continue
        bins = (cov["start"].iloc[idx].to_numpy() // bin_size).astype(np.int64)
        valid = bins < len(masked[chrom])
        if valid.any():
            out[idx[valid]] = masked[chrom][bins[valid]]
    return out


def annotate_centromere_pad(cov: pd.DataFrame, cb: pd.DataFrame,
                            lengths: dict[str, int],
                            pad_kb: float) -> np.ndarray:
    """Per-cov-row bool: bins within ``pad_kb`` of an ``acen`` band.

    Peri-centromeric coverage is inflated by satellite-derived
    multi-mappers; many of those bins still pass the short-read
    accessibility mask, and when the rolling-median window crosses a
    centromere the surviving flanking bins dominate the median and
    produce a spurious smooth-line spike.
    """
    if pad_kb <= 0:
        return np.zeros(len(cov), dtype=bool)
    pad_bp = int(pad_kb * 1000)
    out = np.zeros(len(cov), dtype=bool)
    acen = cb[cb["stain"] == "acen"]
    for chrom, sub in acen.groupby("chrom", observed=True, sort=False):
        if chrom not in lengths:
            continue
        cen_start = max(0, int(sub["start"].min()) - pad_bp)
        cen_end = min(lengths[chrom], int(sub["end"].max()) + pad_bp)
        idx = ((cov["chrom"] == chrom)
               & (cov["start"] < cen_end)
               & (cov["end"] > cen_start)).to_numpy()
        out |= idx
    return out


def annotate_gc(cov: pd.DataFrame, gc_bed: Path,
                gc_bin_size: int = 10_000) -> np.ndarray:
    """Per-cov-row GC% (int 0-100; :data:`GC_MISSING` if unknown).

    Each cov row gets the GC% of the ``gc_bin_size``-bp window
    containing its midpoint. Rows whose midpoint falls outside every
    GC bin (e.g. cov rows past the GC table's chromosome length)
    return :data:`GC_MISSING`.
    """
    gc = pd.read_csv(
        gc_bed, sep="\t", header=None,
        names=["chrom", "start", "end", "gc"],
        dtype={"chrom": str, "start": np.int64, "end": np.int64,
               "gc": np.int16},
    )
    gc = gc[gc["chrom"].isin(CHROM_SET)]
    cov_mid = ((cov["start"] + cov["end"]) // 2).to_numpy()
    cov_bin_start = (cov_mid // gc_bin_size) * gc_bin_size
    lookup = pd.DataFrame({
        "chrom": cov["chrom"].astype(str).to_numpy(),
        "start": cov_bin_start,
    })
    merged = lookup.merge(
        gc[["chrom", "start", "gc"]], on=["chrom", "start"], how="left",
    )
    return merged["gc"].fillna(GC_MISSING).astype(np.int16).to_numpy()


def fit_gc_correction(depth: np.ndarray, gc_pct: np.ndarray,
                      mask_pass: np.ndarray, chroms: pd.Series,
                      min_bin_count: int = 100) -> dict[int, float]:
    """Per-1 % GC bucket median ratio for depth correction.

    Restricted to autosomal, non-masked, GC-defined bins. Buckets
    with fewer than ``min_bin_count`` supporting bins or a
    non-positive median get factor 1.0. Otherwise the factor is
    ``global_autosomal_median / bucket_median`` -- multiplying raw
    depth by ``factor`` flattens the GC bias.
    """
    use = (
        mask_pass
        & chroms.isin(KARY_AUTOSOMES).to_numpy()
        & (gc_pct != GC_MISSING)
    )
    if not use.any():
        return {}
    df = pd.DataFrame({"gc": gc_pct[use].astype(int), "depth": depth[use]})
    grouped = df.groupby("gc")["depth"].agg(["median", "count"])
    global_med = float(np.median(df["depth"]))
    factors: dict[int, float] = {}
    for gc_val, row in grouped.iterrows():
        if row["count"] < min_bin_count or row["median"] <= 0:
            factors[int(gc_val)] = 1.0
        else:
            factors[int(gc_val)] = global_med / float(row["median"])
    return factors


def apply_gc_correction(depth: np.ndarray, gc_pct: np.ndarray,
                        factors: dict[int, float]) -> np.ndarray:
    """Apply :func:`fit_gc_correction` factors via a 256-entry lookup."""
    if not factors:
        return depth.astype(np.float64)
    lookup = np.ones(256, dtype=np.float64)
    for k, v in factors.items():
        if 0 <= k < 256:
            lookup[k] = v
    return depth.astype(np.float64) * lookup[gc_pct.astype(np.int64).clip(0, 255)]


def detect_sex_from_cn(cn: np.ndarray, chroms: pd.Series,
                       mask_pass: np.ndarray) -> str:
    """Call ``"male"`` if median chrY CN > 0.3; else ``"female"``.

    The threshold is a heuristic that works for both hg38 (mappable
    chrY fraction is small) and T2T-CHM13v2 (chrY fully resolved).
    XYY / XXY samples will still be called male, which is fine for
    the expected-CN-line purpose; HTML metadata surfaces this as
    "inferred genomic sex" rather than a clinical sex call.
    """
    is_y = (chroms == "chrY").to_numpy()
    use = is_y & mask_pass & np.isfinite(cn)
    if not use.any():
        return "female"
    return "male" if float(np.median(cn[use])) > 0.3 else "female"


def expected_copy_number(chrom: str, sex: str) -> float:
    """Reference copy number for ``chrom`` under the given ``sex``."""
    if chrom in KARY_AUTOSOMES:
        return 2.0
    if chrom == "chrX":
        return 1.0 if sex == "male" else 2.0
    if chrom == "chrY":
        return 1.0 if sex == "male" else 0.0
    return float("nan")


def rolling_median_per_chrom(df: pd.DataFrame, bin_size: int,
                             window_mb: float,
                             mask_pass: np.ndarray) -> pd.DataFrame:
    """Per-chrom rolling median of ``df['cn']``, ignoring masked bins.

    Requires at least half the window's worth of non-masked support
    to emit a value -- otherwise the smooth would jitter wildly
    across the sliver of unmasked bins next to a big masked stretch.

    Returns a copy of ``df`` with a new ``smooth`` column.
    """
    win = max(1, int(round(window_mb * 1e6 / bin_size)))
    df = df.copy()
    df["_in"] = np.where(mask_pass, df["cn"].to_numpy(), np.nan)
    parts = []
    for _, sub in df.groupby("chrom", observed=True, sort=False):
        sub = sub.copy()
        sub["smooth"] = (
            sub["_in"]
            .rolling(win, center=True, min_periods=max(1, win // 2))
            .median()
        )
        parts.append(sub)
    out = pd.concat(parts, ignore_index=True)
    return out.drop(columns=["_in"])


def aggregate_for_scatter(cov: pd.DataFrame, factor: int) -> pd.DataFrame:
    """Per-chrom median-aggregate every ``factor`` consecutive bins.

    Each output row is the median of ``factor`` non-masked bins.
    Median is robust to outlier bins (e.g. residual segdup
    contamination that survives the mask). Fully-masked groups are
    dropped. ``factor <= 1`` returns the non-masked rows unchanged.
    """
    if factor <= 1:
        out = cov[cov["mask_pass"]].copy()
        return out[["chrom", "xpos", "cn", "mask_pass"]].reset_index(drop=True)

    parts = []
    for chrom, sub in cov.groupby("chrom", observed=True, sort=False):
        sub = sub.reset_index(drop=True)
        n = len(sub)
        if n == 0:
            continue
        groups = np.arange(n) // factor
        cn_masked = np.where(sub["mask_pass"], sub["cn"], np.nan)
        df = pd.DataFrame({
            "chrom": [chrom] * n,
            "xpos": sub["xpos"].to_numpy(),
            "_cn_masked": cn_masked,
            "_grp": groups,
        })
        agg = df.groupby("_grp", sort=False, observed=True).agg(
            chrom=("chrom", "first"),
            xpos=("xpos", "mean"),
            cn=("_cn_masked", "median"),
        )
        agg = agg.dropna(subset=["cn"]).reset_index(drop=True)
        agg["mask_pass"] = True
        parts.append(agg)
    if not parts:
        return cov.iloc[0:0][["chrom", "xpos", "cn", "mask_pass"]].copy()
    return pd.concat(parts, ignore_index=True)


def downsample_systematic(df: pd.DataFrame, max_points: int) -> pd.DataFrame:
    """Keep every Nth row to cap a DataFrame at ``max_points`` rows.

    Used for BAF, where averaging would erase the very 0.33 / 0.67 /
    0 / 1 deviations a karyotype reviewer is looking for.
    """
    if len(df) <= max_points:
        return df
    step = int(np.ceil(len(df) / max_points))
    return df.iloc[::step].copy()


# ---------------------------------------------------------------------------
# Karyotype-mode plotting primitives
# ---------------------------------------------------------------------------
#
# All styling (palette, fonts, line widths) is passed per-call. Do NOT
# call ``plt.rcParams.update()`` from any function in this section:
# rcParams is global module state and any leakage would silently shift
# the SV-mode and compound-het renders that share the same Python
# process.

_KARY_DARK_STAINS: frozenset[str] = frozenset({"gpos75", "gpos100", "acen"})


def _kary_rounded_pill_path(x: float, y: float, w: float, h: float,
                            rx: float) -> MplPath:
    """Path for a rectangle with circular caps at the left and right ends.

    ``rx`` is the horizontal radius in data coords; the vertical
    radius is half the height (so the caps are full half-circles in
    the axes aspect, scaled into ellipses by the data transform).
    """
    ry = h / 2
    cy = y + ry
    verts = [
        (x + rx, y),
        (x + w - rx, y),
        (x + w, y),
        (x + w, cy),
        (x + w, y + h),
        (x + w - rx, y + h),
        (x + rx, y + h),
        (x, y + h),
        (x, cy),
        (x, y),
        (x + rx, y),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.LINETO,
        MplPath.CURVE3, MplPath.CURVE3,
        MplPath.CURVE3, MplPath.CURVE3,
        MplPath.LINETO,
        MplPath.CURVE3, MplPath.CURVE3,
        MplPath.CURVE3, MplPath.CURVE3,
    ]
    return MplPath(verts, codes)


def _draw_kary_cytoband_strip(ax, bands: pd.DataFrame,
                              label_min_mb: float = 5.0) -> None:
    """Karyotype-mode cytoband strip as a rounded-pill silhouette.

    Bands are clipped to the silhouette so the leftmost and rightmost
    stains take the rounded telomere shape. Band names are rendered
    on bands wider than ``label_min_mb`` Mb; ``acen`` / dark stains
    get a paper-colour label, light stains get an INK_2 label.
    """
    if bands.empty:
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.grid(False)
        for s in ("right", "top", "left", "bottom"):
            ax.spines[s].set_visible(False)
        return

    chrom_start = int(bands["start"].min())
    chrom_end = int(bands["end"].max())
    chrom_w = chrom_end - chrom_start
    rx = chrom_w * 0.008
    y0, y1 = 0.15, 0.85
    h = y1 - y0

    silhouette = _kary_rounded_pill_path(chrom_start, y0, chrom_w, h, rx)
    clip_patch = mpatches.PathPatch(
        silhouette, facecolor="none", edgecolor="none",
        transform=ax.transData,
    )
    ax.add_patch(clip_patch)

    for _, row in bands.iterrows():
        rect = mpatches.Rectangle(
            (row["start"], y0), row["end"] - row["start"], h,
            facecolor=KARY_CYTO_STAIN_COLOR.get(row["stain"], "#dddddd"),
            edgecolor="none",
        )
        ax.add_patch(rect)
        rect.set_clip_path(clip_patch)

    border = mpatches.PathPatch(
        silhouette, facecolor="none", edgecolor=KARY_INK_2,
        linewidth=0.7, transform=ax.transData,
    )
    ax.add_patch(border)

    for _, row in bands.iterrows():
        if (row["end"] - row["start"]) > label_min_mb * 1e6:
            text_color = (KARY_PAPER if row["stain"] in _KARY_DARK_STAINS
                          else KARY_INK_2)
            ax.text(
                (row["start"] + row["end"]) / 2, 0.5, row["name"],
                ha="center", va="center", fontsize=9,
                fontfamily=list(KARY_FONT_SANS), color=text_color,
            )

    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.grid(False)
    for s in ("right", "top", "left", "bottom"):
        ax.spines[s].set_visible(False)


def _draw_kary_centromere_ticks(ax, cb: pd.DataFrame,
                                offsets: dict[str, int]) -> None:
    """Small INK_3 tick at each chromosome's p/q boundary on the top edge.

    Acts as an orientation cue inside the genome-wide coverage panel.
    """
    for chrom, sub in cb.groupby("chrom", observed=True, sort=False):
        sub = sub.sort_values("start")
        p_bands = sub[sub["arm"] == "p"]
        if p_bands.empty or chrom not in offsets:
            continue
        cen = offsets[chrom] + int(p_bands["end"].max())
        ax.plot(
            [cen, cen], [0.985, 1.015],
            transform=ax.get_xaxis_transform(),
            color=KARY_INK_3, linewidth=0.5, clip_on=False, zorder=5,
        )


def _kary_build_arm_ticks(cb: pd.DataFrame,
                          offsets: dict[str, int]) -> tuple[list[float], list[str]]:
    """Tick (xs, labels) for the chromosome-arm axis above the coverage panel.

    Acrocentric p-arms (chr13/14/15/21/22, chrY) are dropped from the
    label list — in matplotlib's real font metrics their p-label
    crashes into the q-label of the same chromosome.
    """
    xs: list[float] = []
    labels: list[str] = []
    for chrom, sub in cb.groupby("chrom", observed=True, sort=False):
        sub = sub.sort_values("start")
        if chrom not in offsets:
            continue
        short = chrom.replace("chr", "")
        for arm in ("p", "q"):
            if arm == "p" and chrom in KARY_ACROCENTRIC_FOR_LABELS:
                continue
            bands = sub[sub["arm"] == arm]
            if bands.empty:
                continue
            centre = offsets[chrom] + (int(bands["start"].min())
                                       + int(bands["end"].max())) / 2
            xs.append(centre)
            labels.append(f"{short}{arm}")
    return xs, labels


def _kary_plot_coverage(ax, scatter_df: pd.DataFrame,
                        smooth_df: pd.DataFrame,
                        expected_lines: list[tuple[float, float, float]],
                        ymax: float, x_col: str = "xpos") -> None:
    """CN scatter + per-chrom rolling-median smooth + expected-CN dashes."""
    sc_p = scatter_df[scatter_df["mask_pass"]]
    ax.scatter(
        sc_p[x_col], sc_p["cn"], s=2.4, alpha=0.32,
        c=KARY_SCATTER, linewidths=0, rasterized=True,
    )
    for _chrom, sub in smooth_df.groupby("chrom", observed=True, sort=False):
        ax.plot(
            sub[x_col], sub["smooth"], color=KARY_ROSE,
            lw=2.0, solid_capstyle="round",
        )
    for x0, x1, y in expected_lines:
        ax.plot(
            [x0, x1], [y, y], color=KARY_OXFORD,
            lw=0.9, alpha=0.85, dashes=(4, 2),
        )
    ax.set_ylim(0, ymax)
    ax.set_ylabel("CN")


def _kary_plot_baf(ax, baf_df: pd.DataFrame, x_col: str = "xpos") -> None:
    """BAF scatter with reference grid lines at 0.25 / 0.5 / 0.75."""
    ax.scatter(
        baf_df[x_col], baf_df["baf"], s=2.0, alpha=0.30,
        c=KARY_SCATTER, linewidths=0, rasterized=True,
    )
    for y in (0.25, 0.5, 0.75):
        ax.axhline(y, color=KARY_OXFORD, lw=0.6, dashes=(3, 2), alpha=0.55)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel("BAF")


def _kary_apply_tabular_numerics(*axes) -> None:
    """Switch tick label font to the karyotype mono stack.

    Prevents digit-width jitter between ticks like ``0.5`` and ``0.50``
    on the BAF panel.
    """
    mono = list(KARY_FONT_MONO)
    for ax in axes:
        for tl in ax.get_xticklabels() + ax.get_yticklabels():
            tl.set_fontfamily(mono)


def _kary_align_panel_ylabels(*axes) -> None:
    """Pin each panel's y-axis label to ``KARY_YLABEL_X``.

    Otherwise the BAF panel's 4-char ticks (``0.50``) push its
    y-label further left than the CN panel's; with this pinning,
    ``CN`` and ``BAF`` align vertically.
    """
    for ax in axes:
        ax.yaxis.set_label_coords(KARY_YLABEL_X, 0.5)


def _kary_format_bin_size(bp: float) -> str:
    """Human-readable bin size (Mb / kb / bp) for run-metadata strips."""
    if bp >= 1e6:
        return f"{bp / 1e6:.1f} Mb"
    if bp >= 1e3:
        return f"{bp / 1e3:.0f} kb"
    return f"{bp:.0f} bp"


def _karyotype_meta_chips(args: argparse.Namespace, scatter_bin_label: str,
                          sex: str) -> list[str]:
    """Compose the karyotype run-metadata strip as a list of chips.

    Used by both the in-figure metadata line and the HTML report's
    run-metadata section so the same provenance shows up in both
    places.
    """
    bits: list[str] = [
        args.mosdepth.name,
        args.reference,
        f"sex={sex}",
        f"bin={scatter_bin_label}",
        f"smooth={args.smooth_window_mb} Mb",
    ]
    if getattr(args, "no_mask", False):
        bits.append("mask=off")
    if getattr(args, "no_gc", False):
        bits.append("GC=off")
    return bits


def _kary_attach_xpos(df: pd.DataFrame, pos_col: str,
                      offsets: dict[str, int]) -> pd.DataFrame:
    """Add a genome-wide ``xpos`` column = ``pos_col + offset[chrom]``."""
    df = df.copy()
    df["xpos"] = (
        df[pos_col].astype(np.int64)
        + df["chrom"].map(offsets).astype(np.int64)
    )
    return df


def render_karyotype_genome_png(
    cov: pd.DataFrame, cb: pd.DataFrame, lengths: dict[str, int],
    baf_df: pd.DataFrame | None, sex: str, bin_size: int,
    args: argparse.Namespace,
) -> tuple[bytes, str]:
    """Render the genome-wide karyotype figure to PNG bytes.

    ``cov`` must already carry ``cn`` (autosomal-median-normalised
    copy number), ``mask_pass`` (``True`` to keep), and ``smooth``
    (rolling-median per chrom) columns; :func:`karyotype_main`
    prepares those before calling.

    Returns ``(png_bytes, scatter_bin_label)``. The scatter-bin label
    (e.g. ``"50 kb"``) is also stamped onto the figure's metadata
    strip — returning it lets the HTML report reuse the same string
    without re-deriving the aggregation factor.
    """
    offsets = cum_offsets(lengths)
    cov_xy = _kary_attach_xpos(cov, "start", offsets)

    factor = max(1, int(round(args.scatter_bin_kb * 1000 / bin_size)))
    cap_factor = max(1, int(np.ceil(len(cov_xy) / args.max_points)))
    factor = max(factor, cap_factor)
    scatter = aggregate_for_scatter(cov_xy, factor)
    scatter_bin_label = _kary_format_bin_size(factor * bin_size)
    print(
        f"[info] coverage scatter: {len(scatter):,} points "
        f"(median of {factor} bins, ~{scatter_bin_label})",
    )

    expected_lines: list[tuple[float, float, float]] = []
    for chrom in CHROM_ORDER:
        chrom_len = lengths.get(chrom, 0)
        if chrom_len == 0:
            continue
        x0 = float(offsets[chrom])
        expected_lines.append(
            (x0, x0 + float(chrom_len), expected_copy_number(chrom, sex)),
        )

    baf_plot: pd.DataFrame | None = None
    if baf_df is not None and len(baf_df) > 0:
        baf_xy = _kary_attach_xpos(baf_df, "pos", offsets)
        baf_plot = downsample_systematic(baf_xy, args.max_baf_points)
        print(
            f"[info] BAF scatter: {len(baf_plot):,} points "
            f"(downsampled from {len(baf_df):,})",
        )

    if baf_plot is not None and len(baf_plot) > 0:
        fig, (ax_cov, ax_baf) = plt.subplots(
            2, 1,
            figsize=(KARY_FIG_W, KARY_FIG_H_GENOME_BAF),
            gridspec_kw={"height_ratios": list(KARY_HEIGHT_RATIOS_CN_BAF)},
            sharex=True,
        )
    else:
        fig, ax_cov = plt.subplots(
            figsize=(KARY_FIG_W, KARY_FIG_H_GENOME_ONLY),
        )
        ax_baf = None
    fig.subplots_adjust(
        left=0.045, right=0.996, top=0.880, bottom=0.085, hspace=0.06,
    )
    fig.patch.set_alpha(0)

    _kary_plot_coverage(
        ax_cov, scatter, cov_xy, expected_lines, args.ymax,
    )

    chrom_centres: list[float] = []
    chrom_labels: list[str] = []
    for chrom in CHROM_ORDER:
        chrom_len = lengths.get(chrom, 0)
        if chrom_len == 0:
            continue
        for ax in (ax_cov, ax_baf) if ax_baf else (ax_cov,):
            ax.axvline(offsets[chrom], color=KARY_RULE, lw=0.5, zorder=0.5)
        chrom_centres.append(offsets[chrom] + chrom_len / 2)
        chrom_labels.append(chrom.replace("chr", ""))
    total = offsets[CHROM_ORDER[-1]] + lengths.get(CHROM_ORDER[-1], 0)
    for ax in (ax_cov, ax_baf) if ax_baf else (ax_cov,):
        ax.axvline(total, color=KARY_RULE, lw=0.5, zorder=0.5)
        ax.set_xlim(0, total)

    if ax_baf is not None:
        _kary_plot_baf(ax_baf, baf_plot)
        ax_cov.tick_params(labelbottom=False)

    bottom_ax = ax_baf if ax_baf else ax_cov
    bottom_ax.set_xticks(chrom_centres)
    bottom_ax.set_xticklabels(
        chrom_labels, fontfamily=list(KARY_FONT_MONO),
        color=KARY_INK, fontsize=11.5, fontweight="medium",
    )

    arm_xs, arm_labels = _kary_build_arm_ticks(cb, offsets)
    sec = ax_cov.secondary_xaxis("top")
    sec.set_xticks(arm_xs)
    sec.set_xticklabels(
        arm_labels, fontfamily=list(KARY_FONT_SANS),
        color=KARY_INK_2, fontsize=9.5, rotation=0,
    )
    sec.tick_params(length=0, pad=2)
    for spine in sec.spines.values():
        spine.set_visible(False)

    _draw_kary_centromere_ticks(ax_cov, cb, offsets)

    fig.text(
        0.045, 0.965, "Genome-wide coverage",
        fontfamily=list(KARY_FONT_SANS), fontsize=14, color=KARY_INK,
        ha="left", va="center",
    )
    fig.text(
        0.045, 0.935,
        "  ·  ".join(_karyotype_meta_chips(args, scatter_bin_label, sex)),
        fontfamily=list(KARY_FONT_MONO), fontsize=10.5, color=KARY_INK_2,
        ha="left", va="center",
    )

    if ax_baf is not None:
        _kary_apply_tabular_numerics(ax_cov, ax_baf)
        _kary_align_panel_ylabels(ax_cov, ax_baf)
    else:
        _kary_apply_tabular_numerics(ax_cov)
        _kary_align_panel_ylabels(ax_cov)

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=200)
    plt.close(fig)
    return buf.getvalue(), scatter_bin_label


def _render_kary_region_panels(
    cov: pd.DataFrame, cb: pd.DataFrame, sex: str, bin_size: int,
    chrom: str, start: int, end: int,
    ax_band, ax_cov, args: argparse.Namespace,
    show_xlabel: bool = True,
) -> str:
    """Draw the cytoband strip + CN panel for one chrom region into the
    supplied axes. Returns the scatter-bin label string used."""
    sub_cov = cov[
        (cov["chrom"] == chrom)
        & (cov["start"] < end)
        & (cov["end"] > start)
    ].copy()
    sub_bands = cb[
        (cb["chrom"] == chrom)
        & (cb["start"] < end)
        & (cb["end"] > start)
    ].copy()

    factor = max(1, int(round(args.scatter_bin_kb * 1000 / bin_size)))
    cap_factor = max(1, int(np.ceil(max(len(sub_cov), 1) / args.max_points)))
    factor = max(factor, cap_factor)
    sub_cov_xy = sub_cov.copy()
    sub_cov_xy["xpos"] = sub_cov_xy["start"]
    scatter = aggregate_for_scatter(sub_cov_xy, factor)
    scatter["xpos_local"] = scatter["xpos"]

    expected_lines = [(float(start), float(end),
                       expected_copy_number(chrom, sex))]

    if ax_band is not None:
        _draw_kary_cytoband_strip(ax_band, sub_bands)

    sub_cov["xpos_local"] = sub_cov["start"]
    _kary_plot_coverage(
        ax_cov, scatter, sub_cov, expected_lines, args.ymax,
        x_col="xpos_local",
    )

    ax_cov.set_xlim(start, end)
    if show_xlabel:
        ax_cov.set_xlabel(f"{chrom} position (Mb)")
    ax_cov.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: f"{x / 1e6:.1f}"),
    )

    _kary_apply_tabular_numerics(ax_cov)
    _kary_align_panel_ylabels(ax_cov)

    return _kary_format_bin_size(factor * bin_size)


#: Hardcoded per-chrom layout: 3 columns × 8 rows = 24 cells matches
#: chr1..22 + chrX + chrY. Locked at port time (handoff decision
#: 2026-05-11); no CLI knob.
_KARY_PER_CHROM_COLS: int = 3
_KARY_PER_CHROM_ROWS: int = 8


def render_karyotype_per_chrom_png(
    cov: pd.DataFrame, cb: pd.DataFrame, lengths: dict[str, int],
    sex: str, bin_size: int, args: argparse.Namespace,
) -> tuple[bytes, str]:
    """Render the per-chromosome 3 × 8 A4-portrait karyotype grid to PNG bytes.

    ``cov`` must already carry the ``cn``, ``mask_pass``, ``smooth``
    columns prepared by :func:`karyotype_main`. BAF is intentionally
    omitted from this view -- per-chrom panels are small and CN is
    the primary karyotype signal; the genome-wide panel covers BAF.

    Returns ``(png_bytes, scatter_bin_label)``. The scatter-bin label
    is the one used by the last chromosome rendered; it is included
    in the figure's metadata strip and returned for HTML reuse.
    """
    chroms = [c for c in CHROM_ORDER if lengths.get(c, 0) > 0]
    if not chroms:
        raise ValueError("no chromosomes with non-zero length in cytoband")

    cols = _KARY_PER_CHROM_COLS
    rows = _KARY_PER_CHROM_ROWS
    fig_w, fig_h = 8.27, 11.69  # A4 portrait, inches
    header_h = 0.55
    gs_top = 1.0 - header_h / fig_h
    gs_left, gs_right = 0.060, 0.985
    gs_bottom = 0.030
    hspace = 0.60
    wspace = 0.20
    title_y = 1.0 - 0.18 / fig_h
    meta_y = 1.0 - 0.40 / fig_h
    inner_hspace = 0.05

    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.patch.set_alpha(0)
    outer = fig.add_gridspec(
        rows, cols,
        top=gs_top - 0.005, bottom=gs_bottom,
        left=gs_left, right=gs_right,
        hspace=hspace, wspace=wspace,
    )

    last_scatter_label = "?"
    for i, chrom in enumerate(chroms):
        r, c = i // cols, i % cols
        if r >= rows:
            break
        inner = outer[r, c].subgridspec(
            len(KARY_HEIGHT_RATIOS_BAND_CN), 1,
            height_ratios=list(KARY_HEIGHT_RATIOS_BAND_CN),
            hspace=inner_hspace,
        )
        ax_band = fig.add_subplot(inner[0])
        ax_cov = fig.add_subplot(inner[1], sharex=ax_band)
        chrom_len = lengths[chrom]
        last_scatter_label = _render_kary_region_panels(
            cov, cb, sex, bin_size, chrom, 0, chrom_len,
            ax_band, ax_cov, args, show_xlabel=True,
        )
        ax_band.set_title(
            chrom, fontsize=8, fontweight=500, color=KARY_INK,
            loc="left", pad=2,
        )
        if c != 0:
            ax_cov.set_ylabel("")
            ax_cov.tick_params(labelleft=False)

    fig.text(
        gs_left, title_y, "Per-chromosome coverage",
        fontfamily=list(KARY_FONT_SANS), fontsize=11.0, color=KARY_INK,
        ha="left", va="top",
    )
    fig.text(
        gs_left, meta_y,
        "  ·  ".join(_karyotype_meta_chips(args, last_scatter_label, sex)),
        fontfamily=list(KARY_FONT_MONO), fontsize=7.0, color=KARY_INK_2,
        ha="left", va="top",
    )

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=200)
    plt.close(fig)
    return buf.getvalue(), last_scatter_label


def load_cytobands(path: Path) -> dict[str, list[tuple[int, int, str, str]]]:
    """Load a UCSC ``cytoBand.txt(.gz)`` file.

    Parameters
    ----------
    path : Path
        Cytoband file (5 columns: chrom, start, end, name, gieStain).

    Returns
    -------
    dict[str, list[tuple[int, int, str, str]]]
        ``{chrom: [(start, end, band_name, stain), ...]}``. Only
        chromosomes in :data:`CHROM_SET` are kept.
    """
    out: dict = {}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, start, end, name, stain = f[:5]
            if chrom not in CHROM_SET:
                continue
            out.setdefault(chrom, []).append((int(start), int(end), name, stain))
    return out


def acrocentric_p_arm_ends(cytobands: dict) -> dict[str, int]:
    """Return ``{chrom: p/q-boundary}`` for the five acrocentric chromosomes.

    The p/q boundary is taken as the start coordinate of the first
    ``acen`` band on each chromosome - i.e. the p-arm spans
    ``[0, that_position)``. Acrocentric chromosomes whose first
    ``acen`` band is missing from the cytoband file are omitted.
    """
    out: dict[str, int] = {}
    for chrom in ACROCENTRIC_CHROMS:
        for (start, _end, _name, stain) in cytobands.get(chrom, []):
            if stain == "acen":
                out[chrom] = start
                break
    return out


def in_acrocentric_p(chrom: str, pos: int, p_arm_ends: dict[str, int]) -> bool:
    """True iff ``(chrom, pos)`` falls in an acrocentric short arm."""
    end = p_arm_ends.get(chrom)
    return end is not None and pos < end


def resolve_cytoband(chrom: str, pos: int, cytobands: dict) -> str:
    """Return the cytoband name (e.g. ``"q11.23"``) at ``chrom:pos``.

    Returns ``""`` if the position falls outside every band on
    ``chrom`` or if ``chrom`` isn't in the cytoband file.
    """
    for (start, end, name, _stain) in cytobands.get(chrom, []):
        if start <= pos < end:
            return name
    return ""


def iscn_label(event, cytobands: dict) -> str:
    """Return an ISCN-style nomenclature string for an SV or BND event.

    Forms produced (without the leading karyotype prefix):

    - ``BND``: ``t(7;17)(q11.23;q12)`` - canonical chromosome ordering.
    - ``DEL`` / ``DUP`` / ``INV``: ``del(7)(q11.23q12.2)`` (start band
      to end band; collapsed to a single band when both ends are in the
      same band).
    - ``INS``: ``ins(7)(q11.23)`` - non-strict, since ISCN's full
      insertion form requires the donor coords; a point-form is used
      for completeness.

    Bands that can't be resolved render as ``"?"`` so the output stays
    syntactically intact.
    """
    if isinstance(event, BND):
        i1 = CHROM_ORDER.index(event.chr1) if event.chr1 in CHROM_SET else 99
        i2 = CHROM_ORDER.index(event.chr2) if event.chr2 in CHROM_SET else 99
        if (i1, event.pos1) <= (i2, event.pos2):
            c1, p1, c2, p2 = event.chr1, event.pos1, event.chr2, event.pos2
        else:
            c1, p1, c2, p2 = event.chr2, event.pos2, event.chr1, event.pos1
        n1 = c1.replace("chr", "")
        n2 = c2.replace("chr", "")
        b1 = resolve_cytoband(c1, p1, cytobands) or "?"
        b2 = resolve_cytoband(c2, p2, cytobands) or "?"
        return f"t({n1};{n2})({b1};{b2})"

    n = event.chrom.replace("chr", "")
    b1 = resolve_cytoband(event.chrom, event.start, cytobands) or "?"
    if event.svtype == "INS":
        return f"ins({n})({b1})"
    b2 = resolve_cytoband(event.chrom, event.end, cytobands) or "?"
    bands = b1 if b1 == b2 else b1 + b2
    if event.svtype == "DEL":
        return f"del({n})({bands})"
    if event.svtype == "DUP":
        return f"dup({n})({bands})"
    if event.svtype == "INV":
        return f"inv({n})({bands})"
    return ""


def _parse_cov_ratio(s: str) -> float | str:
    """Accept either a positive number (int or float) or the literal ``"auto"``."""
    if s.lower() == "auto":
        return "auto"
    try:
        v = float(s)
    except ValueError as e:
        raise argparse.ArgumentTypeError(
            f"--cov-ratio must be a number (int or float) or 'auto', "
            f"got {s!r}: {e}",
        )
    if v <= 0:
        raise argparse.ArgumentTypeError(
            f"--cov-ratio must be > 0, got {v}",
        )
    return v


def auto_cov_ratio_threshold(
    bnds: list[BND],
    svs: list[SV],
    median_cov: float,
    *,
    floor: float = 2.0,
    quantile: float = 0.99,
) -> tuple[float, int]:
    """Compute an empirical cov-anomaly threshold from in-memory events.

    Returns ``(threshold, n_events_used)``. The threshold is
    ``max(floor, quantile of max_coverage / median_coverage)`` over all
    PASS events with non-zero max coverage. The floor prevents a
    sample with no real outliers from setting an unreasonably loose
    cutoff.
    """
    if median_cov <= 0:
        return floor, 0
    ratios: list[float] = []
    for b in bnds:
        if not b.is_pass:
            continue
        if b.max_coverage > 0:
            ratios.append(b.max_coverage / median_cov)
    for s in svs:
        if not s.is_pass:
            continue
        if s.max_coverage > 0:
            ratios.append(s.max_coverage / median_cov)
    if not ratios:
        return floor, 0
    ratios.sort()
    # Linear interpolation at `quantile` (matches numpy.quantile default).
    idx = quantile * (len(ratios) - 1)
    lo, hi = int(idx), min(int(idx) + 1, len(ratios) - 1)
    frac = idx - lo
    p = ratios[lo] + frac * (ratios[hi] - ratios[lo])
    return max(floor, p), len(ratios)


def annotate_noise(
    bnds: list[BND],
    cytobands: dict,
    median_cov: float,
    *,
    mark_acrocentric: bool,
    cov_ratio_thr: float,
    cov_vaf_max: float,
) -> None:
    """Tag each BND in-place with :attr:`BND.noise_flags`.

    Two flag types are set:

    - ``"acrocentric"`` if both endpoints fall in an acrocentric
      short arm (chr13/14/15/21/22 p) and ``mark_acrocentric`` is True.
    - ``"cov_anomaly"`` if ``max(COVERAGE) > cov_ratio_thr * median_cov``
      and ``VAF < cov_vaf_max`` - the repeat-collapse signature.
    """
    p_arm_ends = acrocentric_p_arm_ends(cytobands)
    for b in bnds:
        flags: set = set()
        if mark_acrocentric and (
            in_acrocentric_p(b.chr1, b.pos1, p_arm_ends)
            and in_acrocentric_p(b.chr2, b.pos2, p_arm_ends)
        ):
            flags.add("acrocentric")
        if median_cov > 0 and b.max_coverage > 0:
            ratio = b.max_coverage / median_cov
            if ratio > cov_ratio_thr and b.vaf < cov_vaf_max:
                flags.add("cov_anomaly")
        b.noise_flags = flags


def annotate_sv_noise(
    svs: list[SV],
    median_cov: float,
    *,
    cov_ratio_thr: float,
    cov_vaf_max: float,
    types: tuple = SV_COV_FILTER_TYPES,
) -> None:
    """Mark DEL/DUP SVs whose coverage profile suggests repeat collapse.

    Same composite criterion as the BND filter: ``max(COVERAGE) >
    cov_ratio_thr * median_cov`` AND ``VAF < cov_vaf_max``. Real het
    DELs halve coverage and real DUPs sit at 1.5-2x, so the threshold
    spares true events while catching collapse-like outliers.

    Sets :attr:`SV.is_noise = True` in place. SVs of types not in
    ``types`` are left untouched.
    """
    if median_cov <= 0:
        return
    for s in svs:
        if s.svtype not in types:
            continue
        if s.max_coverage <= 0:
            continue
        ratio = s.max_coverage / median_cov
        if ratio > cov_ratio_thr and s.vaf < cov_vaf_max:
            s.noise_flags.add("cov_anomaly")


def sv_noise_breakdown(svs: list[SV]) -> dict:
    """Return a per-type ``{"pass", "clean", "noise", "cov"}`` summary."""
    out = {}
    for t in SV_TYPES:
        n_pass = sum(1 for s in svs if s.svtype == t and s.is_pass)
        n_clean = sum(1 for s in svs
                      if s.svtype == t and s.is_pass and not s.is_noise)
        n_cov = sum(1 for s in svs
                    if s.svtype == t and s.is_pass
                    and "cov_anomaly" in s.noise_flags)
        out[t] = {
            "pass": n_pass, "clean": n_clean,
            "noise": n_pass - n_clean,
            "cov": n_cov,
        }
    return out


def noise_breakdown(bnds: list[BND]) -> dict:
    """Return per-flag BND noise counts."""
    return {
        "n": len(bnds),
        "clean": sum(1 for b in bnds if not b.is_noise),
        "acrocentric": sum(1 for b in bnds if "acrocentric" in b.noise_flags),
        "cov_anomaly": sum(1 for b in bnds if "cov_anomaly" in b.noise_flags),
        "any_noise": sum(1 for b in bnds if b.is_noise),
    }


def title_suffix(bd: dict) -> str:
    """One-line BND noise breakdown for figure titles."""
    return (
        f"{bd['n']} BNDs: "
        f"clean={bd['clean']}, "
        f"acrocentric={bd['acrocentric']}, "
        f"cov-anomaly={bd['cov_anomaly']}"
    )


# ---------------------------------------------------------------------------
# Shared rendering helpers
# ---------------------------------------------------------------------------

def support_to_lw(
    support: float,
    support_max: float,
    lo: float = 0.4,
    hi: float = 2.4,
) -> float:
    """Map read-support to a line-width in ``[lo, hi]`` (linear, clipped)."""
    if support_max <= 0:
        return lo
    return lo + (hi - lo) * min(support / support_max, 1.0)


def vaf_to_color(vaf: float, cmap=VAF_CMAP):
    """Map a VAF in ``[0, 1]`` to an RGBA colour using the plasma cmap."""
    v = (vaf - VAF_VMIN) / (VAF_VMAX - VAF_VMIN)
    return cmap(min(max(v, 0.0), 1.0))


def render_props(b: BND) -> tuple:
    """Return ``(color, alpha, linestyle)`` for a BND in any plot.

    Noise-flagged events are rendered grey/dashed/faint; non-PASS
    events are rendered VAF-coloured but dashed and dimmer; PASS,
    non-noise events are rendered solid in the VAF colour.
    """
    if b.is_noise:
        return NOISE_COLOR, 0.18, (0, (3, 2))
    if not b.is_pass:
        return vaf_to_color(b.vaf), 0.30, (0, (3, 2))
    return vaf_to_color(b.vaf), 0.70, "-"


# ---------------------------------------------------------------------------
# A) Circos plot
# ---------------------------------------------------------------------------

def plot_circos(
    bnds_unique: list[BND],
    contigs: dict[str, int],
    cytoband_path: Path,
    out_path: Path,
    sample: str,
    n_total: int,
    n_pass: int,
    filter_label: str,
    breakdown: dict,
) -> None:
    """Render the circos plot (BND ribbons on hg38 cytoband ideogram).

    Uses pyCirclize. Writes a PNG to ``out_path`` and closes the figure.

    Parameters
    ----------
    bnds_unique : list[BND]
        BNDs to plot, already deduplicated and noise-annotated.
    contigs : dict[str, int]
        Chromosome lengths from the VCF header.
    cytoband_path : Path
        UCSC cytoband file (gzipped .gz is auto-decompressed to a
        temp file because pyCirclize wants a TSV).
    out_path : Path
        PNG output path.
    sample : str
        Sample label used as a temp-file prefix and in the title.
    n_total, n_pass : int
        Total and PASS-count of BND records (for the title).
    filter_label : str
        ``"pass"`` or ``"all"``, displayed in the title.
    breakdown : dict
        Output of :func:`noise_breakdown`.
    """
    import tempfile

    from pycirclize import Circos

    # pyCirclize needs file paths for the chr-BED and (uncompressed)
    # cytoband TSV. We render in-memory but pyCirclize doesn't, so we
    # stage these as ephemeral files in the OS tempdir and clean up
    # at the end of the function.
    tmp_dir = Path(tempfile.mkdtemp(prefix="molamola_circos_"))
    tmp_bed = tmp_dir / f"{sample}.chr.bed"
    bed_lines = ["#chrom\tchromStart\tchromEnd\tname"]
    for c in CHROM_ORDER:
        if c in contigs:
            bed_lines.append(f"{c}\t0\t{contigs[c]}\t{c}")
    tmp_bed.write_text("\n".join(bed_lines) + "\n")

    cyto_for_pycirclize = cytoband_path
    tmp_cyto = None
    if str(cytoband_path).endswith(".gz"):
        tmp_cyto = tmp_dir / f"{sample}.cytoband.tsv"
        with gzip.open(cytoband_path, "rt") as inp, open(tmp_cyto, "w") as outp:
            outp.write("#chrom\tchromStart\tchromEnd\tname\tgieStain\n")
            for line in inp:
                if line.startswith("#"):
                    continue
                outp.write(line)
        cyto_for_pycirclize = tmp_cyto

    circos = Circos.initialize_from_bed(str(tmp_bed), space=2)
    circos.add_cytoband_tracks((95, 100), str(cyto_for_pycirclize))

    for sector in circos.sectors:
        sector.text(sector.name.replace("chr", ""), r=108, size=10)
        track = sector.get_track("cytoband")
        track.xticks_by_interval(
            100_000_000,
            outer=False,
            label_formatter=lambda v: f"{int(v / 1e6)}",
            label_size=6,
            tick_length=1,
            label_orientation="vertical",
        )

    supports = np.array([b.support for b in bnds_unique], dtype=float)
    smax = supports.max() if supports.size else 1.0

    pad = 250_000
    ordered = sorted(bnds_unique, key=lambda b: (not b.is_noise, b.support))
    for b in ordered:
        if b.chr1 not in contigs or b.chr2 not in contigs:
            continue
        L1, L2 = contigs[b.chr1], contigs[b.chr2]
        s1 = (b.chr1, max(0, b.pos1 - pad), min(L1, b.pos1 + pad))
        s2 = (b.chr2, max(0, b.pos2 - pad), min(L2, b.pos2 + pad))
        color, alpha, _ls = render_props(b)
        try:
            circos.link(s1, s2, color=color, alpha=alpha,
                        height_ratio=0.55,
                        linewidth=support_to_lw(b.support, smax, 0.2, 1.0))
        except Exception as e:  # noqa: BLE001
            print(f"  warn: circos link failed for {b.sv_id}: {e}", file=sys.stderr)

    fig = circos.plotfig()
    # Widen the figure so the VAF colorbar sits clearly right of the
    # circos disk instead of overlapping the outer cytoband ring or
    # its position tick labels.
    w, h = fig.get_size_inches()
    fig.set_size_inches(w * 1.30, h)

    cax = fig.add_axes([0.95, 0.30, 0.012, 0.40])
    sm = plt.cm.ScalarMappable(cmap=VAF_CMAP, norm=plt.Normalize(VAF_VMIN, VAF_VMAX))
    cb = fig.colorbar(sm, cax=cax, label="VAF")
    cb.ax.tick_params(labelsize=8)

    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    # Clean up the staged tempdir.
    import shutil
    shutil.rmtree(tmp_dir, ignore_errors=True)


# ---------------------------------------------------------------------------
# B) Genome SV map
# ---------------------------------------------------------------------------

def _bin_sv_counts(
    svs_per_chrom: dict[str, list[SV]],
    contigs: dict[str, int],
    bin_size: int,
) -> tuple[dict, dict[str, float]]:
    """Bin SV counts per chromosome × type.

    Parameters
    ----------
    svs_per_chrom : dict[str, list[SV]]
        SVs grouped by chromosome.
    contigs : dict[str, int]
        Chromosome lengths.
    bin_size : int
        Bin width in bp.

    Returns
    -------
    bins_cache : dict[str, dict[str, np.ndarray]]
        ``{chrom: {svtype: counts_array}}``.
    type_max : dict[str, float]
        Max bin-count seen per SV type, used for colour normalisation.
    """
    bins_cache: dict = {c: {} for c in svs_per_chrom}
    type_max: dict[str, float] = {t: 1.0 for t in SV_TYPES}
    for c, chr_svs in svs_per_chrom.items():
        L = contigs[c]
        n_bins = (L // bin_size) + 1
        for t in SV_TYPES:
            counts = np.zeros(n_bins, dtype=float)
            for s in chr_svs:
                if s.svtype != t:
                    continue
                b = s.start // bin_size
                if 0 <= b < n_bins:
                    counts[b] += 1
            bins_cache[c][t] = counts
            if counts.size and counts.max() > type_max[t]:
                type_max[t] = float(counts.max())
    return bins_cache, type_max


def _draw_chromosome_row(
    ax,
    chrom: str,
    y0: float,
    chr_h: float,
    strip_h: float,
    contig_len: int,
    cytobands: dict,
    bins: dict,
    type_max: dict[str, float],
    max_len: int,
) -> None:
    """Render one chromosome row: cytoband bar + 4 density strips + label."""
    # Cytoband bar
    for (start, end, _name, stain) in cytobands.get(chrom, []):
        color = CYTOBAND_COLORS.get(stain, "#FFFFFF")
        ax.add_patch(mpatches.Rectangle(
            (start, y0), end - start, chr_h,
            facecolor=color, edgecolor="none",
        ))
    ax.add_patch(mpatches.Rectangle(
        (0, y0), contig_len, chr_h,
        facecolor="none", edgecolor="black", linewidth=0.6,
    ))
    # Density strips (alpha-encoded count above each chromosome)
    for k, t in enumerate(SV_TYPES):
        counts = bins[t]
        y_strip = y0 + chr_h + k * strip_h
        base_rgb = np.array(to_rgb(SV_TYPE_COLOR[t]))
        type_max_t = max(type_max[t], 1.0)
        alphas = np.zeros_like(counts)
        mask = counts > 0
        if mask.any():
            alphas[mask] = (
                0.40 + 0.60 * np.sqrt(counts[mask] / type_max_t).clip(0, 1)
            )
        rgba = np.zeros((1, len(counts), 4))
        rgba[0, :, :3] = base_rgb
        rgba[0, :, 3] = alphas
        ax.imshow(
            rgba,
            extent=(0, contig_len, y_strip, y_strip + strip_h),
            aspect="auto", origin="lower",
            interpolation="nearest",
        )
    # Chromosome label, centred vertically over the cytoband + strip stack
    ax.text(
        -max_len * 0.018,
        y0 + (chr_h + len(SV_TYPES) * strip_h) / 2,
        chrom.replace("chr", ""),
        ha="right", va="center", fontsize=10, fontweight="bold",
    )


def _draw_bnd_arcs(
    ax,
    bnds: list[BND],
    chr_index: dict[str, int],
    y_for: dict[str, float],
    chr_h: float,
    strip_h: float,
    row_h: float,
    max_len: int,
) -> None:
    """Render BND arcs above the chromosome stack.

    Each arc is a quadratic Bezier whose apex sits well above the
    density-strip top of both endpoints, so it visibly clears the
    strips. Apex height grows with the row gap and horizontal span.
    """
    supports = np.array([b.support for b in bnds], dtype=float)
    smax = supports.max() if supports.size else 1.0
    strips_top_offset = chr_h + len(SV_TYPES) * strip_h

    ordered = sorted(bnds, key=lambda x: (not x.is_noise, x.support))
    for b in ordered:
        if b.chr1 not in y_for or b.chr2 not in y_for:
            continue
        y_chr_top1 = y_for[b.chr1] + chr_h
        y_chr_top2 = y_for[b.chr2] + chr_h
        strips_top1 = y_for[b.chr1] + strips_top_offset
        strips_top2 = y_for[b.chr2] + strips_top_offset

        x1, x2 = b.pos1, b.pos2
        row_gap_steps = abs(chr_index[b.chr1] - chr_index[b.chr2])
        span_frac = abs(x2 - x1) / max_len
        apex_y = (
            max(strips_top1, strips_top2)
            + 0.6
            + 0.45 * row_gap_steps * row_h
            + 0.6 * span_frac
        )
        cx = (x1 + x2) / 2
        path = MplPath(
            [(x1, y_chr_top1), (cx, apex_y), (x2, y_chr_top2)],
            [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3],
        )
        color, alpha, ls = render_props(b)
        lw = support_to_lw(b.support, smax, 0.4, 2.0)
        ax.add_patch(mpatches.PathPatch(
            path, facecolor="none", edgecolor=color,
            alpha=alpha, linewidth=lw, linestyle=ls,
        ))
        ax.scatter(
            [x1, x2], [y_chr_top1, y_chr_top2], s=8,
            facecolors=[color], edgecolors="black",
            linewidths=0.3, zorder=4, alpha=alpha,
        )


def _add_genome_map_decor(
    fig,
    ax,
    sample: str,
    breakdown: dict,
    type_max: dict[str, float],
    n_svs_by_type: dict[str, int],
    n_total: int,
    n_pass: int,
    filter_label: str,
    bin_size: int,
    min_svlen: int,
    max_len: int,
    n_chr: int,
    row_h: float,
) -> None:
    """Apply axes formatting, title, legends, and colorbar to the genome map."""
    top_y = n_chr * row_h + 4.0
    ax.set_xlim(-max_len * 0.06, max_len * 1.02)
    ax.set_ylim(-0.4, top_y)
    ax.set_yticks([])
    ax.set_xticks(np.arange(0, max_len + 1, 50_000_000))
    ax.set_xticklabels([str(int(x / 1e6)) for x in np.arange(0, max_len + 1, 50_000_000)])
    ax.set_xlabel("Position (Mb)")
    ax.spines[["top", "right", "left"]].set_visible(False)

    bin_mb = bin_size // 1_000_000
    cyto_legend = [
        mpatches.Patch(facecolor=CYTOBAND_COLORS[s], edgecolor="black", label=s)
        for s in ("gneg", "gpos50", "gpos100", "acen", "gvar")
    ]
    leg1 = ax.legend(
        handles=cyto_legend, loc="upper center",
        bbox_to_anchor=(0.5, -0.04),
        fontsize=7, frameon=True, framealpha=0.85,
        edgecolor="#cccccc", ncols=5, title="cytoband",
    )
    ax.add_artist(leg1)

    type_handles = [
        mpatches.Patch(facecolor=SV_TYPE_COLOR[t], edgecolor="black",
                       label=f"{t} (peak {int(type_max[t])})")
        for t in SV_TYPES
    ]
    type_handles.append(plt.Line2D([0], [0], color=NOISE_COLOR,
                                    linestyle=(0, (3, 2)), linewidth=1.4,
                                    label="BND noise"))
    leg2 = ax.legend(
        handles=type_handles, loc="upper center",
        bbox_to_anchor=(0.5, -0.10),
        fontsize=7, frameon=True, framealpha=0.85,
        edgecolor="#cccccc", ncols=5,
        title=f"SV density per {bin_mb} Mb bin (alpha scaled to peak)",
    )
    ax.add_artist(leg2)

    sm = plt.cm.ScalarMappable(cmap=VAF_CMAP, norm=plt.Normalize(VAF_VMIN, VAF_VMAX))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.022, pad=0.02)
    cbar.set_label("BND VAF", fontsize=9)


def plot_genome_sv_map(
    bnds_unique: list[BND],
    svs: list[SV],
    contigs: dict[str, int],
    cytobands: dict,
    out_path: Path,
    sample: str,
    n_total: int,
    n_pass: int,
    filter_label: str,
    breakdown: dict,
    *,
    bin_size: int = 1_000_000,
    min_svlen: int = 50,
) -> None:
    """Render the genome SV map (cytobands + density strips + BND arcs).

    Parameters
    ----------
    bnds_unique : list[BND]
        BNDs (deduplicated, noise-annotated) to overlay as arcs.
    svs : list[SV]
        Non-BND SVs. Already size-filtered upstream (`--min-svlen`);
        here we additionally drop non-PASS and noise-flagged events.
    contigs : dict[str, int]
        Chromosome lengths.
    cytobands : dict
        Loaded by :func:`load_cytobands`.
    out_path : Path
        PNG output path.
    sample : str
        Sample label for the title.
    n_total, n_pass : int
        Total and PASS counts for BNDs (for the title).
    filter_label : str
        ``"pass"`` or ``"all"``.
    breakdown : dict
        Output of :func:`noise_breakdown`.
    bin_size : int, optional
        Density-track bin width in bp (default 1,000,000).
    min_svlen : int, optional
        The current `--min-svlen` value; used in the figure title for
        annotation only (the actual filtering is applied upstream).
    """
    chroms_present = [c for c in CHROM_ORDER if c in contigs]
    n = len(chroms_present)

    svs_filt = [s for s in svs if s.is_pass and not s.is_noise]
    svs_per_chrom = {c: [s for s in svs_filt if s.chrom == c] for c in chroms_present}
    bins_cache, type_max = _bin_sv_counts(svs_per_chrom, contigs, bin_size)

    fig, ax = plt.subplots(figsize=(17, 14))

    chr_h = 0.28
    strip_h = 0.13
    inter_row_gap = 0.06
    row_h = chr_h + len(SV_TYPES) * strip_h + inter_row_gap

    chr_index = {c: i for i, c in enumerate(chroms_present)}
    y_for = {c: (n - i - 1) * row_h for i, c in enumerate(chroms_present)}
    max_len = max(contigs[c] for c in chroms_present)

    for c in chroms_present:
        _draw_chromosome_row(
            ax, c, y_for[c], chr_h, strip_h,
            contigs[c], cytobands, bins_cache[c], type_max, max_len,
        )

    _draw_bnd_arcs(
        ax, bnds_unique, chr_index, y_for,
        chr_h, strip_h, row_h, max_len,
    )

    n_svs_by_type = {t: sum(1 for s in svs_filt if s.svtype == t) for t in SV_TYPES}
    _add_genome_map_decor(
        fig, ax, sample, breakdown, type_max, n_svs_by_type,
        n_total, n_pass, filter_label,
        bin_size, min_svlen, max_len, n, row_h,
    )

    # The two legends sit BELOW the axes via bbox_to_anchor; tight-cropping
    # needs them explicitly listed or it can omit entries that overflow.
    from matplotlib.legend import Legend
    extra = list(ax.findobj(Legend))
    fig.savefig(out_path, dpi=200, bbox_inches="tight",
                 bbox_extra_artists=extra)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Compound-het renderer
# ---------------------------------------------------------------------------

# Build label used in the panel title and HTML report. ``hg38`` and
# ``t2t`` are stored lowercase internally; surface the cytogeneticist-
# friendly capitalised forms for display.
_BUILD_LABEL: dict[str, str] = {
    "hg38": "GRCh38",
    "t2t":  "T2T-CHM13v2",
}


def _build_phase_blocks(
    variants: list[PhasedVariant],
) -> list[PhaseBlock]:
    """Group variants into phase blocks by ``ps``.

    Each block's ``start``/``end`` is the min/max position of the
    variants in that phase set. ``n_phased`` counts variants of any
    consequence. Returns blocks sorted by ``start``.
    """
    by_ps: dict[int, list[PhasedVariant]] = {}
    for v in variants:
        by_ps.setdefault(v.ps, []).append(v)
    blocks: list[PhaseBlock] = []
    for ps, vs in by_ps.items():
        positions = [v.pos for v in vs]
        blocks.append(PhaseBlock(
            ps=ps, start=min(positions), end=max(positions),
            n_phased=len(vs),
        ))
    blocks.sort(key=lambda b: b.start)
    return blocks


def draw_compound_het_panel(
    ax,
    gene: Gene,
    variants: list[PhasedVariant],
    phase_blocks: list[PhaseBlock],
) -> tuple[int, int]:
    """Render one gene's phased-haplotype panel onto ``ax``.

    Visual spec (locked 2026-05-02; do not relitigate) ports verbatim
    from the prototype's ``draw_panel`` (compound-het-plotter
    ``plot_gene.py`` lines 266-415):

    - exon track at top (thin grey line + dark-grey rectangles)
    - two horizontal hap lines H1 / H2
    - mint phase-block rectangles spanning both haps with off-edge
      arrows when the block extends past the gene window
    - missense lollipops on the canonical transcript (stem outward
      from hap line; coloured circle at tip)
    - non-missense ticks on the hap line for synonymous / intronic /
      etc. that landed in the same phase block
    - panel title with build, locus, and ``T trans, C cis`` summary

    The renderer does NOT draw trans/cis arcs. Opposite-hap dots in
    the same block already read as trans, same-hap dots read as cis,
    and arcs were visual clutter (removed in the prototype). Pair
    counts are computed via :func:`classify_pairs` for the title.

    Returns ``(n_trans, n_cis)`` for the report's run-metadata.
    """
    pad = 0.02 * (gene.end - gene.start) if gene.end > gene.start else 0
    xmin = gene.start - pad
    xmax = gene.end + pad
    span = xmax - xmin if xmax > xmin else 1.0
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-1.4, 2.0)
    ax.set_yticks([HAP_Y[0], HAP_Y[1]])
    ax.set_yticklabels(["H1", "H2"])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)

    # Format x-axis as genome positions. matplotlib defaults to a
    # scientific-notation offset (e.g. ``1.427e8``) for the large
    # integers typical of genomic coordinates; that's unreadable for
    # cytogeneticists. Always show in Mb (genome positions are
    # always Mb-scale), with decimals chosen by the gene span so
    # adjacent tick labels stay distinguishable. For tiny synthetic
    # windows (< 1 kb) drop to bp. FuncFormatter replaces the
    # default ScalarFormatter so the offset is gone implicitly.
    if span >= 10e6:
        decimals, unit, divisor = 0, "Mb", 1e6
    elif span >= 1e6:
        decimals, unit, divisor = 1, "Mb", 1e6
    elif span >= 100e3:
        decimals, unit, divisor = 2, "Mb", 1e6
    elif span >= 10e3:
        decimals, unit, divisor = 3, "Mb", 1e6
    elif span >= 1e3:
        decimals, unit, divisor = 1, "kb", 1e3
    else:
        decimals, unit, divisor = 0, "bp", 1
    def _genome_pos_fmt(x, _pos, _d=decimals, _u=unit, _v=divisor):
        if _v == 1:
            return f"{int(x):,} {_u}"
        return f"{x / _v:.{_d}f} {_u}"
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(_genome_pos_fmt))

    block_y_top = 1.2
    block_y_bot = -1.2
    arrow_dx = span * 0.012
    for idx, b in enumerate(phase_blocks):
        x_left_clip = max(b.start, xmin)
        x_right_clip = min(b.end, xmax)
        if x_right_clip <= x_left_clip:
            continue
        ax.add_patch(mpatches.Rectangle(
            (x_left_clip, block_y_bot),
            x_right_clip - x_left_clip,
            block_y_top - block_y_bot,
            facecolor=BLOCK_FILL,
            edgecolor=BLOCK_EDGE, linewidth=BLOCK_EDGE_LW, zorder=0,
        ))
        if b.start < xmin:
            ax.annotate(
                "", xy=(xmin - arrow_dx * 0.4, 0),
                xytext=(xmin + arrow_dx, 0),
                arrowprops=dict(
                    arrowstyle="->,head_width=0.4,head_length=0.5",
                    color="#666", lw=0.9,
                ),
                zorder=2,
            )
        if b.end > xmax:
            ax.annotate(
                "", xy=(xmax + arrow_dx * 0.4, 0),
                xytext=(xmax - arrow_dx, 0),
                arrowprops=dict(
                    arrowstyle="->,head_width=0.4,head_length=0.5",
                    color="#666", lw=0.9,
                ),
                zorder=2,
            )
        label_x = (x_left_clip + x_right_clip) / 2
        block_kb = (b.end - b.start) / 1000.0
        ax.text(
            label_x, block_y_top + 0.05,
            f"block {idx + 1}  ·  PS {b.ps}  ·  {b.n_phased} phased  ·  "
            f"{block_kb:.0f} kb",
            ha="center", va="bottom", fontsize=6.5, color="#666",
            zorder=2,
        )

    if gene.canonical_exons:
        ax.plot(
            [gene.start, gene.end], [EXON_Y, EXON_Y],
            color=EXON_CONNECTOR_COLOR, lw=0.7, zorder=1,
        )
        for s, e in gene.canonical_exons:
            ax.add_patch(mpatches.Rectangle(
                (s, EXON_Y - EXON_H / 2), max(e - s, 1), EXON_H,
                facecolor=EXON_FILL, edgecolor="none", zorder=2,
            ))

    ax.axhline(HAP_Y[0], color="#888", lw=0.7, zorder=1)
    ax.axhline(HAP_Y[1], color="#888", lw=0.7, zorder=1)

    missense = [
        v for v in variants
        if v.is_canonical_transcript and is_missense_consequence(v.consequence)
    ]
    other = [v for v in variants if v not in missense]

    for v in other:
        if v.variant_hap not in (1, 2):
            continue
        y = HAP_Y[v.variant_hap - 1]
        ax.plot(
            v.pos, y, marker=NON_MISSENSE_MARKER,
            color=NON_MISSENSE_COLOR, markersize=NON_MISSENSE_SIZE,
            markeredgecolor=NON_MISSENSE_COLOR,
            markeredgewidth=NON_MISSENSE_EDGE_WIDTH, zorder=3,
        )

    for v in missense:
        if v.variant_hap not in (1, 2):
            continue
        y_line = HAP_Y[v.variant_hap - 1]
        y_dot = LOLLIPOP_TIP_Y[v.variant_hap - 1]
        color = CLNSIG_COLOR.get(v.clnsig, DEFAULT_CLNSIG_COLOR)
        ax.plot(
            [v.pos, v.pos], [y_line, y_dot],
            color="#666", lw=0.6, zorder=3,
        )
        ax.plot(
            v.pos, y_dot, marker="o", color=color, markersize=LOLLIPOP_MS,
            markeredgecolor="black", markeredgewidth=LOLLIPOP_EDGE_LW,
            zorder=5,
        )

    n_trans, n_cis = classify_pairs(variants)
    # Just the gene symbol baked into the figure so a downloaded PNG
    # carries its identity. Locus, build, and stats live in the HTML
    # <h2> + chip row above the embedded image.
    ax.set_title(gene.symbol, fontsize=14, loc="left", fontweight="bold")
    return n_trans, n_cis


def render_compound_het_png(
    gene: Gene,
    variants: list[PhasedVariant],
    phase_blocks: list[PhaseBlock],
    *,
    reference: str,
) -> tuple[bytes, dict]:
    """Render one gene panel + legend to a PNG byte string.

    Returns ``(png_bytes, stats)`` where ``stats`` is a small dict
    of per-gene counts (``n_trans``, ``n_cis``, ``n_missense``,
    ``n_total``, ``n_blocks``) for inclusion in the HTML report's
    run-metadata.

    Width auto-scales with variant count, capped at 24 inches; dpi
    is reduced so the rendered PNG never exceeds
    :data:`MAX_PX_WIDTH` pixels (Anthropic many-image upload cap).
    """
    from matplotlib.lines import Line2D

    n_vars = len(variants)
    width = max(12, min(24, 8 + n_vars * 0.4))
    # No constrained_layout: it does not reserve space for a
    # figure-level legend, which leaves the legend handles overlapping
    # the x-axis tick labels. Manage spacing with subplots_adjust
    # instead, then drop bbox_inches="tight" so savefig respects the
    # reserved bottom margin.
    fig, ax = plt.subplots(figsize=(width, 5.0))
    fig.subplots_adjust(left=0.04, right=0.99, top=0.90, bottom=0.18)
    n_trans, n_cis = draw_compound_het_panel(
        ax, gene, variants, phase_blocks,
    )

    legend_marker_size = 11
    legend_handles = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=CLNSIG_COLOR["p_or_lp"],
               markersize=legend_marker_size, label="P / LP"),
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=CLNSIG_COLOR["vus"],
               markersize=legend_marker_size, label="VUS"),
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=CLNSIG_COLOR["conflicting"],
               markersize=legend_marker_size, label="conflicting"),
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=CLNSIG_COLOR["benign"],
               markersize=legend_marker_size, label="benign"),
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=DEFAULT_CLNSIG_COLOR,
               markersize=legend_marker_size, label="no ClinVar"),
        Line2D([0], [0], marker=NON_MISSENSE_MARKER, color="w",
               markerfacecolor=NON_MISSENSE_COLOR,
               markeredgecolor=NON_MISSENSE_COLOR,
               markersize=max(5, NON_MISSENSE_SIZE) + 2,
               label="synonymous"),
    ]
    fig.legend(
        handles=legend_handles, loc="lower center", ncol=6,
        fontsize=10.5, frameon=False,
        bbox_to_anchor=(0.5, 0.02),
        bbox_transform=fig.transFigure,
    )

    target_dpi = min(150.0, MAX_PX_WIDTH / fig.get_figwidth())
    buf = io.BytesIO()
    fig.savefig(buf, dpi=target_dpi, format="png")
    plt.close(fig)
    stats = {
        "n_trans": n_trans, "n_cis": n_cis,
        "n_missense": len([
            v for v in variants
            if v.is_canonical_transcript
            and is_missense_consequence(v.consequence)
        ]),
        "n_total": len(variants), "n_blocks": len(phase_blocks),
    }
    return buf.getvalue(), stats


# ---------------------------------------------------------------------------
# Self-contained HTML report
# ---------------------------------------------------------------------------

def _esc(s) -> str:
    """Shorthand for html.escape on arbitrary values."""
    return html_mod.escape(str(s))


def make_html_report(
    sample: str,
    circos_png: bytes,
    sv_map_png: bytes,
    out_path: Path,
    *,
    bd: dict,
    sv_bd: dict,
    median_cov: float,
    cov_ratio: float,
    cov_vaf_max: float,
    n_total_bnd: int,
    n_pass_bnd: int,
    filter_label: str,
    caller: str,
    caller_basis: str,
) -> None:
    """Write a single-file, self-contained HTML report to ``out_path``.

    Embeds both PNGs as base64 data URIs (so the file works when
    emailed or copied around) plus a collapsible run-metadata block
    at the bottom (caller, filter, coverage thresholds, noise
    breakdown). No external CSS / JS / font dependencies.
    """
    circos_uri = "data:image/png;base64," + base64.b64encode(circos_png).decode("ascii")
    sv_map_uri = "data:image/png;base64," + base64.b64encode(sv_map_png).decode("ascii")

    sv_summary = " · ".join(
        f"{t}={sv_bd[t]['pass']}"
        + (f" (cov-noise={sv_bd[t]['cov']})"
           if sv_bd[t].get("cov", 0) else "")
        for t in SV_TYPES
    )

    bd_chips = []
    for k, label in [
        ("clean", "clean"),
        ("acrocentric", "acrocentric"),
        ("cov_anomaly", "cov-anomaly"),
    ]:
        v = bd.get(k, 0)
        bd_chips.append(f'<span class="chip">{label} = {v}</span>')

    caller_chip_class = "ok" if caller == "sniffles2" else "warn"
    caller_chip_text = (
        f"caller: {caller}"
        if caller == "sniffles2"
        else f"caller: {caller} (NOTE — parser tuned for sniffles2)"
    )

    sample_e = _esc(sample)
    html = (
        "<!DOCTYPE html>\n<html lang=\"en\"><head>\n<meta charset=\"UTF-8\">\n"
        f"<title>molamola — {sample_e}</title>\n"
        f"<style>{_HTML_REPORT_CSS}</style>\n"
        "</head>\n<body>\n"
        f"<header>\n{_header_svg()}</header>\n"
        f"<section class=\"figure\" id=\"fig-circos\">\n"
        f"  <img src=\"{circos_uri}\" alt=\"circos plot for {sample_e}\">\n"
        f"</section>\n"
        f"<section class=\"figure\" id=\"fig-svmap\">\n"
        f"  <img src=\"{sv_map_uri}\" alt=\"genome SV map for {sample_e}\">\n"
        f"</section>\n"
        f"<details class=\"meta\">\n"
        f"  <summary>run metadata</summary>\n"
        f"  <div class=\"chips\">\n"
        f"    <span class=\"chip {caller_chip_class}\">{_esc(caller_chip_text)}</span>\n"
        f"    <span class=\"chip\">filter: {_esc(filter_label)} ({n_pass_bnd}/{n_total_bnd} BNDs PASS)</span>\n"
        f"    <span class=\"chip\">baseline cov: {median_cov:.1f}×</span>\n"
        f"    <span class=\"chip\">cov-anomaly: max≥{cov_ratio*median_cov:.0f}× AND VAF&lt;{cov_vaf_max}</span>\n"
        f"  </div>\n"
        f"  <div class=\"chips\"><strong>BND noise:</strong> {''.join(bd_chips)}</div>\n"
        f"  <div class=\"chips\"><strong>Non-BND SVs (PASS):</strong> {_esc(sv_summary)}</div>\n"
        f"  <div class=\"basis\">caller fingerprint basis: {_esc(caller_basis)}</div>\n"
        f"</details>\n"
        f"<footer>generated by <code>molamola</code> "
        f"&middot; circos via "
        f"<a href=\"https://github.com/moshi4/pyCirclize\">pyCirclize</a> "
        f"&middot; genome SV map via "
        f"<a href=\"https://matplotlib.org/\">matplotlib</a></footer>\n"
        "</body></html>\n"
    )
    out_path.write_text(html)


def make_compound_het_report(
    sample: str,
    gene_panels: list[tuple[Gene, bytes, dict]],
    out_path: Path,
    *,
    reference: str,
    n_genes_scanned: int,
    n_genes_plotted: int,
    n_genes_capped: int,
    clinvar_source: str,
    canonical_exons_source: str,
    refusal_counts: dict[str, object],
    selection_rule: str,
    empty_genes: list[str],
    strict_symbols: set[str] | None = None,
    is_auto_select: bool = False,
) -> None:
    """Assemble the self-contained compound-het HTML report.

    One section per plotted gene with anchor ``id="gene-<symbol>"``.
    Run-metadata sits in a collapsible ``<details>`` block at the
    bottom: ClinVar source, canonical-exons source, refusal counters,
    selection rule, genes-scanned vs plotted, and a T2T-mode
    disclaimer about HGVSc-based ClinVar matching being best-effort.

    In auto-select mode the panels are split into two clearly
    labelled groups:

    - *strict* (both variants P/LP or VUS — true compound-het)
    - *extended* (anchor P/LP-or-VUS, partner not benign)

    The strict heading appears even when its subset is empty, so
    readers always see the dichotomy. Explicit ``--gene`` mode
    skips the split and renders panels in a single section.

    Empty genes (``--gene FOO`` with 0 phased hets in the gene
    window) get a placeholder section so multi-gene runs don't
    silently drop them.
    """
    sample_e = _esc(sample)
    build_label_e = _esc(_BUILD_LABEL.get(reference, reference))
    strict_set = strict_symbols or set()

    def _render_gene_panel(gene: Gene, png: bytes, stats: dict) -> str:
        png_uri = (
            "data:image/png;base64,"
            + base64.b64encode(png).decode("ascii")
        )
        symbol_e = _esc(gene.symbol)
        locus_e = _esc(
            f"{gene.chrom}:{gene.start:,}-{gene.end:,} ({gene.strand})"
        )
        stats_chips = (
            f'<span class="chip">{stats["n_missense"]} missense / '
            f'{stats["n_total"]} phased hets</span>'
            f'<span class="chip">{stats["n_blocks"]} block(s)</span>'
            f'<span class="chip">{stats["n_trans"]} trans</span>'
            f'<span class="chip">{stats["n_cis"]} cis</span>'
        )
        return (
            f'<section class="figure gene-panel" id="gene-{symbol_e}">\n'
            f'  <h2>{build_label_e} &middot; {locus_e}</h2>\n'
            f'  <div class="chips">{stats_chips}</div>\n'
            f'  <img src="{png_uri}" '
            f'alt="phased haplotype panel for {symbol_e}">\n'
            f'</section>\n'
        )

    def _empty_panel(symbol: str, blurb: str) -> str:
        symbol_e = _esc(symbol)
        return (
            f'<section class="figure gene-panel empty" id="gene-{symbol_e}">\n'
            f'  <h2><span class="sample">{symbol_e}</span></h2>\n'
            f'  <p class="placeholder">{blurb}</p>\n'
            f'</section>\n'
        )

    sections: list[str] = []
    if is_auto_select:
        strict_panels = [
            p for p in gene_panels if p[0].symbol in strict_set
        ]
        extended_panels = [
            p for p in gene_panels if p[0].symbol not in strict_set
        ]
        sections.append(
            '<div class="auto-select-section strict-section">\n'
            '  <h2 class="section-heading">Strict compound-het '
            '<span class="section-blurb">(both variants ClinVar P/LP '
            'or VUS)</span></h2>\n'
        )
        if strict_panels:
            for tup in strict_panels:
                sections.append(_render_gene_panel(*tup))
        else:
            sections.append(
                '<p class="section-empty">No genes meet the strict '
                'rule in this VCF.</p>\n'
            )
        sections.append('</div>\n')

        sections.append(
            '<div class="auto-select-section extended-section">\n'
            '  <h2 class="section-heading">Extended candidates '
            '<span class="section-blurb">(anchor P/LP or VUS, '
            'partner conflicting / no-ClinVar / P/LP / VUS)'
            '</span></h2>\n'
        )
        if extended_panels:
            for tup in extended_panels:
                sections.append(_render_gene_panel(*tup))
        else:
            sections.append(
                '<p class="section-empty">No additional extended '
                'candidates in this VCF.</p>\n'
            )
        sections.append('</div>\n')
    else:
        # Explicit --gene mode: single flat section, no strict/extended split.
        for tup in gene_panels:
            sections.append(_render_gene_panel(*tup))

    for symbol in empty_genes:
        sections.append(_empty_panel(
            symbol,
            f"No phased het variants in {_esc(symbol)} in the input "
            f"VCF; nothing to plot.",
        ))

    refusal_keys = ("kept", "unphased", "non_het", "multi_allelic",
                    "no_canonical_csq")
    refusal_chips = "".join(
        f'<span class="chip">{_esc(k.replace("_", " "))} = '
        f'{_esc(refusal_counts.get(k, 0))}</span>'
        for k in refusal_keys
    )

    cap_chip = (
        f'<span class="chip warn">capped at {n_genes_capped} genes; '
        f'pass --max-genes higher to plot all</span>'
        if n_genes_capped > 0 else ""
    )
    html = (
        "<!DOCTYPE html>\n<html lang=\"en\"><head>\n"
        "<meta charset=\"UTF-8\">\n"
        f"<title>molamola compound-het — {sample_e}</title>\n"
        f"<style>{_HTML_REPORT_CSS}{_COMPOUND_HET_REPORT_CSS}</style>\n"
        "</head>\n<body>\n"
        f"<header>\n{_header_svg()}</header>\n"
        f"<div class=\"chips\">"
        f'<span class="chip">sample: {sample_e}</span>'
        f'<span class="chip">build: {build_label_e}</span>'
        f'<span class="chip">{n_genes_plotted} of '
        f'{n_genes_scanned} candidate gene(s) plotted</span>'
        f"{cap_chip}"
        f"</div>\n"
        + "".join(sections)
        + f"<details class=\"meta\">\n"
        f"  <summary>run metadata</summary>\n"
        f"  <div class=\"chips\"><strong>VCF parse:</strong> {refusal_chips}</div>\n"
        f"  <div class=\"chips\"><strong>selection rule:</strong> "
        f"{_esc(selection_rule)}</div>\n"
        f"  <div class=\"basis\">ClinVar source: "
        f"{_esc(clinvar_source)}</div>\n"
        f"  <div class=\"basis\">canonical exons source: "
        f"{_esc(canonical_exons_source)}</div>\n"
        f"</details>\n"
        f"<footer>generated by <code>molamola</code> "
        f"&middot; figures via "
        f"<a href=\"https://matplotlib.org/\">matplotlib</a></footer>\n"
        "</body></html>\n"
    )
    out_path.write_text(html)


_COMPOUND_HET_REPORT_CSS = """
.gene-panel { margin: 2.2em 0; padding-top: 0.6em;
              border-top: 1px solid var(--border); }
.auto-select-section > .gene-panel:first-of-type { border-top: none;
                                                    padding-top: 0; }
.gene-panel h2 { text-align: left; margin: 0.2em 0 0.4em; }
.gene-panel.empty .placeholder { color: var(--muted); font-style: italic;
                                  text-align: left; margin: 0.4em 0; }
.auto-select-section { margin: 2.4em 0 2em; }
.section-heading { font-size: 1.15em; font-weight: 600;
                    border-bottom: 2px solid var(--border);
                    padding: 0 0 0.3em; margin: 0.2em 0 0.6em;
                    text-align: left; }
.section-blurb { font-weight: 400; color: var(--muted);
                  font-size: 0.85em; margin-left: 0.4em; }
.strict-section .section-heading { border-bottom-color: #c0143c; }
.extended-section .section-heading { border-bottom-color: #f4a013; }
.section-empty { color: var(--muted); font-style: italic;
                  margin: 0.4em 0 1em; }
"""


_HTML_REPORT_CSS = """
:root { --fg: #1a1a1a; --bg: #fff; --muted: #6b6b6b; --border: #ddd;
        --chip-bg: #f0f0f0; --chip-warn: #fff4d6; }
* { box-sizing: border-box; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, system-ui, sans-serif;
       max-width: 1500px; margin: 1.5em auto; color: var(--fg); padding: 0 1.5em; line-height: 1.5; }
header { margin-bottom: 1.5em; }
header svg.banner { display: block; width: 100%; height: auto;
                    margin: 0 auto 0.4em; }
.chips { margin-bottom: 0.5em; font-size: 0.9em; }
.chip { display: inline-block; padding: 0.15em 0.55em; margin: 0.2em 0.3em 0.2em 0;
        background: var(--chip-bg); border-radius: 4px; }
.chip.ok { background: #e6f5e6; }
.chip.warn { background: var(--chip-warn); border: 1px solid #e0c060; }
.basis { color: var(--muted); font-size: 0.85em; margin-top: 0.5em; }
.figure { text-align: center; margin: 1.5em 0; }
.figure h2 { font-weight: 500; font-size: 1.05em; color: var(--muted);
              margin: 0.5em 0; text-align: center; }
.figure h2 .sample { font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
                     color: var(--fg); }
.figure img { max-width: 100%; height: auto; }
details.meta { margin: 3em 0 1em; padding: 0.8em 1em; background: #fafafa;
               border: 1px solid var(--border); border-radius: 6px; }
details.meta > summary { cursor: pointer; font-weight: 500; color: var(--muted);
                          padding: 0.2em 0; user-select: none; }
details.meta[open] > summary { margin-bottom: 0.6em; }
footer { color: var(--muted); font-size: 0.8em; text-align: center; margin: 4em 0 1em; }
footer a { color: inherit; text-decoration: underline; text-decoration-color: #cfcfcf;
           text-underline-offset: 2px; }
footer a:hover { color: var(--fg); text-decoration-color: var(--muted); }
"""


def _header_svg() -> str:
    """Build the inline header SVG, embedding the bundled fish PNG."""
    fish_path = Path(__file__).resolve().parent / "data" / "header_fish.png"
    if fish_path.exists():
        fish_uri = (
            "data:image/png;base64,"
            + base64.b64encode(fish_path.read_bytes()).decode("ascii")
        )
        fish_image = (
            f'    <image x="1340" y="10" width="130" height="130" '
            f'href="{fish_uri}" preserveAspectRatio="xMidYMid meet"/>\n'
        )
    else:
        fish_image = ""
    return (
        '  <svg class="banner" xmlns="http://www.w3.org/2000/svg" '
        'viewBox="0 0 1500 150" preserveAspectRatio="xMidYMid meet" '
        'role="img" aria-label="molamola">\n'
        '    <text x="40" y="88" font-family="-apple-system, '
        'BlinkMacSystemFont, \'Segoe UI\', sans-serif" font-weight="800" '
        'font-size="72" fill="#1e3a5f">molamola</text>\n'
        '    <text x="46" y="118" font-family="-apple-system, '
        'BlinkMacSystemFont, \'Segoe UI\', sans-serif" font-weight="400" '
        'font-size="14" fill="#59636e" font-style="italic">'
        'python tool for visualizing variation</text>\n'
        f'{fish_image}'
        '  </svg>\n'
    )


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _add_common_args(p) -> None:
    """Top-level args shared by all modes."""
    p.add_argument("--vcf", default=None, type=Path,
                   help="input VCF. molamola auto-detects the plot type "
                        "from the header: ##INFO=<ID=SVTYPE> selects the "
                        "SV / cytogenetics report; ##INFO=<ID=CSQ> + "
                        "##FORMAT=<ID=PS> selects the per-gene phased-"
                        "haplotype panels (compound-het). Required "
                        "unless --mosdepth is given; when combined with "
                        "--mosdepth the VCF is consumed only as the BAF "
                        "source for the karyotype panel.")
    p.add_argument("--mosdepth", default=None, type=Path,
                   help="mosdepth regions.bed.gz output (uniform-bin runs "
                        "only). Activates karyotype coverage mode. "
                        "Required unless --vcf is given.")
    p.add_argument("--out", default=None, type=Path,
                   help="output directory (default: parent directory of "
                        "the input file)")
    p.add_argument(
        "--reference",
        choices=list(SUPPORTED_REFERENCES),
        default="hg38",
        help="reference assembly the input VCF was called against. "
             "SV mode supports hg38 + T2T-CHM13v2.0; compound-het "
             "mode is hg38-only (the bundled canonical-exon and "
             "ClinVar refs are hg38-coordinate). Default hg38.",
    )
    p.add_argument("--sample", default=None,
                   help="sample label for the report header "
                        "(default: VCF basename)")
    p.add_argument("--force", action="store_true",
                   help="bypass the safety check that errors out when the "
                        "VCF filename hints at a reference different from "
                        "--reference (e.g. 'sample.t2t.vcf' with "
                        "--reference hg38).")
    p.add_argument("--png", action="store_true",
                   help="alongside the HTML report, also write each "
                        "embedded figure as a standalone PNG file in "
                        "the same output directory (useful for MultiQC "
                        "and other pipeline-report embeds).")


def _add_sv_args(p) -> None:
    """SV-mode flags. Used when the VCF carries ##INFO=<ID=SVTYPE,...>."""
    p.add_argument("--filter", choices=["pass", "all"], default="pass",
                   help="'pass' (default) keeps PASS BNDs only; "
                        "'all' keeps PASS + GT-filtered")
    p.add_argument("--mark-acrocentric", action=argparse.BooleanOptionalAction,
                   default=None,
                   help="grey out BNDs with both ends in acrocentric p-arms "
                        "(chr13/14/15/21/22). Default ON for hg38 (where "
                        "those p-arms are mostly N-padded mappability sinks) "
                        "and OFF for T2T (where the p-arms are real, "
                        "fully-resolved sequence). Set explicitly to override.")
    p.add_argument("--cov-filter", choices=["none", "mark", "drop"],
                   default="mark",
                   help="how to handle BNDs / DEL / DUP whose breakpoint sits "
                        "on a coverage spike (default 'mark' = grey BNDs and "
                        "drop noisy DEL/DUP from density; 'drop' = also remove "
                        "flagged BNDs entirely; 'none' = ignore)")
    p.add_argument("--cov-ratio", type=_parse_cov_ratio, default="auto",
                   help="max(COVERAGE) / genome-median-coverage above which "
                        "an event is suspicious. Default 'auto' computes "
                        "max(2.0, p99 of the in-sample distribution) so "
                        "the threshold adapts to each sample's coverage "
                        "profile. Pass a number to override.")
    p.add_argument("--cov-vaf-max", type=float, default=0.35,
                   help="VAF below which a high-coverage event is treated as "
                        "repeat-collapse noise (default 0.35)")
    p.add_argument("--focus", action="append", type=parse_focus,
                   default=None, metavar="CHR:POS",
                   help="show only BNDs with an endpoint within --focus-window "
                        "of CHR:POS; flag repeatable. Filters BNDs only; non-BND "
                        "SVs are still shown genome-wide on the SV map.")
    p.add_argument("--focus-window", type=int, default=1000,
                   help="+/-bp tolerance for --focus matching (default 1000)")
    p.add_argument("--min-svlen", type=int, default=50,
                   help="hard SVLEN cutoff (bp) for non-BND SVs. "
                        "Events shorter than this are dropped from every "
                        "downstream consumer. BNDs are unaffected. "
                        "Default 50; set 0 to disable.")
    p.add_argument("--bin-size", type=int, default=1_000_000,
                   help="bin width (bp) for genome SV map density tracks "
                        "(default 1,000,000)")
    p.add_argument(
        "--caller",
        choices=["auto", "sniffles2", "sniffles1", "cutesv",
                 "svim", "pbsv", "nanovar"],
        default="auto",
        help="SV caller that produced the VCF. 'auto' (default) runs "
             "the fingerprint detector; falls back to sniffles2 on "
             "no match. Set explicitly to override.",
    )


def _add_compound_het_args(p) -> None:
    """Compound-het flags. Used when the VCF carries CSQ + PS headers."""
    p.add_argument("--gene", action="append", default=None, metavar="SYMBOL",
                   help="plot exactly this gene; repeatable. Plots "
                        "regardless of variant count (with clear "
                        "messaging for empty cases). When omitted, the "
                        "auto-select rule picks candidate genes.")
    p.add_argument("--clinvar", type=Path, default=None,
                   help="override the bundled ClinVar lookup. Accepts "
                        "either molamola's reduced TSV (~13 MB; default "
                        "at data/clinvar.hg38.tsv.xz) or NCBI's raw "
                        "ClinVar VCF (190+ MB). Format auto-detected.")
    p.add_argument("--canonical-exons", type=Path, default=None,
                   help="override the bundled canonical-exon TSV "
                        "(default: data/canonical_exons.hg38.tsv.gz)")
    p.add_argument("--min-pair-count", type=int, default=1,
                   help="auto-select threshold: gene needs >= N trans "
                        "pairs in a single phase set, where one anchor "
                        "is ClinVar P/LP or VUS and the partner is not "
                        "benign. Ignored when --gene is given. Default 1.")
    p.add_argument("--max-genes", type=int, default=50,
                   help="cap number of auto-selected genes. Default 50.")


def build_argparser() -> argparse.ArgumentParser:
    """Construct the molamola CLI parser.

    Single flat parser; no subcommands. The plot type is chosen by
    inspecting the VCF header at run time (see :func:`detect_vcf_mode`),
    so the user only needs ``--vcf``. Mode-specific flag groups are
    visually separated in ``--help``.
    """
    p = argparse.ArgumentParser(
        prog="molamola",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p)
    sv_group = p.add_argument_group(
        "SV-mode flags (used when input VCF has ##INFO=<ID=SVTYPE,...>)"
    )
    _add_sv_args(sv_group)
    ch_group = p.add_argument_group(
        "Compound-het mode flags (used when input VCF has "
        "##INFO=<ID=CSQ,...> + ##FORMAT=<ID=PS,...>)"
    )
    _add_compound_het_args(ch_group)
    return p


def plot_main(args: argparse.Namespace) -> int:
    """Execute the `plot` subcommand."""
    hint = detect_reference_hint(args.vcf.name)
    if hint is not None and hint != args.reference:
        msg = (f"VCF filename '{args.vcf.name}' suggests --reference "
               f"{hint} but --reference {args.reference} was given")
        if args.force:
            print(f"  WARNING: {msg}. Continuing because of --force.",
                  file=sys.stderr)
        else:
            print(f"ERROR: {msg}. Use --force to override.", file=sys.stderr)
            return 2

    out_dir = args.out if args.out is not None else args.vcf.resolve().parent
    out_dir.mkdir(parents=True, exist_ok=True)
    sample = args.sample or args.vcf.name.replace(".vcf.gz", "").replace(".vcf", "")

    detected, basis = detect_caller(args.vcf)
    if args.caller == "auto":
        caller = "sniffles2" if detected == "unknown" else detected
        print(f"Detected caller: {caller}  ({basis}; --caller auto)")
    else:
        caller = args.caller
        if detected != "unknown" and detected != caller:
            print(
                f"  WARNING: --caller {caller} but fingerprint suggests "
                f"{detected} ({basis}). Continuing with the explicit choice.",
                file=sys.stderr,
            )
        print(f"Caller: {caller}  (--caller explicit; fingerprint: {detected})")

    contigs, bnds_all, svs, median_cov = read_vcf(args.vcf, caller=caller)
    if not bnds_all and not svs:
        print("no SV records found", file=sys.stderr)
        return 1

    # Sanity check: VCF contigs must intersect the canonical chr1-22/X/Y set
    # AFTER normalization. If none do, the VCF is on an assembly molamola
    # doesn't support (anything other than hg38 / T2T-CHM13v2.0).
    canonical = [c for c in contigs if c in CHROM_SET]
    if not canonical:
        sample_names = ", ".join(list(contigs)[:6]) or "(none in header)"
        print(
            f"ERROR: no contigs in {args.vcf.name} match the canonical "
            f"chromosome set (chr1-22, chrX, chrY). VCF contigs sampled: "
            f"{sample_names}. molamola supports hg38 and T2T-CHM13v2.0; "
            f"non-UCSC contig naming ('1' instead of 'chr1') is auto-"
            f"normalized, so seeing this message means the underlying "
            f"assembly is unsupported (e.g. GRCh37/hg19).",
            file=sys.stderr,
        )
        return 1

    # Hard size filter on non-BND SVs (--min-svlen). BNDs are unaffected
    # (no SVLEN). The filter applies once here so every downstream
    # consumer (cov-anomaly stats, density tracks, HTML table, stdout
    # summaries) sees the same set.
    if args.min_svlen > 0:
        before = len(svs)
        svs = [s for s in svs if s.svlen >= args.min_svlen]
        dropped = before - len(svs)
        if dropped:
            print(f"--min-svlen {args.min_svlen}: dropped {dropped:,} non-BND "
                  f"SVs shorter than {args.min_svlen} bp "
                  f"({before:,} -> {len(svs):,})")

    n_total = len(bnds_all)
    n_pass = sum(1 for b in bnds_all if b.is_pass)
    keep_bnds = [b for b in bnds_all if b.is_pass] if args.filter == "pass" else bnds_all
    unique_bnds = deduplicate_reciprocal(keep_bnds)

    cytoband_file = find_cytoband_file(args.reference)
    cytobands = load_cytobands(cytoband_file)

    # Assembly sanity check: declared VCF contig lengths should agree
    # with the cytoband's per-chromosome max-end (which is the assembly
    # length). hg19/hg38/T2T differ by tens-to-hundreds of kb on most
    # chromosomes; a 100 kb absolute tolerance accommodates patch
    # versions while still catching wrong-assembly cases.
    cyto_lens = {c: max(end for _s, end, _n, _st in bands)
                 for c, bands in cytobands.items() if c in CHROM_SET}
    mismatched: list[tuple[str, int, int]] = []
    for c in canonical:
        v_len = contigs.get(c)
        ref_len = cyto_lens.get(c)
        if v_len is None or ref_len is None:
            continue
        if abs(v_len - ref_len) > 100_000:
            mismatched.append((c, v_len, ref_len))
    if mismatched and len(mismatched) >= len(canonical) // 2:
        head = "; ".join(
            f"{c}: VCF={v:,} vs {args.reference}={r:,}"
            for c, v, r in mismatched[:3]
        )
        print(
            f"ERROR: contig lengths in {args.vcf.name} disagree with "
            f"the bundled {args.reference} cytoband on "
            f"{len(mismatched)}/{len(canonical)} canonical chromosomes "
            f"({head}{'; ...' if len(mismatched) > 3 else ''}). "
            f"Looks like the VCF was called against a different "
            f"assembly. molamola supports hg38 and T2T-CHM13v2.0 only — "
            f"refusing to plot rather than render misleading coordinates.",
            file=sys.stderr,
        )
        return 1

    # Acrocentric default flips with reference: on for hg38 (where those
    # p-arms are N-padded mappability sinks producing phantom
    # translocations), off for T2T (where the p-arms are real
    # fully-resolved sequence). Explicit --mark-acrocentric / --no-...
    # always wins (default=None means user did not set it).
    mark_acrocentric = args.mark_acrocentric
    if mark_acrocentric is None:
        mark_acrocentric = args.reference == "hg38"

    if args.cov_ratio == "auto":
        cov_ratio_thr, n_used = auto_cov_ratio_threshold(
            unique_bnds, svs, median_cov,
        )
        print(f"--cov-ratio auto: {cov_ratio_thr:.2f}x  "
              f"(p99 of {n_used:,} PASS events; floor 2.0x)")
    else:
        cov_ratio_thr = args.cov_ratio

    annotate_noise(
        unique_bnds, cytobands, median_cov,
        mark_acrocentric=mark_acrocentric,
        cov_ratio_thr=cov_ratio_thr,
        cov_vaf_max=args.cov_vaf_max,
    )
    if args.cov_filter != "none":
        annotate_sv_noise(
            svs, median_cov,
            cov_ratio_thr=cov_ratio_thr,
            cov_vaf_max=args.cov_vaf_max,
        )

    if args.cov_filter == "drop":
        before = len(unique_bnds)
        unique_bnds = [b for b in unique_bnds if "cov_anomaly" not in b.noise_flags]
        print(f"--cov-filter=drop removed {before - len(unique_bnds)} BNDs")
    if args.cov_filter == "none":
        for b in unique_bnds:
            b.noise_flags.discard("cov_anomaly")

    focus_tag = ""
    if args.focus:
        # Validate band foci early — typos should fail fast.
        for fc, key in args.focus:
            if isinstance(key, str) and resolve_band_range(fc, key, cytobands) is None:
                print(f"ERROR: --focus {fc}:{key} did not match any cytoband "
                      f"on {fc}. Check the band name (e.g. 'q11.23').",
                      file=sys.stderr)
                return 2
        before = len(unique_bnds)
        unique_bnds = [b for b in unique_bnds
                       if matches_focus(b, args.focus, args.focus_window, cytobands)]
        print(f"--focus filtered: {before} -> {len(unique_bnds)} BNDs "
              f"(window +/-{args.focus_window} bp for positions; "
              f"band-range overlap for ISCN bands)")
        for b in unique_bnds:
            iscn = iscn_label(b, cytobands)
            print(f"  match: {b.sv_id}  {b.chr1}:{b.pos1}  <->  {b.chr2}:{b.pos2}  "
                  f"ISCN={iscn}  "
                  f"FILTER={b.filter_}  SUPPORT={b.support}  VAF={b.vaf:.3f}  "
                  f"noise={'+'.join(sorted(b.noise_flags)) or 'none'}")
        if not unique_bnds:
            print("no BNDs matched --focus; nothing to plot", file=sys.stderr)
            return 2
        if len(args.focus) == 1:
            c, p_ = args.focus[0]
            focus_tag = f".focus_{c}_{p_}"
        else:
            focus_tag = f".focus_{len(args.focus)}sites"

    bd = noise_breakdown(unique_bnds)
    sv_bd = sv_noise_breakdown(svs)
    print(f"BND records: {n_total} total, {n_pass} PASS")
    print(f"After --filter={args.filter}: {len(keep_bnds)} kept")
    print(f"Unique BND events after reciprocal dedupe: {len(unique_bnds)}")
    def _sv_summary_chunk(t: str) -> str:
        s = sv_bd[t]
        if not s["cov"]:
            return f"{t}={s['pass']}"
        return f"{t}={s['pass']} (clean={s['clean']}, cov-noise={s['cov']})"

    print("Non-BND SVs (PASS):  " + "  ".join(_sv_summary_chunk(t) for t in SV_TYPES))
    print(f"Genome median coverage (from INS+DEL center): {median_cov:.1f}x  "
          f"(cov-anomaly threshold: max>={cov_ratio_thr*median_cov:.0f}x AND "
          f"VAF<{args.cov_vaf_max})")
    print(f"Cytoband file: {cytoband_file}")
    print(f"BND noise breakdown: clean={bd['clean']}, "
          f"acrocentric={bd['acrocentric']}, cov_anomaly={bd['cov_anomaly']}, "
          f"any_noise={bd['any_noise']}")

    # Render both figures into in-memory PNG buffers — no temp files.
    circos_buf = io.BytesIO()
    plot_circos(unique_bnds, contigs, cytoband_file,
                circos_buf, sample, n_total, n_pass, args.filter, bd)
    sv_map_buf = io.BytesIO()
    plot_genome_sv_map(unique_bnds, svs, contigs, cytobands,
                        sv_map_buf, sample, n_total, n_pass, args.filter, bd,
                        bin_size=args.bin_size, min_svlen=args.min_svlen)

    out_html = out_dir / f"{sample}{focus_tag}.report.html"
    make_html_report(
        sample=sample,
        circos_png=circos_buf.getvalue(),
        sv_map_png=sv_map_buf.getvalue(),
        out_path=out_html,
        bd=bd,
        sv_bd=sv_bd,
        median_cov=median_cov,
        cov_ratio=cov_ratio_thr,
        cov_vaf_max=args.cov_vaf_max,
        n_total_bnd=n_total,
        n_pass_bnd=n_pass,
        filter_label=args.filter,
        caller=caller,
        caller_basis=basis,
    )
    print(f"wrote {out_html}")
    if args.png:
        base = out_html.stem
        circos_path = out_dir / f"{base}.circos.png"
        sv_map_path = out_dir / f"{base}.sv_map.png"
        circos_path.write_bytes(circos_buf.getvalue())
        sv_map_path.write_bytes(sv_map_buf.getvalue())
        print(f"wrote {circos_path}")
        print(f"wrote {sv_map_path}")
    return 0


def compound_het_main(args: argparse.Namespace) -> int:
    """Execute the compound-het render path.

    Reads a phased small-variant VCF with VEP CSQ annotation, loads
    bundled canonical-exon and ClinVar references, picks gene(s) to
    plot (explicit ``--gene`` or auto-select), renders one panel per
    gene, and assembles a self-contained HTML report.

    Compound-het mode is hg38-only: ClinVar is hg38-coordinate and
    coordinate-based lookup needs the input VCF on the same build.
    T2T support will land once a non-coordinate matching path is in
    place.

    Returns the process exit code:
    - ``0`` on success (including "0 candidate genes" — empty result
      is a valid finding).
    - ``1`` on hard refusals (T2T input, unphased VCF, missing CSQ,
      unknown gene symbol, missing bundled reference, etc).
    - ``2`` on reference-hint mismatch without ``--force``.
    """
    if args.reference != "hg38":
        print(f"ERROR: compound-het mode is hg38-only "
              f"(--reference {args.reference} not supported). "
              f"ClinVar coordinates are hg38; T2T support will land "
              f"once a non-coordinate matching path is in place.",
              file=sys.stderr)
        return 1

    hint = detect_reference_hint(args.vcf.name)
    if hint is not None and hint != args.reference:
        msg = (
            f"VCF filename {args.vcf.name!r} suggests --reference "
            f"{hint} but --reference {args.reference} was given"
        )
        if args.force:
            print(f"  WARNING: {msg}. Continuing because of --force.",
                  file=sys.stderr)
        else:
            print(f"ERROR: {msg}. Use --force to override.",
                  file=sys.stderr)
            return 2

    out_dir = args.out if args.out is not None else args.vcf.resolve().parent
    out_dir.mkdir(parents=True, exist_ok=True)
    sample = (args.sample
              or args.vcf.name.replace(".vcf.gz", "").replace(".vcf", ""))

    # Load canonical-exon table (bundled or override).
    if args.canonical_exons is not None:
        exon_path = args.canonical_exons
        exon_source = str(exon_path)
    else:
        try:
            exon_path = find_canonical_exon_file()
        except FileNotFoundError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 1
        exon_source = f"(bundled) {exon_path.name}"
    try:
        genes = load_canonical_exons(exon_path)
    except (ValueError, OSError) as e:
        print(f"ERROR: failed to load canonical exons: {e}",
              file=sys.stderr)
        return 1
    print(f"Loaded {len(genes)} canonical-transcript genes from "
          f"{exon_path.name}")

    # Read phased VCF.
    try:
        variants, parse_summary = read_phased_vcf(args.vcf)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    print(f"Phased VCF: kept {parse_summary['kept']} records, "
          f"refused unphased={parse_summary['unphased']}, "
          f"non_het={parse_summary['non_het']}, "
          f"multi_allelic={parse_summary['multi_allelic']}, "
          f"no_canonical_csq={parse_summary['no_canonical_csq']}")
    if not variants:
        print("ERROR: VCF parsed cleanly but yielded no phased het "
              "variants; nothing to plot", file=sys.stderr)
        return 1

    # Load ClinVar (bundled or override) and annotate by chrom+pos+ref+alt.
    # Compound-het mode is hg38-only as of v0.3.1 (T2T support deferred
    # until a non-coordinate-based matching path lands; the reference
    # check above already refused t2t).
    # Loaders return buckets directly; no canon_clnsig pass needed here.
    clinvar_path: Path | None
    if args.clinvar is not None:
        clinvar_path = args.clinvar
        clinvar_source = str(clinvar_path)
    else:
        try:
            clinvar_path = find_clinvar_file()
            clinvar_source = f"(bundled) {clinvar_path.name}"
        except FileNotFoundError:
            clinvar_path = None
            clinvar_source = "(none — all lollipops grey)"

    if clinvar_path is None:
        annotated = list(variants)
    else:
        keys = {(v.chrom, v.pos, v.ref, v.alt) for v in variants}
        try:
            lookup = load_clinvar_lookup(
                clinvar_path, keys_of_interest=keys,
            )
        except (OSError, ValueError) as e:
            print(f"ERROR: failed to load ClinVar source: {e}",
                  file=sys.stderr)
            return 1
        annotated = [
            v.with_clnsig(lookup.get((v.chrom, v.pos, v.ref, v.alt)))
            for v in variants
        ]
    print(f"ClinVar annotation: {sum(1 for v in annotated if v.clnsig)} "
          f"of {len(annotated)} variants matched")

    # Determine gene set to plot.
    empty_genes: list[str] = []
    strict_set: set[str] = set()
    if args.gene:
        unknown = [s for s in args.gene if s not in genes]
        if unknown:
            print(f"ERROR: gene(s) not found in {args.reference} "
                  f"canonical exon table: {', '.join(unknown)}",
                  file=sys.stderr)
            return 1
        plot_symbols = list(args.gene)
        selection_rule = (
            f"explicit --gene: {', '.join(plot_symbols)}"
        )
    else:
        strict, extended_only = find_compound_het_candidates(
            annotated, genes,
            min_pair_count=args.min_pair_count,
        )
        strict_set = set(strict)
        plot_symbols = strict + extended_only
        selection_rule = (
            f"auto-select: >= {args.min_pair_count} trans pair(s) in "
            f"same phase set with anchor ClinVar P/LP or VUS, partner "
            f"not benign. Strict subset: pairs where BOTH variants "
            f"are P/LP or VUS."
        )

    n_genes_scanned = len(plot_symbols)
    n_genes_capped = 0
    if not args.gene and n_genes_scanned > args.max_genes:
        n_genes_capped = n_genes_scanned - args.max_genes
        plot_symbols = plot_symbols[:args.max_genes]
        print(f"  WARNING: {n_genes_scanned} candidate genes; capping "
              f"at {args.max_genes}; pass --max-genes higher to plot all",
              file=sys.stderr)

    if not plot_symbols and not args.gene:
        print(f"0 candidate genes (rule: {selection_rule})",
              file=sys.stderr)
        # Still emit a report so the strict-section heading is visible
        # even on an empty result.

    # Render each gene.
    gene_panels: list[tuple[Gene, bytes, dict]] = []
    for symbol in plot_symbols:
        gene = genes[symbol]
        gene_vars = [
            v for v in annotated
            if v.chrom == gene.chrom and gene.start <= v.pos <= gene.end
        ]
        if not gene_vars:
            print(f"{symbol}: 0 phased hets in {gene.chrom}:"
                  f"{gene.start}-{gene.end}; no plot produced",
                  file=sys.stderr)
            empty_genes.append(symbol)
            continue
        blocks = _build_phase_blocks(gene_vars)
        png, stats = render_compound_het_png(
            gene, gene_vars, blocks, reference=args.reference,
        )
        gene_panels.append((gene, png, stats))

    n_genes_plotted = len(gene_panels)
    out_html = out_dir / f"{sample}.compound_het.report.html"
    make_compound_het_report(
        sample=sample,
        gene_panels=gene_panels,
        out_path=out_html,
        reference=args.reference,
        n_genes_scanned=n_genes_scanned,
        n_genes_plotted=n_genes_plotted,
        n_genes_capped=n_genes_capped,
        clinvar_source=clinvar_source,
        canonical_exons_source=exon_source,
        refusal_counts=parse_summary,
        selection_rule=selection_rule,
        empty_genes=empty_genes,
        strict_symbols=strict_set,
        is_auto_select=not args.gene,
    )
    print(f"wrote {out_html}")
    if args.png and gene_panels:
        base = out_html.stem
        for gene, png, _stats in gene_panels:
            png_path = out_dir / f"{base}.{gene.symbol}.png"
            png_path.write_bytes(png)
            print(f"wrote {png_path}")
    return 0


def karyotype_main(args: argparse.Namespace) -> int:
    """Karyotype coverage mode entry point (mosdepth-driven).

    Activated when --mosdepth is given. When --vcf is also given, the
    VCF is consumed as the BAF source for the karyotype panel rather
    than running through SV / compound-het header dispatch.
    """
    print("[info] karyotype coverage mode")
    if args.vcf is not None:
        print(f"[info] BAF source: {args.vcf}")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Top-level CLI entry point. Returns a process exit code.

    Either ``--vcf`` or ``--mosdepth`` (or both) is required. When
    ``--mosdepth`` is set, karyotype coverage mode runs and any
    accompanying ``--vcf`` is consumed only as the BAF source.
    Otherwise the plot mode is auto-detected from the VCF header
    (SV vs. compound-het). molamola refuses cleanly rather than
    rendering a misleading default.
    """
    args = build_argparser().parse_args(argv)
    if args.mosdepth is None and args.vcf is None:
        print("ERROR: either --vcf or --mosdepth is required",
              file=sys.stderr)
        return 2
    if args.mosdepth is not None:
        return karyotype_main(args)
    try:
        mode = detect_vcf_mode(args.vcf)
    except (ValueError, FileNotFoundError, OSError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    if mode == "compound-het":
        return compound_het_main(args)
    return plot_main(args)


if __name__ == "__main__":
    sys.exit(main())
