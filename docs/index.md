<p align="center"><img src="header_fish.png" alt="molamola" width="128"></p>

# molamola

A Python cytogenetics plotting tool for Oxford Nanopore data.
**One input in, one self-contained HTML report out.**

molamola picks the right report from its input — SV mode is
auto-detected from the VCF header, karyotype mode is selected by
the `--mosdepth` flag:

- **SV / cytogenetics report** for long-read SV VCFs
  (Sniffles2, cuteSV, SVIM, pbsv, NanoVar).
- **Karyotype coverage report** for a mosdepth `regions.bed.gz`
  (via `--mosdepth`) — genome-wide log2 relative-depth scatter + rolling-median
  smooth, with an optional BAF panel beneath when paired with a
  small-variant VCF.

## Install

```sh
pip install molamola
```

Or via conda — note that both bioconda and conda-forge channels are
needed (pycirclize lives on conda-forge):

```sh
conda create -n molamola -c bioconda -c conda-forge molamola
```

## Quick start

```sh
molamola --vcf sample.vcf --out reports/
# or, for a coverage karyotype:
molamola --mosdepth sample.regions.bed.gz --reference hg38 --out reports/
```

The plot type is picked from the input (VCF header, or `--mosdepth`
for karyotype mode). `--out` is required — molamola refuses rather
than silently writing the report next to the input file. Output is
a single self-contained HTML report — figures embedded as base64,
no external assets, opens offline.

## Example output

Figures below come from running molamola's SV mode on sample MH001
(ONT LSK114 library prep, aligned-read N50 10.4 kb, median autosomal
coverage 54x). The HTML report embeds both plots back-to-back; shown
separately here for clarity.

### Circos plot

![SV circos plot](example_sv_circos.png)

22 autosomes plus X and Y arranged around the disc, with greyscale
ISCN cytobands on the rim and a red centromere. Inside the rim sit
four 1-Mb-bin SV density rings -- INS, DEL, DUP, INV, outermost
first -- so located events read on the rings while connections read
across the disc. Each ribbon across the disc is a BND (translocation
or large rearrangement); ribbon colour encodes VAF as one of three
discrete classes -- mosaic (0-33 %), het (33-66 %), hom (66-100 %) --
rather than a continuous ramp, so the bands are separable at a
glance. At-a-glance view for inter-chromosomal events.

### Linear genome map

![SV linear plot](example_sv_linear.png)

One row per chromosome (chr1 at top, chrY at bottom). Cytobands
embedded inside each chromosome track. Above each track sit four
1-Mb-bin density strips — INS = blue, DEL = red, DUP = green,
INV = purple — with alpha encoding per-bin event count on a scale
that saturates at the 99th-percentile bin, identical to the circos
rings. BND arcs hang above the tracks, colour-encoded by
VAF as in the circos. Better for per-chromosome detail and density
hotspots, and for reading a single chromosome end to end.

## Bundled references

molamola ships its own reference data inside `molamola/data/`:

- `cytoBand.txt.gz` (hg38) and `cytoBand.t2t.txt.gz`
  (T2T-CHM13v2.0) — UCSC cytoband annotations, used for the SV
  ideogram tracks.
- `exclusion.hg38.bed.gz`, `exclusion.t2t.bed.gz` — karyotype-mode
  exclusion masks (low-mappability ∪ polymorphic-TR catalog),
  built from 500 bp runs.
- `gc_10kb.hg38.bed.gz`, `gc_10kb.t2t.bed.gz` — 10 kb GC tables
  driving karyotype mode's per-1 % GC-bucket median-ratio
  correction.

Bundled-only by design: molamola does not auto-download or look up
online. The karyotype refs are reproducibly regeneratable from
public sources via `scripts/derive_karyotype_refs.py` in the repo.

> **Compound-het mode was removed after v0.5.1.** The per-gene
> phased-haplotype panels for recessive-disease workup, and the
> bundled ClinVar and MANE Select references they needed, are gone.
> molamola is a cytogenetics tool now. Install `molamola==0.5.1`
> if you need them.

## Documentation

- [CLI reference](CLI.md) — every flag, with defaults and meanings.
- [Examples](EXAMPLES.md) — worked commands for each report type.
- [Filters](FILTERS.md) — noise heuristics, focus windows, mask and GC handling.
- [Output spec](OUTPUTS.md) — what's in the HTML; how figures are encoded.

## Source

- Repo: [github.com/martinandclaude/molamola](https://github.com/martinandclaude/molamola)
- PyPI: [pypi.org/project/molamola](https://pypi.org/project/molamola/)
- Changelog: [CHANGELOG.md](https://github.com/martinandclaude/molamola/blob/main/CHANGELOG.md)
- Issues: [github.com/martinandclaude/molamola/issues](https://github.com/martinandclaude/molamola/issues)
