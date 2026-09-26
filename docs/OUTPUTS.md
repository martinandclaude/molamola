# Outputs

Each run writes a single self-contained HTML artefact into `--out` (required — molamola refuses rather than silently writing next to the input file; the directory is created if absent). The filename and contents differ by mode:

- **SV mode** → `<sample>.report.html` (the circos plot, then one derivative-chromosome panel per candidate rearrangement, embedded as base64 PNG data URIs, plus a collapsible run-metadata block).
- **Karyotype mode** → `<sample>.karyotype.report.html` (one genome-wide CN figure, with a BAF panel beneath when `--vcf` is supplied; same metadata block).

With `--png` (added in v0.2.0 for MultiQC / pipeline embeds), each embedded figure is *also* written as a standalone PNG alongside the HTML — SV mode writes `<sample>.report.circos.png` and one `<sample>.report.rearrangement_<n>.png` per panel, karyotype mode writes `<sample>.karyotype.genome.png`. Without `--png` the figures live only inside the HTML; right-click → save to extract one ad hoc.

## Reference data

All bundled in `molamola/data/`:

- `cytoBand.txt.gz` (hg38) and `cytoBand.t2t.txt.gz` (T2T-CHM13v2.0) — UCSC cytoband annotations for SV and karyotype modes.
- `exclusion.hg38.bed.gz`, `exclusion.t2t.bed.gz` — karyotype-mode exclusion masks (~5.7 MB / ~6 MB; low-mappability ∪ polymorphic-TR catalog), built from 500 bp runs. A bin is dropped when more than `--mask-overlap` of it is masked. Override with `--mask`, disable with `--no-mask`.
- `gc_10kb.hg38.bed.gz`, `gc_10kb.t2t.bed.gz` — karyotype-mode 10 kb GC tables (~1.5 MB each) driving the per-1 % GC-bucket median-ratio correction. Override with `--gc`, disable with `--no-gc`.

- `genes.hg38.bed.gz`, `genes.t2t.bed.gz` — SV-mode gene tables (~275 kB each) for breakpoint and fusion labels: ~20,000 protein-coding genes plus the IG / TR loci (IGH, IGK, IGL, TRA/TRD, TRB, TRG) as one span each. From NCBI RefSeq annotation release RS_2025_08, which annotates GRCh38 and T2T-CHM13v2.0 natively in one run: the same gene symbols on both builds, and T2T coordinates that are the annotation's own rather than a liftover (they match NASVAR's T2T configuration).

molamola is bundled-only by design — no auto-download, no online lookups. The bundled-data total is ~15.5 MB. All bundled refs are reproducibly regeneratable from public sources by running `scripts/derive_karyotype_refs.py` and `scripts/derive_gene_tables.py`.

## Circos plot (SV mode)

The SV report's one figure. From the rim inward:

- **Cytoband ideogram** — molamola's greyscale ISCN ramp, with the centromere in red (black is not separable from `gpos100` on a ring five radial units thick).
- **Four SV density rings**, outermost first: INS = blue, DEL = red, DUP = green, INV = purple. Each is 1 Mb-binned (`--bin-size`); alpha encodes count, sqrt-scaled and saturating at the 99th-percentile occupied bin rather than at the peak. Anchoring on the peak looked principled and was not — SV density is heavily skewed, so on a normal genome the peak bin is a far outlier from the typical one (INS: median 3 events per bin, peak 67), which pushed almost every bin into a 0.15-wide alpha band and made the track a near-uniform wash. The ~1 % of bins above the anchor clip to full alpha; single-event bins stay visible via a 0.20 alpha floor. Only occupied bins are drawn — an empty ring is bare paper, not a tinted band. Noise-flagged and non-PASS events are left out of the rings.
  - The rings are deliberately **not** equal width: INS and DEL are drawn thin and capped below full alpha, DUP and INV get most of the radius. Per-1 Mb-bin INS and DEL density is ~85 % identical between two unrelated people (measured on the GIAB trio, holding everything but relatedness constant), so those tracks describe the species more than the sample and are drawn as background texture.
  - Type is encoded by *radius* first and colour second: the four-colour palette is not separable under protanopia or deuteranopia on its own, and giving each type its own ring keeps colour redundant rather than load-bearing.
- **Rearrangement arcs** across the disc. BND junctions are first grouped into events (see [Rearrangement tiers](FILTERS.md#rearrangement-tiers-sv-mode)), and each event is **one** arc — a reciprocal translocation's two junctions are not drawn twice. Colour is VAF class, three discrete bands: 0-33 %, 33-66 % and 66-100 %. They are named by range only: germline names (mosaic / het / hom) would misread tumour samples, where VAF follows clone fraction. The classes are blue, burnt orange and near-black. The colours are chosen under simulated colour-vision deficiency, judged as drawn on the page, and stay well apart under protanopia, deuteranopia and tritanopia. Lightness decreases as VAF rises, so the ordering survives even total loss of hue discrimination, and the fullest events are the darkest. Width scales with read support. Emphasis follows the tier:
  - **candidate** (both junctions found, outside repeats) — widest, drawn on top, and numbered: a badge at each breakpoint matches a numbered name beside the disc, e.g. `1  t(8;21)(q21.3;q22.12)`. Badges that would collide step inward; the two ends of a short event share one badge.
  - **paired, but in repeats** — thinner than a candidate but just as opaque: the VAF colours hold their contrast and colour-blind separation only at that opacity, and a real fusion demoted by the mask rule lands here, so it has to stay fully readable.
  - **single junction** — faint background. On a normal genome this is ~90 % of all arcs.
  - **noise-flagged** — grey and fainter still, whatever the tier. **non-PASS** (`--filter all`) — half the opacity of its tier. Arcs are filled ribbons, so there are no dashed styles.
  - Large inversions, from `INV` records or a pair of intrachromosomal BNDs, are drawn as arcs within their chromosome rather than as a tick in the INV ring. An insertion's ribbon spans the whole donor segment. With `--plotvaf`, each drawn arc also carries its VAF as a percentage on the rim.

Keys sit to the right of the disc. At the top, side by side: the density rings with their per-type event counts (events on the drawn chromosomes only), and the arc tiers. Below them, a small horizontal VAF bar; below that, the numbered candidate list, which has the rest of the column and names up to 25 candidates ("none" when there are none — the absence is stated, not implied). There is no cytoband key: the greyscale is the ISCN convention its readers already know. The ring key deliberately quotes no number for the alpha scale — the ramp saturates at the 99th-percentile bin, which its title says, and printing the busiest bin beside that only invited the reader to treat the wrong number as the top of the scale.

`--only-sv-chroms` draws only the chromosomes that carry an arc — for targeted or adaptive-sampling runs where most of a genome-wide disc would be empty. The run metadata block also states the rearrangement summary line.

Up to v0.6 the report also carried a linear genome map (one row per chromosome with density strips and BND arcs). It was removed: with ~100 arcs on a normal genome crossing 24 rows, an arc's endpoints could not be read, and the circos carries the same density signal.

## Derivative-chromosome panels (SV mode)

One panel per candidate rearrangement, below the circos and numbered as on it (up to 12; any further candidates are counted). Each shows the chromosomes as a cytogeneticist reads a karyotype, G-banded, pter at the top:

- **The normal chromosomes and the derivatives** the junctions build — `8`, `der(8)`, `21`, `der(21)` for a translocation; `16`, `inv(16)` for an inversion; the recipient and donor chromosomes with their derivatives for an insertion. Derivatives are constructed from which side of each breakpoint every junction keeps, oriented so the centromere-bearing chromosome reads pter to qter, and named by the centromere they carry (`der`, `dic` for two, `ace` for none). A paint bar beside each derivative shows which chromosome each segment came from; an arrow marks an inverted segment.
- **Breakpoints** on the normal chromosomes, labelled with their band and the gene they fall in — or the nearest gene within 500 kb (`near MYC (157 kb)`), since enhancer-hijacking breakpoints often sit well outside the gene they drive. A breakpoint inside an IG / TR locus is named after the locus.
- **Junctions** on the derivatives, labelled with the fusion they make, written 5'::3' (`RUNX1::RUNX1T1` on der(8), `RUNX1T1::RUNX1` on der(21)). The orientation follows from the kept sides and the gene strands; a junction that puts two genes head to head or tail to tail makes no fusion and is not labelled with one. IG / TR loci are written first (`IGH::CRLF2`). Gene names are annotation only: molamola does not say whether a fusion is recurrent or significant.
- **The evidence**: junction count (or "1 INV record" when a single record implies both junctions), read support and VAF per junction, and exact coordinates. A candidate promoted by the gene rule (see [Rearrangement tiers](FILTERS.md#rearrangement-tiers-sv-mode)) says so.

![example panel](example_rearrangement_panel.png)

*Synthetic t(8;21) breakpoints placed in RUNX1T1 and RUNX1 — an illustration, not a sample.*

## Karyotype coverage figure (karyotype mode)

One genome-wide figure, chr1 → chrY left to right:

- **log2 relative-depth scatter** — per-bin `log2(depth / autosomal median)` after masking + GC correction, so CN 2 sits at 0 and a single-copy loss and gain read symmetrically. CN 1 / 2 / 3 reference lines are labelled at the right edge. Adjacent chromosomes alternate between two inks and every other chromosome carries a background band, so a deviation can be attributed without tracing back to the axis; chrX and chrY have their own inks and a legend. Aggregated into ~`--scatter-bin-kb` windows; systematically downsampled above `--max-points`.
- **Rolling-median smooth** — a deep-pink per-chromosome line (`--smooth-window-mb` window, default 10 Mb), the load-bearing signal for arm-scale CN events.
- **Expected-CN dashes** — Oxford-blue horizontal guides at the expected CN for each chromosome given the inferred (or `--sex`-forced) genomic sex.
- **BAF panel** (only when `--vcf` is given) — allele balance at PASS heterozygous SNVs beneath the CN panel, sharing the x-axis. Reference lines sit at 33 / 50 / 67 % — the CN 3 het expectations — so a gain reads as the cloud splitting onto the outer two lines rather than as a vague widening; deletions / LOH splay toward 0 and 100 %.
  - **Phased VCF** (`FORMAT/PS` + `FORMAT/AD` present): haplotype-resolved. Read counts are summed within each phase block over tiling windows of 40 het SNVs — summing reads rather than averaging per-site fractions is the variance-correct estimator — and each window is plotted twice, at *v* and *1−v*, because which haplotype a block labels "1" is arbitrary and flips between blocks. A systematic reference-mapping bias is absorbed by shifting the autosomal median onto 50 %. On a normal ONT genome this tightens the panel about 1.4× versus per-site (SD 0.058 vs 0.079), which is what puts the CN 3 expectation ~2.9 SD clear of balanced.
  - **Unphased VCF**: per-site het allele fraction, as before. Phase is an upgrade, never a requirement; `--no-phased-baf` forces the per-site panel. The mode used is recorded in the report's metadata block.
- **Run metadata** lives in the HTML report's collapsible block, not on the figure (an in-figure strip crowded the arm labels and was removed): mosdepth source, reference, inferred sex, bin sizes, smooth window, mask / GC labels, BAF source, and an **"AS suspected"** chip when the autosomal depth distribution looks bimodal (adaptive-sampling-like — CN normalisation may be biased on such samples; the CLI also prints a one-line warning). A standalone `--png` therefore carries no provenance line.

The per-chromosome 3 × 8 A4-portrait grid that briefly shipped in v0.3 development was dropped before release — the genome-wide figure does the cytogenetics work on its own.

## HTML report

A single self-contained HTML file (no server, no external CSS/JS, opens offline). Figures are embedded as base64 PNG data URIs. A collapsible run-metadata block at the bottom shows mode-specific provenance:

- **SV mode**: caller, filter mode, baseline coverage, cov-anomaly threshold, BND noise breakdown, and the rearrangement summary.
- **Karyotype mode**: mosdepth source, reference, inferred genomic sex, mosdepth bin + scatter bin + smooth window, mask / GC source labels, BAF source + het-site count (when `--vcf` given), and the AS-suspected flag.

Designed as a portable review surface — emailable, archivable, no install needed for the reader.
