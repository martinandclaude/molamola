# Outputs

Each run writes a single self-contained HTML artefact into `--out` (required — molamola refuses rather than silently writing next to the input file; the directory is created if absent). The filename and contents differ by mode:

- **SV mode** → `<sample>.report.html` (circos plot + linear cytoband SV map embedded as base64 PNG data URIs, plus a collapsible run-metadata block).
- **Compound-het mode** → `<sample>.compound_het.report.html` (one phased-haplotype panel per gene, embedded the same way; auto-select runs split into a `strict` and `extended` section).
- **Karyotype mode** → `<sample>.karyotype.report.html` (one genome-wide CN figure, with a BAF panel beneath when `--vcf` is supplied; same metadata block).

With `--png` (added in v0.2.0 for MultiQC / pipeline embeds), each embedded figure is *also* written as a standalone PNG alongside the HTML — SV mode writes the circos + SV-map PNGs, compound-het writes one PNG per gene (gene symbol baked into the title, top-left, bold, so a saved image is self-describing), karyotype mode writes `<sample>.karyotype.genome.png`. Without `--png` the figures live only inside the HTML; right-click → save to extract one ad hoc.

## Reference data

All bundled in `molamola/data/`:

- `cytoBand.txt.gz` (hg38) and `cytoBand.t2t.txt.gz` (T2T-CHM13v2.0) — UCSC cytoband annotations for SV and karyotype modes.
- `canonical_exons.hg38.tsv.gz` — MANE Select v1.5 canonical-transcript span + per-exon coordinates (~19,200 protein-coding genes).
- `clinvar.hg38.tsv.xz` — molamola's reduced ClinVar TSV (~12 MB; xz-compressed; release date logged in each report's run-metadata). The `--clinvar` flag accepts either this TSV or NCBI's raw ClinVar VCF (auto-detected by extension).
- `exclusion.hg38.bed.gz`, `exclusion.t2t.bed.gz` — karyotype-mode exclusion masks (~5.7 MB / ~6 MB; low-mappability ∪ polymorphic-TR catalog), built from 500 bp runs. A bin is dropped when more than `--mask-overlap` of it is masked. Override with `--mask`, disable with `--no-mask`.
- `gc_10kb.hg38.bed.gz`, `gc_10kb.t2t.bed.gz` — karyotype-mode 10 kb GC tables (~1.5 MB each) driving the per-1 % GC-bucket median-ratio correction. Override with `--gc`, disable with `--no-gc`.

molamola is bundled-only by design — no auto-download, no online lookups. The bundled-data total is ~30 MB. All bundled refs are reproducibly regeneratable from public sources by running `scripts/derive_canonical_exons.py`, `scripts/derive_clinvar_for_molamola.py`, and `scripts/derive_karyotype_refs.py`.

## Circos plot (SV mode)

Outer ring: cytoband ideogram (default pyCirclize colours). Inner ribbons: each PASS BND, line-thickness scaled by `SUPPORT`, colour by VAF class (three discrete bands: mosaic 0-33 %, het 33-66 %, hom 66-100 %). Noise-flagged BNDs render grey/dashed at low alpha so the eye goes to the candidate signal. With `--plotvaf`, each drawn arc also carries its VAF as a percentage on the rim.

## Genome SV map (SV mode)

One row per chromosome, top-down chr1 → chrY:

- **Chromosome bar** (cytoband ideogram, greyscale: white → very-dark-grey, centromere = black).
- **Four density strips** above the bar, one per non-BND SV type (INS = blue, DEL = red, DUP = green, INV = purple). Each strip is 1 Mb-binned; alpha encodes count (sqrt-scaled, single-event bins still visible).
- **BND arcs** above the strips, with apex height proportional to row gap and span. Same VAF / SUPPORT encoding as the circos ribbons.

## Phased-haplotype panel (compound-het mode)

One panel per plotted gene. Top-down:

- **Title** — gene symbol baked into the PNG (top-left, bold) so saved images are self-describing.
- **Canonical-transcript exon track** — IGV-style blue rectangles connected by a thin grey line; coordinates from the bundled MANE Select TSV.
- **Mint phase block** — one rectangle per WhatsHap PS group, spanning both haplotypes; off-edge arrows when the block stretches past the gene window.
- **H1 / H2 hap lines** — grey horizontal rules; H1 above H2.
- **Missense lollipops** — one filled circle per phased het missense canonical-transcript variant, hanging downward from its hap line. Colour by ClinVar bucket (P/LP `#c0143c`, VUS `#f4a013`, conflicting `#d9c200`, benign `#5fa860`, no-ClinVar grey).
- **Synonymous-variant ticks** — `x` markers on the hap line for canonical-transcript variants that aren't missense, included as context.
- **Chip row above the figure** — missense / total / blocks / trans / cis stats.
- **Legend** — at the bottom; ClinVar buckets + the synonymous tick.

## Karyotype coverage figure (karyotype mode)

One genome-wide figure, chr1 → chrY left to right:

- **log2 relative-depth scatter** — per-bin `log2(depth / autosomal median)` after masking + GC correction, so CN 2 sits at 0 and a single-copy loss and gain read symmetrically. CN 1 / 2 / 3 reference lines are labelled at the right edge. Adjacent chromosomes alternate between two inks and every other chromosome carries a background band, so a deviation can be attributed without tracing back to the axis; chrX and chrY have their own inks and a legend. Aggregated into ~`--scatter-bin-kb` windows; systematically downsampled above `--max-points`.
- **Rolling-median smooth** — a deep-pink per-chromosome line (`--smooth-window-mb` window, default 10 Mb), the load-bearing signal for arm-scale CN events.
- **Expected-CN dashes** — Oxford-blue horizontal guides at the expected CN for each chromosome given the inferred (or `--sex`-forced) genomic sex.
- **BAF panel** (only when `--vcf` is given) — allele balance at PASS heterozygous SNVs beneath the CN panel, sharing the x-axis. Reference lines sit at 33 / 50 / 67 % — the CN 3 het expectations — so a gain reads as the cloud splitting onto the outer two lines rather than as a vague widening; deletions / LOH splay toward 0 and 100 %.
  - **Phased VCF** (`FORMAT/PS` + `FORMAT/AD` present): haplotype-resolved. Read counts are summed within each phase block over tiling windows of 40 het SNVs — summing reads rather than averaging per-site fractions is the variance-correct estimator — and each window is plotted twice, at *v* and *1−v*, because which haplotype a block labels "1" is arbitrary and flips between blocks. A systematic reference-mapping bias is absorbed by shifting the autosomal median onto 50 %. On a normal ONT genome this tightens the panel about 1.4× versus per-site (SD 0.058 vs 0.079), which is what puts the CN 3 expectation ~2.9 SD clear of balanced.
  - **Unphased VCF**: per-site het allele fraction, as before. Phase is an upgrade, never a requirement; `--no-phased-baf` forces the per-site panel. The mode used is recorded in the report's metadata block.
- **Metadata strip** across the top: mosdepth source, reference, inferred sex, bin sizes, smooth window, mask / GC labels, BAF source, and an **"AS suspected"** chip when the autosomal depth distribution looks bimodal (adaptive-sampling-like — CN normalisation may be biased on such samples; the CLI also prints a one-line warning).

The per-chromosome 3 × 8 A4-portrait grid that briefly shipped in v0.3 development was dropped before release — the genome-wide figure does the cytogenetics work on its own.

## HTML report

A single self-contained HTML file (no server, no external CSS/JS, opens offline). Figures are embedded as base64 PNG data URIs. A collapsible run-metadata block at the bottom shows mode-specific provenance:

- **SV mode**: caller, filter mode, baseline coverage, cov-anomaly threshold, BND noise breakdown.
- **Compound-het mode**: ClinVar source + release date, canonical-exons source, VCF parse refusal counts, selection rule, plotted-vs-scanned gene counts.
- **Karyotype mode**: mosdepth source, reference, inferred genomic sex, mosdepth bin + scatter bin + smooth window, mask / GC source labels, BAF source + het-site count (when `--vcf` given), and the AS-suspected flag.

Designed as a portable review surface — emailable, archivable, no install needed for the reader.
