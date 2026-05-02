# Outputs

Each run writes a single self-contained HTML artefact into `--out` (default: the parent directory of `--vcf`). The filename and contents differ by the auto-detected mode:

- **SV mode** → `<sample>.report.html` (circos plot + linear cytoband SV map embedded as base64 PNG data URIs, plus a collapsible run-metadata block).
- **Compound-het mode** → `<sample>.compound_het.report.html` (one phased-haplotype panel per gene, embedded the same way; auto-select runs split into a `strict` and `extended` section).

No separate PNG files are written; figures live inside the HTML. Right-click → save to extract a PNG locally if you need one — each compound-het panel's PNG carries the gene symbol baked in (top-left, bold) so a saved image is self-describing.

## Reference data

All bundled in `molamola/data/`:

- `cytoBand.txt.gz` (hg38) and `cytoBand.t2t.txt.gz` (T2T-CHM13v2.0) — UCSC cytoband annotations for SV mode.
- `canonical_exons.hg38.tsv.gz` — MANE Select v1.5 canonical-transcript span + per-exon coordinates (~19,200 protein-coding genes).
- `clinvar.hg38.tsv.xz` — molamola's reduced ClinVar TSV (~13 MB; xz-compressed; release date logged in each report's run-metadata). The `--clinvar` flag accepts either this TSV or NCBI's raw ClinVar VCF (auto-detected by extension).

molamola is bundled-only by design — no auto-download, no online lookups. All bundled refs are reproducibly regeneratable from public sources by running `scripts/derive_canonical_exons.py` and `scripts/derive_clinvar_for_molamola.py`.

## Circos plot (SV mode)

Outer ring: cytoband ideogram (default pyCirclize colours). Inner ribbons: each PASS BND, line-thickness scaled by `SUPPORT`, colour by VAF (plasma, saturating above 0.7). Noise-flagged BNDs render grey/dashed at low alpha so the eye goes to the candidate signal.

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

## HTML report

A single self-contained HTML file (no server, no external CSS/JS, opens offline). Figures are embedded as base64 PNG data URIs. A collapsible run-metadata block at the bottom shows mode-specific provenance:

- **SV mode**: caller, filter mode, baseline coverage, cov-anomaly threshold, BND noise breakdown.
- **Compound-het mode**: ClinVar source + release date, canonical-exons source, VCF parse refusal counts, selection rule, plotted-vs-scanned gene counts.

Designed as a portable review surface — emailable, archivable, no install needed for the reader.
