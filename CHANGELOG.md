# Changelog

All notable changes to molamola are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project follows [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.2.0] — 2026-05-05

### Added

- `--png` flag — when set, writes each embedded figure as a standalone PNG file in the same directory as the HTML report. Default behaviour unchanged (HTML remains the headline output and is always written). SV mode produces `<basename>.report.circos.png` and `<basename>.report.sv_map.png`; compound-het mode produces one `<basename>.compound_het.report.<GENE_SYMBOL>.png` per plotted gene. Addresses [#1](https://github.com/martinandclaude/molamola/issues/1) — useful for embedding molamola figures in MultiQC and other pipeline reports.

## [0.1.0] — 2026-05-03

Initial public release.

### Plotting modes

- **SV / cytogenetics report** — circos plot (pyCirclize) plus linear genome SV map. Embeds two figures into a single self-contained HTML report. Reads long-read SV VCFs from Sniffles2, cuteSV, SVIM, pbsv, and NanoVar; auto-detects the caller from `##source=` and INFO-field fingerprinting. Supports hg38 and T2T-CHM13v2.0 via bundled UCSC cytobands. ISCN-style nomenclature on BND events. Two on-by-default noise heuristics: acrocentric short-arm BNDs (chr13/14/15/21/22 p-arms — off by default on T2T where those regions are properly resolved) and coverage anomalies (`max(COVERAGE) >= --cov-ratio × genome-median-coverage AND VAF < --cov-vaf-max`). Adaptive `--cov-ratio auto` (default) computes `max(2.0, p99 of in-sample distribution)` per sample. ISCN-band focus filter `--focus chr7:q11.23` accepts both genomic positions and cytoband names.

- **Per-gene phased-haplotype panels** (compound-het) — one panel per candidate gene: IGV-style blue canonical-transcript exon track, H1 / H2 hap lines, mint phase blocks across both haps with off-edge arrows when a block stretches past the gene window, ClinVar-coloured missense lollipops hanging downward, synonymous-variant `x` markers on the hap line for context. Reads phased + VEP-annotated small-variant VCFs (WhatsHap / HiPhase). hg38-only.

  Auto-select sweep (no `--gene`): gene qualifies iff at least one trans pair has one variant in ClinVar P/LP or VUS and the partner is not benign. Report splits results into a `strict` section (both variants P/LP or VUS) and an `extended` section (anchor P/LP-or-VUS, partner conflicting / no-ClinVar / P/LP / VUS). The strict heading is shown even when its subset is empty.

  Pair colours: P/LP `#c0143c` (merged), VUS `#f4a013`, conflicting `#d9c200`, benign `#5fa860`, no-ClinVar `#888`. Phase blocks: mint fill `#e2f0e3`, dark-green border `#3f6e44`. Exon track: IGV blue `#1E5BA8`. Gene symbol baked into the PNG title (top-left, bold) so saved images are self-describing.

### CLI

- Single flat parser; the plot type is auto-detected from the VCF header: `##INFO=<ID=SVTYPE,...>` → SV mode; `##INFO=<ID=CSQ,...>` AND `##FORMAT=<ID=PS,...>` → compound-het mode. VCFs that match neither shape are refused with a clear error.
- `--vcf PATH` is the only required argument. SV-mode and compound-het flags are visually grouped in `--help`.
- Console script: `molamola --vcf foo.vcf` after `pip install molamola`.

### Bundled references

- `data/cytoBand.txt.gz` (hg38) and `data/cytoBand.t2t.txt.gz` (T2T-CHM13v2.0) — UCSC cytoband annotations.
- `data/canonical_exons.hg38.tsv.gz` — MANE Select v1.5 (~19,200 protein-coding genes, schema `gene_symbol\tchrom\tstart\tend\tstrand\ttranscript_id\texon_starts\texon_ends`).
- `data/clinvar.hg38.tsv.xz` — molamola's reduced ClinVar TSV (`chrom\tpos\tref\talt\tbucket`; xz-compressed; ~13 MB; release date logged in each report's run-metadata). The `--clinvar` flag accepts either this TSV or NCBI's raw ClinVar VCF (auto-detected by extension).
- `scripts/derive_canonical_exons.py` and `scripts/derive_clinvar_for_molamola.py` reproducibly regenerate the bundled refs from public sources.

### Output

- One self-contained HTML report per run, embedded base64 PNG figures, no external CSS/JS, no separate image files. SV mode writes `<sample>.report.html`; compound-het mode writes `<sample>.compound_het.report.html`. Run-metadata in a collapsible `<details>` block at the bottom (caller, filter mode, baseline coverage, ClinVar release date, selection rule, refusal counts as appropriate).

### Hard refusals

- Wrong-assembly VCFs (contig lengths disagree with the chosen reference) → exit 1 with a clear error.
- Compound-het mode on `--reference t2t` → exit 1 (hg38-only in v0.1).
- Phased-VCF without `##INFO=<ID=CSQ,...>` (no VEP annotation) → exit 1.
- Phased-VCF without `##FORMAT=<ID=PS,...>` (unphased) → exit 1.
- `--gene FOO` for a symbol absent from the canonical-exon table → exit 1.

### Dependencies

- Runtime: matplotlib, numpy, pandas, biopython, pyCirclize.
- Dev: pytest, ruff.
- `derive` extra: pyliftover (only needed by `scripts/derive_canonical_exons.py --lift-to <chain>`).
