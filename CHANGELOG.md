# Changelog

All notable changes to molamola are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project follows [Semantic Versioning](https://semver.org/).

## [0.4.0] — 2026-08-21

### Added

- **`--plotvaf`** (SV mode) — prints each BND's VAF as a percentage next to its arc on the circos plot. Off by default, since on a WGS call set the labels overplot; intended for targeted / panel runs where reading the exact VAF off the plot beats reading the class colour. Labels sit on the rim just outside the ideogram, staggered across four rings because breakpoints cluster and a single ring is unreadable. Where a cluster has more breakpoints than rings the labels are drawn anyway rather than silently dropped, and the run reports how many overlap. Noise-flagged BNDs are not labelled — they are deliberately de-emphasised and a label would undo that.

### Changed

- **BND VAF is now drawn as three discrete classes instead of a continuous ramp.** Boundaries are even thirds: `0-33 %` mosaic, `33-66 %` het, `66-100 %` hom. A linear `[0, 1]` plasma ramp compressed the bulk of a typical ONT call set into one narrow purple-to-magenta stretch that was not separable at 1 px linewidth, and its yellow end reached only 1.6-2.0:1 contrast against the page, making the rarest and most interesting arcs (high VAF) the hardest to see. Even thirds keep the scale readable without memorising boundaries and land close enough to the biology to be useful. Every class colour clears 3:1 contrast. Both colorbars are stepped, tick at the class boundaries **as percentages**, and carry the class name inside each band.
- **Figure and page background is now the off-white `PAPER_BG` (`#FAF8F4`)** across all three modes, giving thin BND arcs, pale cytobands, and the karyotype scatter something to sit against. Applied to the SV circos, the linear genome map, the karyotype genome figure, and the compound-het gene panels, and to the shared HTML report CSS so the embedded figures sit flush with the page instead of showing as off-white boxes on white. This colour was already specified as `KARY_PAPER` in the karyotype palette but had never been referenced by any figure; `KARY_PAPER` is now an alias of `PAPER_BG`.
- `vaf_to_color()` and both colorbars are driven by a single shared `VAF_NORM`, so the arcs and the scale they are read against cannot drift apart.

### Fixed

- **BND records with a real reference base in the ALT were silently dropped.** `parse_alt_for_mate()` matched the base-then-bracket forms against a literal `"N"` (`N[chr:pos[`, `N]chr:pos]`), but per the VCF spec that leading character is the actual reference base. Callers that write the real base — Sniffles2 >= 2.8 and DRAGEN_SV among them — produce ALTs like `G]chr16:12345]`, which raised `ValueError`; `_build_event()` catches that and returns `None`, so the record vanished from the report with no warning. Only the bracket-first forms survived. On real ONT call sets this dropped roughly half of every BND set (e.g. 169 of 349, and 2059 of 4184 on a DRAGEN run). The parser now accepts any replacement sequence on either side, matches case-insensitively, and takes `.` for the single-breakend forms. Orientation mapping for the previously-working forms is unchanged.
- **BND ALTs carrying inserted sequence at the breakpoint were dropped in both orientations.** The old fixed-offset slicing (`s[2:-1]` / `s[1:-2]`) assumed a one-character replacement string, so multi-base ALTs such as `GTTTT[chr2:123[` or `]chr2:123]GTTTT` failed to parse. DRAGEN emits these routinely (up to 916 bases in one test call set). Parsing is now delimiter-based rather than offset-based.
- Contig names containing colons (e.g. HLA ALT contigs, `G]HLA-A*01:01:01:01:12345]`) now parse — only the final `:` is treated as the contig/position separator.

## [0.3.0] — 2026-05-12

### Added

- **Karyotype coverage mode** — third plot type, dispatched by a new top-level `--mosdepth PATH` flag. Produces a single self-contained HTML report with one embedded figure: a genome-wide CN scatter + rolling-median smooth, plus an optional BAF panel beneath when `--vcf` is also given. Mosdepth `regions.bed.gz` is the primary input; uniform-bin runs only (`mosdepth --by <int>`). CN is computed as `2.0 × depth / autosomal_non_masked_median`. Sex is auto-detected from chrY CN unless `--sex {male,female}` overrides; HTML metadata surfaces this as "inferred genomic sex" — the threshold is a heuristic, not a clinical sex call. If the autosomal depth distribution looks bimodal (off-target + on-target modes), an "AS suspected" chip is added to the run-metadata block and the CLI prints a one-line warning. Supports hg38 and T2T-CHM13v2.0.
- **Bundled exclusion masks** — `data/exclusion.hg38.bed.gz` (~5.7 MB) and `data/exclusion.t2t.bed.gz` (~6 MB). Bins overlapping these intervals (low mappability ∪ polymorphic-TR catalog) drop out of the CN scatter, smooth, and normalisation. Override with `--mask PATH` or disable with `--no-mask`. Best fit for ~500 bp mosdepth runs.
- **Bundled coarse (10 kb) GC tables** — `data/gc_10kb.hg38.bed.gz` and `data/gc_10kb.t2t.bed.gz` (~1.5 MB each). Per-1 % GC bucket median-ratio correction applied to raw depth before the autosomal-median normalisation. Override with `--gc PATH` or disable with `--no-gc`.
- **Plain-text BAF parser** — `read_baf_vcf()` streams het allele fractions from a small-variant VCF without shelling out to bcftools. Reads `FORMAT/AF` when present; falls back to `FORMAT/AD` parsed as `ref,alt` and computed as `alt / (ref + alt)`. Filters: `FILTER == "PASS"`, biallelic heterozygous GT, **REF and every ALT a single nucleotide** (SNV-only — indel-het AF is noisier on ONT data, alignment ambiguity around the breakpoint inflates the estimate), no symbolic ALT (`<DEL>` etc.), `FORMAT/DP >= --min-baf-dp`, `FORMAT/GQ >= --min-baf-gq` (applied only when GQ is present in the VCF), canonical chroms only. Hard-refuses upfront when the VCF header carries `##INFO=<ID=SVTYPE,...>` (i.e. an SV / CNV / BND VCF was supplied by mistake) with an error message pointing at the small-variant call set.
- **Karyotype-mode flags** — `--mosdepth`, `--mask`, `--no-mask`, `--gc`, `--no-gc`, `--centromere-pad-kb` (default 1000), `--scatter-bin-kb` (default 50), `--smooth-window-mb` (default 0.5), `--max-points` (default 200000), `--max-baf-points` (default 80000), `--min-baf-dp` (default 10), `--min-baf-gq` (default 20), `--ymax` (default 5.0), `--sex {male,female,auto}` (default auto).
- **`scripts/derive_karyotype_refs.py`** — one-shot derive script that copies the upstream exclusion BED verbatim into the molamola data tree and aggregates the upstream 500 bp GC table down to 10 kb bins (mean of non-sentinel values per window, all-N runs stay flagged as 255). Re-run on upstream refresh.

### Changed

- **`--vcf` is no longer required** when `--mosdepth` is given. The CLI now needs either `--vcf` or `--mosdepth`; both may be combined to add a BAF panel to the karyotype figure. When `--mosdepth` is set, karyotype mode runs and any accompanying `--vcf` is consumed only as the BAF source — VCF-header-based dispatch is skipped.
- **`--out <directory>` is now required** for every mode (SV, compound-het, karyotype). Previously molamola silently defaulted the output directory to the parent of the input file; the new behaviour is to exit 2 with a clear usage error if `--out` is omitted, consistent with molamola's "refuse rather than render misleading data" rule. The directory is created if it does not exist.
- `build_argparser()` gains a third `add_argument_group` ("Karyotype-mode flags") so `--help` shows SV / compound-het / karyotype flag groups side-by-side.
- Bundled-data total grows from ~14 MB to ~28 MB (the new exclusion masks + 10 kb GC tables).

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
