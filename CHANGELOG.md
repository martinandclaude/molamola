# Changelog

All notable changes to molamola are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project follows [Semantic Versioning](https://semver.org/).

## [0.6.0] — 2026-08-30

### Removed

- **Compound-het mode is gone, and molamola is a cytogenetics tool now.** The per-gene phased-haplotype panels for recessive-disease workup — canonical-transcript exon track, H1/H2 hap lines, phase blocks, ClinVar-coloured missense lollipops — have been removed along with everything that only served them: ~1,300 lines, the `--gene` / `--clinvar` / `--min-pair-count` / `--max-genes` flags, 74 tests, and both derivation scripts. molamola covers SV / cytogenetics reports and karyotype coverage; those two share an audience and a workflow, and the third mode did not.
- **The bundled ClinVar and MANE Select references went with it.** `clinvar.hg38.tsv.xz` (12 MB) and `canonical_exons.hg38.tsv.gz` (1.8 MB) are no longer shipped, so the bundled-data total falls from **28.5 MB to 14.7 MB**. That also retires the "is this much bundled data acceptable in a noarch package" question every bioconda submission has raised. The `derive` optional-dependency extra (`pyliftover`) is dropped too — it existed only for the canonical-exon derivation.

### Added

- **SV density rings on the circos.** The circos now carries INS / DEL / DUP / INV as four alpha-encoded rings inside the cytoband ring, outermost first, matching the linear map's strip order. Same bins, same per-type normalisation and same alpha ramp as the linear map — the computation was hoisted into one shared helper so an identical bin cannot render at two different alphas across the two panels of one report. Before this the circos showed translocations against bare cytobands and positional SVs existed only on the linear map. Type is deliberately encoded by **radius** first and colour second: measured under simulated colour-vision deficiency, `SV_TYPE_COLOR` collapses on its own (INS/INV ΔE 4.0 under protanopia, DEL/DUP ΔE 5.9 under deuteranopia), and it is safe on the linear map only because each type owns a strip row there. Giving each type its own ring keeps colour redundant rather than load-bearing, so the locked palette is unchanged. Arcs now start at the inner edge of the ring stack, so connections occupy the disc and located events occupy the rings and the two never overprint.
- **A legend on the circos.** It was the only molamola figure without one, which left the new rings unlabelled. One box now sits right of the disc, keying the density rings with per-type event counts and per-bin peaks, with the VAF colorbar beneath it. The cytoband greyscale and the grey/dashed noise-arc style are keyed on the facing linear map rather than repeated here — the two figures are embedded back to back in one report, so a second copy cost disc space without telling the reader anything new.

- **The circos ring stack is weighted by how much each type actually says about the sample.** INS and DEL get 0.60 of the per-type radius share and an alpha ceiling of 0.62; DUP and INV get 1.40 and the full range. Drawn at equal weight, INS and DEL were the two thickest, densest bands on the figure — and they are the two that say least about the individual. Measured across the GIAB trio, which holds lab, prep, pipeline, caller, coverage and population constant and varies only relatedness: two **unrelated** founders already agree at INS R² 89.1 % / DEL 84.3 % per 1 Mb bin, and extrapolating the IBD 0 → 0.5 slope out to identical genomes adds only 7 pp / 10 pp. The per-bin INS/DEL landscape is a property of the species rather than of the person, so it is now drawn as background texture. This is also why the alpha ramp alone could not fix it: Spearman and Kendall are invariant under any monotone count→alpha mapping, so no ramp changes the ordering, and empirically the darkest 1–50 % of bins are 51–71 % *the same bins* in two strangers. Note this bounds background genome-to-genome variation only — whole-arm and aneuploidy-scale changes still register, and a real multi-Mb deletion is emitted as its own DEL record regardless.

### Changed

- **A phased, VEP-annotated VCF is now recognised and refused by name.** `detect_vcf_mode` still detects the `##INFO=<ID=CSQ,...>` + `##FORMAT=<ID=PS,...>` shape and returns `"compound-het-removed"`; the CLI exits 1 saying the mode was removed after v0.5.1 and pointing at that release. Such a VCF carries no `SVTYPE`, so without this it would have fallen through to the generic "unsupported shape" error and a returning user could not tell a removed capability from a malformed file. `SVTYPE` is now checked **first**, so a VEP-annotated phased SV VCF — rare, but plottable — is plotted rather than refused.

### Fixed

- **The SV density scale wasted most of its range, in both SV figures.** `sv_density_alpha` normalised to the peak bin, and SV density is heavily skewed enough that the peak is a far outlier from the typical bin — on MH001, INS has a median of 3 events per 1 Mb bin against a peak of 67. Combined with a 0.40 alpha floor eating 40 % of the range before any data was drawn, the p10–p90 alpha spread was **0.15 for INS and 0.12 for DEL**: nearly every bin rendered at the same ~0.5 and the track encoded almost nothing, reading as a uniform wash rather than as density. The ramp now saturates at the 99th-percentile occupied bin (genome-wide per type, never per chromosome) with the floor lowered to 0.20, giving spreads of **0.27 and 0.24** while leaving the faintest bin at 0.34 — still plainly visible. The ~1 % of bins above the anchor clip to full alpha, which costs nothing a reader needs: a hotspot that far out is already the darkest thing in the track. Legends now report both the real peak per type and the fact that the scale saturates at p99. This affects the circos rings and the linear map's strips identically, since they share the helper.
- **Empty SV density rings are no longer tinted.** The first version of the rings drew a faint full-length baseline so an empty ring still read as a ring. On a normal genome DUP and INV are two orders of magnitude rarer than INS and DEL, so those rings were almost entirely baseline, and a continuous pale band reads as low-level density everywhere — the opposite of what it meant. Only occupied bins get ink; the legend carries the ring order and the per-type counts.
- **The circos and the linear map drew cytobands with different palettes.** `plot_circos` never passed a `cytoband_cmap`, so it silently took pyCirclize's own colormap while the linear map used molamola's — the two figures in a single report disagreed on the entire grey ramp (`gpos50` was `#C8C8C8` on one and `#9C9C9C` on the other), not just on the centromere. The circos now uses molamola's greys, with one deliberate carve-out: the centromere stays red, because black `acen` is not separable from `gpos100` `#2A2A2A` on a ring five radial units thick, and it is the landmark that orients a circle with no axis. The linear map has room to render it black and continues to.
- **Chromosome labels could be clipped by the figure edge.** pyCirclize builds the figure with `tight_layout=True`, which re-lays-out the polar axes at draw time; the disc is now placed explicitly with a fixed margin, which both stops `--plotvaf`'s sector names at r=116 running off the canvas and silences the "Axes not compatible with tight_layout" warning the colorbar axes provoked.

## [0.5.1] — 2026-08-29

### Changed

- **VAF class colours recoloured for colour-vision deficiency.** `mosaic` / `het` / `hom` are now `#035AF3` / `#634980` / `#781B00`. The previous palette put `het` and `hom` at nearly the same lightness, so they differed almost only in hue and merged under tritanopia (ΔE 11.3 as drawn — below the level at which two thin arcs can be told apart). Colours are now selected under simulated CVD (Machado et al. 2009) and, importantly, evaluated on the colours **as drawn**: arcs render at alpha 0.70 over the page, which costs about a third of the nominal contrast, so judging raw hex values overstates how distinct they actually are. Worst-case class separation across normal, protan, deutan and tritan vision is now ΔE 37.2, every class still clears 3:1 against the page as composited, no class lands within ΔE 25 of the noise grey, and lightness is monotone with VAF so the ordering survives even total loss of hue discrimination. `tests/test_vaf_colors.py` pins all four of these.

## [0.5.0] — 2026-08-21

### Added

- **Haplotype-resolved BAF panel, used automatically when the VCF is phased.** If the small-variant VCF carries `FORMAT/PS` and `FORMAT/AD`, het read counts are summed within each `(chrom, PS)` phase block over tiling windows of 40 het SNVs (also split on a >1 Mb positional gap), giving a depth-weighted haplotype-1 fraction per window. Summing reads rather than averaging per-site fractions is the variance-correct estimator. Each window is plotted twice, at *v* and *1−v*, because which haplotype a block labels "1" is arbitrary and flips between blocks, so the track is invariant to that flip. A systematic reference-mapping bias is absorbed by shifting the autosomal median onto 50 %. On a normal ONT genome this tightens the panel ~1.4× versus per-site (SD 0.058 vs 0.079), putting the CN 3 expectation ~2.9 SD clear of balanced. Phase is **never required** — molamola is a plotting tool and takes whatever VCF it is handed, so an unphased VCF still gets the per-site panel, and `--no-phased-baf` forces per-site regardless. The mode used is printed at run time and recorded in the report metadata.

### Fixed

- **The exclusion mask over-excluded on any mosdepth bin coarser than 500 bp.** `annotate_mask()` flagged a bin if it touched a mask interval at all. The bundled masks are built from 500 bp runs (median interval exactly 500 bp), so on a 1 kb mosdepth run a half-bin hit removed the whole bin: a mask covering **39.8 %** of hg38 excluded **57.8 %** of bins. Those bins are dropped from the CN normalisation anchor as well as from the plot, so this moved the numbers, not just the picture. Overlap is now measured in base pairs and a bin is dropped only when **more than** `--mask-overlap` (default 0.5) of it is masked. On a 1 kb run the excluded fraction falls from 57.8 % to 23.6 %. The cut is empirical, not just geometric: half-masked 1 kb bins measure like clean sequence (SD 0.230 log2 vs 0.225 for clean bins), while fully-masked bins have SD 1.072. `--mask-overlap 0` restores the old any-overlap behaviour. At 500 bp bins the two rules agree exactly.

### Changed

- **Karyotype CN panel is now log2 relative depth instead of linear CN.** A linear 0-5 CN axis spent more than half its height above CN 2.5 where there is no data, compressing the deviations that matter. The panel plots `log2(depth / autosomal median)`, so a single-copy loss and a single-copy gain read symmetrically, with CN 1 / 2 / 3 reference lines labelled at the right edge. **`--ymax` is now in log2 units** (default 1.5, was 5.0 in CN units) and a matching `--ymin` (default -2.0) is added — a value that made sense as a CN limit will not make sense here.
- **Chromosomes are now visually separable.** Adjacent chromosomes alternate between two scatter inks and every other chromosome gets a background band, so a deviation can be attributed without tracing back to the axis. chrX and chrY get their own inks with a legend, making a single-copy sex-chromosome pattern legible at a glance.
- **BAF reference lines moved from 25 / 50 / 75 % to 33 / 50 / 67 %** — the CN 3 het expectations, so a trisomy reads as the cloud splitting onto two lines rather than as a vague widening. Tick labels are percentages, matching the rest of molamola's user-facing output.
- **`--smooth-window-mb` default raised from 0.5 to 10.** The 0.5 Mb default was tuned for the per-chromosome grid dropped in v0.3; across a single 3.1 Gb axis it produces roughly 6,200 sub-pixel wiggles that render as a solid band over the scatter rather than as a trend line.
- **The run-metadata strip is no longer drawn on the figure.** Sample / build / sex / bin / smooth crowded the arm ticks; the same information (including the "AS suspected" chip) is already in the HTML report's metadata block, which is the headline output. Note this does mean a standalone `--png` no longer carries its own provenance line.
- **Arm labels are rotated upright** so the short arms on chr17-22 stop overprinting each other, and the panels are taller (genome-only 4.8 -> 6.4 in, with BAF 6.4 -> 8.6 in) to give the point cloud room to read as a distribution.

### Notes

- GC correction was reviewed against PIKA's LOESS approach and deliberately left alone: molamola's per-1 %-bucket median-ratio correction already removes GC bias completely (per-bucket median spread 0.085 log2 raw -> 0.000 corrected). LOESS is a smoother fit, not a better one.

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
