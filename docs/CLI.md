# CLI reference

```sh
molamola --vcf VCF --out DIR [--reference hg38|t2t] [...]
molamola --mosdepth REGIONS.bed.gz --out DIR [--vcf VCF] [...]
```

Single flat parser. The plot type is picked from the input:

- **VCF inputs are header-driven.** `##INFO=<ID=SVTYPE,...>` selects the SV / cytogenetics report; `##INFO=<ID=CSQ,...>` AND `##FORMAT=<ID=PS,...>` selects the per-gene compound-het panels.
- **`--mosdepth` selects karyotype mode.** When `--mosdepth` is given, karyotype coverage mode runs regardless of any `--vcf`; an accompanying `--vcf` is consumed only as the BAF source (its header shape is not used for dispatch).

Either `--vcf` or `--mosdepth` is required. Mode-specific flags are silently ignored when they don't apply to the active mode (e.g. passing `--gene` against an SV VCF — the SV mode just doesn't read it).

VCFs that match neither shape are refused with a clear error rather than rendering a misleading default.

See [`FILTERS.md`](FILTERS.md) for what every threshold does. See [`OUTPUTS.md`](OUTPUTS.md) for the output format.

## Common flags

| flag | default | description |
|---|---|---|
| `--vcf PATH` | — | input VCF (gzipped OK). Required unless `--mosdepth` is given. Selects SV or compound-het mode by header; consumed only as the BAF source when `--mosdepth` is also given. |
| `--mosdepth PATH` | — | mosdepth `regions.bed.gz` (uniform-bin runs only — `mosdepth --by <int>`). Activates karyotype coverage mode. Required unless `--vcf` is given. |
| `--out DIR` | **required** | output directory (created if absent). molamola exits 2 with a usage error if omitted, rather than silently writing next to the input. Filenames: `<sample>.report.html` (SV), `<sample>.compound_het.report.html` (compound-het), `<sample>.karyotype.report.html` (karyotype). |
| `--reference {hg38,t2t}` | `hg38` | reference assembly the input was called against. SV and karyotype modes support both; compound-het mode is hg38-only (the bundled canonical-exon and ClinVar refs are hg38-coordinate). |
| `--sample NAME` | input basename | sample label shown in the report header. |
| `--png` | off | also write each embedded figure as a standalone PNG alongside the HTML (for MultiQC / pipeline embeds). Karyotype mode writes `<sample>.karyotype.genome.png`. |
| `--force` | off | bypass the safety check that errors out when the input filename hints at a reference different from `--reference` (e.g. `sample.t2t.vcf` with `--reference hg38`). |

## SV-mode flags

Active when the input VCF carries `##INFO=<ID=SVTYPE,...>` (Sniffles2 / cuteSV / SVIM / pbsv / NanoVar).

| flag | default | description |
|---|---|---|
| `--filter {pass,all}` | `pass` | keep PASS BNDs only, or include GT-filtered events. |
| `--caller {auto,sniffles2,sniffles1,cutesv,svim,pbsv,nanovar}` | `auto` | SV caller; `auto` runs an INFO-fingerprint detector and falls back to sniffles2 on no match. Override useful for bcftools-merged or re-headered VCFs. |
| `--mark-acrocentric` / `--no-mark-acrocentric` | on for hg38, off for t2t | grey out BNDs with both ends in chr13/14/15/21/22 p-arms (mostly mapping artefacts on hg38; real sequence on T2T-CHM13v2.0). |
| `--cov-filter {none,mark,drop}` | `mark` | how to handle coverage-spike BNDs and DEL/DUP. `mark` greys flagged BNDs and drops noisy DEL/DUP from density strips; `drop` also removes flagged BNDs entirely; `none` ignores. |
| `--cov-ratio R` | `auto` | `max(COVERAGE) / baseline` threshold above which an event is suspicious. Default `auto` = `max(2.0, p99 of in-sample distribution)` per sample. The chosen value is printed at run start. |
| `--cov-vaf-max V` | `0.35` | VAF below which a high-coverage event is treated as repeat-collapse noise. |
| `--focus CHR:POS` | none | show only BNDs with an endpoint within `--focus-window` of `CHR:POS`. The second part can be either a position (`chr7:57716411`) or an ISCN cytoband (`chr7:q11.23`). Repeatable. |
| `--focus-window N` | `1000` | +/-bp tolerance for `--focus` matching. |
| `--min-svlen N` | `50` | hard SVLEN cutoff (bp) for non-BND SVs. Set `0` to disable. BNDs are unaffected. |
| `--bin-size N` | `1,000,000` | density-track bin width in bp. |

## Compound-het mode flags

Active when the input VCF carries `##INFO=<ID=CSQ,...>` AND `##FORMAT=<ID=PS,...>` (phased + VEP-annotated small-variant VCF).

| flag | default | description |
|---|---|---|
| `--gene SYMBOL` | none | plot exactly this gene; repeatable. Plots regardless of variant count (with a clear placeholder for empty cases). When omitted, the auto-select rule picks candidate genes. |
| `--clinvar PATH` | bundled | override the bundled ClinVar lookup. Accepts either molamola's reduced TSV (default at `data/clinvar.hg38.tsv.xz`) or NCBI's raw ClinVar VCF. Format auto-detected by `.tsv` vs `.vcf` extension. |
| `--canonical-exons PATH` | bundled | override the bundled canonical-exon TSV (default: `data/canonical_exons.hg38.tsv.gz`). |
| `--min-pair-count N` | `1` | auto-select threshold: gene qualifies iff at least one phase set has `>= N` trans pairs where one anchor is ClinVar `P/LP` or `VUS` and the partner is not benign. The HTML splits results into a **strict** section (both P/LP or VUS) and an **extended** section (anchor P/LP-or-VUS, partner conflicting / no-ClinVar / P/LP / VUS). Ignored when `--gene` is given. |
| `--max-genes N` | `50` | cap on the number of auto-selected genes; capped runs emit a stderr warning. |

## Karyotype-mode flags

Active when `--mosdepth PATH` is given. Produces one HTML with a genome-wide CN scatter + rolling-median smooth; an optional BAF panel is added beneath when `--vcf` is also supplied.

| flag | default | description |
|---|---|---|
| `--mask PATH` | bundled | override the bundled exclusion mask (`data/exclusion.{hg38,t2t}.bed.gz`). Bins overlapping the mask drop out of the CN scatter, the smooth line, and the autosomal-median normalisation anchor. |
| `--no-mask` | off | disable masking entirely (use every bin). |
| `--gc PATH` | bundled | override the bundled 10 kb GC table (`data/gc_10kb.{hg38,t2t}.bed.gz`). Drives a per-1 % GC-bucket median-ratio correction applied to depth before normalisation. |
| `--no-gc` | off | disable GC correction. |
| `--centromere-pad-kb N` | `1000` | extend the exclusion mask by N kb either side of each `acen` band. Centromeric depth is unreliable on ONT even after the polymorphic-TR mask. |
| `--scatter-bin-kb N` | `50` | aggregation window (kb) for the CN scatter cloud. |
| `--smooth-window-mb N` | `0.5` | rolling-median window (Mb) for the deep-pink smooth line. |
| `--max-points N` | `200000` | hard cap on scatter points (systematic downsample above this). |
| `--max-baf-points N` | `80000` | hard cap on BAF points. |
| `--min-baf-dp N` | `10` | minimum `FORMAT/DP` for a het site to enter the BAF panel. |
| `--min-baf-gq N` | `20` | minimum `FORMAT/GQ` for a het site to enter the BAF panel. Applied only when GQ is present in the VCF (some callers don't emit it). |
| `--ymax N` | `5.0` | upper limit of the CN axis; higher events clip to the top edge by design. |
| `--sex {male,female,auto}` | `auto` | genomic sex for the expected-CN dashes. `auto` calls male iff chrY median CN > 0.3. The HTML wording is "inferred genomic sex" — a heuristic, not a clinical call. |
