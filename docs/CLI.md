# CLI reference

```sh
molamola --vcf VCF [--out DIR] [--reference hg38|t2t] [...]
```

Single flat parser: the plot type is auto-detected from the VCF header. `##INFO=<ID=SVTYPE,...>` selects the SV / cytogenetics report; `##INFO=<ID=CSQ,...>` AND `##FORMAT=<ID=PS,...>` selects the per-gene compound-het panels. Mode-specific flags are silently ignored when they don't apply to the detected mode (e.g. passing `--gene` against an SV VCF — the SV mode just doesn't read it).

VCFs that match neither shape are refused with a clear error rather than rendering a misleading default.

See [`FILTERS.md`](FILTERS.md) for what every threshold does. See [`OUTPUTS.md`](OUTPUTS.md) for the output format.

## Common flags

| flag | default | description |
|---|---|---|
| `--vcf PATH` | required | input VCF (gzipped OK). |
| `--out DIR` | parent of `--vcf` | output directory. The report `<sample>.report.html` (SV) or `<sample>.compound_het.report.html` (compound-het) is written here. |
| `--reference {hg38,t2t}` | `hg38` | reference assembly the input VCF was called against. SV mode supports both; compound-het mode is hg38-only (the bundled canonical-exon and ClinVar refs are hg38-coordinate). |
| `--sample NAME` | VCF basename | sample label shown in the report header. |
| `--force` | off | bypass the safety check that errors out when the VCF filename hints at a reference different from `--reference` (e.g. `sample.t2t.vcf` with `--reference hg38`). |

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
