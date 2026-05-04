```
                 _                       _              .--.
 _ __ ___   ___ | | __ _ _ __ ___   ___ | | __ _      _/    \___
| '_ ` _ \ / _ \| |/ _` | '_ ` _ \ / _ \| |/ _` |    ( o        )
| | | | | | (_) | | (_| | | | | | | (_) | | (_| |     \___..___/
|_| |_| |_|\___/|_|\__,_|_| |_| |_|\___/|_|\__,_|         ||
```

A Python plotting tool for Oxford Nanopore variation data. **One VCF in, one self-contained HTML report out.** molamola inspects the VCF header and picks the right plot type automatically — no flags or subcommands to remember:

- **SV / cytogenetics report** for long-read SV VCFs (Sniffles2 / cuteSV / SVIM / pbsv / NanoVar). Cytoband-ideogram circos plot plus a linear genome SV map with per-type density tracks (INS / DEL / DUP / INV) and BND arcs.
- **Per-gene phased-haplotype panels** for phased + VEP-annotated small-variant VCFs (WhatsHap / HiPhase). One panel per candidate gene: canonical-transcript exon track, H1 / H2 hap lines, mint phase blocks across both haps, ClinVar-coloured missense lollipops and synonymous-variant ticks for context.

Both produce one self-contained HTML report — figures embedded as base64 PNGs, no external assets, opens offline.

## Install

```sh
pip install molamola
```

Or via conda from the bioconda channel:

```sh
conda install -c bioconda molamola
```

Or for development from a clone:

```sh
git clone https://github.com/martinandclaude/molamola.git
cd molamola
pip install -e .[dev]
pytest -v
```

**Documentation:** <https://martinandclaude.github.io/molamola/>

## Quick start

```sh
# Long-read SV VCF (Sniffles2 etc.) → cytogenetics report
molamola --vcf sample.sniffles.vcf
open path/to/sample.report.html

# Phased + VEP-annotated VCF → compound-het workup, all candidate genes
molamola --vcf sample.phased.vep.vcf.gz
open path/to/sample.compound_het.report.html

# Just one gene from a phased + VEP VCF
molamola --vcf sample.phased.vep.vcf.gz --gene NEB
```

The plot type is auto-detected from the VCF header: `##INFO=<ID=SVTYPE>` selects SV mode; `##INFO=<ID=CSQ>` + `##FORMAT=<ID=PS>` selects compound-het mode. VCFs that match neither shape are refused with a clear error.

## Preparing a phased VCF for compound-het mode

Compound-het mode needs both phasing (`PS` FORMAT field) and VEP annotation (`CSQ` INFO field). A raw phased small-variant VCF — e.g. straight Clair3 output, even with `--use_whatshap_for_final_output_phasing` — has the first but not the second, and molamola will refuse it. Annotate with [Ensembl VEP](https://github.com/Ensembl/ensembl-vep) first.

**1. Download the matching VEP cache once** (one-time, ~20 GB), on a machine with internet:

```sh
wget https://ftp.ensembl.org/pub/release-105/variation/indexed_vep_cache/homo_sapiens_vep_105_GRCh38.tar.gz
```

If your VEP isn't 105, swap `105` for your release number (it appears twice in the URL); the cache release must match the VEP release exactly. Transfer to wherever you run VEP if that's a different machine.

**2. Unpack into a stable cache directory:**

```sh
mkdir -p VEP_cache && cd VEP_cache
tar -xzf ../homo_sapiens_vep_105_GRCh38.tar.gz
# creates VEP_cache/homo_sapiens/105_GRCh38/
```

**3. Run VEP fully offline**, with `--canonical --symbol --pick` so the CSQ shape matches what molamola consumes:

```sh
vep --input_file sample.phased.vcf \
    --output_file sample.phased.vep.vcf \
    --vcf --offline \
    --cache --dir_cache /path/to/VEP_cache \
    --assembly GRCh38 \
    --fasta /path/to/hg38.fa \
    --canonical --symbol --pick \
    --force_overwrite
```

Then `molamola --vcf sample.phased.vep.vcf` picks it up as compound-het mode.

**Notes on VEP.** VEP is third-party software (Ensembl); molamola does not bundle or wrap it. The cache release and VEP binary release must match exactly — a mismatch leads to silent mis-annotation rather than a clean error. Compound-het mode reads VEP's `Consequence`, `SYMBOL`, and `CANONICAL` fields as-is; any quirks of a particular VEP build are inherited. `--pick` is what reduces multi-transcript CSQ entries to one per variant; without it molamola will pick `CANONICAL == YES` itself, which is fine but slower on large VCFs.

## What it produces

### SV / cytogenetics report

- **Circos plot** (pyCirclize) — cytoband ideogram with BND ribbons; line thickness scaled by `SUPPORT`, colour by VAF.
- **Linear genome SV map** — chr1 → chrY, one row each. Greyscale ISCN-style cytobands, four per-type density strips (INS = blue, DEL = red, DUP = green, INV = purple) at 1 Mb bins, BND arcs above. Annotated with ISCN nomenclature like `t(7;17)(q11.23;q12)`.
- Two noise heuristics that work on the VCF alone — no external reference data needed: acrocentric short-arm BNDs (chr13/14/15/21/22 p-arms; on by default for hg38, off for T2T) and coverage-anomaly BNDs (`max(COVERAGE) >= --cov-ratio × baseline AND VAF < --cov-vaf-max`).
- Supports hg38 and T2T-CHM13v2.0 via bundled cytobands.

### Compound-het panels

- One panel per gene: IGV-style blue canonical-transcript exon track on top, two horizontal H1 / H2 hap lines, mint phase-block rectangles spanning both haps (with off-edge arrows when a block stretches past the gene window), ClinVar-coloured missense lollipops hanging downward, synonymous-variant `x` markers on the hap line for context.
- Auto-select sweep when `--gene` is omitted: gene qualifies iff at least one trans pair has one variant in ClinVar P/LP or VUS and the partner is not benign. The report splits results into a `strict` section (both variants P/LP or VUS — true compound-het) and an `extended` section (anchor P/LP-or-VUS, partner conflicting / no-ClinVar / P/LP / VUS). The strict heading is shown even when its subset is empty so the dichotomy is always visible.
- Use `--gene SYMBOL` to plot a specific gene regardless of the auto-select rule (useful for manual review of P/LP + benign or no-ClinVar + no-ClinVar pairs).
- Tunable via `--min-pair-count` (raise for stricter sweeps) and `--max-genes` (default 50).
- hg38-only: ClinVar coordinates are hg38, and coordinate-based lookup is the matching path.

ClinVar usage in compound-het is purely a colour key on data points — no clinical interpretation is performed or implied.

## Bundled references

All in `molamola/data/`:

- `cytoBand.txt.gz` (hg38), `cytoBand.t2t.txt.gz` (T2T-CHM13v2.0) — UCSC cytoband annotations for SV mode.
- `canonical_exons.hg38.tsv.gz` — MANE Select v1.x canonical transcripts and exon coordinates.
- `clinvar.hg38.tsv.xz` — molamola's reduced ClinVar TSV (chrom, pos, ref, alt, significance bucket; xz-compressed). Release date logged in each report's run-metadata.

Bundled-only by design: molamola does not auto-download or look up online. Use `--clinvar PATH` or `--canonical-exons PATH` to override.

The reduced TSVs are reproducibly regeneratable from public sources via `scripts/derive_canonical_exons.py` and `scripts/derive_clinvar_for_molamola.py`.

## CLI

```sh
molamola --vcf VCF [--out DIR] [--reference hg38|t2t] [...]
```

Full flag list: [`docs/CLI.md`](docs/CLI.md). Filter explanations: [`docs/FILTERS.md`](docs/FILTERS.md). Output formats: [`docs/OUTPUTS.md`](docs/OUTPUTS.md). Worked examples: [`docs/EXAMPLES.md`](docs/EXAMPLES.md). Per-release changes: [`CHANGELOG.md`](CHANGELOG.md).

## Development

```sh
pytest -v
ruff check .
```

CI runs lint + pytest on Python 3.10 / 3.11 / 3.12 on every push and PR.

## Acknowledgements

- [Sniffles2](https://github.com/fritzsedlazeck/Sniffles), [cuteSV](https://github.com/tjiangHIT/cuteSV), [SVIM](https://github.com/eldariont/svim), [pbsv](https://github.com/PacificBiosciences/pbsv), [NanoVar](https://github.com/cytham/nanovar) — long-read SV callers.
- [WhatsHap](https://github.com/whatshap/whatshap), [HiPhase](https://github.com/PacificBiosciences/HiPhase) — long-read phasing.
- [VEP](https://github.com/Ensembl/ensembl-vep), [MANE Select](https://www.ncbi.nlm.nih.gov/refseq/MANE/), [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) — variant annotation and significance.
- [pyCirclize](https://github.com/moshi4/pyCirclize) — circos plot.
- [matplotlib](https://github.com/matplotlib/matplotlib), [numpy](https://github.com/numpy/numpy).
- [bcftools / samtools / htslib](https://github.com/samtools/bcftools) — VCF pre-processing helpers.
- [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/) — hg38 and T2T-CHM13v2.0 cytobands.
- [iconsdb.com](https://www.iconsdb.com/) — header fish icon (deep-pink, mirrored).

## License

[MIT](LICENSE).
