# Examples

molamola auto-detects the plot type from the VCF header, so the same command shape works for both SV and compound-het VCFs:

```sh
molamola --vcf path/to/sample.vcf --out reports/
open reports/sample.report.html        # SV mode
open reports/sample.compound_het.report.html  # compound-het mode
```

`--out DIR` is required — molamola refuses rather than silently writing the report next to the input file. The directory is created if it does not exist.

## SV / cytogenetics

### Genome-wide view (defaults)

```sh
molamola --vcf sample.sniffles.vcf --out reports/
```

### T2T-CHM13v2.0 reference

```sh
molamola --vcf sample.t2t.sniffles.vcf --reference t2t --out reports/
```

The bundled T2T cytoband (`molamola/data/cytoBand.t2t.txt.gz`) is selected automatically. The acrocentric noise flag flips off for T2T (those p-arms are real sequence on T2T-CHM13v2.0).

### Zoom in on a specific BND endpoint

```sh
molamola --vcf sample.sniffles.vcf --out reports/ \
    --focus chr7:57716411 --focus-window 10000
```

Output filenames get a `.focus_<chr>_<pos>` tag so focused renders don't overwrite genome-wide ones. ISCN labels (e.g. `t(7;17)(q11.23;q12)`) print to stdout for each matched event.

### Filter by ISCN cytoband

```sh
molamola --vcf sample.sniffles.vcf --focus chr7:q11.23 --out reports/
```

The second part of `--focus` can be either a position (`chr7:57716411`) or an ISCN cytoband (`chr7:q11.23`). Bands can also be specified as a prefix (`chr1:p36`) to match every sub-band on the chromosome whose name starts with that prefix.

### Override the adaptive coverage threshold

`--cov-ratio` defaults to `auto`, computed as `max(2.0, p99 of the in-sample max-coverage / median-coverage distribution)`. Pass a number (int or float) to fix the threshold instead — useful for cross-sample consistency:

```sh
molamola --vcf sample.sniffles.vcf --cov-ratio 2.5 --cov-vaf-max 0.45 \
    --out reports/
```

### Coarser density tracks

```sh
molamola --vcf sample.sniffles.vcf --bin-size 5000000 --out reports/
```

5 Mb bins instead of 1 Mb — smoother visual, less detail.

### Run a non-Sniffles2 caller

`--caller auto` (default) detects the caller via INFO-field fingerprinting. Override when the fingerprint is ambiguous (bcftools-merged or re-headered VCFs):

```sh
molamola --vcf sample.cutesv.vcf --caller cutesv --out reports/
```

## Compound-het

### Prepare a phased VCF for compound-het mode

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

Then `molamola --vcf sample.phased.vep.vcf --out reports/` picks it up as compound-het mode.

**Notes on VEP.** VEP is third-party software (Ensembl); molamola does not bundle or wrap it. The cache release and VEP binary release must match exactly — a mismatch leads to silent mis-annotation rather than a clean error. Compound-het mode reads VEP's `Consequence`, `SYMBOL`, and `CANONICAL` fields as-is; any quirks of a particular VEP build are inherited. `--pick` reduces multi-transcript CSQ entries to one per variant; without it molamola will pick `CANONICAL == YES` itself, which is fine but slower on large VCFs.

### One gene, explicit

```sh
molamola --vcf sample.phased.vep.vcf.gz --gene NEB --out reports/
```

The panel always plots regardless of variant count. If `--gene FOO` lands in a gene with zero phased het variants, the report still gets a placeholder section so multi-gene runs (`--gene A --gene B --gene C`) don't silently lose entries.

### Multi-gene focus

```sh
molamola --vcf sample.phased.vep.vcf.gz \
    --gene NEB --gene LDLR --gene CFTR --out reports/
```

### Auto-select sweep

Drop `--gene` to let molamola surface candidate genes:

```sh
molamola --vcf sample.phased.vep.vcf.gz --out reports/
```

The HTML splits the result into two clearly-labelled sections:

- **strict** — both variants P/LP or VUS (true compound-het).
- **extended** — anchor P/LP-or-VUS, partner conflicting / no-ClinVar / P/LP / VUS.

The strict heading is shown even when its subset is empty so the dichotomy is always visible.

Tune via `--min-pair-count N` (raise for stricter sweeps) and `--max-genes N` (default 50). Conflicting+conflicting and no-ClinVar+no-ClinVar pairs are excluded from the sweep — use explicit `--gene` to plot those one-off.

### Override ClinVar with a fresher snapshot

The bundled ClinVar (`molamola/data/clinvar.hg38.tsv.xz`, ~13 MB) is dated; pass a fresher copy with `--clinvar`:

```sh
# Reduced TSV (smallest)
molamola --vcf sample.phased.vep.vcf.gz --out reports/ \
    --clinvar /path/to/clinvar.hg38.tsv.xz

# Or NCBI's raw VCF (190+ MB) — auto-detected by extension
molamola --vcf sample.phased.vep.vcf.gz --out reports/ \
    --clinvar /path/to/clinvar.vcf.gz
```

The release date of whichever file is loaded is logged in the HTML report's run-metadata.

### Override the canonical-exon table

```sh
molamola --vcf sample.phased.vep.vcf.gz --out reports/ \
    --canonical-exons /path/to/custom_exons.tsv.gz
```

Schema: `gene_symbol\tchrom\tstart\tend\tstrand\ttranscript_id\texon_starts\texon_ends` (gzipped TSV; comma-separated exon coords; 0-based half-open).

## Verify a VCF before plotting

```sh
./summarize_bnds.sh path/to/sample.sniffles.vcf
```

Reports SV-type counts, BND filter distribution, intra/inter-chromosomal split, and the top 10 BND chromosome pairs.
