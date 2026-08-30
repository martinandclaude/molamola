# Examples

molamola selects SV mode from the VCF header and karyotype mode from the `--mosdepth` flag:

```sh
molamola --vcf path/to/sample.sniffles.vcf --out reports/
open reports/sample.report.html                 # SV mode

molamola --mosdepth path/to/sample.regions.bed.gz --out reports/
open reports/sample.karyotype.report.html       # karyotype mode
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

## Karyotype coverage

Karyotype mode is selected by `--mosdepth`, not by a VCF header. The input is a mosdepth `regions.bed.gz` from a **uniform-bin** run — generate it with a single fixed bin size:

```sh
mosdepth --by 1000 --no-per-base --fast-mode sample sample.bam
# → sample.regions.bed.gz
```

### Genome-wide CN (no BAF)

```sh
molamola --mosdepth sample.regions.bed.gz --reference hg38 --out reports/
```

### Add a BAF panel from a small-variant VCF

Pair with a phased / unphased **small-variant** VCF (Clair3 / DeepVariant — *not* an SV or CNV VCF; those are refused):

```sh
molamola --mosdepth sample.regions.bed.gz \
    --vcf sample.clair3.vcf.gz \
    --reference hg38 --out reports/
```

The VCF is consumed only as the BAF source — its header shape is not used for dispatch. A BAF panel renders beneath the depth panel, sharing the x-axis — haplotype-resolved when the VCF carries `FORMAT/PS` + `FORMAT/AD`, per-site otherwise.

### T2T-CHM13v2.0

```sh
molamola --mosdepth sample.t2t.regions.bed.gz --reference t2t --out reports/
```

### Loosen / tighten the BAF QC

```sh
# keep more (noisier) het sites
molamola --mosdepth sample.regions.bed.gz --vcf sample.clair3.vcf.gz \
    --min-baf-dp 5 --min-baf-gq 0 --reference hg38 --out reports/
```

### Disable the mask / GC correction

```sh
molamola --mosdepth sample.regions.bed.gz \
    --no-mask --no-gc --reference hg38 --out reports/
```

Useful as a diagnostic when you want to see the raw coverage shape without molamola's exclusion mask or GC correction applied.

### Heavier smoothing, coarser scatter

```sh
molamola --mosdepth sample.regions.bed.gz --reference hg38 --out reports/ \
    --smooth-window-mb 1.0 --scatter-bin-kb 100 --ymax 4
```

If the run-metadata shows an **"AS suspected"** chip, the sample looks like adaptive sampling and the CN scale is biased (anchored on the off-target background) — interpret the depth axis with care.

## Verify a VCF before plotting

```sh
./summarize_bnds.sh path/to/sample.sniffles.vcf
```

Reports SV-type counts, BND filter distribution, intra/inter-chromosomal split, and the top 10 BND chromosome pairs.
