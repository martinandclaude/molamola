```
                 _                       _              .--.
 _ __ ___   ___ | | __ _ _ __ ___   ___ | | __ _      _/    \___
| '_ ` _ \ / _ \| |/ _` | '_ ` _ \ / _ \| |/ _` |    ( o        )
| | | | | | (_) | | (_| | | | | | | (_) | | (_| |     \___..___/
|_| |_| |_|\___/|_|\__,_|_| |_| |_|\___/|_|\__,_|         ||
```

A Python cytogenetics plotting tool for Oxford Nanopore data. **One input in, one self-contained HTML report out.** molamola picks between two report types based on its input:

- a VCF with `##INFO=<ID=SVTYPE,...>` → SV / cytogenetics report (long-read SV VCFs from Sniffles2 / cuteSV / SVIM / pbsv / NanoVar);
- a mosdepth `regions.bed.gz` (via `--mosdepth`) → karyotype coverage report (genome-wide log2 relative-depth scatter + rolling-median smooth, with an optional BAF panel beneath when paired with a small-variant VCF — haplotype-resolved if that VCF is phased).

Figures embedded as base64 PNGs — no external assets, opens offline.

**Full documentation:** <https://martinandclaude.github.io/molamola/>

## Install

```sh
pip install molamola
```

Or via conda — note that both bioconda and conda-forge channels are needed (pycirclize lives on conda-forge):

```sh
conda create -n molamola -c bioconda -c conda-forge molamola
```

Or for development from a clone:

```sh
git clone https://github.com/martinandclaude/molamola.git
cd molamola
pip install -e .[dev]
pytest -v
```

## Quick start

```sh
# Long-read SV VCF (Sniffles2 etc.) → cytogenetics report
molamola --vcf sample.sniffles.vcf --out reports/
open reports/sample.report.html

# Mosdepth output → karyotype coverage report (genome-wide log2 depth)
molamola --mosdepth sample.regions.bed.gz --reference hg38 --out reports/
open reports/sample.karyotype.report.html

# Add a BAF panel beneath the genome-wide depth panel
molamola --mosdepth sample.regions.bed.gz --vcf sample.phased.vcf.gz \
         --reference hg38 --out reports/
```

Either `--vcf` or `--mosdepth` is required; `--out <directory>` is required too (molamola refuses rather than silently writing the report next to the input file). With `--vcf`, SV mode is selected from the VCF header (`##INFO=<ID=SVTYPE>`); anything else is refused rather than plotted on a guess. With `--mosdepth`, karyotype coverage mode runs and any accompanying `--vcf` is used as the BAF source. Karyotype mode expects uniform-bin mosdepth runs (`mosdepth --by <int>`).

See the [docs](https://martinandclaude.github.io/molamola/) for example output, the full CLI reference, filter explanations, and worked examples.

> **Compound-het mode was removed after v0.5.1.** The per-gene phased-haplotype panels for recessive-disease workup, and the bundled ClinVar and MANE Select references they needed, are gone; molamola is now a cytogenetics tool only. Install `molamola==0.5.1` if you need them.

## Acknowledgements

- [Sniffles2](https://github.com/fritzsedlazeck/Sniffles), [cuteSV](https://github.com/tjiangHIT/cuteSV), [SVIM](https://github.com/eldariont/svim), [pbsv](https://github.com/PacificBiosciences/pbsv), [NanoVar](https://github.com/cytham/nanovar) — long-read SV callers.
- [WhatsHap](https://github.com/whatshap/whatshap), [HiPhase](https://github.com/PacificBiosciences/HiPhase) — long-read phasing, which the BAF panel uses opportunistically when present.
- [mosdepth](https://github.com/brentp/mosdepth) — per-bin coverage for karyotype mode.
- [pyCirclize](https://github.com/moshi4/pyCirclize) — circos plot.
- [matplotlib](https://github.com/matplotlib/matplotlib), [numpy](https://github.com/numpy/numpy).
- [bcftools / samtools / htslib](https://github.com/samtools/bcftools) — VCF pre-processing helpers.
- [mosdepth](https://github.com/brentp/mosdepth) — per-bin coverage input to karyotype mode.
- [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/) — hg38 and T2T-CHM13v2.0 cytobands.
- [iconsdb.com](https://www.iconsdb.com/) — header fish icon (deep-pink, mirrored).

## License

[MIT](LICENSE).
