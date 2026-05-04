```
                 _                       _              .--.
 _ __ ___   ___ | | __ _ _ __ ___   ___ | | __ _      _/    \___
| '_ ` _ \ / _ \| |/ _` | '_ ` _ \ / _ \| |/ _` |    ( o        )
| | | | | | (_) | | (_| | | | | | | (_) | | (_| |     \___..___/
|_| |_| |_|\___/|_|\__,_|_| |_| |_|\___/|_|\__,_|         ||
```

A Python plotting tool for Oxford Nanopore variation data. **One VCF in, one self-contained HTML report out.** molamola inspects the VCF header and picks between an SV / cytogenetics report (long-read SV VCFs from Sniffles2 / cuteSV / SVIM / pbsv / NanoVar) and per-gene phased-haplotype panels (phased + VEP-annotated small-variant VCFs from WhatsHap / HiPhase). Figures embedded as base64 PNGs — no external assets, opens offline.

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
molamola --vcf sample.sniffles.vcf
open path/to/sample.report.html

# Phased + VEP-annotated VCF → compound-het workup, all candidate genes
molamola --vcf sample.phased.vep.vcf.gz
open path/to/sample.compound_het.report.html

# Just one gene from a phased + VEP VCF
molamola --vcf sample.phased.vep.vcf.gz --gene NEB
```

The plot type is auto-detected from the VCF header: `##INFO=<ID=SVTYPE>` selects SV mode; `##INFO=<ID=CSQ>` + `##FORMAT=<ID=PS>` selects compound-het mode. VCFs that match neither shape are refused with a clear error.

See the [docs](https://martinandclaude.github.io/molamola/) for example output, VEP annotation prep, the full CLI reference, filter explanations, and worked examples.

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
